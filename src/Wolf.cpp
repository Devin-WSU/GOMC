/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Wolf.h"
#include "Ewald.h"
#include "CalculateEnergy.h"
#include "EnergyTypes.h"            //Energy structs
#include "EnsemblePreprocessor.h"   //Flags
#include "BasicTypes.h"             //uint
#include "System.h"                 //For init
#include "StaticVals.h"             //For init
#include "Forcefield.h"
#include "MoleculeKind.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include "NumLib.h"
#include <cassert>
#ifdef GOMC_CUDA
#include "CalculateEwaldCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif
#include "GOMCEventsProfile.h"


Wolf::Wolf(StaticVals & stat, System & sys) :
  Ewald(stat, sys) {
    for(uint b = 0 ; b < BOX_TOTAL; b++) {
        wolfAlpha[b] = ff.wolfAlpha[b];
        wolfFactor1[b] = ff.wolfFactor1[b];
        rCutCoulomb[b] = ff.rCutCoulomb[b];
    }
}

void Wolf::Init() {
    for(uint m = 0; m < mols.count; ++m) {
        const MoleculeKind& molKind = mols.GetKind(m);
        for(uint a = 0; a < molKind.NumAtoms(); ++a) {
          particleKind.push_back(molKind.AtomKind(a));
          particleMol.push_back(m);
          particleCharge.push_back(molKind.AtomCharge(a));
          if(std::abs(molKind.AtomCharge(a)) < 0.000000001) {
              particleHasNoCharge.push_back(true);
          } else {
              particleHasNoCharge.push_back(false);
          }
        }
    }

    // initialize starting index and length index of each molecule
    startMol.resize(currentCoords.Count());
    lengthMol.resize(currentCoords.Count());

    for(int atom = 0; atom < (int) currentCoords.Count(); atom++) {
        startMol[atom] = mols.MolStart(particleMol[atom]);
        lengthMol[atom] = mols.MolLength(particleMol[atom]);
    }

    double molSelfEnergy;
    uint i, j, length;
    molSelfEnergies.resize(mols.kindsCount);

// Each thread calculates a kind
#ifdef _OPENMP
    #pragma omp parallel for default(none) private(molSelfEnergy, i, j, length) \
    shared(mols, molLookup, lambdaRef, molSelfEnergies)
#endif
    for (i = 0; i < mols.GetKindsCount(); i++) {
      MoleculeKind const& thisKind = mols.kinds[i];
      length = thisKind.NumAtoms();
      molSelfEnergy = 0.0;
      for (j = 0; j < length; j++) {
          molSelfEnergy += (thisKind.AtomCharge(j) * thisKind.AtomCharge(j));
      }
      // Store this for quick access in ChangeSelf
      molSelfEnergies[i] = molSelfEnergy;
    }

    oneThree = ff.OneThree;
    oneFour = ff.OneFour;
    scaling_14 = ff.scaling_14;
    switch(ff.wolfKind) {
      //WOLF_HYBRID_KIND
      case 0:
        isVlugtWolf = false;
        isGrossWolf = false;
        isHybridWolf = true;
        break;
      // WOLF_VLUGT_KIND
      case 1:
        isVlugtWolf = true;
        isGrossWolf = false;
        isHybridWolf = false;
        break;
      // WOLF_GROSS_KIND
      case 2:
        isVlugtWolf = false;
        isGrossWolf = true;
        isHybridWolf = false;
        break;
      default:
        std::cout << "Error ff.WolfKind has invalid value!  Check WolfKind in Config File!" << std::endl;
        exit(1);
    }
}

void Wolf::AllocMem()
{
  return;
}

void Wolf::RecipInit(uint box, BoxDimensions const& boxAxes)
{
  return;
}

//compute reciprocal term for a box with a new volume
void Wolf::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
  return;
}

//compute reciprocal term for a box when not testing a volume change
void Wolf::BoxReciprocalSums(uint box, XYZArray const& molCoords)
{
  return;
}

//calculate reciprocal term for a box
double Wolf::BoxReciprocal(uint box, bool isNewVolume) const
{
  return 0.0;
}

//calculate reciprocal force term for a box with molCoords
void Wolf::BoxForceReciprocal(XYZArray const& molCoords,
                                 XYZArray& atomForceRec, XYZArray& molForceRec,
                                 uint box)
{
  return;
}

//calculate reciprocal force term for a box
Virial Wolf::VirialReciprocal(Virial& virial, uint box) const
{
  return virial;
}


//calculate reciprocal term for displacement and rotation move
double Wolf::MolReciprocal(XYZArray const& molCoords,
                              const uint molIndex, const uint box)
{
  return 0.0;
}


//calculate reciprocal term for lambdaNew and Old with same coordinates
double Wolf::ChangeLambdaRecip(XYZArray const& molCoords, const double lambdaOld,
                                  const double lambdaNew, const uint molIndex,
                                  const uint box)
{
  return 0.0;
}


//calculate self term for a box
double Wolf::BoxSelf(uint box) const
{
    if (box >= BOXES_WITH_U_NB){
        return 0.0;
    } else {
        GOMC_EVENT_START(1, GomcProfileEvent::SELF_BOX);
        double self = 0.0;
        double molSelfEnergy;
        uint i, j, length, molNum;
        double lambdaCoef = 1.0;

// Each thread calculates a kind
#ifdef _OPENMP
    #pragma omp parallel for default(none) private(molSelfEnergy, i, j, length, molNum, lambdaCoef) \
    shared(box, mols, molLookup, lambdaRef, molSelfEnergies) \
    reduction(+:self)
#endif
        for (i = 0; i < mols.GetKindsCount(); i++) {
          MoleculeKind const& thisKind = mols.kinds[i];
          molNum = molLookup.NumKindInBox(i, box);
          molSelfEnergy = 0.0;
          if(lambdaRef.KindIsFractional(i, box)) {
              //If a molecule is fractional, we subtract the fractional molecule and
              // add it later
              --molNum;
              //returns lambda and not sqrt(lambda)
              lambdaCoef = lambdaRef.GetLambdaCoulomb(i, box);
          }
          // Store this for quick access in ChangeSelf
          molSelfEnergy = molSelfEnergies[i];
          self += (molSelfEnergy * molNum);
          if(lambdaRef.KindIsFractional(i, box)) {
              //Add the fractional molecule part
              self += (molSelfEnergy * lambdaCoef);
          }
        }
        // M_2_SQRTPI is 2/sqrt(PI), so need to multiply by 0.5 to get sqrt(PI)
        if (isVlugtWolf){
          self *= ((ff.wolfAlpha[box] * M_2_SQRTPI * 0.5) + ff.wolfFactor1[box]);
        } else {
          // we eliminate the alpha/root(pi) using Wolf,mod from Gross et al
          self *= -1.0 * ff.wolfFactor1[box] * num::qqFact;
        }

        GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_BOX);
        return self;
    }
}


//calculate correction term for a molecule
double Wolf::MolCorrection(uint molIndex, uint box) const
{
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_MOL);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0, undampenedCorr = 0.0;
  XYZ virComponents;

  MoleculeKind& thisKind = mols.kinds[mols.kIndex[molIndex]];
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);

  for (uint i = 0; i < atomSize; i++) {
    if(particleHasNoCharge[start + i]) {
      continue;
    }
    // This term only needs to be calculated once, if the molecules are rigid
    // Or of the maximum radius of gyration < RCutCoulomb
    // Otherwise, parts of the molecule may extend out of range of each other
    // For now assume the latter.
    for (uint j = i + 1; j < atomSize; j++) {
      // Need to check for cutoff for all kinds
      if(currentAxes.InRcut(distSq, virComponents, currentCoords,
                         start + i, start + j, box) && 
        distSq < rCutCoulomb[box]){
          dampenedCorr = 0.0;
          // Gross scales intra dampened by sf
          if (isGrossWolf || isVlugtWolf){
            dist = sqrt(distSq);
            // 1-3
            if (i + 2 == j && oneThree) {
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;   
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            }
            // 1-4
            if(i + 3 == j && oneFour){
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            } 
            // 1-N
            if (i + 3 < j){
              // Unscaled
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;
            }
          } 
          // Hybrid is only the negative wolfFactor1
          dampenedCorr -= wolfFactor1[box];
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_MOL);
  return num::qqFact * correction * lambdaCoef * lambdaCoef;
}


//calculate reciprocal term in destination box for swap move
double Wolf::SwapDestRecip(const cbmc::TrialMol &newMol,
                              const uint box,
                              const int molIndex)
{
  return 0.0;
}


//calculate reciprocal term in source box for swap move
double Wolf::SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                const uint box, const int molIndex)
{
  return 0.0;
}


//calculate reciprocal term for inserting some molecules (kindA) in destination
// box and removing a molecule (kindB) from destination box
double Wolf::MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                                      const std::vector<cbmc::TrialMol> &oldMol,
                                      const std::vector<uint> &molIndexNew,
                                      const std::vector<uint> &molIndexold,
                                      bool first_call)
{
  return 0.0;
}


//calculate self term after swap move
double Wolf::SwapSelf(const cbmc::TrialMol& trialMol) const
{
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::SELF_SWAP);
  MoleculeKind const& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  double en_self = 0.0;

  en_self = molSelfEnergies[thisKind.kindIndex];

  GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_SWAP);
  if (isVlugtWolf){
    return (en_self *= -1.0 * ((ff.wolfAlpha[box] * M_2_SQRTPI * 0.5) + ff.wolfFactor1[box]));
  } else {
    // we eliminate the alpha/root(pi) using Wolf,mod from Gross et al
    return (en_self *= -1.0 * ff.wolfFactor1[box] * num::qqFact);
  }
}

//calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol) const
{
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0;
  XYZ virComponents;
  const MoleculeKind& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();

  for (uint i = 0; i < atomSize; i++) {
    for (uint j = i + 1; j < atomSize; j++) {
      if(currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, j, box) && 
        distSq < rCutCoulomb[box]){
          dampenedCorr = 0.0;
          // Gross scales intra dampened by sf
          if (isGrossWolf || isVlugtWolf){
            dist = sqrt(distSq);
            // 1-3
            if (i + 2 == j && oneThree) {
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;   
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            }
            // 1-4
            if(i + 3 == j && oneFour){
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            } 
            // 1-N
            if (i + 3 < j){
              // Unscaled
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;
            }
          } 
          // Hybrid is only the negative wolfFactor1
          dampenedCorr -= wolfFactor1[box];
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return -1.0 * num::qqFact * correction;
}
//calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol,
                               const uint molIndex) const
{
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0;
  XYZ virComponents;
  const MoleculeKind& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);

  for (uint i = 0; i < atomSize; i++) {
    if(particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      if(currentAxes.InRcut(distSq, virComponents, currentCoords,
                         start + i, start + j, box) && 
        distSq < rCutCoulomb[box]){
          dampenedCorr = 0.0;
          // Gross scales intra dampened by sf
          if (isGrossWolf || isVlugtWolf){
            dist = sqrt(distSq);
            // 1-3
            if (i + 2 == j && oneThree) {
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;   
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            }
            // 1-4
            if(i + 3 == j && oneFour){
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            } 
            // 1-N
            if (i + 3 < j){
              // Unscaled
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;
            }
          } 
          // Hybrid is only the negative wolfFactor1
          dampenedCorr -= wolfFactor1[box];
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
                        dampenedCorr;
      }
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return -1.0 * num::qqFact * correction * lambdaCoef * lambdaCoef;
}

//It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void Wolf::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                         const std::vector<double> &lambda_Coul,
                         const uint iState, const uint molIndex,
                         const uint box) const
{
    uint lambdaSize = lambda_Coul.size();
    double coefDiff, en_self = 0.0;
    //Load the self energy with lambda = 1
    en_self = molSelfEnergies[molIndex];
    
    // M_2_SQRTPI is 2/sqrt(PI), so need to multiply by 0.5 to get sqrt(PI)
    //Vlugt
    //en_self *= ((ff.wolfAlpha[box] * M_2_SQRTPI * 0.5) +  ff.wolfFactor1[box] );
    // We eliminate the alpha/root(pi) using Wolf,mod
    if (isVlugtWolf){
      en_self *= -1.0 * ((ff.wolfAlpha[box] * M_2_SQRTPI * 0.5) + ff.wolfFactor1[box]);
    } else {
      // we eliminate the alpha/root(pi) using Wolf,mod from Gross et al
      en_self *= -1.0 * ff.wolfFactor1[box] * num::qqFact;
    }
    //Calculate the energy difference for each lambda state
    for (uint s = 0; s < lambdaSize; s++) {
      coefDiff = lambda_Coul[s] - lambda_Coul[iState];
      energyDiff[s].self += coefDiff * en_self;
    }
    //Calculate du/dl of self for current state, for linear scaling
    dUdL_Coul.self += en_self;
}

//It's called in free energy calculation to calculate the change in
// correction energy in all lambda states
void Wolf::ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                               const std::vector<double> &lambda_Coul,
                               const uint iState, const uint molIndex,
                               const uint box) const
{
  uint atomSize = mols.GetKind(molIndex).NumAtoms();
  uint start = mols.MolStart(molIndex);
  uint lambdaSize = lambda_Coul.size();
  double coefDiff, distSq, dist, correction = 0.0, dampenedCorr;
  XYZ virComponents;

  //Calculate the correction energy with lambda = 1
  for (uint i = 0; i < atomSize; i++) {
    if(particleHasNoCharge[start + i]) {
      continue;
    }

    for (uint j = i + 1; j < atomSize; j++) {
      distSq = 0.0;
      if(currentAxes.InRcut(distSq, virComponents, currentCoords,
                         start + i, start + j, box) && 
        distSq < rCutCoulomb[box]){
          dampenedCorr = 0.0;
          // Gross scales intra dampened by sf
          if (isGrossWolf || isVlugtWolf){
            dist = sqrt(distSq);
            // 1-3
            if (i + 2 == j && oneThree) {
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;   
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            }
            // 1-4
            if(i + 3 == j && oneFour){
              // Eq (5) Rahbari 2019, 2nd term
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
              // Vlugt doesnt scale
              if(isGrossWolf)       
                dampenedCorr *= scaling_14;
            } 
            // 1-N
            if (i + 3 < j){
              // Unscaled
              dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;
            }
          } 
          // Hybrid is only the negative wolfFactor1
          dampenedCorr -= wolfFactor1[box];
          correction += particleCharge[i + start] * particleCharge[j + start] * dampenedCorr;
      }
    }
  }
  correction *= -1.0 * num::qqFact;
  //Calculate the energy difference for each lambda state
  for (uint s = 0; s < lambdaSize; s++) {
    coefDiff = lambda_Coul[s] - lambda_Coul[iState];
    energyDiff[s].correction += coefDiff * correction;
  }
  //Calculate du/dl of correction for current state, for linear scaling
  dUdL_Coul.correction += correction;
}

//It's called in free energy calculation to calculate the change in
// reciprocal energy in all lambda states
void Wolf::ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const
{
  return;
}


//back up reciprocal value to Ref (will be called during initialization)
void Wolf::SetRecipRef(uint box)
{
  return;
}

//update reciprocal values
void Wolf::UpdateRecip(uint box)
{
  return;
}

//copy reciprocal values from ref to new
void Wolf::CopyRecip(uint box)
{
  return;
}

//update the kx, ky, kz, hsqr and prefact
void Wolf::UpdateRecipVec(uint box)
{
  return;
}

//restore cosMol and sinMol
void Wolf::RestoreMol(int molIndex)
{
  return;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Wolf::exgMolCache()
{
  return;
}

//backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Wolf::backupMolCache()
{
  return;
}

void Wolf::UpdateVectorsAndRecipTerms(bool output)
{
  return;
}
