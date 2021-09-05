/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Wolf.h"
#include "GOMCEventsProfile.h"


Wolf::Wolf(StaticVals & stat, System & sys) :
  Ewald(stat, sys) {}

void Wolf::Init() {}

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

        for (i = 0; i < mols.GetKindsCount(); i++) {
        MoleculeKind const& thisKind = mols.kinds[i];
        length = thisKind.NumAtoms();
        molNum = molLookup.NumKindInBox(i, box);
        molSelfEnergy = 0.0;
        if(lambdaRef.KindIsFractional(i, box)) {
            //If a molecule is fractional, we subtract the fractional molecule and
            // add it later
            --molNum;
            //returns lambda and not sqrt(lambda)
            lambdaCoef = lambdaRef.GetLambdaCoulomb(i, box);
        }

        for (j = 0; j < length; j++) {
            molSelfEnergy += (thisKind.AtomCharge(j) * thisKind.AtomCharge(j));
        }
        self += (molSelfEnergy * molNum);
        if(lambdaRef.KindIsFractional(i, box)) {
            //Add the fractional molecule part
            self += (molSelfEnergy * lambdaCoef);
        }
        }
        // M_2_SQRTPI is 2/sqrt(PI), so need to multiply by 0.5 to get sqrt(PI)
        self *= ((ff.wolfAlpha[box] * M_2_SQRTPI * 0.5) +  0.5 * ff.wolfFactor1[box]);
        self *= -1.0 * num::qqFact;

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
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, currentCoords,
                         start + i, start + j, box);
      dist = sqrt(distSq);
      // Eq (5) Gezelter 2006, 2nd term
      dampenedCorr = erfc(forcefield.wolfAlpha[box] * dist)/dist;
      dampenedCorr -= forcefield.wolfFactor1[box];
      // Eq (5) Gezelter 2006, 3rd term
      undampenedCorr = 1.0/dist;
      correction += (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * (dampenedCorr - undampenedCorr);
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
  return 0.0;
}

//calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol) const
{
  return 0.0;
}
//calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol,
                               const uint molIndex) const
{
  return 0.0;
}

//It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void Wolf::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                         const std::vector<double> &lambda_Coul,
                         const uint iState, const uint molIndex,
                         const uint box) const
{
    uint atomSize = mols.GetKind(molIndex).NumAtoms();
    uint start = mols.MolStart(molIndex);
    uint lambdaSize = lambda_Coul.size();
    double coefDiff, en_self = 0.0;
    //Calculate the self energy with lambda = 1
    for (uint i = 0; i < atomSize; i++) {
      en_self += (particleCharge[i + start] * particleCharge[i + start]);
    }
    // M_2_SQRTPI is 2/sqrt(PI), so need to multiply by 0.5 to get sqrt(PI)
    en_self *= ((ff.wolfAlpha[box] * M_2_SQRTPI * 0.5) +  0.5 * ff.wolfFactor1[box]);
    en_self *= -1.0 * num::qqFact;

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
  return;
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
