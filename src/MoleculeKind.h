/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_MOLECULE_H
#define FF_MOLECULE_H

#include "BasicTypes.h"          //For uint
#include "EnsemblePreprocessor.h"
#include "SubdividedArray.h"
#include "Geometry.h"            //members
#include "CBMC.h"

#include <string>
#include <vector>

#include <cassert>

/* For checkpointing serialization */
// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

BOOST_CLASS_VERSION(MoleculeKind, 1);


namespace mol_setup
{
class Atom;
class MolKind;
}
namespace ff_setup
{
class Bond;
class FFBase;
}
namespace cbmc
{
class TrialMol;
}

class FFSetup;
class PRNG;
struct MolPick;
class System;
class Forcefield;
class Setup;

class MoleculeKind
{
public:

  MoleculeKind();
  ~MoleculeKind();

  uint NumAtoms() const
  {
    return numAtoms;
  }
  uint NumBonds() const
  {
    return bondList.count;
  }
  uint NumAngles() const
  {
    return angles.Count();
  }
  uint NumDihs() const
  {
    return dihedrals.Count();
  }
  uint AtomKind(const uint a) const
  {
    return atomKind[a];
  }
  double AtomCharge(const uint a) const
  {
    return atomCharge[a];
  }

  //Initialize this kind
  //Exits program if param and psf files don't match
  void Init(std::string const& l_name,
            Setup const& setup,
            Forcefield const& forcefield,
            System & sys);

  //Invoke CBMC, oldMol and newMol will be modified
  void Build(cbmc::TrialMol& oldMol, cbmc::TrialMol& newMol,
             const uint molIndex)
  {
    builder->Build(oldMol, newMol, molIndex);
  }

  //CBMC for regrowth move
  void Regrowth(cbmc::TrialMol& oldMol, cbmc::TrialMol& newMol,
                const uint molIndex)
  {
    builder->Regrowth(oldMol, newMol, molIndex);
  }

  //Crank Shaft move
  void CrankShaft(cbmc::TrialMol& oldMol, cbmc::TrialMol& newMol,
                  const uint molIndex)
  {
    builder->CrankShaft(oldMol, newMol, molIndex);
  }

  // Targeted Swap move
  void BuildGrowInCav(cbmc::TrialMol& oldMol, cbmc::TrialMol& newMol,
                      const uint molIndex)
  {
    builder->BuildGrowInCav(oldMol, newMol, molIndex);
  }

  //Used in MEMC move
  void BuildIDNew(cbmc::TrialMol& newMol, const uint molIndex)
  {
    builder->BuildIDNew(newMol, molIndex);
  }

  void BuildIDOld(cbmc::TrialMol& oldMol, const uint molIndex)
  {
    builder->BuildIDOld(oldMol, molIndex);
  }

  void BuildNew(cbmc::TrialMol& newMol, const uint molIndex)
  {
    builder->BuildNew(newMol, molIndex);
  }

  void BuildOld(cbmc::TrialMol& oldMol, const uint molIndex)
  {
    builder->BuildOld(oldMol, molIndex);
  }

  void BuildGrowNew(cbmc::TrialMol& newMol, const uint molIndex)
  {
    builder->BuildGrowNew(newMol, molIndex);
  }

  void BuildGrowOld(cbmc::TrialMol& oldMol, const uint molIndex)
  {
    builder->BuildGrowOld(oldMol, molIndex);
  }

  double GetMoleculeCharge();

  bool MoleculeHasCharge();

  SortedNonbond sortedNB, sortedNB_1_4, sortedNB_1_3, sortedEwaldNB;


  //these are used for total energy calculations, see Geometry.h/cpp
  Nonbond nonBonded;
  Nonbond_1_4 nonBonded_1_4;
  Nonbond_1_3 nonBonded_1_3;
  EwaldNonbond nonEwaldBonded;

  BondList bondList;
  GeomFeature angles;
  GeomFeature dihedrals;

// Straightforward to serialize 

  bool oneThree, oneFour;

  std::string name;
  std::vector<std::string> atomNames, atomTypeNames, resNames;
  double molMass;

  double * atomMass;

  bool isMultiResidue;
  std::vector<uint> intraMoleculeResIDs;

#if ENSEMBLE == GCMC
  double chemPot;
#endif
  cbmc::CBMC* builder;

private:
/* Boost supports serializing polymorphic pointers without any special handling by the user */

  uint numAtoms;
  uint * atomKind;
  double * atomCharge;

  void InitAtoms(mol_setup::MolKind const& molData);

  //uses buildBonds to check if molecule is branched
  //bool CheckBranches();
  void InitCBMC(System& sys, Forcefield& ff,
                Setup& set);

  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {

    ar & numAtoms;
    ar & molMass;
    ar & name;
    ar & sortedNB;
    ar & sortedNB_1_4; 
    ar & sortedNB_1_3; 
    ar & sortedEwaldNB;
    ar & nonBonded;
    ar & nonBonded_1_4;
    ar & nonBonded_1_3;
    ar & nonEwaldBonded;
    ar & bondList;
    ar & angles;
    ar & dihedrals;

    ar & oneThree;
    ar & oneFour;
    ar & name;
    ar & atomNames;
    ar & atomTypeNames;
    ar & resNames;
    if (Archive::is_loading::value)
    {
        assert(atomMass == nullptr);
        atomMass = new double[numAtoms];
        assert(atomKind == nullptr);
        atomKind = new uint[numAtoms];
        assert(atomCharge == nullptr);
        atomCharge = new double[numAtoms];        
    }
    ar & boost::serialization::make_array<double>(atomMass, numAtoms);  
    ar & boost::serialization::make_array<uint>(atomKind, numAtoms);  
    ar & boost::serialization::make_array<double>(atomCharge, numAtoms);  
    
    
  }
};






#endif /*FF_MOLECULE_H*/
