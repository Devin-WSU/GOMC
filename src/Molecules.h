/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOLECULES_H
#define MOLECULES_H

#include "BasicTypes.h" //For uint
#include "MolSetup.h"
#include <map>
#include <string>

namespace pdb_setup
{
class Atoms;
}
class FFSetup;
class Forcefield;
class System;

#include "MoleculeKind.h" //For member var.

//Note: This info is static and will never change in current ensembles
class Molecules
{
public:
  Molecules();
  ~Molecules();

  const MoleculeKind& GetKind(const uint molIndex) const
  {
    return kinds[kIndex[molIndex]];
  }

  uint GetMolKind(const uint molIndex) const
  {
    return kIndex[molIndex];
  }

  void Init(Setup& setup, Forcefield& forcefield,
            System& sys);

  uint NumAtomsByMol(const uint m) const
  {
    return start[m + 1] - start[m];
  }
  uint NumAtoms(const uint mk) const
  {
    return kinds[mk].NumAtoms();
  }

  int MolStart(const uint molIndex) const
  {
    return start[molIndex];
  }

  int MolEnd(const uint molIndex) const
  {
    return start[molIndex + 1];
  }

  int MolLength(const uint molIndex) const
  {
    return MolEnd(molIndex) - MolStart(molIndex);
  }

  void GetRange(uint & _start, uint & stop, uint & len, const uint m) const
  {
    _start = start[m];
    stop = start[m + 1];
    len = stop - _start;
  }

  void GetRangeStartStop(uint & _start, uint & stop, const uint m) const
  {
    _start = start[m];
    stop = start[m + 1];
  }
  void GetRangeStartLength(uint & _start, uint & len, const uint m) const
  {
    _start = start[m];
    len = start[m + 1] - _start;
  }

  uint GetKindsCount() const
  {
    return kindsCount;
  }

  void PrintLJInfo(std::vector<uint> &totAtomKind,
                   std::vector<std::string> &names,
                   Forcefield & forcefield);

  //private:
  //Kind index of each molecule and start in master particle array
  //Plus counts
  uint* start;
  uint count;
  uint* kIndex;
  uint kIndexCount;
  uint* countByKind;
  uint chainCount;
  char* chain;

  MoleculeKind * kinds;

  uint kindsCount;
  uint fractionKind, lambdaSize;
  double* pairEnCorrections;
  double* pairVirCorrections;

  bool printFlag;

private:
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & count;
    ar & chainCount;
    ar & kIndexCount;
    ar & kindsCount;
    ar & fractionKind;
    ar & lambdaSize; 
    ar & printFlag; 

    if (Archive::is_loading::value)
    {
        assert(start == nullptr);
        start = new uint [count + 1];
        assert(kIndex == nullptr);
        kIndex = new uint[count];
        assert(countByKind == nullptr);
        countByKind = new uint[kindsCount];     
        assert(chain == nullptr);
        chain = new char[chainCount];   
        assert(pairEnCorrections == nullptr);
        pairEnCorrections = new double[kindsCount * kindsCount];    
        assert(pairVirCorrections == nullptr);
        pairVirCorrections = new double[kindsCount * kindsCount];   
    }
    ar & boost::serialization::make_array<uint>(start, count + 1);  
    ar & boost::serialization::make_array<uint>(kIndex, count);  
    ar & boost::serialization::make_array<uint>(countByKind, kindsCount);  
    ar & boost::serialization::make_array<char>(chain, chainCount);  
    ar & boost::serialization::make_array<double>(pairEnCorrections, kindsCount * kindsCount);  
    ar & boost::serialization::make_array<double>(pairVirCorrections, kindsCount * kindsCount);  
  }

};


#endif /*MOLECULES_H*/
