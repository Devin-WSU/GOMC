/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "BasicTypes.h"
#include "SubdividedArray.h"
#include <assert.h>

#include <vector>

#include <boost/serialization/array.hpp>

namespace mol_setup
{
class Atom;
class Bond;
class Angle;
class Dihedral;
class MolKind;
}
namespace ff_setup
{
class Bond;
class FFBase;
}
class FFSetup;

// for 1-5 and more interaction
struct Nonbond {
  uint* part1;
  uint* part2;
  uint count;

  virtual void Init(const mol_setup::MolKind& molData);
  Nonbond();
  ~Nonbond();

  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & count;
      if (Archive::is_loading::value)
      {
          assert(part1 == nullptr);
          part1 = new uint[count];
          assert(part2 == nullptr);
          part2 = new uint[count];
      }
      ar & boost::serialization::make_array<uint>(part1, count);  
      ar & boost::serialization::make_array<uint>(part2, count);  
  }
};
BOOST_CLASS_VERSION(Nonbond, 1);

// for 1-4  interaction
struct Nonbond_1_4 : public Nonbond {
  virtual void Init(const mol_setup::MolKind& molData);
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & boost::serialization::base_object<Nonbond>(*this);
      /*
      ar & count;
      if (Archive::is_loading::value)
      {
          assert(part1 == nullptr);
          part1 = new uint[count];
          assert(part2 == nullptr);
          part2 = new uint[count];
      }
      ar & boost::serialization::make_array<uint>(part1, count);  
      ar & boost::serialization::make_array<uint>(part2, count);  
      */
  }
};
BOOST_CLASS_VERSION(Nonbond_1_4, 1);

// for 1-3  interaction, used for Martini ForceField
struct Nonbond_1_3 : public Nonbond {
  virtual void Init(const mol_setup::MolKind& molData);
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & boost::serialization::base_object<Nonbond>(*this);
      /*
      ar & count;
      if (Archive::is_loading::value)
      {
          assert(part1 == nullptr);
          part1 = new uint[count];
          assert(part2 == nullptr);
          part2 = new uint[count];
      }
      ar & boost::serialization::make_array<uint>(part1, count);  
      ar & boost::serialization::make_array<uint>(part2, count);  
      */
  }
};
BOOST_CLASS_VERSION(Nonbond_1_3, 1);


// for ewald correction energy calculation
struct EwaldNonbond : public Nonbond {
  virtual void Init(const mol_setup::MolKind& molData);
  private:
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & boost::serialization::base_object<Nonbond>(*this);
      /*
      ar & count;
      if (Archive::is_loading::value)
      {
          assert(part1 == nullptr);
          part1 = new uint[count];
          assert(part2 == nullptr);
          part2 = new uint[count];
      }
      ar & boost::serialization::make_array<uint>(part1, count);  
      ar & boost::serialization::make_array<uint>(part2, count);
      */  
  }
};
BOOST_CLASS_VERSION(EwaldNonbond, 1);


//!List of all pairs of particles in bonds.
struct BondList {
  uint* part1;
  uint* part2;
  uint* kinds;
  uint count;

  void Init(const std::vector<mol_setup::Bond>& bonds);
  bool IsBonded(const uint &i, const uint &j);
  std::vector<uint> GetBondIndices(uint atom) const;

  BondList();
  ~BondList();

  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & count;
      if (Archive::is_loading::value)
      {
          assert(part1 == nullptr);
          part1 = new uint[count];
          assert(part2 == nullptr);
          part2 = new uint[count];
          assert(kinds == nullptr);
          kinds = new uint[count];
      }
      ar & boost::serialization::make_array<uint>(part1, count);  
      ar & boost::serialization::make_array<uint>(part2, count);  
      ar & boost::serialization::make_array<uint>(kinds, count);  
  }

};
BOOST_CLASS_VERSION(BondList, 1);


//!List describing a geometry feature of a molecule by the bonds it contains.
/*!Contains bond indices as used in bondlist, bonds greater than bondList.count
 * represent bonds that need to be flipped to describe the feature in the same
 * order as the PSF file (head to tail). Bonds are used to avoid recalcuating
 * vectors
 */
class GeomFeature
{
public:

  //!Prepares topology data for a geometry feature
  void Init(const std::vector<mol_setup::Angle>& angles, const BondList& bList);
  void Init(const std::vector<mol_setup::Dihedral>& dihs, const BondList& bList);
  explicit GeomFeature(uint atomsPer);
  ~GeomFeature();

  std::vector<uint> GetMidAnglesIndices(uint focus) const;
  std::vector<uint> GetDihsOnBondIndices(uint atom, uint partner) const;


  //Return number of features in molecule
  uint Count() const
  {
    return count;
  }


  //!return bondIndex of a bond in a feature
  //!/param feature index of angle/dihedral; in [0, Count)
  //!/param b index of bond within feature; in [0, bondsPer)
  uint GetBond(uint feature, uint b) const
  {
    assert((feature * bondsPer + b) < (count * bondsPer));
    return bondIndices[feature * bondsPer + b];
  }

  //!returns kind index of feature i
  uint GetKind(uint i) const
  {
    return kinds[i];
  }


private:
  //number of bonds per feature e.g. dihedral = 3, angle = 2 etc
  uint bondsPer;
  //bond indices, size = count * bondsPer
  //bonds needing reversal are represented as index + bondCount
  uint* bondIndices;
  //feature kind indices, size = count
  uint* kinds;
  //number of features. If count == 0, do not attempt to deref the arrays
  uint count;
  uint bondIndicesCount;

  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & bondsPer;
    ar & count;
    ar & bondIndicesCount;
    if (Archive::is_loading::value)
    {
        assert(bondIndices == nullptr);
        bondIndices = new uint[bondIndicesCount];
        assert(kinds == nullptr);
        kinds = new uint[count];
    }
    ar & boost::serialization::make_array<uint>(bondIndices, bondIndicesCount);  
    ar & boost::serialization::make_array<uint>(kinds, count);  
  }
};
BOOST_CLASS_VERSION(GeomFeature, 1);


//List of all pairs that interact via nonbonded LJ potentials, i.e. all pairs
//that are at least 4 bonds apart.

//List of all pairs that interact via nonbonded LJ potentials, i.e. all pairs
//that are at least 4 bonds apart.
//Sorted (redundantly) into a series of lists of all partners of each particle
class SortedNonbond
{
public:
  //!Returns iterable pointer to first pair containing atom
  const uint* Begin(uint atom) const
  {
    return partners + subdiv.Begin(atom);
  }
  //!Returns pointer to position after the last pair containing atom
  const uint* End(uint atom) const
  {
    return partners + subdiv.End(atom);
  }

  //!Initialize from a nonbond structure
  void Init(const Nonbond& nb, uint numAtoms);

  SortedNonbond() : partners(NULL) {}
  ~SortedNonbond()
  {
    delete[] partners;
  }

private:
  uint* partners;
  SubdividedArray subdiv;
  uint partnersSize;
  
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & partnersSize;
      ar & subdiv;
      if (Archive::is_loading::value)
      {
          assert(partners == nullptr);
          partners = new uint[partnersSize];
      }
      ar & boost::serialization::make_array<uint>(partners, partnersSize);  
  }
};
BOOST_CLASS_VERSION(SortedNonbond, 1)


#endif
