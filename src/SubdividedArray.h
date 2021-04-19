/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SUBDIV_ARRAY
#define SUBDIV_ARRAY

#include <cstddef>
#include "BasicTypes.h"

/* For checkpointing serialization */
// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/array.hpp>
//Common class used for dihedral, sorted kind array, topology arrays, etc.

class SubdividedArray
{
public:
  SubdividedArray(): start(NULL),  subdivCount(0) {}
  SubdividedArray(SubdividedArray const& other)
  {
    subdivCount = other.subdivCount;
    startCount = other.subdivCount + 1;
    start = new uint[startCount];
    for(uint i = 0; i <= subdivCount; ++i)
      start[i] = other.start[i];
  }
  void Init(const uint subdiv)
  {
    if (start != NULL) Cleanup();
    subdivCount = subdiv;
    startCount = subdiv + 1;
    start = new uint[startCount];
  }
  void Set(const uint div, const uint first, const uint len)
  {
    start[div] = first;
    start[div + 1] = first + len;
  }

  ~SubdividedArray(void)
  {
    Cleanup();
  }
  void Cleanup(void)
  {
    delete[] start;
    start = NULL;
  }

  //returns index of the offset-th element of kind
  uint Index(const uint kind, const uint offset) const
  {
    return start[kind] + offset;
  }

  //Return first el, last el, or length of subdiv
  uint Begin(const uint kind) const
  {
    return start[kind];
  }
  uint End(const uint kind) const
  {
    return start[kind + 1];
  }
  uint Length(const uint kind) const
  {
    return End(kind) - Begin(kind);
  }

  SubdividedArray& operator=(SubdividedArray other)
  {
    subdivCount = other.subdivCount;
    startCount = other.startCount;
    unsigned int* tmp = other.start;
    other.start = start;
    start = tmp;
    return *this;
  }

private:
  //note start is one longer than subdivCount to account for length of last
  //actual element
  uint * start;
  uint startCount, subdivCount;

  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
      ar & startCount;
      ar & subdivCount;
      if (Archive::is_loading::value)
      {
          assert(start == nullptr);
          start = new uint[startCount];
      }
      ar & boost::serialization::make_array<uint>(start, startCount);  
  }
};

#endif /*SUBDIV_ARRAY*/
