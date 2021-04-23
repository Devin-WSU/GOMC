/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCDATA_H
#define DCDATA_H
#include "BasicTypes.h"
#include "XYZArray.h"
#include "Setup.h"
#include "System.h"
#include "CBMC.h"
#include <vector>
#include <algorithm>

/* For checkpointing serialization */
// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

class Forcefield;


namespace cbmc
{
//Class to avoid reallocating arrays for CBMC
//Could be refactored into an object pool. This would be easier if I had bothered to write accessors
class DCData
{
public:
  explicit  DCData(System& sys, const Forcefield& forcefield,
                   const Setup& set);
  ~DCData();

  const CalculateEnergy& calc;

  const Ewald  *calcEwald;

  const Forcefield& ff;
  const BoxDimensions& axes;
  PRNG& prng;

  const uint nAngleTrials;
  const uint nDihTrials;
  const uint nLJTrialsFirst;
  const uint nLJTrialsNth;

  XYZArray& positions;     //candidate positions for inclusion (alias for multiPositions[0])

  uint maxLJTrials;
  /* Arrays of size maxLJTrials */
  double* inter;          //intermolecule energies, reused for new and old
  double* real;           //short range coulomb interaction
  double* bonded;
  double* oneFour;
  double* nonbonded;      //calculated nonbonded 1_N LJ and coulomb energies
  double* ljWeights;
  bool* overlap;      //For detecting overlap for each LJ trial
  /* Arrays of size maxLJTrials */


  uint totalTrials;
  /* Arrays of size totalTrials */
  double* interT;     //For DCRotateCOM, we have combined first and Nth trial
  double* realT;      //For DCRotateCOM, we have combined first and Nth trial
  double* ljWeightsT; //For DCRotateCOM, we have combined first and Nth trial
  bool* overlapT;     //For detecting overlap for each LJ trial. Used in DCRotateCOM
  /* Arrays of size totalTrials */

  uint trialMax;
  /* Arrays of size trialMax */
  //used for both angles and dihedrals
  double* angles;
  double* angleWeights;
  double* angleEnergy;
  double* nonbonded_1_4;  //calculated nonbonded 1_4 LJ and coulomb energies
  double* nonbonded_1_3;  //calculated nonbonded 1_3 LJ and coulomb energies
  /* Arrays of size trialMax */

  
  XYZArray multiPositions[MAX_BONDS];

private:
  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar &   nAngleTrials;
    ar &   nDihTrials;
    ar &   nLJTrialsFirst;
    ar &   nLJTrialsNth;

    ar &   maxLJTrials;
    ar &   totalTrials;
    ar &   trialMax;
    if (Archive::is_loading::value)
    {
        assert(inter == nullptr);
        inter = new double[maxLJTrials];
        assert(real == nullptr);
        real = new double[maxLJTrials];
        assert(bonded == nullptr);
        bonded = new double[maxLJTrials];        
        assert(oneFour == nullptr);
        oneFour = new double[maxLJTrials];
        assert(nonbonded == nullptr);
        nonbonded = new double[maxLJTrials];
        assert(ljWeights == nullptr);
        ljWeights = new double[maxLJTrials]; 
        assert(overlap == nullptr);
        overlap = new bool[maxLJTrials];

        assert(interT == nullptr);
        interT = new double[totalTrials];
        assert(realT == nullptr);
        realT = new double[totalTrials];        
        assert(ljWeightsT == nullptr);
        ljWeightsT = new double[totalTrials];
        assert(overlapT == nullptr);
        overlapT = new bool[totalTrials];

        assert(angles == nullptr);
        angles = new double[trialMax];
        assert(angleWeights == nullptr);
        angleWeights = new double[trialMax];
        assert(angleEnergy == nullptr);
        angleEnergy = new double[trialMax];        
        assert(nonbonded_1_4 == nullptr);
        nonbonded_1_4 = new double[trialMax];
        assert(nonbonded_1_3 == nullptr);
        nonbonded_1_3 = new double[trialMax];

        assert(multiPositions == nullptr);
        multiPositions = new XYZArray[MAX_BONDS];        
    }
    ar & boost::serialization::make_array<double>(inter, maxLJTrials);  
    ar & boost::serialization::make_array<double>(real, maxLJTrials);  
    ar & boost::serialization::make_array<double>(bonded, maxLJTrials);
    ar & boost::serialization::make_array<double>(oneFour, maxLJTrials);  
    ar & boost::serialization::make_array<double>(nonbonded, maxLJTrials);  
    ar & boost::serialization::make_array<double>(ljWeights, maxLJTrials);
    ar & boost::serialization::make_array<bool>(overlap, maxLJTrials); 

    ar & boost::serialization::make_array<double>(interT, totalTrials);  
    ar & boost::serialization::make_array<double>(realT, totalTrials);
    ar & boost::serialization::make_array<double>(ljWeightsT, totalTrials);  
    ar & boost::serialization::make_array<bool>(overlapT, totalTrials);

    ar & boost::serialization::make_array<double>(angles, trialMax);
    ar & boost::serialization::make_array<double>(angleWeights, trialMax);  
    ar & boost::serialization::make_array<double>(angleEnergy, trialMax);  
    ar & boost::serialization::make_array<double>(nonbonded_1_4, trialMax);
    ar & boost::serialization::make_array<double>(nonbonded_1_3, trialMax);

    ar & boost::serialization::make_array<XYZArray>(multiPositions, MAX_BONDS);
  }

};

inline DCData::DCData(System& sys, const Forcefield& forcefield, const Setup& set):

  calc(sys.calcEnergy), ff(forcefield),
  axes(sys.boxDimRef), prng(sys.prng),
  nAngleTrials(set.config.sys.cbmcTrials.bonded.ang),
  nDihTrials(set.config.sys.cbmcTrials.bonded.dih),
  nLJTrialsFirst(set.config.sys.cbmcTrials.nonbonded.first),
  nLJTrialsNth(set.config.sys.cbmcTrials.nonbonded.nth),
  positions(*multiPositions)
{
  calcEwald = sys.GetEwald();
  uint maxLJTrials = nLJTrialsFirst;
  if ( nLJTrialsNth > nLJTrialsFirst )
    maxLJTrials = nLJTrialsNth;

  totalTrials = nLJTrialsFirst * nLJTrialsNth;
  if(totalTrials == 0)
    totalTrials = maxLJTrials;

  for(uint i = 0; i < MAX_BONDS; ++i) {
    multiPositions[i] = XYZArray(maxLJTrials);
  }
  inter = new double[maxLJTrials];
  real = new double[maxLJTrials];
  bonded = new double[maxLJTrials];
  oneFour = new double[maxLJTrials];
  nonbonded = new double[maxLJTrials];
  ljWeights = new double[maxLJTrials];
  overlap = new bool[maxLJTrials];

  interT = new double[totalTrials];
  realT = new double[totalTrials];
  ljWeightsT = new double[totalTrials];
  overlapT = new bool[totalTrials];

  trialMax = std::max(nAngleTrials, nDihTrials);
  angleEnergy = new double[trialMax];
  angleWeights = new double[trialMax];
  angles = new double[trialMax];
  nonbonded_1_3 = new double[trialMax];
  nonbonded_1_4 = new double[trialMax];
}

inline DCData::~DCData()
{
  delete[] inter;
  delete[] real;
  delete[] bonded;
  delete[] oneFour;
  delete[] nonbonded;
  delete[] nonbonded_1_4;
  delete[] nonbonded_1_3;
  delete[] ljWeights;
  delete[] angles;
  delete[] angleWeights;
  delete[] angleEnergy;
  delete[] interT;
  delete[] realT;
  delete[] ljWeightsT;
  delete[] overlap;
  delete[] overlapT;
}

}

#endif
