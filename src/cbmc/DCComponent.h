/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCCOMPONENT_H
#define DCCOMPONENT_H
#include "BasicTypes.h"
#include "MoleculeKind.h"
//Class for Deferred Coupling CBMC components

namespace cbmc
{
class TrialMol;

class DCComponent
{
public:
  //Perform Decoupled portions of CBMC
  virtual void PrepareNew(TrialMol& newMol, uint molIndex) {}
  virtual void PrepareOld(TrialMol& oldMol, uint molIndex) {}

  //Perform Coupled final build
  virtual void BuildOld(TrialMol& oldMol, uint molIndex) = 0;
  virtual void BuildNew(TrialMol& newMol, uint molIndex) = 0;

  virtual void UpdateAcceptance(const TrialMol& mol) {}
  virtual ~DCComponent() {};

  //List of bonds with atom at one end, atom first
  uint GetBondKind(const MoleculeKind& molKind, uint atom, uint focus){
    uint result;
    for (int i = 0; i < molKind.bondList.count; ++i){
      /* H - O - C - C */
      /* a ^ f         */
      if ((molKind.bondList.part1[i] == atom )|| (molKind.bondList.part2[i] == atom)){
        if((molKind.bondList.part1[i] == focus )|| (molKind.bondList.part2[i] == focus)){
          result = molKind.bondList.kinds[i];
          break;
        }
      }
    }
    return result;
  }

};
}

#endif
