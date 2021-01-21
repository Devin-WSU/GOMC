/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "CBMC.h"
#include "MoleculeKind.h"
#include "MolSetup.h"
#include "DCLinear.h"
#include "DCGraph.h"
#include "DCCyclic.h"
#include <vector>


namespace cbmc
{

CBMC* MakeCBMC(System& sys, const Forcefield& ff,
               const MoleculeKind& kind, const Setup& set)
{
// We used to do boundCount here for no apparent reason unless scope somehow
// extended to these subclasses

  bool cyclic = (kind.NumBonds() > kind.NumAtoms() - 1) ? true : false;

  if(cyclic) {
    return new DCCyclic(sys, ff, kind, set);
  } else if (kind.NumAtoms() > 2) {
    //Any molecule woth 3 atoms and more will be built in DCGraph
    return new DCGraph(sys, ff, kind, set);
  } else {
    return new DCLinear(sys, ff, kind, set);
  }
}

}
