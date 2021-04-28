/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PSFOUTPUT_H
#define PSFOUTPUT_H

#include <vector>
#include <string>

#include "BasicTypes.h"
#include "MolSetup.h"
#include "System.h"
#include "MoleculeLookup.h"
#include "OutputAbstracts.h"
#include "Writer.h"

class Molecules;

class PSFOutput : public OutputableBase {

public:
  PSFOutput(const Molecules& molecules, const System &sys,
            Setup & set);

  //PSF does not need to sample on every step, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const& atoms, config_setup::Output const& output);

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);
  //virtual void Output(const ulong step);

private:
  const Molecules* molecules;
  const MoleculeLookup & molLookRef;
  std::vector<mol_setup::MolKind> molKinds;
  std::vector<std::string> molNames;
  uint totalAtoms;
  uint totalBonds;
  uint totalAngles;
  uint totalDihs;

  uint boxAtoms[BOX_TOTAL];
  uint boxBonds[BOX_TOTAL];
  uint boxAngles[BOX_TOTAL];
  uint boxDihs[BOX_TOTAL];

  void PrintRemarks(FILE* outfile, const std::vector<std::string>& remarks) const;
  void PrintAtoms(FILE* outfile) const;
  void PrintBonds(FILE* outfile) const;
  void PrintAngles(FILE* outfile) const;
  void PrintDihedrals(FILE* outfile) const;

  void PrintRemarksInBox(FILE* outfile, uint b) const;
  void PrintAtomsInBox(FILE* outfile, uint b) const;
  void PrintBondsInBox(FILE* outfile, uint b) const;
  void PrintAnglesInBox(FILE* outfile, uint b) const;
  void PrintDihedralsInBox(FILE* outfile, uint b) const;

  void PrintNAMDCompliantSuffix(FILE* outfile) const;
  void PrintNAMDCompliantSuffixInBox(FILE* outfile) const;

  void CountMolecules();
  void CountMoleculesInBoxes();

  //NEW_RESTART_CODE
  FILE * outRebuildRestart[BOX_TOTAL];
  std::string outFName;
  std::string outRebuildRestartFName[BOX_TOTAL];
};

#endif
