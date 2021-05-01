/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PDB_SETUP_H
#define PDB_SETUP_H

#include <vector>
#include <map> //for function lookup table.

#include "InputAbstracts.h" //For FWReadableBase
#include "BasicTypes.h" //For uint
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "PDBConst.h" //For fields positions, etc.
#include "XYZArray.h" //For box dimensions.

// Recalc Traj from Binary
#include "DCDlib.h" // for Error output

namespace config_setup
{
struct RestartSettings;
}
struct FWReadableBase;


namespace pdb_setup
{
struct Remarks : FWReadableBase {
  uint currBox;
  double disp[BOX_TOTAL], rotate[BOX_TOTAL], vol[BOX_TOTAL];
  ulong step[BOX_TOTAL];
  uint frameNumber[BOX_TOTAL], targetFrame[BOX_TOTAL];
  std::vector<ulong> frameSteps;
  bool restart, reached[BOX_TOTAL], recalcTrajectory;
  bool restartFromXSC, restartFromBinary;
  Remarks () {
    currBox = 0;
    restart = recalcTrajectory = false;
    restartFromXSC = restartFromBinary = false;
    for (int b = 0; b < BOX_TOTAL; b++) {
      reached[b] = false;
      frameNumber[b] = targetFrame[b] = 0;
      //set to -1, so in Movesetting, we initialize it
      disp[b] = rotate[b] = vol[b] = -1.0; 
      step[b] = 0;
    }
  }

  void SetRestart(config_setup::RestartSettings const& r);
  void Read(FixedWidthReader & pdb);
  void SetBox(const uint b)
  {
    currBox = b;
  }
  void SetFrameNumber(const uint b, const uint frameNum)
  {
    targetFrame[b] = frameNum;
  }
  void Clear();

private:
  void CheckGOMC(std::string const& varName);
};

struct Cryst1 : FWReadableBase {
  //box dimensions
  uint currBox;
  bool hasVolume[BOX_TOTAL];    //If reads cellBasis info from PDB
  bool hasCellBasis[BOX_TOTAL]; //If reads cellBasis info from XSC
  XYZArray axis, cellBasis[BOX_TOTAL];
  double cellAngle[BOX_TOTAL][3];
  Cryst1(void) : currBox(0), axis(BOX_TOTAL) 
  {
    for(int b = 0; b < BOX_TOTAL; b++) {
      hasCellBasis[b] = false;
      hasVolume[b] = false;
      cellBasis[b] = XYZArray(3);
    }
  }
  void SetBox(const uint b)
  {
    currBox = b;
  }
  void Read(FixedWidthReader & pdb);
};

class Atoms : public FWReadableBase
{
public:
  //Set the current residue to something other than 1
  Atoms(void) : restart(false), currBox(0), count(0)
    {
      for (int b = 0; b < BOX_TOTAL; b++) {
        numAtomsInBox[b] = 0;
      }
    }
  void SetRestart(config_setup::RestartSettings const& r);
  void SetBox(const uint b)
  {
    currBox = b;
  }
  void Assign(std::string const& resName,
              const char l_chain,
              const double l_x,
              const double l_y,
              const double l_z,
              const double l_beta);

  void Read(FixedWidthReader & file);
  void Clear();

  //private:
  //member data
  std::vector<char> chainLetter; //chain ids of each atom respectively
  std::vector<double> x, y, z; //coordinates of each particle
  std::vector<double> beta;  //beta value of each molecule
  std::vector<uint> box;
  std::vector<std::string> resNames;
  bool restart, firstResInFile, recalcTrajectory, restartFromBinary;
  //CurrRes is used to store res vals, currBox is used to
  //determine box either via the file (new) or the occupancy
  //(restart), count allows overwriting of coordinates during
  //second box read (restart only)
  uint currBox, count;
  uint numAtomsInBox[BOX_TOTAL]; // number of atom in each box
};

}

struct PDBSetup {
  pdb_setup::Atoms atoms;
  pdb_setup::Cryst1 cryst;
  pdb_setup::Remarks remarks;
  FixedWidthReader pdb[BOX_TOTAL];
  PDBSetup(void) : dataKinds(SetReadFunctions()) {}
  void Init(config_setup::RestartSettings const& restart,
            std::string const*const name, uint frameNumber = 1);
  std::vector<ulong> GetFrameSteps(std::string const*const name);
  std::vector<ulong> GetFrameStepsFromBinary(std::string const*const name, uint * numAtomsInBox);

private:
  //Map variable names to functions
  std::map<std::string, FWReadableBase *>  SetReadFunctions(void)
  {
    std::map<std::string, FWReadableBase *> funct;
    funct[pdb_entry::label::REMARK] = &remarks;
    funct[pdb_entry::label::CRYST1] = &cryst;
    funct[pdb_entry::label::ATOM] = &atoms;
    funct[pdb_entry::label::HETATM] = &atoms;
    return funct;
  }
  const std::map<std::string, FWReadableBase *> dataKinds;
  static const std::string pdbAlias[];
};

#endif /*PDB_SETUP_H*/
