#ifndef DCFREEHEDRON_H
#define DCFREEHEDRON_H
#include "DCComponent.h"
#include "DCSingle.h"
#include "DCHedron.h"
#include "../CBMC.h"

namespace mol_setup { struct MolKind; }

namespace cbmc {
   class DCData;

   class DCFreeHedron : public DCComponent
   {  
   public:
      DCFreeHedron(DCData* data, const mol_setup::MolKind& kind,
		   uint focus, uint prev);
      void PrepareNew(TrialMol& newMol, uint molIndex);
      void PrepareOld(TrialMol& oldMol, uint molIndex);
      void BuildOld(TrialMol& oldMol, uint molIndex);
      void BuildNew(TrialMol& newMol, uint molIndex);
      DCComponent* Clone() { return new DCFreeHedron(*this); };

   private:
      DCData* data;
      DCSingle seed;
      DCHedron hed;
      double anchorBond;
      //anchor bond energy of old molecule
      double oldBondEnergy;
      uint anchorBondKind;
   };
}

#endif
