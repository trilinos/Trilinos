#include "Petra_LocalMap.h"


//============================================================================
Petra_LocalMap::Petra_LocalMap(int NumMyElements, int IndexBase, 
			       const Petra_Comm& Comm)
: Petra_Map(NumMyElements, NumMyElements, IndexBase, Comm) {
  DistributedGlobal_ = false;
  if (CheckInput()!=0) {
    cout << "Replicated Local Map not the same size on all PEs" << endl;
    abort();
  }
}
//============================================================================
Petra_LocalMap::Petra_LocalMap(const Petra_LocalMap& map)
: Petra_Map(map) {
  DistributedGlobal_ = false;
  if (CheckInput()!=0) {
    cout << "Replicated Local Map not the same size on all PEs" << endl;
    abort();
  }
}
 
//============================================================================
int Petra_LocalMap::CheckInput() {
  DistributedGlobal_ = false;
  int * tmp = new int[4];
  tmp[0] = NumMyElements_;
  tmp[1] = - NumMyElements_;
  Comm().MaxAll(tmp, tmp+2, 2);

  int tmp1 = tmp[2]; // Max of all NumMyElements across all processors
  int tmp2 = - tmp[3]; // Min of all ...
  delete [] tmp;

  if (tmp1==tmp2) return(0);
  else return(-1);
}
//=========================================================================
Petra_LocalMap::~Petra_LocalMap(){}
