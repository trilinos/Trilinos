#include "checkmap.h"
int checkmap(Petra_BlockMap & Map, int NumGlobalElements, int NumMyElements, 
	     int *MyGlobalElements, int ElementSize, int * ElementSizeList,
	     int NumGlobalEquations, int NumMyEquations,
	     int IndexBase, Petra_Comm& Comm,
	     bool DistributedGlobal)
{
  int i;

  if (ElementSizeList==0)
    {
      if (!Map.ConstantElementSize()) return(-1);
    }
  else
    if (Map.ConstantElementSize()) return(-1);
  
  if (DistributedGlobal!=Map.DistributedGlobal()) return(-3);

      int *MyElementSizeList;

  if (ElementSizeList==0)
    {
      if (Map.ElementSize()!=ElementSize) return(-4);
      
      MyElementSizeList = new int[NumMyElements];
      
      if (Map.ElementSizeList(MyElementSizeList)!=0) return(-5);
      for (i=0; i<NumMyElements; i++) 
        if (MyElementSizeList[i]!=ElementSize) return(-5);

      delete [] MyElementSizeList;
    }
  else
    {
      MyElementSizeList = new int[NumMyElements];
      if (Map.ElementSizeList(MyElementSizeList)!=0) return(-5);

      for (i=0; i<NumMyElements; i++) 
        if (MyElementSizeList[i]!=ElementSizeList[i]) return(-5);

      delete [] MyElementSizeList;
    }

  const Petra_Comm & Comm1 = Map.Comm();

  if (Comm1.NumProc()!=Comm.NumProc()) return(-6);

  if (Comm1.MyPID()!=Comm.MyPID()) return(-7);

  if (Map.IndexBase()!=IndexBase) return(-8);

  if (!Map.LinearMap() && MyGlobalElements==0) return(-9);

  if (Map.LinearMap() && MyGlobalElements!=0) return(-9);

  if (Map.MaxAllGID()!=NumGlobalElements-1+IndexBase) return(-10);

  if (Map.MaxElementSize()!=ElementSize) return(-11);

  if (Map.MaxLID()!=NumMyElements-1+IndexBase) return(-12);

  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  if (!DistributedGlobal) MaxMyGID = NumMyElements-1+IndexBase;
  if (Map.MaxMyGID()!=MaxMyGID) return(-13);

  if (Map.MinAllGID()!=IndexBase) return(-14);

  if (ElementSizeList==0)
    {
      if (Map.MinElementSize()!=ElementSize) return(-15);
    }
  else if (Map.MinElementSize()!=2) return(-15);

  if (Map.MinLID()!=IndexBase) return(-16);

  int MinMyGID = Comm.MyPID()*NumMyElements+IndexBase;
  if (Comm.MyPID()>2) MinMyGID+=3;
  if (!DistributedGlobal) MinMyGID = 0;
  if (Map.MinMyGID()!=MinMyGID) return(-17);
  
  int * MyGlobalElements1 = new int[NumMyElements];
  if (Map.MyGlobalElements(MyGlobalElements1)!=0) return(-18);
  
  if (MyGlobalElements==0)
    {
      for (i=0; i<NumMyElements; i++) 
	if (MyGlobalElements1[i]!=MinMyGID+i) return(-18);
    }
  else
    for (i=0; i<NumMyElements; i++)
      if (MyGlobalElements[i]!=MyGlobalElements1[i]) return(-18);
  
  delete [] MyGlobalElements1;

  if (Map.NumGlobalElements()!=NumGlobalElements) return(-19);
  
  if (Map.NumGlobalEquations()!=NumGlobalEquations) return(-20);
  
  if (Map.NumMyElements()!=NumMyElements) return(-21);  

  if (Map.NumMyEquations()!=NumMyEquations) return(-22);

  // Check RemoteIDList function (assumes all maps are linear, even if not stored that way)

  if (Map.LinearMap()) {

    int * GIDList = new int[3];
    int * PIDList = new int[3];
    int * LIDList = new int[3];
    int MyPID = Map.Comm().MyPID();
  
    int NumIDs = 0;
    //GIDList[NumIDs++] = Map.MaxAllGID()+1; // Should return -1 for both PID and LID
    if (Map.MinMyGID()-1>=Map.MinAllGID()) GIDList[NumIDs++] = Map.MinMyGID()-1;
    if (Map.MaxMyGID()+1<=Map.MaxAllGID()) GIDList[NumIDs++] = Map.MaxMyGID()+1;

    Map.RemoteIDList(NumIDs, GIDList, PIDList, LIDList);

    NumIDs = 0;

    //assert(PIDList[NumIDs]==-1);
    //assert(LIDList[NumIDs++]==-1);

    if (Map.MinMyGID()-1>=Map.MinAllGID()) assert(PIDList[NumIDs++]==MyPID-1);
    if (Map.MaxMyGID()+1<=Map.MaxAllGID()) assert(PIDList[NumIDs]==MyPID+1);
    if (Map.MaxMyGID()+1<=Map.MaxAllGID()) assert(LIDList[NumIDs++]==0);

    delete [] GIDList;
    delete [] PIDList;
    delete [] LIDList;


  }
  return(0);
}
