#include "checkmap.h"
int checkmap(Epetra_BlockMap & Map, int NumGlobalElements, int NumMyElements, 
	     int *MyGlobalElements, int ElementSize, int * ElementSizeList,
	     int NumGlobalPoints, int NumMyPoints,
	     int IndexBase, Epetra_Comm& Comm,
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

      if (Map.MaxMyElementSize() != ElementSize) return(-54);
      if (Map.MinMyElementSize() != ElementSize) return(-56);
      delete [] MyElementSizeList;
    }
  else
    {
      MyElementSizeList = new int[NumMyElements];
      if (Map.ElementSizeList(MyElementSizeList)!=0) return(-5);
      int MaxSize = MyElementSizeList[0];
      int MinSize = MyElementSizeList[0];
      for (i=0; i<NumMyElements; i++) {
        if (MyElementSizeList[i]!=ElementSizeList[i]) return(-5);
	if (MyElementSizeList[i] > MaxSize)
	  MaxSize = MyElementSizeList[i];
	if (MyElementSizeList[i] < MinSize)
	  MinSize = MyElementSizeList[i];

	// Test ElementSize(int LID) method	

	if (Map.ElementSize(Map.LID(MyGlobalElements[i])) != ElementSizeList[i])
	  return(-51);
      }
      if (MaxSize !=Map.MaxMyElementSize()) return(-53);
      if (MinSize !=Map.MinMyElementSize()) return(-55);
      delete [] MyElementSizeList;
    }

  const Epetra_Comm & Comm1 = Map.Comm();

  if (Comm1.NumProc()!=Comm.NumProc()) return(-6);

  if (Comm1.MyPID()!=Comm.MyPID()) return(-7);

  if (Map.IndexBase()!=IndexBase) return(-8);

  if (!Map.LinearMap() && MyGlobalElements==0) return(-9);

  if (Map.LinearMap() && MyGlobalElements!=0) return(-9);

  if (Map.MaxAllGID()!=NumGlobalElements-1+IndexBase) return(-10);

  if (Map.MaxElementSize()!=ElementSize) return(-11);

  int MaxLID = Map.MaxLID();
  if (MaxLID!=NumMyElements-1+IndexBase) return(-12);

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

  int MinLID = Map.MinLID();
  if (MinLID!=IndexBase) return(-16);

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
  
  if (Map.NumGlobalPoints()!=NumGlobalPoints) return(-20);
  
  if (Map.NumMyElements()!=NumMyElements) return(-21);  

  if (Map.NumMyPoints()!=NumMyPoints) return(-22);

  int MaxMyGID2 = Map.GID(Map.LID(MaxMyGID));
  if (MaxMyGID2 != MaxMyGID) return (-23);
  int MaxLID2 = Map.LID(Map.GID(MaxLID));
  if (MaxLID2 != MaxLID) return(-24);

  if (Map.GID(MaxLID+1) != -1) return(-25);// MaxLID+1 doesn't exist
  if (Map.LID(MaxMyGID+1) != -1) return(-26);// MaxMyGID+1 doesn't exist or is on a different processor

  if (!Map.MyGID(MaxMyGID)) return (-27);
  if (Map.MyGID(MaxMyGID+1)) return (-28);

  if (!Map.MyLID(MaxLID)) return (-29);
  if (Map.MyLID(MaxLID+1)) return (-30);

  if (!Map.MyGID(Map.GID(MaxLID))) return(-31);
  if (Map.MyGID(Map.GID(MaxLID+1))) return(-32);

  if (!Map.MyLID(Map.LID(MaxMyGID))) return(-33);
  if (Map.MyLID(Map.LID(MaxMyGID+1))) return(-34);

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
