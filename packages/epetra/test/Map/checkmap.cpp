#include "checkmap.h"
int checkmap(Epetra_Map & Map, int NumGlobalElements, int NumMyElements, 
	     int *MyGlobalElements, int IndexBase, Epetra_Comm& Comm,
	     bool DistributedGlobal)
{
  int i;

  if (!Map.ConstantElementSize()) return(-1);

  if (DistributedGlobal!=Map.DistributedGlobal()) return(-3);


  if (Map.ElementSize()!=1) return(-4);
  int *MyElementSizeList = new int[NumMyElements];

  if (Map.ElementSizeList(MyElementSizeList)!=0) return(-5);
  for (i=0; i<NumMyElements; i++) if (MyElementSizeList[i]!=1) return(-5);

  delete [] MyElementSizeList;

  const Epetra_Comm & Comm1 = Map.Comm();

  if (Comm1.NumProc()!=Comm.NumProc()) return(-6);

  if (Comm1.MyPID()!=Comm.MyPID()) return(-7);

  if (Map.IndexBase()!=IndexBase) return(-8);

  if (!Map.LinearMap() && MyGlobalElements==0) return(-9);

  if (Map.LinearMap() && MyGlobalElements!=0) return(-9);

  if (Map.MaxAllGID()!=NumGlobalElements-1+IndexBase) return(-10);

  if (Map.MaxElementSize()!=1) return(-11);

  int MaxLID = Map.MaxLID();
  if (MaxLID!=NumMyElements-1) return(-12);

  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  if (!DistributedGlobal) MaxMyGID = NumMyElements-1+IndexBase;
  if (Map.MaxMyGID()!=MaxMyGID) return(-13);

  if (Map.MinAllGID()!=IndexBase) return(-14);

  if (Map.MinElementSize()!=1) return(-15);

  if (Map.MinLID()!=0) return(-16);

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

  if (Map.NumGlobalElements()!=NumGlobalElements) return(-19);
  
  if (Map.NumGlobalPoints()!=NumGlobalElements) return(-20);
  
  if (Map.NumMyElements()!=NumMyElements) return(-21);  

  if (Map.NumMyPoints()!=NumMyElements) return(-22);

  int MaxMyGID2 = Map.GID(Map.LID(MaxMyGID));
  if (MaxMyGID2 != MaxMyGID) return (-23);
  int MaxLID2 = Map.LID(Map.GID(MaxLID));
  if (MaxLID2 != MaxLID) return(-24);

  if (Map.GID(MaxLID+1) != IndexBase-1) return(-25);// MaxLID+1 doesn't exist
  if (Map.LID(MaxMyGID+1) != -1) return(-26);// MaxMyGID+1 doesn't exist or is on a different processor

  if (!Map.MyGID(MaxMyGID)) return (-27);
  if (Map.MyGID(MaxMyGID+1)) return (-28);

  if (!Map.MyLID(MaxLID)) return (-29);
  if (Map.MyLID(MaxLID+1)) return (-30);

  if (!Map.MyGID(Map.GID(MaxLID))) return(-31);
  if (Map.MyGID(Map.GID(MaxLID+1))) return(-32);

  if (!Map.MyLID(Map.LID(MaxMyGID))) return(-33);
  if (Map.MyLID(Map.LID(MaxMyGID+1))) return(-34);

  // Check RemoteIDList function
  // Get some GIDs off of each processor to test
  int TotalNumEle, NumElePerProc, NumProc = Comm.NumProc();
  int MinNumEleOnProc;
  int NumMyEle;
  Comm.MinAll(&NumMyEle,&MinNumEleOnProc,1);
  if (MinNumEleOnProc > 5) NumElePerProc = 6;
  else NumElePerProc = MinNumEleOnProc;
  if (NumElePerProc > 0) {
    TotalNumEle = NumElePerProc*NumProc;
    int * MyGIDlist = new int[NumElePerProc];
    int * GIDlist = new int[TotalNumEle];
    int * PIDlist = new int[TotalNumEle];
    int * LIDlist = new int[TotalNumEle];
    for (i=0; i<NumElePerProc; i++)
	  MyGIDlist[i] = MyGlobalElements1[i];
    Comm.GatherAll(MyGIDlist,GIDlist,NumElePerProc);// Get a few values from each proc
    Map.RemoteIDList(TotalNumEle, GIDlist, PIDlist, LIDlist);
    int MyPID= Comm.MyPID();
    for (i=0; i<TotalNumEle; i++) {
      if (Map.MyGID(GIDlist[i])) {
	if (PIDlist[i] != MyPID) return(-44);
	if (!Map.MyLID(Map.LID(GIDlist[i])) || Map.LID(GIDlist[i]) != LIDlist[i] || Map.GID(LIDlist[i]) != GIDlist[i]) return(-45);
      }
      else {
	if (PIDlist[i] == MyPID) return(-47); // If MyGID comes back false, the PID listed should be that of another proc
      }
    }

    delete [] MyGIDlist;
    delete [] GIDlist;
    delete [] PIDlist;
    delete [] LIDlist;
  }

  delete [] MyGlobalElements1;

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

 
