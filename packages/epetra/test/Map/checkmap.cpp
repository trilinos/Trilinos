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

  if (Map.MaxLID()!=NumMyElements-1+IndexBase) return(-12);

  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  if (!DistributedGlobal) MaxMyGID = NumMyElements-1+IndexBase;
  if (Map.MaxMyGID()!=MaxMyGID) return(-13);

  if (Map.MinAllGID()!=IndexBase) return(-14);

  if (Map.MinElementSize()!=1) return(-15);

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
  
  if (Map.NumGlobalPoints()!=NumGlobalElements) return(-20);
  
  if (Map.NumMyElements()!=NumMyElements) return(-21);  

  if (Map.NumMyPoints()!=NumMyElements) return(-22);

  return(0);
}
