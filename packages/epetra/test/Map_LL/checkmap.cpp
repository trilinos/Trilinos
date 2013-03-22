//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER


#include "checkmap.h"
int checkmap(Epetra_Map & Map, long long NumGlobalElements, int NumMyElements, 
	     long long *MyGlobalElements, long long IndexBase, Epetra_Comm& Comm,
	     bool DistributedGlobal)
{
  int i, ierr=0, forierr = 0;

  EPETRA_TEST_ERR(!Map.ConstantElementSize(),ierr);

  EPETRA_TEST_ERR(DistributedGlobal!=Map.DistributedGlobal(),ierr);


  EPETRA_TEST_ERR(Map.ElementSize()!=1,ierr);
  int *MyElementSizeList = new int[NumMyElements];

  EPETRA_TEST_ERR(Map.ElementSizeList(MyElementSizeList)!=0,ierr);

  forierr = 0;
  for (i=0; i<NumMyElements; i++) forierr += MyElementSizeList[i]!=1;
  EPETRA_TEST_ERR(forierr,ierr);

  delete [] MyElementSizeList;

  const Epetra_Comm & Comm1 = Map.Comm();

  EPETRA_TEST_ERR(Comm1.NumProc()!=Comm.NumProc(),ierr);

  EPETRA_TEST_ERR(Comm1.MyPID()!=Comm.MyPID(),ierr);

  EPETRA_TEST_ERR(Map.IndexBase64()!=IndexBase,ierr);

  EPETRA_TEST_ERR(!Map.LinearMap() && MyGlobalElements==0,ierr);

  EPETRA_TEST_ERR(Map.LinearMap() && MyGlobalElements!=0,ierr);

  EPETRA_TEST_ERR(Map.MaxAllGID64()!=NumGlobalElements-1+IndexBase,ierr);

  EPETRA_TEST_ERR(Map.MaxElementSize()!=1,ierr);

  int MaxLID = Map.MaxLID();
  EPETRA_TEST_ERR(MaxLID!=NumMyElements-1,ierr);

  long long MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  if (!DistributedGlobal) MaxMyGID = NumMyElements-1+IndexBase;
  EPETRA_TEST_ERR(Map.MaxMyGID64()!=MaxMyGID,ierr);

  EPETRA_TEST_ERR(Map.MinAllGID64()!=IndexBase,ierr);

  EPETRA_TEST_ERR(Map.MinElementSize()!=1,ierr);

  EPETRA_TEST_ERR(Map.MinLID()!=0,ierr);

  long long MinMyGID = Comm.MyPID()*NumMyElements+IndexBase;
  if (Comm.MyPID()>2) MinMyGID+=3;
  if (!DistributedGlobal) MinMyGID = 0;
  EPETRA_TEST_ERR(Map.MinMyGID64()!=MinMyGID,ierr);
  
  long long * MyGlobalElements1 = new long long[NumMyElements];
  EPETRA_TEST_ERR(Map.MyGlobalElements(MyGlobalElements1)!=0,ierr);

  forierr = 0;
  if (MyGlobalElements==0)
    {
      for (i=0; i<NumMyElements; i++) 
	forierr += MyGlobalElements1[i]!=MinMyGID+i;
      EPETRA_TEST_ERR(forierr,ierr);
    }
  else {
    for (i=0; i<NumMyElements; i++)
      forierr += MyGlobalElements[i]!=MyGlobalElements1[i];
    EPETRA_TEST_ERR(forierr,ierr);
  }
  EPETRA_TEST_ERR(Map.NumGlobalElements64()!=NumGlobalElements,ierr);
  
  EPETRA_TEST_ERR(Map.NumGlobalPoints64()!=NumGlobalElements,ierr);
  
  EPETRA_TEST_ERR(Map.NumMyElements()!=NumMyElements,ierr);  

  EPETRA_TEST_ERR(Map.NumMyPoints()!=NumMyElements,ierr);

  long long MaxMyGID2 = Map.GID64(Map.LID(MaxMyGID));
  EPETRA_TEST_ERR(MaxMyGID2 != MaxMyGID,ierr);
  int MaxLID2 = Map.LID(Map.GID64(MaxLID));
  EPETRA_TEST_ERR(MaxLID2 != MaxLID,ierr);

  EPETRA_TEST_ERR(Map.GID64(MaxLID+1) != IndexBase-1,ierr);// MaxLID+1 doesn't exist
  EPETRA_TEST_ERR(Map.LID(MaxMyGID+1) != -1,ierr);// MaxMyGID+1 doesn't exist or is on a different processor

  EPETRA_TEST_ERR(!Map.MyGID(MaxMyGID),ierr);
  EPETRA_TEST_ERR(Map.MyGID(MaxMyGID+1),ierr);

  EPETRA_TEST_ERR(!Map.MyLID(MaxLID),ierr);
  EPETRA_TEST_ERR(Map.MyLID(MaxLID+1),ierr);

  EPETRA_TEST_ERR(!Map.MyGID(Map.GID64(MaxLID)),ierr);
  EPETRA_TEST_ERR(Map.MyGID(Map.GID64(MaxLID+1)),ierr);

  EPETRA_TEST_ERR(!Map.MyLID(Map.LID(MaxMyGID)),ierr);
  EPETRA_TEST_ERR(Map.MyLID(Map.LID(MaxMyGID+1)),ierr);

  // Check RemoteIDList function
  // Get some GIDs off of each processor to test
  int TotalNumEle, NumElePerProc, NumProc = Comm.NumProc();
  int MinNumEleOnProc;
  int NumMyEle=Map.NumMyElements();
  Comm.MinAll(&NumMyEle,&MinNumEleOnProc,1);
  if (MinNumEleOnProc > 5) NumElePerProc = 6;
  else NumElePerProc = MinNumEleOnProc;
  if (NumElePerProc > 0) {
    TotalNumEle = NumElePerProc*NumProc;
    long long * MyGIDlist = new long long[NumElePerProc];
    long long * GIDlist = new long long[TotalNumEle];
    int * PIDlist = new int[TotalNumEle];
    int * LIDlist = new int[TotalNumEle];
    for (i=0; i<NumElePerProc; i++)
	  MyGIDlist[i] = MyGlobalElements1[i];
    Comm.GatherAll(MyGIDlist,GIDlist,NumElePerProc);// Get a few values from each proc
    Map.RemoteIDList(TotalNumEle, GIDlist, PIDlist, LIDlist);
    int MyPID= Comm.MyPID();

    forierr = 0;
    for (i=0; i<TotalNumEle; i++) {
      if (Map.MyGID(GIDlist[i])) {
	forierr += PIDlist[i] != MyPID;
	forierr += !Map.MyLID(Map.LID(GIDlist[i])) || Map.LID(GIDlist[i]) != LIDlist[i] || Map.GID64(LIDlist[i]) != GIDlist[i];
      }
      else {
	forierr += PIDlist[i] == MyPID; // If MyGID comes back false, the PID listed should be that of another proc
      }
    }
    EPETRA_TEST_ERR(forierr,ierr);

    delete [] MyGIDlist;
    delete [] GIDlist;
    delete [] PIDlist;
    delete [] LIDlist;
  }

  delete [] MyGlobalElements1;

  // Check RemoteIDList function (assumes all maps are linear, even if not stored that way)

  if (Map.LinearMap()) {

    long long * GIDList = new long long[3];
    int * PIDList = new int[3];
    int * LIDList = new int[3];
    int MyPID = Map.Comm().MyPID();
  
    int NumIDs = 0;
    //GIDList[NumIDs++] = Map.MaxAllGID64()+1; // Should return -1 for both PID and LID
    if (Map.MinMyGID64()-1>=Map.MinAllGID64()) GIDList[NumIDs++] = Map.MinMyGID64()-1;
    if (Map.MaxMyGID64()+1<=Map.MaxAllGID64()) GIDList[NumIDs++] = Map.MaxMyGID64()+1;

    Map.RemoteIDList(NumIDs, GIDList, PIDList, LIDList);

    NumIDs = 0;

    //EPETRA_TEST_ERR(!(PIDList[NumIDs]==-1),ierr);
    //EPETRA_TEST_ERR(!(LIDList[NumIDs++]==-1),ierr);

    if (Map.MinMyGID64()-1>=Map.MinAllGID64()) EPETRA_TEST_ERR(!(PIDList[NumIDs++]==MyPID-1),ierr);
    if (Map.MaxMyGID64()+1<=Map.MaxAllGID64()) EPETRA_TEST_ERR(!(PIDList[NumIDs]==MyPID+1),ierr);
    if (Map.MaxMyGID64()+1<=Map.MaxAllGID64()) EPETRA_TEST_ERR(!(LIDList[NumIDs++]==0),ierr);

    delete [] GIDList;
    delete [] PIDList;
    delete [] LIDList;

  }
  
  return (ierr);
}

 
