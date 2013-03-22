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
int checkmap(Epetra_BlockMap & Map, long long NumGlobalElements, int NumMyElements, 
	     long long *MyGlobalElements, int ElementSize, int * ElementSizeList,
	     long long NumGlobalPoints, int NumMyPoints,
	     long long IndexBase, Epetra_Comm& Comm,
	     bool DistributedGlobal,
	     bool IsOneToOne)
{

  int i, ierr=0, forierr=0;// forierr is used in for loops, then is tested
  // after for loop completes to see if it is non zero - potentially prevents
  // thousands of error messages

  if (ElementSizeList==0)
    {
      EPETRA_TEST_ERR(!Map.ConstantElementSize(),ierr);
    }
  else
    EPETRA_TEST_ERR(Map.ConstantElementSize(),ierr);
  
  EPETRA_TEST_ERR(DistributedGlobal!=Map.DistributedGlobal(),ierr);

  EPETRA_TEST_ERR(IsOneToOne!=Map.IsOneToOne(),ierr);

  int *MyElementSizeList;

  if (ElementSizeList==0)
    {
      EPETRA_TEST_ERR(Map.ElementSize()!=ElementSize,ierr);
      
      MyElementSizeList = new int[NumMyElements];
      
      EPETRA_TEST_ERR(Map.ElementSizeList(MyElementSizeList)!=0,ierr);
      forierr = 0;
      for (i=0; i<NumMyElements; i++) 
        forierr += MyElementSizeList[i]!=ElementSize;
      EPETRA_TEST_ERR(forierr,ierr);

      EPETRA_TEST_ERR(Map.MaxMyElementSize() != ElementSize,ierr);
      EPETRA_TEST_ERR(Map.MinMyElementSize() != ElementSize,ierr);
    }
  else
    {
      MyElementSizeList = new int[NumMyElements];
      EPETRA_TEST_ERR(Map.ElementSizeList(MyElementSizeList)!=0,ierr);
      int MaxSize = MyElementSizeList[0];
      int MinSize = MyElementSizeList[0];
      forierr=0;
      for (i=0; i<NumMyElements; i++) {
        forierr += MyElementSizeList[i]!=ElementSizeList[i];
	if (MyElementSizeList[i] > MaxSize)
	  MaxSize = MyElementSizeList[i];
	if (MyElementSizeList[i] < MinSize)
	  MinSize = MyElementSizeList[i];

	// Test ElementSize(int LID) method	

	forierr += Map.ElementSize(Map.LID(MyGlobalElements[i])) != ElementSizeList[i];
      }
      EPETRA_TEST_ERR(forierr,ierr);
   
      EPETRA_TEST_ERR(MaxSize !=Map.MaxMyElementSize(),ierr);
      EPETRA_TEST_ERR(MinSize !=Map.MinMyElementSize(),ierr);
    }

  const Epetra_Comm & Comm1 = Map.Comm();

  EPETRA_TEST_ERR(Comm1.NumProc()!=Comm.NumProc(),ierr);

  EPETRA_TEST_ERR(Comm1.MyPID()!=Comm.MyPID(),ierr);

  EPETRA_TEST_ERR(Map.IndexBase64()!=IndexBase,ierr);

  EPETRA_TEST_ERR(!Map.LinearMap() && MyGlobalElements==0,ierr);

  EPETRA_TEST_ERR(Map.LinearMap() && MyGlobalElements!=0,ierr);

  EPETRA_TEST_ERR(Map.MaxAllGID64()!=NumGlobalElements-1+IndexBase,ierr);

  EPETRA_TEST_ERR(Map.MaxElementSize()!=ElementSize,ierr);

  int MaxLID = Map.MaxLID();
  EPETRA_TEST_ERR(MaxLID!=NumMyElements-1,ierr);

  long long MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  if (!DistributedGlobal) MaxMyGID = NumMyElements-1+IndexBase;
  EPETRA_TEST_ERR(Map.MaxMyGID64()!=MaxMyGID,ierr);

  EPETRA_TEST_ERR(Map.MinAllGID64()!=IndexBase,ierr);

  if (ElementSizeList==0)
    {
      EPETRA_TEST_ERR(Map.MinElementSize()!=ElementSize,ierr);
    }
  else EPETRA_TEST_ERR(Map.MinElementSize()!=2,ierr);

  int MinLID = Map.MinLID();
  EPETRA_TEST_ERR(MinLID!=0,ierr);

  long long MinMyGID = Comm.MyPID()*NumMyElements+IndexBase;
  if (Comm.MyPID()>2) MinMyGID+=3;
  if (!DistributedGlobal) MinMyGID = IndexBase; // Not really needed
  EPETRA_TEST_ERR(Map.MinMyGID64()!=MinMyGID,ierr);

  long long * MyGlobalElements1 = new long long[NumMyElements];
  EPETRA_TEST_ERR(Map.MyGlobalElements(MyGlobalElements1)!=0,ierr);
  
  forierr = 0;
  if (MyGlobalElements==0) {
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
  
  EPETRA_TEST_ERR(Map.NumGlobalPoints64()!=NumGlobalPoints,ierr);
  
  EPETRA_TEST_ERR(Map.NumMyElements()!=NumMyElements,ierr);  

  EPETRA_TEST_ERR(Map.NumMyPoints()!=NumMyPoints,ierr);

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

  // Test the FirstPointInElementList methods, begin by testing that they produce identical results
  int * FirstPointInElementList = new int[NumMyElements+1];
  Map.FirstPointInElementList(FirstPointInElementList);
  int * FirstPointInElementList1 = Map.FirstPointInElementList();
  forierr = 0;
  for (i=0; i<=NumMyElements; i++)
    forierr += FirstPointInElementList[i]!=FirstPointInElementList1[i];
  EPETRA_TEST_ERR(forierr,ierr);
  // Now make sure values are correct
  forierr = 0;
  if (Map.ConstantElementSize()) {
    for (i=0; i<=NumMyElements; i++)
      forierr += FirstPointInElementList1[i]!=(i*ElementSize);// NOTE:FirstPointInElement[NumMyElements] is not the first point of an element
    EPETRA_TEST_ERR(forierr,ierr);
  }
  else {
    int FirstPoint = 0;
    for (i=0; i<NumMyElements; i++) {
      forierr += FirstPointInElementList1[i]!=FirstPoint;
      FirstPoint += ElementSizeList[i];
    }
    EPETRA_TEST_ERR(forierr,ierr);
    EPETRA_TEST_ERR(FirstPointInElementList[NumMyElements] != NumMyPoints,ierr);// The last entry in the array = the total number of Points on the proc
  }
  delete [] FirstPointInElementList;

  // Declare some variables for the FindLocalElementID test
  int ElementID, Offset;
  // Test the PointToElementList methods, begin by testing that they produce identical results
  int * PointToElementList = new int[NumMyPoints];
  Map.PointToElementList(PointToElementList);
  int * PointToElementList1 = Map.PointToElementList();
  forierr = 0;
  for (i=0; i<NumMyPoints; i++)
    forierr += PointToElementList1[i] != PointToElementList[i];
  EPETRA_TEST_ERR(forierr,ierr);
  //Now make sure values are correct
  forierr=0;
  if (Map.ConstantElementSize()) {
    for (i=0; i<NumMyElements; i++)
      for (int j=0; j<ElementSize; j++) {
	forierr += PointToElementList[i*ElementSize+j] != i;
	// Test FindLocalElementID method
	Map.FindLocalElementID(i*ElementSize+j,ElementID,Offset);
	forierr += ElementID != i || Offset != j;
      }
    EPETRA_TEST_ERR(forierr,ierr);
  }
  else {
    int MyPointTot = 0; // Keep track of total number of points in all previously completely checked elements
    for (i=0; i<NumMyElements; i++) {
      for (int j=0; j<ElementSizeList[i]; j++) {
	forierr += PointToElementList[MyPointTot+j] != i;
	// Test FindLocalElementID method
	Map.FindLocalElementID(MyPointTot+j,ElementID,Offset);
	forierr += ElementID != i || Offset != j;
      }
      MyPointTot += ElementSizeList[i];
    }
    EPETRA_TEST_ERR(forierr,ierr);
  }
  delete [] PointToElementList;

  // Check RemoteIDList function that includes a parameter for size
  // Get some GIDs off of each processor to test
  int TotalNumEle, NumElePerProc, NumProc = Comm.NumProc();
  int MinNumEleOnProc;
  int NumMyEle = Map.NumMyElements();
  Comm.MinAll(&NumMyEle,&MinNumEleOnProc,1);
  if (MinNumEleOnProc > 5) NumElePerProc = 6;
  else NumElePerProc = MinNumEleOnProc;
  if (NumElePerProc > 0) {
    TotalNumEle = NumElePerProc*NumProc;
    long long * MyGIDlist = new long long[NumElePerProc];
    long long * GIDlist = new long long[TotalNumEle];
    int * PIDlist = new int[TotalNumEle];
    int * LIDlist = new int[TotalNumEle];
    int * SizeList = new int[TotalNumEle];
    for (i=0; i<NumElePerProc; i++)
	  MyGIDlist[i] = MyGlobalElements1[i];
    Comm.GatherAll(MyGIDlist,GIDlist,NumElePerProc);// Get a few values from each proc
    Map.RemoteIDList(TotalNumEle, GIDlist, PIDlist, LIDlist, SizeList);
    int MyPID= Comm.MyPID();
    forierr = 0;
    for (i=0; i<TotalNumEle; i++) {
      if (Map.MyGID(GIDlist[i])) {
	forierr += PIDlist[i] != MyPID;
	forierr += !Map.MyLID(Map.LID(GIDlist[i])) || Map.LID(GIDlist[i]) != LIDlist[i] || Map.GID64(LIDlist[i]) != GIDlist[i];
	forierr += SizeList[i] != Map.ElementSize(LIDlist[i]);
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
    delete [] SizeList;
  }

  delete [] MyGlobalElements1;
  delete [] MyElementSizeList;

  // Check RemoteIDList function (assumes all maps are linear, even if not stored that way)

  if (Map.LinearMap()) {

    long long * GIDList = new long long[3];
    int * PIDList = new int[3];
    int * LIDList = new int[3];
    int MyPID = Map.Comm().MyPID();
  
    int NumIDs = 0;
    //GIDList[NumIDs++] = Map.MaxAllGID()+1; // Should return -1 for both PID and LID
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
