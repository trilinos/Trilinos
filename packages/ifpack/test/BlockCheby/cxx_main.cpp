/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AZTECOO) && defined(HAVE_IFPACK_EPETRAEXT)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "AztecOO.h"
#include "EpetraExt_BlockDiagMatrix.h"
#include "EpetraExt_PointToBlockDiagPermute.h"
#include "Ifpack_Chebyshev.h"

using namespace Trilinos_Util;
/*
(unused; commented out to avoid build warnings)
static bool verbose = false;
static bool SymmetricGallery = false;
static bool Solver = AZ_gmres;
*/

//=============================================
// Test the BlockDiagMatrix
bool TestBlockDiagMatrix(const Epetra_Comm& Comm){
  const int NUM_BLOCKS=30;
  bool TestPassed=true;
  int my_blockgids[NUM_BLOCKS];
  int my_blocksizes[NUM_BLOCKS];

  for(int i=0;i<NUM_BLOCKS;i++){
    my_blockgids[i]=i+NUM_BLOCKS*Comm.MyPID();
    if(i<NUM_BLOCKS/3)
      my_blocksizes[i]=1;
    else if(i<2*NUM_BLOCKS/3)
      my_blocksizes[i]=2;
    else
      my_blocksizes[i]=3;
   }


  // Build me a map and a DBM to go with it...
  Epetra_BlockMap BDMap(-1,NUM_BLOCKS,my_blockgids,my_blocksizes,0,Comm);
  EpetraExt_BlockDiagMatrix BMAT(BDMap,true);

  // Fill the matrix - This tests all three block size cases in the code, 1x1, 2x2 and larger.
  for(int i=0;i<BMAT.NumMyBlocks();i++){
    double*start=BMAT[i];
    if(BMAT.BlockSize(i)==1)
      *start=2.0;
    else if(BMAT.BlockSize(i)==2){
      start[0]=2.0;
      start[1]=-1.0;
      start[2]=-1.0;
      start[3]=2.0;
    }
    else if(BMAT.BlockSize(i)==3){
      start[0]=4.0;
      start[1]=-1.0;
      start[2]=-1.0;
      start[3]=-1.0;
      start[4]=4.0;
      start[5]=-1.0;
      start[2]=-1.0;
      start[7]=-1.0;
      start[8]=4.0;
    }
    else
      exit(1);
  }


  // Allocations for tests
  double norm2;
  Epetra_MultiVector X(BDMap,1);
  Epetra_MultiVector Y(BDMap,1);
  Epetra_MultiVector Z(BDMap,1);
  EpetraExt_BlockDiagMatrix BMAT_forward(BMAT);
  EpetraExt_BlockDiagMatrix BMAT_factor(BMAT);
  Teuchos::ParameterList List;

  //***************************
  // Test the Forward/Invert
  //***************************
  List.set("apply mode","invert");
  BMAT.SetParameters(List);
  X.SetSeed(24601); X.Random();
  BMAT_forward.ApplyInverse(X,Y);
  BMAT.Compute();
  BMAT.ApplyInverse(Y,Z);
  X.Update(1.0,Z,-1.0);
  X.Norm2(&norm2);
  if(!Comm.MyPID()) cout<<"Forward/Invert Error = "<<norm2<<endl;
  if(norm2 > 1e-12) TestPassed=false;

  //***************************
  // Test the Forward/Factor
  //***************************
  List.set("apply mode","factor");
  BMAT_factor.SetParameters(List);
  X.SetSeed(24601); X.Random();
  BMAT_forward.ApplyInverse(X,Y);
  BMAT_factor.Compute();
  BMAT_factor.ApplyInverse(Y,Z);
  X.Update(1.0,Z,-1.0);
  X.Norm2(&norm2);
  if(!Comm.MyPID()) cout<<"Forward/Factor Error = "<<norm2<<endl;
  if(norm2 > 1e-12) TestPassed=false;

  return TestPassed;
}

//=============================================
//=============================================
//=============================================
void Build_Local_Contiguous_Size2_BlockMatrix(const Epetra_Comm & Comm, int NUM_ROWS,int &NumBlocks,
                                        int *&Blockstart_, int *&Blockids_,Epetra_CrsMatrix *&MAT){
  double values[3]={-1.0,2.0,-1.0};
  int indices[3];

  // Build me a CrsMatrix
  Epetra_Map Map(-1,NUM_ROWS,0,Comm);
  MAT=new Epetra_CrsMatrix(Copy,Map,2);
  assert(MAT->NumMyRows()%2==0);

  NumBlocks=MAT->NumMyRows()/2;
  Blockstart_=new int [NumBlocks+1];
  Blockids_=new int [MAT->NumMyRows()];
  Blockstart_[0]=0;
  int curr_block=0;

  for(int i=0;i<MAT->NumMyRows();i++){
    // local contiguous blocks of constant size 2
    int row_in_block=i%2;
    if(row_in_block==0){
      Blockstart_[curr_block]=i;
      indices[0]=i; indices[1]=i+1;
      Blockids_[i]=i;
      MAT->InsertGlobalValues(Map.GID(i),2,&values[1],&indices[0]);
    }
    else if(row_in_block==1){
      indices[0]=i-1; indices[1]=i;
      MAT->InsertGlobalValues(Map.GID(i),2,&values[0],&indices[0]);
      Blockids_[i]=i;
      curr_block++;
    }
  }
  Blockstart_[NumBlocks]=MAT->NumMyRows();

  MAT->FillComplete();
}

//=============================================
void Build_Local_NonContiguous_Size2_BlockMatrix(const Epetra_Comm & Comm, int NUM_ROWS,int &NumBlocks,
                                                 int *&Blockstart_, int *&Blockids_,Epetra_CrsMatrix *&MAT){
  double values[3]={-1.0,2.0,-1.0};
  int indices[3];

  // Build me a CrsMatrix
  Epetra_Map Map(-1,NUM_ROWS,0,Comm);
  MAT=new Epetra_CrsMatrix(Copy,Map,2);
  assert(MAT->NumMyRows()%2==0);

  NumBlocks=MAT->NumMyRows()/2;
  Blockstart_=new int [NumBlocks+1];
  Blockids_=new int [MAT->NumMyRows()];
  Blockstart_[0]=0;
  int curr_block=0;

  for(int i=0;i<MAT->NumMyRows();i++){
    int row_in_block=(i%2)?1:0;
    if(row_in_block==0){
      int idx=i/2;
      Blockstart_[curr_block]=i;
      indices[0]=Map.GID(idx); indices[1]=Map.GID(idx+NumBlocks);
      Blockids_[i]=idx;
      MAT->InsertGlobalValues(Map.GID(idx),2,&values[1],&indices[0]);

    }
    else if(row_in_block==1){
      int idx=(i-1)/2+NumBlocks;
      indices[0]=Map.GID(idx-NumBlocks); indices[1]=Map.GID(idx);
      MAT->InsertGlobalValues(Map.GID(idx),2,&values[0],&indices[0]);
      Blockids_[i]=idx;
      curr_block++;
    }
  }
  Blockstart_[NumBlocks]=MAT->NumMyRows();

  MAT->FillComplete();
}
//=============================================
void Build_NonLocal_BlockMatrix(const Epetra_Comm & Comm, int NUM_ROWS,int &NumBlocks,
                                int *&Blockstart_, int *&Blockids_,Epetra_CrsMatrix *&MAT){
  double values[3]={-1.0,2.0,-1.0};
  int indices[3];
  int NumProcs=Comm.NumProc();

  // Build me a CrsMatrix
  Epetra_Map Map(-1,NUM_ROWS,0,Comm);
  MAT=new Epetra_CrsMatrix(Copy,Map,2);
  int MyPID=Comm.MyPID();

  for(int i=0;i<MAT->NumMyRows();i++){
    int GID=Map.GID(i);
    indices[0]=GID;
    if(i==0) values[1]=2.0;
    else values[1]=NumProcs+1;

    MAT->InsertGlobalValues(GID,1,&values[1],&indices[0]);

    // A little off-diagonal for good measure
    if(NumProcs>1 && MyPID==0 && i>0){
      indices[1]=GID;
      for(int j=1;j<NumProcs;j++){
        indices[1]+=NUM_ROWS;//HAQ
        MAT->InsertGlobalValues(GID,1,&values[0],&indices[1]);
      }
    }
    else if(NumProcs>1 && MyPID!=0 && i>0){
      indices[1]=GID%NUM_ROWS;//HAQ
      MAT->InsertGlobalValues(GID,1,&values[0],&indices[1]);
    }
  }
  MAT->FillComplete();

  // Build me some block structure
  // The first row on each proc is a solo block.  All others get blocked to
  // PID=0's second block
  if(MyPID==0) NumBlocks=NUM_ROWS;
  else NumBlocks=1;
  Blockstart_=new int [NumBlocks+1];
  Blockstart_[0]=0;
  int curr_idx,curr_block;

  if(MyPID){
    // PID > 0
    Blockids_=new int[1];
    Blockstart_[0]=0; Blockstart_[1]=1;
    Blockids_[0]=Map.GID(0);
  }
  else{
    // PID 0
    int nnz=NumProcs*NumBlocks;
    Blockids_=new int[nnz+1];
    Blockstart_[0]=0;
    Blockids_[0]=Map.GID(0);
    curr_idx=1; curr_block=1;
    for(int j=1;j<NUM_ROWS;j++){
      Blockstart_[curr_block]=curr_idx;
      for(int i=0;i<NumProcs;i++){
        Blockids_[curr_idx]=NUM_ROWS*i+j;//FIX: THIS IS A HACK
        curr_idx++;
      }
      curr_block++;
    }
    Blockstart_[curr_block]=curr_idx;
  }
}

//=============================================
void Build_DiagonalStructure(const Epetra_Map &Map,int &NumBlocks,int *&Blockstart_, int *&Blockids_,bool local_ids){
  NumBlocks=Map.NumMyElements();
  Blockstart_=new int[NumBlocks+1];
  Blockids_=new int[NumBlocks];
  for(int i=0;i<NumBlocks;i++){
    Blockstart_[i]=i;
    if(local_ids) Blockids_[i]=i;
    else Blockids_[i]=Map.GID(i);
  }
  Blockstart_[NumBlocks]=NumBlocks;
}



//=============================================
//=============================================
//=============================================
double Test_PTBDP(const Epetra_CrsMatrix& MAT, int NumBlocks,int* Blockstart_,int* Blockids_,bool is_lid){
  // Build the block lists
  Teuchos::ParameterList List,Sublist;
  List.set("number of local blocks",NumBlocks);
  List.set("block start index",Blockstart_);
  if(is_lid) List.set("block entry lids",Blockids_);
  else List.set("block entry gids",Blockids_);

  Sublist.set("apply mode","invert");
  //Sublist.set("apply mode","multiply");
  List.set("blockdiagmatrix: list",Sublist);

  EpetraExt_PointToBlockDiagPermute Perm(MAT);
  Perm.SetParameters(List);

  Perm.Compute();
  Epetra_MultiVector X(MAT.RowMap(),1);
  Epetra_MultiVector Y(MAT.RowMap(),1);
  Epetra_MultiVector Z(MAT.RowMap(),1);
  X.SetSeed(24601); X.Random();

  double norm2;
  Perm.ApplyInverse(X,Y);
  MAT.Apply(Y,Z);
  X.Update(1.0,Z,-1.0);
  X.Norm2(&norm2);
  return norm2;
}


//=============================================
double Test_PTBDP_C(const Epetra_CrsMatrix& MAT,int BlockSize){
  // Build the block lists
  Teuchos::ParameterList List,Sublist;
  List.set("contiguous block size",BlockSize);

  Sublist.set("apply mode","invert");
  //Sublist.set("apply mode","multiply");
  List.set("blockdiagmatrix: list",Sublist);

  EpetraExt_PointToBlockDiagPermute Perm(MAT);
  Perm.SetParameters(List);

  Perm.Compute();
  Epetra_MultiVector X(MAT.RowMap(),1);
  Epetra_MultiVector Y(MAT.RowMap(),1);
  Epetra_MultiVector Z(MAT.RowMap(),1);
  X.SetSeed(24601); X.Random();

  double norm2;
  Perm.ApplyInverse(X,Y);
  MAT.Apply(Y,Z);
  X.Update(1.0,Z,-1.0);
  X.Norm2(&norm2);
  return norm2;
}

//=============================================
bool TestPointToBlockDiagPermute(const Epetra_Comm & Comm){
  const int NUM_ROWS=64;

  bool TestPassed=true;
  Epetra_CrsMatrix *MAT;
  int NumBlocks, *Blockstart_,*Blockids_;
  double norm2;

  // TEST #1 - Local, Contiguous
  Build_Local_Contiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_PTBDP(*MAT,NumBlocks,Blockstart_,Blockids_,true);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"P2BDP LCMAT    Error = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // TEST #2 - Local, Non-Contiguous
  Build_Local_NonContiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_PTBDP(*MAT,NumBlocks,Blockstart_,Blockids_,true);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"P2BDP LNCMat   Error = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // TEST #3 - Non-Local
  Build_NonLocal_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_PTBDP(*MAT,NumBlocks,Blockstart_,Blockids_,false);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"P2BDP NLMat    Error = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // TEST #4 - Local, Contiguous in ContiguousMode
  Build_Local_Contiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_PTBDP_C(*MAT,2);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"P2BDP LCMAT-C  Error = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;


  return TestPassed;
}


//=============================================
double Test_Cheby(const Epetra_CrsMatrix& MAT, int NumBlocks,int* Blockstart_,int* Blockids_,int maxits,bool is_lid){
  double norm2,norm0;
  // Build the block lists
  Teuchos::ParameterList ChebyList,List,Sublist;
  List.set("number of local blocks",NumBlocks);
  List.set("block start index",Blockstart_);
  if(is_lid) List.set("block entry lids",Blockids_);
  else List.set("block entry gids",Blockids_);

  Sublist.set("apply mode","invert");
  List.set("blockdiagmatrix: list",Sublist);

  ChebyList.set("chebyshev: use block mode",true);
  ChebyList.set("chebyshev: block list",List);
  ChebyList.set("chebyshev: eigenvalue autocompute ratio",30.0);//HAQ
  ChebyList.set("chebyshev: degree",maxits);

  // Build a Chebyshev
  Ifpack_Chebyshev Cheby(&MAT);
  Cheby.SetParameters(ChebyList);
  Cheby.Compute();

  Epetra_MultiVector X(MAT.RowMap(),1);
  Epetra_MultiVector Y(MAT.RowMap(),1);
  Epetra_MultiVector Z(MAT.RowMap(),1);
  X.SetSeed(24601); X.Random();
  MAT.Apply(X,Y);
  Y.Norm2(&norm0);

  Cheby.ApplyInverse(Y,Z);
  X.Update(1.0,Z,-1.0);
  X.Norm2(&norm2);
  return norm2 / norm0;
}

//=============================================
double Test_Cheby_C(const Epetra_CrsMatrix& MAT, int BlockSize,int maxits){
  double norm2,norm0;
  // Build the block lists
  Teuchos::ParameterList ChebyList,List,Sublist;
  List.set("contiguous block size",BlockSize);
  Sublist.set("apply mode","invert");
  List.set("blockdiagmatrix: list",Sublist);

  ChebyList.set("chebyshev: use block mode",true);
  ChebyList.set("chebyshev: block list",List);
  ChebyList.set("chebyshev: eigenvalue autocompute ratio",30.0);//HAQ
  ChebyList.set("chebyshev: degree",maxits);

  // Build a Chebyshev
  Ifpack_Chebyshev Cheby(&MAT);
  Cheby.SetParameters(ChebyList);
  Cheby.Compute();

  Epetra_MultiVector X(MAT.RowMap(),1);
  Epetra_MultiVector Y(MAT.RowMap(),1);
  Epetra_MultiVector Z(MAT.RowMap(),1);
  X.SetSeed(24601); X.Random();
  MAT.Apply(X,Y);
  Y.Norm2(&norm0);

  Cheby.ApplyInverse(Y,Z);
  X.Update(1.0,Z,-1.0);
  X.Norm2(&norm2);
  return norm2 / norm0;
}


//=============================================
bool TestBlockChebyshev(const Epetra_Comm & Comm){
  const int NUM_ROWS=100;

  bool TestPassed=true;
  Epetra_CrsMatrix *MAT;
  int NumBlocks, *Blockstart_,*Blockids_;
  double norm2;

  // Test #1 - Local, Contiguous matrix w/ diagonal precond
  Build_Local_Contiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  delete [] Blockstart_; delete [] Blockids_;
  Build_DiagonalStructure(MAT->RowMap(),NumBlocks,Blockstart_,Blockids_,true);
  norm2=Test_Cheby(*MAT,NumBlocks,Blockstart_,Blockids_,100,true);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby LC-D   nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // Test #2 - Local, Non-Contiguous matrix w/ diagonal precond
  Build_Local_NonContiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  delete [] Blockstart_; delete [] Blockids_;
  Build_DiagonalStructure(MAT->RowMap(),NumBlocks,Blockstart_,Blockids_,true);
  norm2=Test_Cheby(*MAT,NumBlocks,Blockstart_,Blockids_,100,true);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby LNC-D  nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // TEST #3 - Non-Local matrix w/ diagonal precond
  Build_NonLocal_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  delete [] Blockstart_; delete [] Blockids_;
  Build_DiagonalStructure(MAT->RowMap(),NumBlocks,Blockstart_,Blockids_,false);
  norm2=Test_Cheby(*MAT,NumBlocks,Blockstart_,Blockids_,100,false);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby NL-D   nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // Test #4 - Local, Contiguous matrix w/ exact precond
  Build_Local_Contiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_Cheby(*MAT,NumBlocks,Blockstart_,Blockids_,1,true);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby LC-E   nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // Test #5 - Local, Non-Contiguous matrix w/ exact precond
  Build_Local_NonContiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_Cheby(*MAT,NumBlocks,Blockstart_,Blockids_,1,true);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby LNC-E  nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // TEST #6 - Non-Local matrix w/ exact precond
  Build_NonLocal_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_Cheby(*MAT,NumBlocks,Blockstart_,Blockids_,1,false);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby NL-E   nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // Test #7 - Local, Contiguous matrix w/ diagonal precond (contiguous mode)
  Build_Local_Contiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_Cheby_C(*MAT,1,100);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby LC-Dc  nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  // Test #8 - Local, Contiguous matrix w/ exact precond (contiguous mode)
  Build_Local_Contiguous_Size2_BlockMatrix(Comm,NUM_ROWS,NumBlocks,Blockstart_,Blockids_,MAT);
  norm2=Test_Cheby_C(*MAT,2,1);
  if(norm2 > 1e-12) TestPassed=false;
  if(!Comm.MyPID()) cout<<"Cheby LC-Ec  nrm-red = "<<norm2<<endl;
  delete MAT; delete [] Blockstart_; delete [] Blockids_;

  return TestPassed;
}

//=============================================
//=============================================
//=============================================
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  bool TestPassed=true;

  // BlockDiagMatrix test
  TestPassed=TestPassed && TestBlockDiagMatrix(Comm);

  // PointToBlockDiagPermute tests
  TestPassed=TestPassed && TestPointToBlockDiagPermute(Comm);

  // Block Chebyshev Tests
  TestPassed=TestPassed && TestBlockChebyshev(Comm);

  // ============ //
  // final output //
  // ============ //

  if (!TestPassed) {
    if(!Comm.MyPID()) cout << "Test `BlockCheby.exe' FAILED!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  if(!Comm.MyPID()) cout << "Test `BlockCheby.exe' passed!" << endl;
  exit(EXIT_SUCCESS);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  puts("please configure IFPACK with --eanble-aztecoo --enable-teuchos --enable-epetraext");
  puts("to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
