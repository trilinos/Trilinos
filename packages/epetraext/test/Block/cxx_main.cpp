//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Flops.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "EpetraExt_PointToBlockDiagPermute.h"
#include "Teuchos_ParameterList.hpp"
#include "EpetraExt_MatrixMatrix.h"


int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif
  double total_norm=0;

  int blocksize=3;
  int num_local_blocks=2;

  // Generate the rowmap
  Epetra_Map Map(num_local_blocks*blocksize*Comm.NumProc(),0,Comm);
  int Nrows=Map.NumMyElements();

  // Generate a non-symmetric blockdiagonal matrix, where the blocks are neatly processor-aligned
  Epetra_CrsMatrix Matrix(Copy,Map,0);
  for(int i=0;i<Nrows;i++){
    int gid=Map.GID(i);
    int gidm1=gid-1;
    int gidp1=gid+1;
    double diag=10+gid;
    double left=-1;
    double right=-2;
    if(i%blocksize > 0 ) Matrix.InsertGlobalValues(gid,1,&left,&gidm1);
    Matrix.InsertGlobalValues(gid,1,&diag,&gid);
    if(i%blocksize < 2 ) Matrix.InsertGlobalValues(gid,1,&right,&gidp1);
  }    
  Matrix.FillComplete();

  // *********************************************************************
  // Test #1:  Blocks respect initial ordering, no proc boundaries crossed
  // *********************************************************************
  {
    // Build the block diagonalizer
    Teuchos::ParameterList List;
    List.set("contiguous block size",blocksize);
    List.set("number of local blocks",num_local_blocks);

    EpetraExt_PointToBlockDiagPermute Permute(Matrix);
    Permute.SetParameters(List);
    Permute.Compute();
    
    Epetra_FECrsMatrix* Pmat=Permute.CreateFECrsMatrix();
    
    // Multiply matrices, compute difference
    Epetra_CrsMatrix Res(Copy,Map,0);
    EpetraExt::MatrixMatrix::Multiply(*Pmat,false,Matrix,false,Res);
    EpetraExt::MatrixMatrix::Add(Matrix,false,1.0,*Pmat,-1.0);
    total_norm+=Pmat->NormInf();
    
    // Cleanup
    delete Pmat;
  }

  // *********************************************************************
  // Test #2:  Blocks do not respect initial ordering, lids
  // *********************************************************************
  {
    // Build alternative list - just have each block reversed in place
    int* block_lids=new int [Nrows];
    int* block_starts=new int[num_local_blocks+1];
    for(int i=0;i<num_local_blocks;i++){
      block_starts[i]=i*blocksize;
      for(int j=0;j<blocksize;j++){
	block_lids[i*blocksize+j] = i*blocksize+(blocksize-j-1);
      }
      
    }
    block_starts[num_local_blocks]=Nrows;
    
    // Build the block diagonalizer
    Teuchos::ParameterList List;
    List.set("number of local blocks",num_local_blocks);
    List.set("block start index",block_starts);
    List.set("block entry lids",block_lids);
    
    EpetraExt_PointToBlockDiagPermute Permute(Matrix);
    Permute.SetParameters(List);
    Permute.Compute();

    Epetra_FECrsMatrix* Pmat=Permute.CreateFECrsMatrix();

    // Multiply matrices, compute difference
    Epetra_CrsMatrix Res(Copy,Map,0);
    EpetraExt::MatrixMatrix::Multiply(*Pmat,false,Matrix,false,Res);
    EpetraExt::MatrixMatrix::Add(Matrix,false,1.0,*Pmat,-1.0);
    total_norm+=Pmat->NormInf();
    
    // Cleanup
    delete Pmat;
    delete [] block_lids;
    delete [] block_starts;
  }


  // *********************************************************************
  // Test #3:  Blocks do not respect initial ordering, gids
  // *********************************************************************
  {
    // Build alternative list - just have each block reversed in place
    int* block_gids=new int [Nrows];
    int* block_starts=new int[num_local_blocks+1];
    for(int i=0;i<num_local_blocks;i++){
      block_starts[i]=i*blocksize;
      for(int j=0;j<blocksize;j++){
	block_gids[i*blocksize+j] = Map.GID(i*blocksize+(blocksize-j-1));
      }
      
    }
    block_starts[num_local_blocks]=Nrows;
    
    // Build the block diagonalizer
    Teuchos::ParameterList List;
    List.set("number of local blocks",num_local_blocks);
    List.set("block start index",block_starts);
    List.set("block entry gids",block_gids);
    
    EpetraExt_PointToBlockDiagPermute Permute(Matrix);
    Permute.SetParameters(List);
    Permute.Compute();

    Epetra_FECrsMatrix* Pmat=Permute.CreateFECrsMatrix();

    // Multiply matrices, compute difference
    Epetra_CrsMatrix Res(Copy,Map,0);
    EpetraExt::MatrixMatrix::Multiply(*Pmat,false,Matrix,false,Res);
    EpetraExt::MatrixMatrix::Add(Matrix,false,1.0,*Pmat,-1.0);
    total_norm+=Pmat->NormInf();
    
    // Cleanup
    delete Pmat;
    delete [] block_gids;
    delete [] block_starts;
  }



  // passing check
  if(total_norm > 1e-15){
    if (Comm.MyPID()==0) cout << "EpetraExt:: PointToBlockDiagPermute tests FAILED (||res||="<<total_norm<<")." << endl;
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }
  else{
    if (Comm.MyPID()==0) cout << "EpetraExt:: PointToBlockDiagPermute tests passed (||res||="<<total_norm<<")." << endl;
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
  }

  return 0;
}


