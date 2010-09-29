//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
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

  int blocksize=3;
  int num_local_blocks=2;

  // Generate the rowmap
  Epetra_Map Map(num_local_blocks*blocksize*Comm.NumProc(),0,Comm);

  // Generate a non-symmetric blockdiagonal matrix, where the blocks are neatly processor-aligned
  Epetra_CrsMatrix Matrix(Copy,Map,0);
  for(int i=0;i<Map.NumMyElements();i++){
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
  double norminf=Pmat->NormInf();

  // Cleanup
  delete Pmat;

  // passing check
  if(norminf > 1e-15){
    if (Comm.MyPID()==0) cout << "EpetraExt:: PointToBlockDiagPermute tests FAILED." << endl;
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }
  else{
    if (Comm.MyPID()==0) cout << "EpetraExt:: PointToBlockDiagPermute tests passed." << endl;
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
  }

  return 0;
}


