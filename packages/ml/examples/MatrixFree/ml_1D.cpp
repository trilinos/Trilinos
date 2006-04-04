//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#include "ml_include.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "AztecOO.h"
#include "ml_MatrixFreePreconditioner.h"
#include "ml_ElementByElement_SingleElement.h"
#include "ml_epetra_utils.h"

#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
using namespace Galeri;



using namespace ML_Epetra;
using namespace Teuchos;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumMyFEs = 4;
  int NumGlobalFEs = NumMyFEs * Comm.NumProc() + 1;
  int NumVerticesPerFE = 2;
  int NumPDEEqns = 2;
  vector<int> MyFEs(NumMyFEs * NumVerticesPerFE);

  // populate the distributed 1D finite elements. Each
  // processor has 4 elements; MyFEs must contain the
  // global IDs of the vertices.

  int offset = Comm.MyPID() * NumMyFEs;

  for (int ie = 0; ie < NumMyFEs; ++ie)
  {
    MyFEs[ie * NumVerticesPerFE    ] = offset + ie;
    MyFEs[ie * NumVerticesPerFE + 1] = offset + ie + 1;
  }

  int NumMyBoundaryRows = 0;
  vector<int>    MyBoundaryRows;
  vector<double> MyBoundaryValues;
  
  if (Comm.MyPID() == 0)
  {
    NumMyBoundaryRows += 2;
    MyBoundaryRows.push_back(0);
    MyBoundaryRows.push_back(1);
    MyBoundaryValues.push_back(0.0);
    MyBoundaryValues.push_back(0.0);
  }

  if (Comm.MyPID() == Comm.NumProc() - 1)
  {
    NumMyBoundaryRows += 2;
    MyBoundaryRows.push_back(NumPDEEqns * NumGlobalFEs - 1);
    MyBoundaryRows.push_back(NumPDEEqns * NumGlobalFEs - 2);
    MyBoundaryValues.push_back(0.0);
    MyBoundaryValues.push_back(0.0);
  }
    
  // this is the local finite element matrix
  
  Epetra_SerialDenseMatrix FEMatrix(NumPDEEqns * NumVerticesPerFE, NumPDEEqns * NumVerticesPerFE);
  FEMatrix(0,0) =  1.0; FEMatrix(0,1) =  0.0; FEMatrix(0,2) = -1.0; FEMatrix(0,3) =  0.0;
  FEMatrix(1,0) =  0.0; FEMatrix(1,1) =  1.0; FEMatrix(1,2) =  0.0; FEMatrix(1,3) = -1.0;
  FEMatrix(2,0) = -1.0; FEMatrix(2,1) =  0.0; FEMatrix(2,2) =  1.0; FEMatrix(2,3) =  0.0;
  FEMatrix(3,0) =  0.0; FEMatrix(3,1) = -1.0; FEMatrix(3,2) =  0.0; FEMatrix(3,3) =  1.0;

  // still need to specify the map for the graph. This requires
  // knowledge of the global number of vertices
  
  Epetra_Map GraphMap(NumGlobalFEs, 0, Comm);

  // we can now build the matrix

  ElementByElement_SingleElement EBE(Comm, NumMyFEs, NumVerticesPerFE, &MyFEs[0],
                                     NumPDEEqns, &FEMatrix, NumMyBoundaryRows, 
                                     &MyBoundaryRows[0], &MyBoundaryValues[0],
                                     GraphMap, 0);

  Epetra_Vector X(EBE.OperatorDomainMap());
  Epetra_Vector Y(EBE.OperatorDomainMap());
  Epetra_Vector Y2(EBE.OperatorDomainMap());

  X.Random();
  EBE.SetMyBoundaryRows(X);

  EBE.Apply(X, Y);
  EBE.ResetMyBoundaryRows(Y);

  ParameterList GaleriList;
  GaleriList.set("n", NumGlobalFEs);
  GaleriList.set("m", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Linear", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Laplace1D", Map, GaleriList);
  Epetra_VbrMatrix* A = CreateVbrMatrix(CrsA, NumPDEEqns);

  A->Apply(X, Y2);
  EBE.ResetMyBoundaryRows(Y2);

  double norm;
  Y.Update(1.0, Y2, -1.0);
  Y.Norm2(&norm);

  if (Comm.MyPID() == 0)
    cout << "||Y_EBE - Y_VBR||_2 = " << norm << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
