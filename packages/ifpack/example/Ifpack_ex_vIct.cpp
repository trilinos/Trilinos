// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AZTECOO) && defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS) && defined(HAVE_IFPACK_TRIUTILS)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_vIct.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[])
{

  // initialize MPI and Epetra communicator
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // size of the global matrix (must be a square number)
  const int NumPoints = 90000;

  // build the matrix corresponding to a 2D Laplacian on a
  // structured grid.
  CrsMatrixGallery Gallery("laplace_2d_bc", Comm);
  Gallery.Set("problem_size", NumPoints);
  // for simplicity, linear map.
  Gallery.Set("map_type", "linear");

  // get the pointer to the linear system matrix
  Epetra_RowMatrix* A = Gallery.GetMatrix();

  // =============================================================== //
  // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
  // =============================================================== //

  Teuchos::ParameterList List;

  // builds an Ifpack_AdditiveSchwarz. This is templated with
  // the local solvers, in this case Ifpack_Ict. Note that any
  // other Ifpack_Preconditioner-derived class can be used
  // instead of Ifpack_Ict.

  // In this example the overlap is zero. Use
  // Prec(A,OverlapLevel) for the general case.
  Ifpack_AdditiveSchwarz<Ifpack_vIct> Prec(A);

  List.set("fact: level-of-fill", 5);
  List.set("fact: absolute threshold", 0.0);
  List.set("fact: relative threshold", 0.0);

  // sets the parameters
  IFPACK_CHK_ERR(Prec.SetParameters(List));

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  // At this call, Amesos will perform the symbolic factorization.
  IFPACK_CHK_ERR(Prec.Initialize());

  // Builds the preconditioners, by looking for the values of 
  // the matrix. At this call, Amesos will perform the
  // numeric factorization.
  IFPACK_CHK_ERR(Prec.Compute());

  // =================================================== //
  // E N D   O F   I F P A C K   C O N S T R U C T I O N //
  // =================================================== //

  // At this point, we need some additional objects
  // to define and solve the linear system.

  // defines LHS and RHS
  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorDomainMap());

  LHS.PutScalar(0.0);
  RHS.Random();

  // need an Epetra_LinearProblem to define AztecOO solver
  Epetra_LinearProblem Problem(A,&LHS,&RHS);

  // now we can allocate the AztecOO solver
  AztecOO Solver(Problem);

  // specify solver
  Solver.SetAztecOption(AZ_solver,AZ_cg);
  Solver.SetAztecOption(AZ_output,32);

  // HERE WE SET THE IFPACK PRECONDITIONER
  Solver.SetPrecOperator(&Prec);

  // .. and here we solve
  // NOTE: with one process, the solver must converge in
  // one iteration.
  Solver.Iterate(1550,1e-8);

  // Prints out several information about the preconditioner
  // Note that methods NumInitialize(), NumCompute(),
  // NumApplyInverse(), InitializeTime(), ComputeTime(),
  // ApplyInverseTime() can be used to extract information
  // from the Prec object.
  cout << *(Prec.Inverse());

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

    return(EXIT_SUCCESS);
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

  puts("please configure IFPACK with --eanble-aztecoo --enable-teuchos");
  puts("--enable-amesos --enable-triutils to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
