//@HEADER
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

#include "Ifpack_ConfigDefs.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "AztecOO.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_ICT.h"

int main(int argc, char *argv[])
{

  // initialize MPI and Epetra communicator
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList GaleriList;

  // The problem is defined on a 2D grid, global size is nx * nx.
  int nx = 30;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Cartesian2D", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_RowMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );

  // =============================================================== //
  // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
  // =============================================================== //

  Teuchos::ParameterList List;

  // builds an Ifpack_AdditiveSchwarz. This is templated with
  // the local solvers, in this case Ifpack_ICT. Note that any
  // other Ifpack_Preconditioner-derived class can be used
  // instead of Ifpack_ICT.

  // In this example the overlap is zero. Use
  // Prec(A,OverlapLevel) for the general case.
  Ifpack_AdditiveSchwarz<Ifpack_ICT> Prec(&*A);

  // `1.0' means that the factorization should approximatively
  // keep the same number of nonzeros per row of the original matrix.
  List.set("fact: ict level-of-fill", 1.0);
  // no modifications on the diagonal
  List.set("fact: absolute threshold", 0.0);
  List.set("fact: relative threshold", 1.0);
  List.set("fact: relaxation value", 0.0);
  // matrix `laplace_2d_bc' is not symmetric because of the way
  // boundary conditions are imposed. We can filter the singletons,
  // (that is, Dirichlet nodes) and end up with a symmetric
  // matrix (as ICT requires).
  List.set("schwarz: filter singletons", true);

  // sets the parameters
  IFPACK_CHK_ERR(Prec.SetParameters(List));

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  IFPACK_CHK_ERR(Prec.Initialize());

  // Builds the preconditioners, by looking for the values of
  // the matrix.
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
  Epetra_LinearProblem Problem(&*A,&LHS,&RHS);

  // now we can allocate the AztecOO solver
  AztecOO Solver(Problem);

  // specify solver
  Solver.SetAztecOption(AZ_solver,AZ_cg_condnum);
  Solver.SetAztecOption(AZ_output,32);

  // HERE WE SET THE IFPACK PRECONDITIONER
  Solver.SetPrecOperator(&Prec);

  // .. and here we solve
  // NOTE: with one process, the solver must converge in
  // one iteration.
  Solver.Iterate(1550,1e-5);

  // Prints out some information about the preconditioner
  std::cout << Prec;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}
