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

#include "Ifpack.h"
#include "Ifpack_Hypre.h"

#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"

#include "Epetra_MultiVector.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_InvOperator.h"

#include "BelosEpetraAdapter.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "hypre_Helpers.hpp"

#include <string>
#include <stdio.h>
#include <map>


TEUCHOS_UNIT_TEST(Belos_Hypre, Laplace2D){
  const double tol = 1E-7;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  typedef Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>  LinearProblem;

  //
  // Create Laplace2D
  //
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm();
#endif
  Teuchos::ParameterList GaleriList;
  int nx = 10 * Comm.NumProc();
  int ny = 10 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  Epetra_Map Map(nx*ny,0,Comm);
  RCP<Epetra_CrsMatrix> Crs_Matrix   = rcp(Galeri::CreateCrsMatrix("Laplace2D", &Map, GaleriList));
  int NumProc = Crs_Matrix->Comm().NumProc();

  //
  // Create the hypre preconditioner
  //
  RCP<Ifpack_Hypre> preconditioner = rcp(new Ifpack_Hypre(Crs_Matrix.get()));
  TEST_EQUALITY(preconditioner->Initialize(),0);
  TEST_EQUALITY(preconditioner->SetParameter(Preconditioner, ParaSails),0); // Use a Euclid Preconditioner (but not really used)
  TEST_EQUALITY(preconditioner->SetParameter(Preconditioner),0); // Solve the problem
  TEST_EQUALITY(preconditioner->Compute(),0);

  //
  // Create the solution vector and rhs
  //
  int numVec = 1;
  RCP<Epetra_MultiVector> X = rcp(new Epetra_MultiVector(Crs_Matrix->OperatorDomainMap(), numVec));
  RCP<Epetra_MultiVector> KnownX = rcp(new Epetra_MultiVector(Crs_Matrix->OperatorDomainMap(), numVec));
  KnownX->Random();
  RCP<Epetra_MultiVector> B = rcp(new Epetra_MultiVector(Crs_Matrix->OperatorRangeMap(), numVec));
  Crs_Matrix->Apply(*KnownX, *B);

  //
  // Test the EpetraExt wrapper
  // amk November 24, 2015: Should we deprecate this?
  //
//  RCP<ParameterList> pl = rcp(new ParameterList());
//  TEST_EQUALITY(X->PutScalar(0.0),0);
//  HYPRE_IJMatrix hypre_mat = preconditioner->HypreMatrix();
//  RCP<EpetraExt_HypreIJMatrix> Hyp_Matrix = rcp(new EpetraExt_HypreIJMatrix(hypre_mat));
//  TEST_EQUALITY(Hyp_Matrix->SetParameter(Preconditioner, ParaSails),0);
//  TEST_EQUALITY(Hyp_Matrix->SetParameter(Preconditioner),0);
//  TEST_EQUALITY(EquivalentMatrices(*Hyp_Matrix, *Crs_Matrix, tol), true);
//  RCP<LinearProblem> problem1 = rcp(new LinearProblem(Crs_Matrix,X,B));
//  problem1->setLeftPrec(Hyp_Matrix);
//  TEST_EQUALITY(problem1->setProblem(),true);
//  Belos::PseudoBlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator> solMgr1(problem1,pl);
//  Belos::ReturnType rv1 = solMgr1.solve(); // TEST_EQUALITY(solMgr2.solve(),Belos::Converged);
//  TEST_EQUALITY(rv1,Belos::Converged);
//  TEST_EQUALITY(EquivalentVectors(*X, *KnownX, tol*10*pow(10.0,NumProc)), true);

  //
  // Test the Ifpack hypre interface
  //
  RCP<ParameterList> pl2 = rcp(new ParameterList());
  RCP<Epetra_Operator> invOp = rcp(new Epetra_InvOperator(preconditioner.get()));
  TEST_EQUALITY(X->PutScalar(0.0),0);
  RCP<LinearProblem> problem2 = rcp(new LinearProblem(Crs_Matrix,X,B));
  problem2->setLeftPrec(invOp);
  TEST_EQUALITY(problem2->setProblem(),true);
  Belos::PseudoBlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator> solMgr2(problem2,pl2);
  Belos::ReturnType rv2 = solMgr2.solve(); // TEST_EQUALITY(solMgr2.solve(),Belos::Converged);
  TEST_EQUALITY(rv2,Belos::Converged);
  TEST_EQUALITY(EquivalentVectors(*X, *KnownX, tol*10*pow(10.0,NumProc)), true);
}
