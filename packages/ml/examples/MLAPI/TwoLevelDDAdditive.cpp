
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

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "ml_include.h"
#include "MLAPI.h"

using namespace Teuchos;
using namespace Trilinos_Util;
using namespace MLAPI;

// ======================================= //
// 2-level additive Schwarz preconditioner //
// ======================================= //

class TwoLevelDDAdditive : public Preconditioner {

public:
  // Constructor assumes that all operators and inverse operators are already
  // filled.
  TwoLevelDDAdditive(const Operator FineMatrix, 
                     const InverseOperator FineSolver, 
                     const InverseOperator CoarseSolver, 
                     const Operator R, 
                     const Operator P) :
    FineMatrix_(FineMatrix),
    R_(R),
    P_(P),
    FineSolver_(FineSolver),
    CoarseSolver_(CoarseSolver)
  {}
      
  int Solve(const DoubleVector& r_f, DoubleVector& x_f) const
  {
    
    DoubleVector r_c(FineSolver_.DomainSpace());

    // apply fine level preconditioner
    x_f = FineSolver_ * r_f;
    // restrict to coarse
    r_c = R_ * r_f;
    // solve coarse problem
    r_c = CoarseSolver_ * r_c;
    // prolongate back and add to solution
    x_f = x_f + P_ * r_c;

    return(0);
  }

  const Space DomainSpace() const {
    return(FineMatrix_.DomainSpace());
  }

  const Space RangeSpace() const {
    return(FineMatrix_.RangeSpace());
  }

  ostream& Print(ostream& os, const bool verbose = true) const 
  {
    if (MyPID() == 0) {
      os << "*** MLAPI::TwoLevelAdditiveSA ***" << endl;
    }

    return(os);
  }

private:
  const Operator FineMatrix_;
  const Operator R_;
  const Operator P_;
  const InverseOperator FineSolver_;
  const InverseOperator CoarseSolver_;

}; // TwoLevelDDAdditive

#include "Trilinos_Util_CrsMatrixGallery.h"
// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_Comm Comm(MPI_COMM_WORLD);
#endif

  // Initialize the workspace and set the output level
  Init();

  try {

    int NumGlobalElements = 10000;
    // define the space for fine level vectors and operators.
    Space FineSpace(NumGlobalElements);

    // define the linear system matrix, solution and RHS
    //Operator FineMatrix = Gallery("laplace_2d", FineSpace);

  Trilinos_Util::CrsMatrixGallery Gallery("laplace_2d", GetEpetra_Comm());
  Gallery.Set("problem_size", NumGlobalElements);
  Operator FineMatrix(FineSpace,FineSpace,*(Gallery.GetMatrix()));

    DoubleVector LHS(FineSpace);
    DoubleVector RHS(FineSpace);

    LHS = 0.0;
    RHS.Random();

    Teuchos::ParameterList MLList;

    // CoarseMatrix will contain the coarse level matrix,
    // while FineSolver and CoarseSolver will contain
    // the fine level smoother and the coarse level solver,
    // respectively. We will use symmetric Gauss-Seidel
    // for the fine level, and Amesos (LU) for the coarse level.
    Operator CoarseMatrix;
    InverseOperator FineSolver, CoarseSolver;

    // Now we define restriction (R) and prolongator (P) from the fine space
    // to the coarse space using non-smoothed aggregation.
    // The coarse-level matrix will be defined via a triple
    // matrix-matrix product.
    Operator R, P;

#ifdef HAVE_ML_METIS
    MLList.set("aggregation: type","METIS");
    MLList.set("aggregation: nodes per aggregate",64);
#else
    MLList.set("aggregation: type","Uncoupled");
#endif

    AggregationDataBase  ADB(MLList);
    SmootherDataBase     SDB(MLList);
    CoarseSolverDataBase CDB(MLList);

    BuildPtent(FineMatrix,ADB,P);
    R = Transpose(P);

    CoarseMatrix = RAP(R,FineMatrix,P);
    FineSolver.Reshape(FineMatrix,SDB);
    CoarseSolver.Reshape(CoarseMatrix,CDB);

    // We can now construct a Preconditioner-derived object, that
    // implements the 2-level hybrid domain decomposition preconditioner.
    // Preconditioner `TwoLevelDDHybrid' can be replaced by
    // `TwoLevelDDAdditive' to define an purely additive preconditioner.
    TwoLevelDDAdditive MLAPIPrec(FineMatrix,FineSolver,CoarseSolver,R,P);

    MLList.set("max iterations", 1550);
    MLList.set("tolerance", 1e-9);
    MLList.set("solver",  "gmres");
    MLList.set("output", 16);

    KrylovDataBase KDB; // set default options
    Krylov(FineMatrix, LHS, RHS, MLAPIPrec, KDB);
    
  }
  catch (const char e[]) {
    cerr << "Caught exception: " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

#ifdef HAVE_MPI
  // finalize the MLAPI workspace
  Finalize();
  MPI_Finalize() ;
#endif

  return 0 ;
  
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
