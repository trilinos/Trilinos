
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_common.h"
#if defined(HAVE_ML_MLAPI) && defined(HAVE_ML_GALERI)

#include "Teuchos_ParameterList.hpp"
#include "ml_include.h"
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_Krylov.h"

using namespace Teuchos;
using namespace MLAPI;

// ======================================= //
// 2-level additive Schwarz preconditioner //
// ======================================= //

class TwoLevelDDAdditive : public BaseOperator {

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
      
  TwoLevelDDAdditive(const TwoLevelDDAdditive& rhs) :
    FineMatrix_(rhs.GetFineMatrix()),
    R_(rhs.GetR()),
    P_(rhs.GetP()),
    FineSolver_(rhs.GetFineSolver()),
    CoarseSolver_(rhs.GetCoarseSolver())
  {}
      
  TwoLevelDDAdditive& operator=(const TwoLevelDDAdditive& rhs)
  {
    if (this != &rhs) {
      FineMatrix_ = rhs.GetFineMatrix();
      R_ = rhs.GetR();
      P_ = rhs.GetP();
      FineSolver_ = rhs.GetFineSolver();
      CoarseSolver_ = rhs.GetCoarseSolver();
    }
    return(*this);
  }
    
  int Apply(const MultiVector& r_f, MultiVector& x_f) const
  {
    
    MultiVector r_c(FineSolver_.GetDomainSpace());

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

  const Space GetOperatorDomainSpace() const {
    return(FineMatrix_.GetDomainSpace());
  }

  const Space GetOperatorRangeSpace() const {
    return(FineMatrix_.GetRangeSpace());
  }

  const Space GetDomainSpace() const {
    return(FineMatrix_.GetDomainSpace());
  }

  const Space GetRangeSpace() const {
    return(FineMatrix_.GetRangeSpace());
  }

  const Operator GetFineMatrix() const 
  {
    return(FineMatrix_);
  }

  const Operator GetR() const 
  {
    return(R_);
  }

  const Operator GetP() const 
  {
    return(P_);
  }

  const InverseOperator GetFineSolver() const 
  {
    return(FineSolver_);
  }

  const InverseOperator GetCoarseSolver() const 
  {
    return(CoarseSolver_);
  }

  ostream& Print(ostream& os, const bool verbose = true) const 
  {
    if (GetMyPID() == 0) {
      os << "*** MLAPI::TwoLevelDDAdditive ***" << endl;
    }

    return(os);
  }

private:
  Operator FineMatrix_;
  Operator R_;
  Operator P_;
  InverseOperator FineSolver_;
  InverseOperator CoarseSolver_;

}; // TwoLevelDDAdditive

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Initialize the workspace and set the output level
  Init();

  try {

    int NumGlobalElements = 10000;
    // define the space for fine level vectors and operators.
    Space FineSpace(NumGlobalElements);

    // define the linear system matrix, solution and RHS
    Operator FineMatrix = Gallery("Laplace2D", FineSpace);

    MultiVector LHS(FineSpace);
    MultiVector RHS(FineSpace);

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

    GetPtent(FineMatrix,MLList,P);
    R = GetTranspose(P);

    CoarseMatrix = GetRAP(R,FineMatrix,P);
    FineSolver.Reshape(FineMatrix,"symmetric Gauss-Seidel", MLList);
    CoarseSolver.Reshape(CoarseMatrix,"Amesos-KLU", MLList);

    // We can now construct a Preconditioner-derived object, that
    // implements the 2-level hybrid domain decomposition preconditioner.
    // Preconditioner `TwoLevelDDHybrid' can be replaced by
    // `TwoLevelDDAdditive' to define an purely additive preconditioner.
    TwoLevelDDAdditive MLAPIPrec(FineMatrix,FineSolver,CoarseSolver,R,P);

    MLList.set("krylov: max iterations", 1550);
    MLList.set("krylov: tolerance", 1e-9);
    MLList.set("krylov: type",  "gmres");
    MLList.set("krylov: output", 16);

    Krylov(FineMatrix, LHS, RHS, MLAPIPrec, MLList);
    
  }
  catch (const int e) {
    cerr << "Caught integer exception, code = " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

#ifdef HAVE_MPI
  // finalize the MLAPI workspace
  Finalize();
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
  
}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("This MLAPI example requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("\t--enable-galeri");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif // #if defined(HAVE_ML_MLAPI)
