// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SHIFTEDLAPLACIAN_DEF_HPP
#define MUELU_SHIFTEDLAPLACIAN_DEF_HPP

#include "MueLu_ShiftedLaplacian_decl.hpp"

#if defined(HAVE_MUELU_IFPACK2) and defined(HAVE_MUELU_TPETRA)

#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_CoupledRBMFactory.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Ifpack2Smoother.hpp>
#include <MueLu_PFactory.hpp>
#include <MueLu_PgPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_RAPShiftFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_ShiftedLaplacian.hpp>
#include <MueLu_ShiftedLaplacianOperator.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_SmootherPrototype.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLu {

// Destructor
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~ShiftedLaplacian() {}

// Input
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameters(Teuchos::RCP< Teuchos::ParameterList > paramList) {

  // Parameters
  coarseGridSize_      = paramList->get("MueLu: coarse size", 1000);
  numLevels_           = paramList->get("MueLu: levels",         3);
  int stype            = paramList->get("MueLu: smoother",       8);
  if(stype==1)         { Smoother_="jacobi";                      }
  else if(stype==2)    { Smoother_="gauss-seidel";                }
  else if(stype==3)    { Smoother_="symmetric gauss-seidel";      }
  else if(stype==4)    { Smoother_="chebyshev";                   }
  else if(stype==5)    { Smoother_="krylov";                      }
  else if(stype==6)    { Smoother_="ilut";                        }
  else if(stype==7)    { Smoother_="riluk";                       }
  else if(stype==8)    { Smoother_="schwarz";                     }
  else if(stype==9)    { Smoother_="superilu";                    }
  else if(stype==10)   { Smoother_="superlu";                     }
  else                 { Smoother_="schwarz";                     }
  smoother_sweeps_     = paramList->get("MueLu: sweeps",         5);
  smoother_damping_    = paramList->get("MueLu: relax val",    1.0);
  ncycles_             = paramList->get("MueLu: cycles",         1);
  iters_               = paramList->get("MueLu: iterations",   500);
  solverType_          = paramList->get("MueLu: solver type",    1);
  restart_size_        = paramList->get("MueLu: restart size", 100);
  recycle_size_        = paramList->get("MueLu: recycle size",  25);
  isSymmetric_         = paramList->get("MueLu: symmetric",   true);
  ilu_leveloffill_     = paramList->get("MueLu: level-of-fill",  5);
  ilu_abs_thresh_      = paramList->get("MueLu: abs thresh",   0.0);
  ilu_rel_thresh_      = paramList->get("MueLu: rel thresh",   1.0);
  ilu_diagpivotthresh_ = paramList->get("MueLu: piv thresh",   0.1);
  ilu_drop_tol_        = paramList->get("MueLu: drop tol",    0.01);
  ilu_fill_tol_        = paramList->get("MueLu: fill tol",    0.01);
  schwarz_overlap_     = paramList->get("MueLu: overlap",        0);
  schwarz_usereorder_  = paramList->get("MueLu: use reorder", true);
  int combinemode      = paramList->get("MueLu: combine mode",   1);
  if(combinemode==0)   { schwarz_combinemode_ = Tpetra::ZERO;     }
  else                 { schwarz_combinemode_ = Tpetra::ADD;      }
  tol_                 = paramList->get("MueLu: tolerance",  0.001);

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setProblemMatrix(RCP<Matrix>& A) {

  A_=A;
  if(A_!=Teuchos::null)
    TpetraA_ = Utilities::Op2NonConstTpetraCrs(A_);
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  if(LinearProblem_!=Teuchos::null)
    LinearProblem_ -> setOperator ( TpetraA_ );
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "ShiftedLaplacian only available with Tpetra and GO=int enabled.");
#endif

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setProblemMatrix(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO> >& TpetraA) {

  TpetraA_=TpetraA;
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  if(LinearProblem_!=Teuchos::null)
    LinearProblem_ -> setOperator ( TpetraA_ );
#endif

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setPreconditioningMatrix(RCP<Matrix>& P) {

  P_=P;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setPreconditioningMatrix(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO> >& TpetraP) {

  RCP< Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Atmp
    = rcp( new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraP) );
  P_= rcp( new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Atmp) );
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setstiff(RCP<Matrix>& K) {

  K_=K;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setstiff(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO> >& TpetraK) {

  RCP< Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Atmp
    = rcp( new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraK) );
  K_= rcp( new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Atmp) );

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setmass(RCP<Matrix>& M) {

  M_=M;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setmass(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO> >& TpetraM) {

  RCP< Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Atmp
    = rcp( new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraM) );
  M_= rcp( new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Atmp) );

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setcoords(RCP<MultiVector>& Coords) {

  Coords_=Coords;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setNullSpace(RCP<MultiVector> NullSpace) {

  NullSpace_=NullSpace;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setLevelShifts(std::vector<Scalar> levelshifts) {

  levelshifts_=levelshifts;
  numLevels_=levelshifts_.size();

}

// initialize
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::initialize() {

  TentPfact_     = rcp( new TentativePFactory           );
  Pfact_         = rcp( new SaPFactory                  );
  PgPfact_       = rcp( new PgPFactory                  );
  TransPfact_    = rcp( new TransPFactory               );
  Rfact_         = rcp( new GenericRFactory             );
  Acfact_        = rcp( new RAPFactory                  );
  Acshift_       = rcp( new RAPShiftFactory             );
  Dropfact_      = rcp( new CoalesceDropFactory         );
  Aggfact_       = rcp( new CoupledAggregationFactory   );
  UCaggfact_     = rcp( new UncoupledAggregationFactory );
  CoarseMapfact_ = rcp( new CoarseMapFactory            );
  Manager_       = rcp( new FactoryManager              );
  if(isSymmetric_==true) {
    Manager_   -> SetFactory("P", Pfact_);
    Manager_   -> SetFactory("R", TransPfact_);
  }
  else {
    Manager_   -> SetFactory("P", PgPfact_);
    Manager_   -> SetFactory("R", Rfact_);
    solverType_ = 10;
  }
  Manager_   -> SetFactory("Ptent", TentPfact_);
  Teuchos::ParameterList params;
  params.set("lightweight wrap",true);
  params.set("aggregation: drop scheme","classical");
  Dropfact_  -> SetParameterList(params);
  Manager_   -> SetFactory("Graph", Dropfact_);
  if(Aggregation_=="coupled") {
    Manager_   -> SetFactory("Aggregates", Aggfact_   );
  }
  else {
    Manager_   -> SetFactory("Aggregates", UCaggfact_ );
  }
  Manager_     -> SetFactory("CoarseMap", CoarseMapfact_);

  // choose smoother
  if(Smoother_=="jacobi") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Jacobi");
    precList_.set("relaxation: sweeps", smoother_sweeps_);
    precList_.set("relaxation: damping factor", smoother_damping_);
  }
  else if(Smoother_=="gauss-seidel") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Gauss-Seidel");
    precList_.set("relaxation: sweeps", smoother_sweeps_);
    precList_.set("relaxation: damping factor", smoother_damping_);
  }
  else if(Smoother_=="symmetric gauss-seidel") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Symmetric Gauss-Seidel");
    precList_.set("relaxation: sweeps", smoother_sweeps_);
    precList_.set("relaxation: damping factor", smoother_damping_);
  }
  else if(Smoother_=="chebyshev") {
    precType_ = "CHEBYSHEV";
  }
  else if(Smoother_=="krylov") {
    precType_ = "KRYLOV";
    precList_.set("krylov: iteration type", krylov_type_);
    precList_.set("krylov: number of iterations", krylov_iterations_);
    precList_.set("krylov: residual tolerance",1.0e-8);
    precList_.set("krylov: block size",1);
    precList_.set("krylov: preconditioner type", krylov_preconditioner_);
    precList_.set("relaxation: sweeps",1);
    solverType_=10;
  }
  else if(Smoother_=="ilut") {
    precType_ = "ILUT";
    precList_.set("fact: ilut level-of-fill", ilu_leveloffill_);
    precList_.set("fact: absolute threshold", ilu_abs_thresh_);
    precList_.set("fact: relative threshold", ilu_rel_thresh_);
    precList_.set("fact: drop tolerance",     ilu_drop_tol_);
    precList_.set("fact: relax value",        ilu_relax_val_);
  }
  else if(Smoother_=="riluk") {
    precType_ = "RILUK";
    precList_.set("fact: iluk level-of-fill", ilu_leveloffill_);
    precList_.set("fact: absolute threshold", ilu_abs_thresh_);
    precList_.set("fact: relative threshold", ilu_rel_thresh_);
    precList_.set("fact: drop tolerance",     ilu_drop_tol_);
    precList_.set("fact: relax value",        ilu_relax_val_);
  }
  else if(Smoother_=="schwarz") {
    precType_ = "SCHWARZ";
    precList_.set("schwarz: overlap level", schwarz_overlap_);
    precList_.set("schwarz: combine mode", schwarz_combinemode_);
    precList_.set("schwarz: use reordering", schwarz_usereorder_);
    //    precList_.set("schwarz: filter singletons", true); // Disabled due to issues w/ Ifpack2/Zoltan2 w.r.t. Issue #560 - CMS 8/26/16
    precList_.set("order_method",schwarz_ordermethod_);
    precList_.sublist("schwarz: reordering list").set("order_method",schwarz_ordermethod_);
    precList_.sublist("schwarz: subdomain solver parameters").set("fact: ilut level-of-fill", ilu_leveloffill_);
    precList_.sublist("schwarz: subdomain solver parameters").set("fact: absolute threshold", ilu_abs_thresh_);
    precList_.sublist("schwarz: subdomain solver parameters").set("fact: relative threshold", ilu_rel_thresh_);
    precList_.sublist("schwarz: subdomain solver parameters").set("fact: drop tolerance",     ilu_drop_tol_);
    precList_.sublist("schwarz: subdomain solver parameters").set("fact: relax value",        ilu_relax_val_);
  }
  else if(Smoother_=="superilu") {
    precType_ = "superlu";
    precList_.set("RowPerm", ilu_rowperm_);
    precList_.set("ColPerm", ilu_colperm_);
    precList_.set("DiagPivotThresh", ilu_diagpivotthresh_);
    precList_.set("ILU_DropRule",ilu_drop_rule_);
    precList_.set("ILU_DropTol",ilu_drop_tol_);
    precList_.set("ILU_FillFactor",ilu_leveloffill_);
    precList_.set("ILU_Norm",ilu_normtype_);
    precList_.set("ILU_MILU",ilu_milutype_);
    precList_.set("ILU_FillTol",ilu_fill_tol_);
    precList_.set("ILU_Flag",true);
  }
  else if(Smoother_=="superlu") {
    precType_ = "superlu";
    precList_.set("ColPerm", ilu_colperm_);
    precList_.set("DiagPivotThresh", ilu_diagpivotthresh_);
  }
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  // construct smoother
  smooProto_ = rcp( new Ifpack2Smoother(precType_,precList_) );
  smooFact_  = rcp( new SmootherFactory(smooProto_) );
#if defined(HAVE_MUELU_AMESOS2) and defined(HAVE_AMESOS2_SUPERLU)
  coarsestSmooProto_ = rcp( new DirectSolver("Superlu",coarsestSmooList_) );
#elif defined(HAVE_MUELU_AMESOS2) and defined(HAVE_AMESOS2_KLU2)
  coarsestSmooProto_ = rcp( new DirectSolver("Klu",coarsestSmooList_) );
#elif defined(HAVE_MUELU_AMESOS2) and defined(HAVE_AMESOS2_SUPERLUDIST)
  coarsestSmooProto_ = rcp( new DirectSolver("Superludist",coarsestSmooList_) );
#else
  coarsestSmooProto_ = rcp( new Ifpack2Smoother(precType_,precList_) );
#endif
  coarsestSmooFact_  = rcp( new SmootherFactory(coarsestSmooProto_, Teuchos::null) );

  // For setupSlowRAP and setupFastRAP, the prolongation/restriction matrices
  // are constructed with the stiffness matrix. These matrices are kept for future
  // setup calls; this is achieved by calling Hierarchy->Keep(). It is particularly
  // useful for multiple frequency problems - when the frequency/preconditioner
  // changes, you only compute coarse grids (RAPs) and setup level smoothers when
  // you call Hierarchy->Setup().
  if(K_!=Teuchos::null) {
    Manager_ -> SetFactory("Smoother", Teuchos::null);
    Manager_ -> SetFactory("CoarseSolver", Teuchos::null);
    Hierarchy_ = rcp( new Hierarchy(K_)  );
    if(NullSpace_!=Teuchos::null)
      Hierarchy_ -> GetLevel(0) -> Set("Nullspace", NullSpace_);
    if(isSymmetric_==true) {
      Hierarchy_ -> Keep("P", Pfact_.get());
      Hierarchy_ -> Keep("R", TransPfact_.get());
      Hierarchy_ -> SetImplicitTranspose(true);
    }
    else {
      Hierarchy_ -> Keep("P", PgPfact_.get());
      Hierarchy_ -> Keep("R", Rfact_.get());
    }
    Hierarchy_ -> Keep("Ptent", TentPfact_.get());
    Hierarchy_ -> SetMaxCoarseSize( coarseGridSize_ );
    Hierarchy_ -> Setup(*Manager_, 0, numLevels_);
    GridTransfersExist_=true;
  }
  // Use preconditioning matrix to setup prolongation/restriction operators
  else {
    Manager_ -> SetFactory("Smoother", smooFact_);
    Manager_ -> SetFactory("CoarseSolver", coarsestSmooFact_);
    Hierarchy_ = rcp( new Hierarchy(P_)  );
    if(NullSpace_!=Teuchos::null)
      Hierarchy_ -> GetLevel(0) -> Set("Nullspace", NullSpace_);
    if(isSymmetric_==true)
      Hierarchy_ -> SetImplicitTranspose(true);
    Hierarchy_ -> SetMaxCoarseSize( coarseGridSize_ );
    Hierarchy_ -> Setup(*Manager_, 0, numLevels_);
    GridTransfersExist_=true;
  }

  // Belos Linear Problem and Solver Manager
  BelosList_ = rcp( new Teuchos::ParameterList("GMRES") );
  BelosList_ -> set("Maximum Iterations",iters_ );
  BelosList_ -> set("Convergence Tolerance",tol_ );
  BelosList_ -> set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  BelosList_ -> set("Output Frequency",1);
  BelosList_ -> set("Output Style",Belos::Brief);
  BelosList_ -> set("Num Blocks",restart_size_);
  BelosList_ -> set("Num Recycled Blocks",recycle_size_);
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "ShiftedLaplacian only available with Tpetra and GO=int enabled.");
#endif
}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setupFastRAP() {

  int numLevels = Hierarchy_ -> GetNumLevels();
  Manager_ -> SetFactory("Smoother", smooFact_);
  Manager_ -> SetFactory("CoarseSolver", coarsestSmooFact_);
  Hierarchy_ -> GetLevel(0) -> Set("A", P_);
  Hierarchy_ -> Setup(*Manager_, 0, numLevels);
  setupSolver();

}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setupSlowRAP() {

  int numLevels = Hierarchy_ -> GetNumLevels();
  Acshift_->SetShifts(levelshifts_);
  Manager_ -> SetFactory("Smoother", smooFact_);
  Manager_ -> SetFactory("CoarseSolver", coarsestSmooFact_);
  Manager_ -> SetFactory("A", Acshift_);
  Manager_ -> SetFactory("K", Acshift_);
  Manager_ -> SetFactory("M", Acshift_);
  Hierarchy_ -> GetLevel(0) -> Set("A", P_);
  Hierarchy_ -> GetLevel(0) -> Set("K", K_);
  Hierarchy_ -> GetLevel(0) -> Set("M", M_);
  Hierarchy_ -> Setup(*Manager_, 0, numLevels);
  setupSolver();

}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setupNormalRAP() {

  // Only setup hierarchy again if preconditioning matrix has changed
  if( GridTransfersExist_ == false ) {
    Hierarchy_ = rcp( new Hierarchy(P_)  );
    if(NullSpace_!=Teuchos::null)
      Hierarchy_ -> GetLevel(0) -> Set("Nullspace", NullSpace_);
    if(isSymmetric_==true)
      Hierarchy_ -> SetImplicitTranspose(true);
    Hierarchy_ -> SetMaxCoarseSize( coarseGridSize_ );
    Hierarchy_ -> Setup(*Manager_, 0, numLevels_);
    GridTransfersExist_=true;
  }
  setupSolver();

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setupSolver() {

#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  // Define Preconditioner and Operator
  MueLuOp_ = rcp( new MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO>
                  (Hierarchy_, A_, ncycles_, subiters_, option_, tol_) );
  // Belos Linear Problem
  if(LinearProblem_==Teuchos::null)
    LinearProblem_ = rcp( new LinearProblem );
  LinearProblem_ -> setOperator (  TpetraA_  );
  LinearProblem_ -> setRightPrec(  MueLuOp_  );
  if(SolverManager_==Teuchos::null) {
    std::string solverName;
    SolverFactory_= rcp( new SolverFactory() );
    if(solverType_==1)      { solverName="Block GMRES";      }
    else if(solverType_==2) { solverName="Recycling GMRES";  }
    else                    { solverName="Flexible GMRES";   }
    SolverManager_ = SolverFactory_->create( solverName, BelosList_ );
    SolverManager_ -> setProblem( LinearProblem_ );
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "ShiftedLaplacian only available with Tpetra and GO=int enabled.");
#endif
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::resetLinearProblem()
{
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  LinearProblem_ -> setOperator (  TpetraA_  );
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "ShiftedLaplacian only available with Tpetra and GO=int enabled.");
#endif
}

// Solve phase
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::solve(const RCP<TMV> B, RCP<TMV>& X)
{
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  // Set left and right hand sides for Belos
  LinearProblem_ -> setProblem(X, B);
  // iterative solve
  SolverManager_ -> solve();
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "ShiftedLaplacian only available with Tpetra and GO=int enabled.");
#endif
  return 0;
}

// Solve phase
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::multigrid_apply(const RCP<MultiVector> B,
                                                                               RCP<MultiVector>& X)
{
  // Set left and right hand sides for Belos
  Hierarchy_ -> Iterate(*B, *X, 1, true, 0);
}

// Solve phase
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::multigrid_apply(const RCP<Tpetra::MultiVector<SC,LO,GO,NO> > B,
                                                                               RCP<Tpetra::MultiVector<SC,LO,GO,NO> >& X)
{
  Teuchos::RCP< Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > XpetraX
    = Teuchos::rcp( new Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(X) );
  Teuchos::RCP< Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > XpetraB
    = Teuchos::rcp( new Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(B) );
  // Set left and right hand sides for Belos
  Hierarchy_ -> Iterate(*XpetraB, *XpetraX, 1, true, 0);
}

// Get most recent iteration count
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetIterations()
{
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  int numiters = SolverManager_ -> getNumIters();
  return numiters;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "ShiftedLaplacian only available with Tpetra and GO=int enabled.");
  return -1;
#endif
}

// Get most recent solver tolerance achieved
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>::GetResidual()
{
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  double residual = SolverManager_ -> achievedTol();
  return residual;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "ShiftedLaplacian only available with Tpetra and GO=int enabled.");
  return -1.0;
#endif
}

}

#define MUELU_SHIFTEDLAPLACIAN_SHORT

#endif //if defined(HAVE_MUELU_IFPACK2) and defined(HAVE_MUELU_TPETRA)
#endif // MUELU_SHIFTEDLAPLACIAN_DEF_HPP
