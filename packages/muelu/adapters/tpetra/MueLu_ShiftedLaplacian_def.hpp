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

namespace MueLu {

// Destructor
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::~ShiftedLaplacian() {}

// Input
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setParameters(const Teuchos::ParameterList List) {

  Problem_        =  List.get<std::string>      (  "muelu: problem type"          );
  Smoother_       =  List.get<std::string>      (  "muelu: smoother"              );
  Aggregation_    =  List.get<std::string>      (  "muelu: aggregation"           );
  Nullspace_      =  List.get<std::string>      (  "muelu: nullspace"             );
  UseLaplacian_   =  List.get<bool>             (  "muelu: use laplacian"         );
  VariableShift_  =  List.get<bool>             (  "muelu: variable shift"        );
  numPDEs_        =  List.get<int>              (  "muelu: dofs per node"         );
  numLevels_      =  List.get<int>              (  "muelu: number of levels"      );
  coarseGridSize_ =  List.get<int>              (  "muelu: coarse grid size"      );
  iters_          =  List.get<int>              (  "muelu: number of iterations"  );
  blksize_        =  List.get<int>              (  "muelu: block size"            );
  FGMRESoption_   =  List.get<bool>             (  "muelu: fgmres on/off"         );
  tol_            =  List.get<double>           (  "muelu: residual tolerance"    );  
  omega_          =  List.get<double>           (  "muelu: omega"                 );

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setLaplacian(RCP<Matrix>& L) {

  L_=L;
  L = null;
  LaplaceOperatorSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setProblemMatrix(RCP<Matrix>& A) {

  A_=A;
  A = null;
  ProblemMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setPreconditioningMatrix(RCP<Matrix>& P) {

  P_=P;
  P = null;
  PreconditioningMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setstiff(RCP<Matrix>& K) {

  K_=K;
  K = null;
  StiffMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setmass(RCP<Matrix>& M) {

  M_=M;
  M = null;
  MassMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setdamp(RCP<Matrix>& C) {

  C_=C;
  C = null;
  DampMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setcoords(RCP<MultiVector>& Coords) {

  Coords_=Coords;
  Coords = null;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setProblemShifts(Scalar ashift1, Scalar ashift2) {

  ashift1_=ashift1;
  ashift2_=ashift2;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setPreconditioningShifts(Scalar pshift1, Scalar pshift2) {

  pshift1_=pshift1;
  pshift2_=pshift2;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setLevelShifts(std::vector<Scalar> levelshifts) {

  levelshifts_=levelshifts;
  numLevels_=levelshifts_.size();
  LevelShiftsSet_=true;
  VariableShift_=true;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setAggregation(int stype) {

  if(stype==1) { Aggregation_="coupled";   }
  else         { Aggregation_="uncoupled"; }

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setSmoother(int stype) {

  if(stype==1)      { Smoother_="gmres";      }
  else if(stype==2) { Smoother_="ilu";        }
  else if(stype==3) { Smoother_="schwarz";    }
  else if(stype==4) { Smoother_="relaxation"; }
  else              { Smoother_="gmres";      }

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setSolver(int stype) {

  if(stype==1) { FGMRESoption_="true";   }
  else         { FGMRESoption_="false";  }

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setSweeps(int nsweeps) {
  
  nsweeps_=nsweeps;

}
  
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setCycles(int ncycles) {
  
  ncycles_=ncycles;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setIterations(int iters) {

  iters_=iters;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setTolerance(double tol) {

  tol_=tol;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setCoarseGridSize(int coarsegridsize) {

  coarseGridSize_=coarsegridsize;

}

// initialize
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::initialize() {

}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setupFastRAP() {

  TentPfact_ = rcp( new TentativePFactory           );
  Pfact_     = rcp( new SaPFactory                  );
  Rfact_     = rcp( new TransPFactory               );
  Acfact_    = rcp( new RAPFactory                  );
  Acshift_   = rcp( new RAPShiftFactory             );
  Aggfact_   = rcp( new CoupledAggregationFactory   );
  UCaggfact_ = rcp( new UncoupledAggregationFactory );
  Manager_   = rcp( new FactoryManager              );
  Manager_   -> SetFactory("P", Pfact_);
  Manager_   -> SetFactory("R", Rfact_);
  Manager_   -> SetFactory("Ptent", TentPfact_);
  Manager_   -> SetFactory("Smoother", Teuchos::null);
  Manager_   -> SetFactory("CoarseSolver", Teuchos::null);
  if(Aggregation_=="coupled") {
    Manager_   -> SetFactory("Aggregates", Aggfact_   );
  }
  else {
    Manager_   -> SetFactory("Aggregates", UCaggfact_ );
  }

  // choose smoother
  if(Smoother_=="gmres") {
    // Krylov smoother
    ifpack2Type_ = "KRYLOV";
    ifpack2List_.set("krylov: iteration type",1);
    ifpack2List_.set("krylov: number of iterations", nsweeps_);
    ifpack2List_.set("krylov: residual tolerance",1e-6);
    ifpack2List_.set("krylov: block size",1);
    ifpack2List_.set("krylov: zero starting solution",true);
    ifpack2List_.set("krylov: preconditioner type",0);
    // must use FGMRES for GMRES smoothing
    FGMRESoption_=true;
  }
  else if(Smoother_=="schwarz") {
    // Additive Schwarz smoother
    ifpack2Type_ = "SCHWARZ";
    ifpack2List_.set("fact: ilut level-of-fill", (double)5.0);
    ifpack2List_.set("fact: drop tolerance", (double) 0.01);
    ifpack2List_.set("schwarz: compute condest", false);
    ifpack2List_.set("schwarz: combine mode", "Add");
    ifpack2List_.set("schwarz: use reordering", true);
    ifpack2List_.set("schwarz: filter singletons", false);
    ifpack2List_.set("schwarz: overlap level", 0);
    ifpack2List_.set("order_method","rcm");
    ifpack2List_.sublist("schwarz: reordering list").set("order_method","rcm");
  }
  else if(Smoother_=="ilu") {
    // ILU smoother
    ifpack2Type_ = "RILUK";
    ifpack2List_.set("fact: iluk level-of-fill", (double)1.0);
    ifpack2List_.set("fact: absolute threshold", (double)1.0);
    ifpack2List_.set("fact: relative threshold", (double)1.0);
    ifpack2List_.set("fact: relax value", (double)0.0);
  }
  else if(Smoother_=="relaxation") {
    // Jacobi smoother
    ifpack2Type_ = "RELAXATION";
    ifpack2List_.set("relaxation: type", "Jacobi");
    ifpack2List_.set("relaxation: sweeps", nsweeps_);
    ifpack2List_.set("relaxation: damping factor", (SC) 0.5);
    ifpack2List_.set("relaxation: zero starting solution", true);
  }
  smooProto_ = rcp( new Ifpack2Smoother(ifpack2Type_,ifpack2List_) );
  smooFact_  = rcp( new SmootherFactory(smooProto_) );
  coarsestSmooProto_ = rcp( new DirectSolver("Superlu",coarsestSmooList_) );
  coarsestSmooFact_  = rcp( new SmootherFactory(coarsestSmooProto_, Teuchos::null) );
  
  // Use stiffness matrix to setup prolongation/restriction operators
  Hierarchy_ = rcp( new Hierarchy(K_)  );
  Hierarchy_ -> Keep("P", Pfact_.get());
  Hierarchy_ -> Keep("R", Rfact_.get());
  Hierarchy_ -> Keep("Ptent", TentPfact_.get());
  Hierarchy_ -> SetImplicitTranspose(true);
  Hierarchy_ -> SetMaxCoarseSize( coarseGridSize_ );
  Hierarchy_ -> Setup(*Manager_, 0, numLevels_);
  GridTransfersExist_=true;

  int numLevels = Hierarchy_ -> GetNumLevels();

  Manager_ -> SetFactory("Smoother", smooFact_);
  Manager_ -> SetFactory("CoarseSolver", coarsestSmooFact_);
  Hierarchy_ -> GetLevel(0) -> Set("A", P_);
  Hierarchy_ -> Setup(*Manager_, 0, numLevels);

  // Define Operator and Preconditioner
  MueLuOp_ = rcp( new MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO>(Hierarchy_, A_, ncycles_, subiters_, option_, tol_) );
  TpetraA_ = Utils::Op2NonConstTpetraCrs(A_);

  // Belos Linear Problem and Solver Manager
  BelosList_ = rcp( new Teuchos::ParameterList("GMRES") );
  BelosList_ -> set("Maximum Iterations",iters_ );
  BelosList_ -> set("Convergence Tolerance",tol_ );
  BelosList_ -> set("Flexible Gmres", FGMRESoption_ );
  BelosList_ -> set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  BelosList_ -> set("Output Frequency",1);
  BelosList_ -> set("Output Style",Belos::Brief);
  // Belos Linear Problem and Solver Manager
  BelosLinearProblem_ = rcp( new BelosLinearProblem );
  BelosLinearProblem_ -> setOperator (  TpetraA_  );
  BelosLinearProblem_ -> setRightPrec(  MueLuOp_  );
  BelosSolverManager_ = rcp( new BelosGMRES(BelosLinearProblem_, BelosList_) );

}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setupSlowRAP() {

  TentPfact_ = rcp( new TentativePFactory           );
  Pfact_     = rcp( new SaPFactory                  );
  Rfact_     = rcp( new TransPFactory               );
  Acfact_    = rcp( new RAPFactory                  );
  Acshift_   = rcp( new RAPShiftFactory             );
  Aggfact_   = rcp( new CoupledAggregationFactory   );
  UCaggfact_ = rcp( new UncoupledAggregationFactory );
  Manager_   = rcp( new FactoryManager );
  Manager_   -> SetFactory("P", Pfact_);
  Manager_   -> SetFactory("R", Rfact_);
  Manager_   -> SetFactory("Ptent", TentPfact_);
  Manager_   -> SetFactory("Smoother", Teuchos::null);
  Manager_   -> SetFactory("CoarseSolver", Teuchos::null);
  if(Aggregation_=="coupled") {
    Manager_   -> SetFactory("Aggregates", Aggfact_   );
  }
  else {
    Manager_   -> SetFactory("Aggregates", UCaggfact_ );
  }

  // choose smoother
  if(Smoother_=="gmres") {
    // Krylov smoother
    ifpack2Type_ = "KRYLOV";
    ifpack2List_.set("krylov: iteration type",1);
    ifpack2List_.set("krylov: number of iterations", nsweeps_);
    ifpack2List_.set("krylov: residual tolerance",1e-6);
    ifpack2List_.set("krylov: block size",1);
    ifpack2List_.set("krylov: zero starting solution",true);
    ifpack2List_.set("krylov: preconditioner type",0);
    // must use FGMRES for GMRES smoothing
    FGMRESoption_=true;
  }
  else if(Smoother_=="schwarz") {
    // Additive Schwarz smoother
    ifpack2Type_ = "SCHWARZ";
    ifpack2List_.set("fact: ilut level-of-fill", (double)5.0);
    ifpack2List_.set("fact: drop tolerance", (double) 0.01);
    ifpack2List_.set("schwarz: compute condest", false);
    ifpack2List_.set("schwarz: combine mode", "Add");
    ifpack2List_.set("schwarz: use reordering", true);
    ifpack2List_.set("schwarz: filter singletons", false);
    ifpack2List_.set("schwarz: overlap level", 0);
    ifpack2List_.set("order_method","rcm");
    ifpack2List_.sublist("schwarz: reordering list").set("order_method","rcm");
  }
  else if(Smoother_=="ilu") {
    // ILU smoother
    ifpack2Type_ = "RILUK";
    ifpack2List_.set("fact: iluk level-of-fill", (double)1.0);
    ifpack2List_.set("fact: absolute threshold", (double)1.0);
    ifpack2List_.set("fact: relative threshold", (double)1.0);
    ifpack2List_.set("fact: relax value", (double)0.0);
  }
  else if(Smoother_=="relaxation") {
    // Jacobi smoother
    ifpack2Type_ = "RELAXATION";
    ifpack2List_.set("relaxation: type", "Jacobi");
    ifpack2List_.set("relaxation: sweeps", nsweeps_);
    ifpack2List_.set("relaxation: damping factor", (SC) 0.5);
    ifpack2List_.set("relaxation: zero starting solution", true);
  }
  smooProto_ = rcp( new Ifpack2Smoother(ifpack2Type_,ifpack2List_) );
  smooFact_  = rcp( new SmootherFactory(smooProto_) );
  coarsestSmooProto_ = rcp( new DirectSolver("Superlu",coarsestSmooList_) );
  coarsestSmooFact_  = rcp( new SmootherFactory(coarsestSmooProto_, Teuchos::null) );
  
  // Use stiffness matrix to setup prolongation/restriction operators
  Hierarchy_ = rcp( new Hierarchy(K_)  );
  Hierarchy_ -> Keep("P", Pfact_.get());
  Hierarchy_ -> Keep("R", Rfact_.get());
  Hierarchy_ -> Keep("Ptent", TentPfact_.get());
  Hierarchy_ -> SetImplicitTranspose(true);
  Hierarchy_ -> SetMaxCoarseSize( coarseGridSize_ );
  Hierarchy_ -> Setup(*Manager_, 0, numLevels_);
  GridTransfersExist_=true;

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
  
  // Define Operator and Preconditioner
  MueLuOp_ = rcp( new MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO>(Hierarchy_, A_, ncycles_, subiters_, option_, tol_) );
  TpetraA_ = Utils::Op2NonConstTpetraCrs(A_);

  // Belos Linear Problem and Solver Manager
  BelosList_ = rcp( new Teuchos::ParameterList("GMRES") );
  BelosList_ -> set("Maximum Iterations",iters_ );
  BelosList_ -> set("Convergence Tolerance",tol_ );
  BelosList_ -> set("Flexible Gmres", FGMRESoption_ );
  BelosList_ -> set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  BelosList_ -> set("Output Frequency",1);
  BelosList_ -> set("Output Style",Belos::Brief);
  BelosLinearProblem_ = rcp( new BelosLinearProblem );
  BelosLinearProblem_ -> setOperator (  TpetraA_  );
  BelosLinearProblem_ -> setRightPrec(  MueLuOp_  );
  BelosSolverManager_ = rcp( new BelosGMRES(BelosLinearProblem_, BelosList_) );

}
  
// Solve phase
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::solve(const RCP<TMV> B, RCP<TMV>& X)
{

  // Set left and right hand sides for Belos
  BelosLinearProblem_ -> setProblem(X, B);
  // iterative solve
  BelosSolverManager_ -> solve();

}

}

#define MUELU_SHIFTEDLAPLACIAN_SHORT
#endif // MUELU_SHIFTEDLAPLACIAN_DEF_HPP
