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
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setParameters(Teuchos::RCP< Teuchos::ParameterList > paramList) {

  // Smoother parameters
  overlap_level_        = paramList->get("schwarz: overlap level",2);
  relaxation_sweeps_    = paramList->get("relaxation: sweeps",4);
  relaxation_damping_   = paramList->get("relaxation: damping factor",(SC)1.0);
  krylov_type_          = paramList->get("krylov: iteration type",1);
  krylov_iterations_    = paramList->get("krylov: number of iterations",5);
  ilut_leveloffill_     = paramList->get("fact: ilut level-of-fill",5.0);
  iluk_leveloffill_     = paramList->get("fact: iluk level-of-fill",1.0);
  silu_leveloffill_     = paramList->get("ILU_FillFactor",10.0);
  ilu_abs_thresh_       = paramList->get("fact: absolute threshold",0.0);
  ilu_rel_thresh_       = paramList->get("fact: relative threshold",1.0);
  ilu_diagpivotthresh_  = paramList->get("DiagPivotThresh",0.1);
  ilu_drop_tol_         = paramList->get("fact: drop tolerance",0.0001);
  ilu_fill_tol_         = paramList->get("ILU_FillTol",0.01);
  ilu_relax_val_        = paramList->get("fact: relax value",1.0);
  ilu_rowperm_          = paramList->get("RowPerm","LargeDiag");
  ilu_colperm_          = paramList->get("ColPerm","COLAMD");
  ilu_drop_rule_        = paramList->get("ILU_DropRule","DROP_BASIC");
  ilu_normtype_         = paramList->get("ILU_Norm","INF_NORM");
  ilu_milutype_         = paramList->get("ILU_MILU","SILU");
  schwarz_usereorder_   = paramList->get("schwarz: use reordering",true);
  schwarz_combinemode_  = paramList->get("schwarz: combine mode",Tpetra::ZERO);
  schwarz_ordermethod_  = paramList->get("order_method","rcm");

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setLaplacian(RCP<Matrix>& L) {

  L_=L;
  LaplaceOperatorSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setProblemMatrix(RCP<Matrix>& A) {

  A_=A;
  if(A_!=Teuchos::null)
    TpetraA_ = Utils::Op2NonConstTpetraCrs(A_);  
  ProblemMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setProblemMatrix(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >& TpetraA) {

  TpetraA_=TpetraA;
  ProblemMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setPreconditioningMatrix(RCP<Matrix>& P) {

  P_=P;
  PreconditioningMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setstiff(RCP<Matrix>& K) {

  K_=K;
  StiffMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setmass(RCP<Matrix>& M) {

  M_=M;
  MassMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setdamp(RCP<Matrix>& C) {

  C_=C;
  DampMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setcoords(RCP<MultiVector>& Coords) {

  Coords_=Coords;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setNullSpace(RCP<MultiVector> NullSpace) {

  NullSpace_=NullSpace;

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

  if(stype==1)       { Smoother_="jacobi";                 }
  else if(stype==2)  { Smoother_="gauss-seidel";           }
  else if(stype==3)  { Smoother_="symmetric gauss-seidel"; }
  else if(stype==4)  { Smoother_="chebyshev";              }
  else if(stype==5)  { Smoother_="krylov";                 }
  else if(stype==6)  { Smoother_="ilut";                   }
  else if(stype==7)  { Smoother_="riluk";                  }
  else if(stype==8)  { Smoother_="schwarz";                }
  else if(stype==9)  { Smoother_="superilu";               }
  else if(stype==10) { Smoother_="superlu";                }
  else               { Smoother_="schwarz";                }

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setSolver(int stype) {

  if(stype==1) { FGMRESoption_=true;   }
  else         { FGMRESoption_=false;  }
  if(FGMRESoption_==true) {
    solverType_=1;
  }

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setSolverType(int stype) {
  
  solverType_=stype;

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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setNumLevels(int numlevels) {

  numLevels_=numlevels;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setSymmetric(bool isSymmetric) {

  isSymmetric_=isSymmetric;
    
}

// initialize
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::initialize() {

  TentPfact_ = rcp( new TentativePFactory           );
  Pfact_     = rcp( new SaPFactory                  );
  PgPfact_   = rcp( new PgPFactory                  );
  TransPfact_= rcp( new TransPFactory               );
  Rfact_     = rcp( new GenericRFactory             );
  Acfact_    = rcp( new RAPFactory                  );
  Acshift_   = rcp( new RAPShiftFactory             );
  Dropfact_  = rcp( new CoalesceDropFactory         );
  Aggfact_   = rcp( new CoupledAggregationFactory   );
  UCaggfact_ = rcp( new UncoupledAggregationFactory );
  Manager_   = rcp( new FactoryManager              );
  if(isSymmetric_==true) {
    Manager_   -> SetFactory("P", Pfact_);
    Manager_   -> SetFactory("R", TransPfact_);
  }
  else {
    Manager_   -> SetFactory("P", PgPfact_);
    Manager_   -> SetFactory("R", Rfact_);
    solverType_ = 1;
  }
  Manager_   -> SetFactory("Ptent", TentPfact_);
  Teuchos::ParameterList params;
  params.set("lightweight wrap",true);
  params.set("algorithm","original");
  Dropfact_  -> SetParameterList(params);
  Manager_   -> SetFactory("Graph", Dropfact_);
  Manager_   -> SetFactory("Smoother", Teuchos::null);
  Manager_   -> SetFactory("CoarseSolver", Teuchos::null);
  if(Aggregation_=="coupled") {
    Manager_   -> SetFactory("Aggregates", Aggfact_   );
  }
  else {
    Manager_   -> SetFactory("Aggregates", UCaggfact_ );
  }

  // choose smoother
  if(Smoother_=="jacobi") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Jacobi");
    precList_.set("relaxation: sweeps", relaxation_sweeps_);
    precList_.set("relaxation: damping factor", relaxation_damping_);
  }
  else if(Smoother_=="gauss-seidel") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Gauss-Seidel");
    precList_.set("relaxation: sweeps", relaxation_sweeps_);
    precList_.set("relaxation: damping factor", relaxation_damping_);
  }
  else if(Smoother_=="symmetric gauss-seidel") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Symmetric Gauss-Seidel");
    precList_.set("relaxation: sweeps", relaxation_sweeps_);
    precList_.set("relaxation: damping factor", relaxation_damping_);
  }
  else if(Smoother_=="chebyshev") {
    precType_ = "CHEBYSHEV";
  }
  else if(Smoother_=="krylov") {
    precType_ = "KRYLOV";
    precList_.set("krylov: iteration type",krylov_type_);
    precList_.set("krylov: number of iterations", krylov_iterations_);
    precList_.set("krylov: residual tolerance",1.0e-8);
    precList_.set("krylov: block size",1);
    precList_.set("krylov: preconditioner type",1);
    precList_.set("relaxation: sweeps",1);
    FGMRESoption_=true;
  }
  else if(Smoother_=="ilut") {
    precType_ = "ILUT";
    precList_.set("fact: ilut level-of-fill", ilut_leveloffill_);
    precList_.set("fact: absolute threshold", ilu_abs_thresh_);
    precList_.set("fact: relative threshold", ilu_rel_thresh_);
    precList_.set("fact: drop tolerance",     ilu_drop_tol_);
    precList_.set("fact: relax value",        ilu_relax_val_);
  }
  else if(Smoother_=="riluk") {
    precType_ = "RILUK";
    precList_.set("fact: iluk level-of-fill", iluk_leveloffill_);
    precList_.set("fact: absolute threshold", ilu_abs_thresh_);
    precList_.set("fact: relative threshold", ilu_rel_thresh_);
    precList_.set("fact: drop tolerance",     ilu_drop_tol_);
    precList_.set("fact: relax value",        ilu_relax_val_);
  }
  else if(Smoother_=="schwarz") {
    precType_ = "SCHWARZ";
    precList_.set("fact: ilut level-of-fill", ilut_leveloffill_);
    precList_.set("fact: drop tolerance", ilu_drop_tol_);
    precList_.set("schwarz: compute condest", false);
    precList_.set("schwarz: combine mode", schwarz_combinemode_);
    precList_.set("schwarz: use reordering", schwarz_usereorder_);
    precList_.set("schwarz: filter singletons", true);
    precList_.set("order_method",schwarz_ordermethod_);
    precList_.sublist("schwarz: reordering list").set("order_method",schwarz_ordermethod_);
  }
  else if(Smoother_=="superilu") {
    precType_ = "superilu";
    precList_.set("RowPerm", ilu_rowperm_);
    precList_.set("ColPerm", ilu_colperm_);
    precList_.set("DiagPivotThresh", ilu_diagpivotthresh_);
    precList_.set("ILU_DropRule",ilu_drop_rule_);
    precList_.set("ILU_DropTol",ilu_drop_tol_);
    precList_.set("ILU_FillFactor",silu_leveloffill_);
    precList_.set("ILU_Norm",ilu_normtype_);
    precList_.set("ILU_MILU",ilu_milutype_);
    precList_.set("ILU_FillTol",ilu_fill_tol_);
  }
  else if(Smoother_=="superlu") {
    precType_ = "superlu";
    precList_.set("ColPerm", ilu_colperm_);
    precList_.set("DiagPivotThresh", ilu_diagpivotthresh_);
  }
  // construct smoother
  if(Smoother_=="schwarz") {
    smooProto_ = rcp( new Ifpack2Smoother(precType_,precList_) );
  }
  else {
    smooProto_ = rcp( new SchwarzSmoother(precType_,precList_,overlap_level_) );
  }
  smooFact_  = rcp( new SmootherFactory(smooProto_) );
  coarsestSmooProto_ = rcp( new DirectSolver("Superlu",coarsestSmooList_) );
  coarsestSmooFact_  = rcp( new SmootherFactory(coarsestSmooProto_, Teuchos::null) );
  
  // Use stiffness matrix to setup prolongation/restriction operators
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

  // Belos Linear Problem and Solver Manager
  BelosList_ = rcp( new Teuchos::ParameterList("GMRES") );
  BelosList_ -> set("Maximum Iterations",iters_ );
  BelosList_ -> set("Convergence Tolerance",tol_ );
  BelosList_ -> set("Flexible Gmres", FGMRESoption_ );
  BelosList_ -> set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  BelosList_ -> set("Output Frequency",1);
  BelosList_ -> set("Output Style",Belos::Brief);

}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setupFastRAP() {

  int numLevels = Hierarchy_ -> GetNumLevels();

  Manager_ -> SetFactory("Smoother", smooFact_);
  Manager_ -> SetFactory("CoarseSolver", coarsestSmooFact_);
  Hierarchy_ -> GetLevel(0) -> Set("A", P_);
  Hierarchy_ -> Setup(*Manager_, 0, numLevels);

  // Define Preconditioner and Operator
  MueLuOp_ = rcp( new MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO>(Hierarchy_, A_, ncycles_, subiters_, option_, tol_) );
  // Belos Linear Problem
  BelosLinearProblem_ = rcp( new BelosLinearProblem );
  BelosLinearProblem_ -> setOperator (  TpetraA_  );
  BelosLinearProblem_ -> setRightPrec(  MueLuOp_  );
  if(solverType_==0) {
    BelosSolverManager_ = rcp( new BelosCG(BelosLinearProblem_, BelosList_) );
  }
  else {
    BelosSolverManager_ = rcp( new BelosGMRES(BelosLinearProblem_, BelosList_) );
  }

}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setupSlowRAP() {

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

  // Define Preconditioner and Operator
  MueLuOp_ = rcp( new MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO>(Hierarchy_, A_, ncycles_, subiters_, option_, tol_) );
  // Belos Linear Problem
  BelosLinearProblem_ = rcp( new BelosLinearProblem );
  BelosLinearProblem_ -> setOperator (  TpetraA_  );
  BelosLinearProblem_ -> setRightPrec(  MueLuOp_  );
  if(solverType_==0) {
    BelosSolverManager_ = rcp( new BelosCG(BelosLinearProblem_, BelosList_) );
  }
  else {
    BelosSolverManager_ = rcp( new BelosGMRES(BelosLinearProblem_, BelosList_) );
  }

}

// setup coarse grids for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setupNormalRAP() {

  TentPfact_ = rcp( new TentativePFactory           );
  Pfact_     = rcp( new SaPFactory                  );
  PgPfact_   = rcp( new PgPFactory                  );
  TransPfact_= rcp( new TransPFactory               );
  Rfact_     = rcp( new GenericRFactory             );
  Acfact_    = rcp( new RAPFactory                  );
  Acshift_   = rcp( new RAPShiftFactory             );
  Dropfact_  = rcp( new CoalesceDropFactory         );
  Aggfact_   = rcp( new CoupledAggregationFactory   );
  UCaggfact_ = rcp( new UncoupledAggregationFactory );
  Manager_   = rcp( new FactoryManager              );
  if(isSymmetric_==true) {
    Manager_   -> SetFactory("P", Pfact_);
    Manager_   -> SetFactory("R", TransPfact_);
  }
  else {
    Manager_   -> SetFactory("P", PgPfact_);
    Manager_   -> SetFactory("R", Rfact_);
    solverType_ = 1;
  }
  Manager_   -> SetFactory("Ptent", TentPfact_);
  Teuchos::ParameterList params;
  params.set("lightweight wrap",true);
  params.set("algorithm","original");
  Dropfact_  -> SetParameterList(params);
  Manager_   -> SetFactory("Graph", Dropfact_);
  if(Aggregation_=="coupled") {
    Manager_   -> SetFactory("Aggregates", Aggfact_   );
  }
  else {
    Manager_   -> SetFactory("Aggregates", UCaggfact_ );
  }

  // choose smoother
  if(Smoother_=="jacobi") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Jacobi");
    precList_.set("relaxation: sweeps", relaxation_sweeps_);
    precList_.set("relaxation: damping factor", relaxation_damping_);
  }
  else if(Smoother_=="gauss-seidel") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Gauss-Seidel");
    precList_.set("relaxation: sweeps", relaxation_sweeps_);
    precList_.set("relaxation: damping factor", relaxation_damping_);
  }
  else if(Smoother_=="symmetric gauss-seidel") {
    precType_ = "RELAXATION";
    precList_.set("relaxation: type", "Symmetric Gauss-Seidel");
    precList_.set("relaxation: sweeps", relaxation_sweeps_);
    precList_.set("relaxation: damping factor", relaxation_damping_);
  }
  else if(Smoother_=="chebyshev") {
    precType_ = "CHEBYSHEV";
  }
  else if(Smoother_=="krylov") {
    precType_ = "KRYLOV";
    precList_.set("krylov: iteration type",krylov_type_);
    precList_.set("krylov: number of iterations", krylov_iterations_);
    precList_.set("krylov: residual tolerance",1.0e-8);
    precList_.set("krylov: block size",1);
    precList_.set("krylov: preconditioner type",1);
    precList_.set("relaxation: sweeps",1);
    FGMRESoption_=true;
  }
  else if(Smoother_=="ilut") {
    precType_ = "ILUT";
    precList_.set("fact: ilut level-of-fill", ilut_leveloffill_);
    precList_.set("fact: absolute threshold", ilu_abs_thresh_);
    precList_.set("fact: relative threshold", ilu_rel_thresh_);
    precList_.set("fact: drop tolerance",     ilu_drop_tol_);
    precList_.set("fact: relax value",        ilu_relax_val_);
  }
  else if(Smoother_=="riluk") {
    precType_ = "RILUK";
    precList_.set("fact: iluk level-of-fill", iluk_leveloffill_);
    precList_.set("fact: absolute threshold", ilu_abs_thresh_);
    precList_.set("fact: relative threshold", ilu_rel_thresh_);
    precList_.set("fact: drop tolerance",     ilu_drop_tol_);
    precList_.set("fact: relax value",        ilu_relax_val_);
  }
  else if(Smoother_=="schwarz") {
    precType_ = "SCHWARZ";
    precList_.set("fact: ilut level-of-fill", ilut_leveloffill_);
    precList_.set("fact: drop tolerance", ilu_drop_tol_);
    precList_.set("schwarz: compute condest", false);
    precList_.set("schwarz: combine mode", schwarz_combinemode_);
    precList_.set("schwarz: use reordering", schwarz_usereorder_);
    precList_.set("schwarz: filter singletons", true);
    precList_.set("order_method",schwarz_ordermethod_);
    precList_.sublist("schwarz: reordering list").set("order_method",schwarz_ordermethod_);
  }
  else if(Smoother_=="superilu") {
    precType_ = "superilu";
    precList_.set("RowPerm", ilu_rowperm_);
    precList_.set("ColPerm", ilu_colperm_);
    precList_.set("DiagPivotThresh", ilu_diagpivotthresh_);
    precList_.set("ILU_DropRule",ilu_drop_rule_);
    precList_.set("ILU_DropTol",ilu_drop_tol_);
    precList_.set("ILU_FillFactor",silu_leveloffill_);
    precList_.set("ILU_Norm",ilu_normtype_);
    precList_.set("ILU_MILU",ilu_milutype_);
    precList_.set("ILU_FillTol",ilu_fill_tol_);
  }
  else if(Smoother_=="superlu") {
    precType_ = "superlu";
    precList_.set("ColPerm", ilu_colperm_);
    precList_.set("DiagPivotThresh", ilu_diagpivotthresh_);
  }
  // construct smoother
  if(Smoother_=="schwarz") {
    smooProto_ = rcp( new Ifpack2Smoother(precType_,precList_) );
  }
  else {
    smooProto_ = rcp( new SchwarzSmoother(precType_,precList_,overlap_level_) );
  }
  smooFact_  = rcp( new SmootherFactory(smooProto_) );
  coarsestSmooProto_ = rcp( new DirectSolver("Superlu",coarsestSmooList_) );
  coarsestSmooFact_  = rcp( new SmootherFactory(coarsestSmooProto_, Teuchos::null) );
  Manager_ -> SetFactory("Smoother", smooFact_);
  Manager_ -> SetFactory("CoarseSolver", coarsestSmooFact_);
  
  // Normal setup
  Hierarchy_ = rcp( new Hierarchy(P_)  );
  if(NullSpace_!=Teuchos::null)
    Hierarchy_ -> GetLevel(0) -> Set("Nullspace", NullSpace_);
  if(isSymmetric_==true)
    Hierarchy_ -> SetImplicitTranspose(true);
  Hierarchy_ -> SetMaxCoarseSize( coarseGridSize_ );
  Hierarchy_ -> Setup(*Manager_, 0, numLevels_);
  GridTransfersExist_=true;

  // Define Operator and Preconditioner
  MueLuOp_ = rcp( new MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO>(Hierarchy_, A_, ncycles_, subiters_, option_, tol_) );

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
  if(solverType_==0) {
    BelosSolverManager_ = rcp( new BelosCG(BelosLinearProblem_, BelosList_) );
  }
  else {
    BelosSolverManager_ = rcp( new BelosGMRES(BelosLinearProblem_, BelosList_) );
  }

}
 
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::resetLinearProblem()
{
  
  BelosLinearProblem_ -> setOperator (  TpetraA_  );

}
 
// Solve phase
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
int ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::solve(const RCP<TMV> B, RCP<TMV>& X)
{

  // Set left and right hand sides for Belos
  BelosLinearProblem_ -> setProblem(X, B);
  // iterative solve
  //Belos::ReturnType convergenceStatus = BelosSolverManager_ -> solve();
  BelosSolverManager_ -> solve();
  /*if(convergenceStatus == Belos::Converged) {
    return 0;
  }
  else {
    return 1;
    }*/
  return 0;

}

// Solve phase
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::multigrid_apply(const RCP<MultiVector> B, RCP<MultiVector>& X)
{

  // Set left and right hand sides for Belos
  Hierarchy_ -> Iterate(*B, 1, *X, true, 0);

}

// Get most recent iteration count
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
int ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::GetIterations()
{

  int numiters = BelosSolverManager_ -> getNumIters();
  return numiters;

}

}

#define MUELU_SHIFTEDLAPLACIAN_SHORT
#endif // MUELU_SHIFTEDLAPLACIAN_DEF_HPP
