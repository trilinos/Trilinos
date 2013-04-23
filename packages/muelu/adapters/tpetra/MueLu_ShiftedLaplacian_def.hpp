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
  alpha_          =  List.get<double>           (  "muelu: alpha"                 );
  beta_           =  List.get<double>           (  "muelu: beta"                  );
  omega_          =  List.get<double>           (  "muelu: omega"                 );
  rshift_         =  List.get<double>           (  "muelu: real shift"            );
  ishift_         =  List.get<double>           (  "muelu: imag shift"            );

}

// Preprocessing
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillmatrices() {

  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setLaplacian(RCP<Matrix> L) {

  L_=L;
  LaplaceOperatorSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setHelmholtz(RCP<Matrix> A) {

  A_=A;
  HelmholtzOperatorSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setShiftedLaplacian(RCP<Matrix> S) {

  S_=S;
  ShiftedLaplacianSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setstiff(RCP<Matrix> K) {

  K_=K;
  StiffMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setmass(RCP<Matrix> M) {

  M_=M;
  MassMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setdamp(RCP<Matrix> C) {

  C_=C;
  DampMatrixSet_=true;
  GridTransfersExist_=false;

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setcoords(RCP<MultiVector> Coords) {

  Coords_=Coords;

}

// setup for new frequency
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void ShiftedLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setup(const double omega) {

  if(GridTransfersExist_==false) {

    TentPfact_ = rcp( new TentativePFactory         );
    Pfact_     = rcp( new SaPFactory                );
    Rfact_     = rcp( new TransPFactory             );
    Acfact_    = rcp( new RAPFactory                );
    Acshift_   = rcp( new RAPShiftFactory           );
    Aggfact_   = rcp( new CoupledAggregationFactory );

    // choose smoother
    if(Smoother_=="gmres") {
      // Krylov smoother with Schwarz inner preconditioning
      ifpack2Type_ = "KRYLOV";
      ifpack2List_.set("krylov: iteration type",1);
      ifpack2List_.set("krylov: number of iterations",5);
      ifpack2List_.set("krylov: residual tolerance",1e-6);
      ifpack2List_.set("krylov: block size",1);
      ifpack2List_.set("krylov: zero starting solution",true);
      ifpack2List_.set("krylov: preconditioner type",3);
    }
    else if(Smoother_=="schwarz") {
      // Additive Schwarz smoother
      ifpack2Type_ = "SCHWARZ";
      ifpack2List_.set("schwarz: compute condest", false);
      ifpack2List_.set("schwarz: combine mode", "Add");
      ifpack2List_.set("schwarz: reordering type", "none");
      ifpack2List_.set("schwarz: filter singletons", false);
      ifpack2List_.set("schwarz: overlap level", 0);
    }
    else if(Smoother_=="ilut") {
      // ILUT smoother
      ifpack2Type_ = "ILUT";
      ifpack2List_.set("fact: ilut level-of-fill", (double)1.0);
      ifpack2List_.set("fact: absolute threshold", (double)0.0);
      ifpack2List_.set("fact: relative threshold", (double)1.0);
      ifpack2List_.set("fact: relax value", (double)0.0);
    }
    else if(Smoother_=="relaxation") {
      // Gauss-Seidel smoother
      ifpack2Type_ = "RELAXATION";
      ifpack2List_.set("relaxation: sweeps", (LO) 4);
      ifpack2List_.set("relaxation: damping factor", (SC) 1.0);
      ifpack2List_.set("relaxation: type", "Gauss-Seidel");
    }
    smooProto_ = rcp( new Ifpack2Smoother(ifpack2Type_,ifpack2List_) );
    smooFact_  = rcp( new SmootherFactory(smooProto_) );    
    coarsestSmooProto_ = rcp( new DirectSolver("Superlu",coarsestSmooList_) );
    coarsestSmooFact_  = rcp( new SmootherFactory(coarsestSmooProto_, Teuchos::null) );

    // construct nullspace
    RCP<MultiVector> nullspace;
    if(Nullspace_=="constant") {
      //RCP<const Map> xmap=K_->getDomainMap();
      //nullspace = MultiVectorFactory::Build(xmap,1);
      //Scalar one=((Scalar) 1.0);
      //nullspace -> putScalar(one);      
    }
    else if(Nullspace_=="plane waves") {
      //RBMFactory_ -> BuildPlaneWaves( K_, Coords_, nullspace );
    }
    else if(Nullspace_=="rigid body modes") {
      //RBMFactory_ -> BuildRBM( K_, Coords_, nullspace );
    }

    // Set factory managers for prolongation/restriction
    if(UseLaplacian_==true && LaplaceOperatorSet_==true) {
      Hierarchy_ = rcp( new Hierarchy(L_) );
    }
    else if(StiffMatrixSet_==true) {
      Hierarchy_ = rcp( new Hierarchy(K_) );
    }
    Manager_ = rcp( new FactoryManager );

    Manager_   -> SetFactory("P", Pfact_);
    Manager_   -> SetFactory("R", Rfact_);
    Manager_   -> SetFactory("Ptent", TentPfact_);
    Manager_   -> SetFactory("Aggregates", Aggfact_);
    Hierarchy_ -> Keep("P", Pfact_.get());
    Hierarchy_ -> Keep("R", Rfact_.get());
    Hierarchy_ -> Keep("Ptent", TentPfact_.get());
    Hierarchy_ -> SetImplicitTranspose(true);
    Hierarchy_ -> Setup(*Manager_, 0, numLevels_);
    Hierarchy_ -> Delete("Smoother");
    Hierarchy_ -> Delete("CoarseSolver");    

    GridTransfersExist_=true;

  }

  // determine shifts for RAPShiftFactory
  double omega2 = omega*omega;
  shifts_.clear();
  SC shift1(alpha_,beta_);
  shifts_.push_back(-shift1*omega2);
  for(int i=0; i<numLevels_; i++) {
    double newalpha=alpha_+((double) i)*rshift_;
    double newbeta=beta_+((double) i)*ishift_;
    SC newshift(newalpha,newbeta);
    shifts_.push_back(-newshift*omega2);
  }
  Acshift_->SetShifts(shifts_);

  // Add operators together to make Helmholtz and shifted Laplace operators
  SC ii(0.0,1.0);
  if(DampMatrixSet_==true) {
    MueLu::Utils2<SC,LO,GO,NO,LMO>::TwoMatrixAdd(K_, false, (SC) 1.0, C_, (SC) ii*omega);
  }
  if(HelmholtzOperatorSet_==false) {
    MueLu::Utils2<SC,LO,GO,NO,LMO>::TwoMatrixAdd(K_, false, (SC) 1.0, M_, false, (SC) -omega2,        A_ );
    A_->fillComplete();
  }
  if(ShiftedLaplacianSet_==false) {
    MueLu::Utils2<SC,LO,GO,NO,LMO>::TwoMatrixAdd(K_, false, (SC) 1.0, M_, false, (SC) -shift1*omega2, S_ );
    S_->fillComplete();
  }

  if(VariableShift_==true) {
    // Set factories for smoothing and coarse grid
    Manager_   -> SetFactory("Smoother", smooFact_);
    Manager_   -> SetFactory("CoarseSolver", coarsestSmooFact_);
    Manager_   -> SetFactory("A", Acshift_);
    Manager_   -> SetFactory("K", Acshift_);
    Manager_   -> SetFactory("M", Acshift_);
    Hierarchy_ -> GetLevel(0)->Set("A", S_);
    Hierarchy_ -> GetLevel(0)->Set("K", K_);
    Hierarchy_ -> GetLevel(0)->Set("M", M_);
    Hierarchy_ -> Setup(*Manager_, 0, Hierarchy_ -> GetNumLevels());
  }
  else {
    // Set factories for smoothing and coarse grid
    Manager_   -> SetFactory("Smoother", smooFact_);
    Manager_   -> SetFactory("CoarseSolver", coarsestSmooFact_);
    Manager_   -> SetFactory("A", Acfact_);
    Hierarchy_ -> GetLevel(0)->Set("A", S_);
    Hierarchy_ -> Setup(*Manager_, 0, Hierarchy_ -> GetNumLevels());    
  }

  // Define Operator and Preconditioner
  RCP<TOP> BelosOper = rcp(new Belos::XpetraOp<SC,LO,GO,NO,LMO> ( A_         ) );
  RCP<TOP> BelosPrec = rcp(new Belos::MueLuOp<SC,LO,GO,NO,LMO>  ( Hierarchy_ ) );

  // Belos parameter list
  BelosList_ = rcp( new Teuchos::ParameterList("GMRES") );
  BelosList_ -> set("Maximum Iterations",iters_ );
  BelosList_ -> set("Convergence Tolerance",tol_ );
  BelosList_ -> set("Flexible Gmres", FGMRESoption_ );
  BelosList_ -> set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  BelosList_ -> set("Output Frequency",1);
  BelosList_ -> set("Output Style",Belos::Brief);

  // Belos Linear Problem and Solver Manager
  BelosLinearProblem_ = rcp( new BelosLinearProblem );
  BelosLinearProblem_ -> setOperator (  BelosOper  );
  BelosLinearProblem_ -> setRightPrec(  BelosPrec  );
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
