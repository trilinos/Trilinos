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
#ifndef MUELU_SHIFTEDLAPLACIAN_DECL_HPP
#define MUELU_SHIFTEDLAPLACIAN_DECL_HPP

// Xpetra
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

// MueLu
#include "MueLu.hpp"
#include "MueLu_ConfigDefs.hpp"
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities_fwd.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_CoupledRBMFactory.hpp>
#include <MueLu_RAPShiftFactory.hpp>
#include <MueLu_ShiftedLaplacian_fwd.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_ShiftedLaplacianOperator.hpp>

// Belos
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>

namespace MueLu {

  //! @brief Shifted Laplacian Helmholtz solver
  /*!
    This class provides a black box solver for indefinite Helmholtz problems.
    An AMG-Shifted Laplacian is used as a preconditioner for Krylov iterative
    solvers in Belos.
  */

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class ShiftedLaplacian : public BaseClass {

#undef MUELU_SHIFTEDLAPLACIAN_SHORT
#include "MueLu_UseShortNames.hpp"
    
    typedef Tpetra::Vector<SC,LO,GO,NO>                  TVEC;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>             TMV;
    typedef Tpetra::Operator<SC,LO,GO,NO>                OP;
    typedef Belos::LinearProblem<SC,TMV,OP>              BelosLinearProblem;
    typedef Belos::SolverManager<SC,TMV,OP>              BelosSolverManager;
    typedef Belos::BlockGmresSolMgr<SC,TMV,OP>           BelosGMRES;
    
  public:

    //! Constructors
    ShiftedLaplacian()
      : Problem_("acoustic"), numPDEs_(1), Smoother_("schwarz"), Aggregation_("coupled"), Nullspace_("constant"), numLevels_(5), coarseGridSize_(100),
	omega_(2.0*M_PI), ashift1_((SC) 0.0), ashift2_((SC) -1.0), pshift1_((SC) 0.0), pshift2_((SC) -1.0), iters_(500), blksize_(1), tol_(1.0e-4),
	nsweeps_(5), ncycles_(1), FGMRESoption_(false), cycles_(8), subiters_(10), option_(1), nproblems_(0),
	GridTransfersExist_(false), UseLaplacian_(true), VariableShift_(false),
	LaplaceOperatorSet_(false), ProblemMatrixSet_(false), PreconditioningMatrixSet_(false),
	StiffMatrixSet_(false), MassMatrixSet_(false), DampMatrixSet_(false),
	LevelShiftsSet_(false)
    { }

    // Destructor
    virtual ~ShiftedLaplacian();
    
    // Input
    void setParameters(const Teuchos::ParameterList List);

    // Set matrices
    void setLaplacian(RCP<Matrix>& L);
    void setProblemMatrix(RCP<Matrix>& A);
    void setPreconditioningMatrix(RCP<Matrix>& P);
    void setstiff(RCP<Matrix>& K);
    void setmass(RCP<Matrix>& M);
    void setdamp(RCP<Matrix>& C);

    // set parameters
    void setcoords(RCP<MultiVector>& Coords);
    void setProblemShifts(Scalar ashift1, Scalar ashift2);
    void setPreconditioningShifts(Scalar pshift1, Scalar pshift2);
    void setLevelShifts(std::vector<Scalar> levelshifts);
    void setAggregation(int stype);
    void setSmoother(int stype);
    void setSolver(int stype);
    void setSweeps(int nsweeps);
    void setCycles(int ncycles);
    void setIterations(int iters);
    void setTolerance(double tol);
    void setCoarseGridSize(int coarsegridsize);

    // various initialization/setup functions
    void initialize();
    void setupFastRAP();
    void setupSlowRAP();

    // Solve phase
    void solve(const RCP<TMV> B, RCP<TMV>& X);

  private:

    // Problem options
    // Problem  -> acoustic, elastic, acoustic-elastic
    // numPDEs_ -> number of DOFs at each node
    
    std::string Problem_; 
    int numPDEs_;
    int numSetups_;

    // Multigrid options
    // numLevels_      -> number of Multigrid levels
    // coarseGridSize_ -> size of coarsest grid (if current level has less DOFs, stop coarsening)

    std::string Smoother_, Aggregation_, Nullspace_;
    int numLevels_, coarseGridSize_;

    // Shifted Laplacian parameters
    // To be compatible with both real and complex scalar types,
    // problem and preconditioning matrices are constructed in the following way:
    //    A = K + ashift1*omega*C + (ashift2*omega^2)*M
    //    P = K + pshift1*omega*C + (pshift2*omega^2)*M
    // where K, C, and M are the stiffness, damping, and mass matrices, and
    // ashift1, ashift2, pshift1, pshift2 are user-defined scalar values.

    double     omega_;
    SC         ashift1_, ashift2_;
    SC         pshift1_, pshift2_;
    std::vector<SC> levelshifts_;

    // Krylov solver inputs
    // iters  -> max number of iterations
    // tol    -> residual tolerance
    // FMGRES -> if true, FGMRES is chosen as solver

    int    iters_;
    int    blksize_;
    double tol_;
    int    nsweeps_;
    int    ncycles_;
    bool   FGMRESoption_;
    int    cycles_;
    int    subiters_;
    int    option_;
    int    nproblems_;

    // flags for setup
    bool GridTransfersExist_;
    bool UseLaplacian_, VariableShift_;
    bool LaplaceOperatorSet_, ProblemMatrixSet_, PreconditioningMatrixSet_;
    bool StiffMatrixSet_, MassMatrixSet_, DampMatrixSet_, LevelShiftsSet_;

    // Xpetra matrices
    // K_ -> stiffness matrix
    // C_ -> damping matrix
    // M_ -> mass matrix
    // L_ -> Laplacian
    // A_ -> Problem matrix
    // P_ -> Preconditioning matrix
    RCP<Matrix>                       K_, C_, M_, L_, A_, P_;
    RCP<MultiVector>                  Coords_;

    // Multigrid Hierarchy and Factory Manager
    RCP<Hierarchy>                    Hierarchy_;
    RCP<FactoryManager>               Manager_;

    // Factories and prototypes
    RCP<TentativePFactory>            TentPfact_;
    RCP<SaPFactory>                   Pfact_;
    RCP<TransPFactory>                Rfact_;
    RCP<RAPFactory>                   Acfact_;
    RCP<RAPShiftFactory>              Acshift_;
    RCP<CoupledAggregationFactory>    Aggfact_;
    RCP<UncoupledAggregationFactory>  UCaggfact_;
    RCP<SmootherPrototype>            smooProto_, coarsestSmooProto_;
    RCP<SmootherFactory>              smooFact_,  coarsestSmooFact_;
    Teuchos::ParameterList            coarsestSmooList_;
    std::string                       ifpack2Type_;
    Teuchos::ParameterList            ifpack2List_;
    
    // Operator and Preconditioner
    RCP< MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO> > MueLuOp_;
    RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >           TpetraA_;

    // Belos Linear Problem and Solver
    RCP<BelosLinearProblem>           BelosLinearProblem_;
    RCP<BelosSolverManager>           BelosSolverManager_;
    RCP<Teuchos::ParameterList>       BelosList_;

  };

}

#define MUELU_SHIFTEDLAPLACIAN_SHORT
#endif // MUELU_SHIFTEDLAPLACIAN_DECL_HPP
