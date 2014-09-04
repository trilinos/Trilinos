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
#include <Xpetra_TpetraMultiVector.hpp>

// MueLu
#include "MueLu.hpp"
#include "MueLu_ConfigDefs.hpp"
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities_fwd.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_CoupledRBMFactory.hpp>
#include <MueLu_RAPShiftFactory.hpp>
#include <MueLu_PgPFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_SchwarzSmoother.hpp>
#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_ShiftedLaplacian_fwd.hpp>
#include <MueLu_ShiftedLaplacianOperator.hpp>

// Belos
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>

namespace MueLu {

  //! @brief Shifted Laplacian Helmholtz solver
  /*!
    This class provides a black box solver for indefinite Helmholtz problems.
    An AMG-Shifted Laplacian is used as a preconditioner for Krylov iterative
    solvers in Belos.
  */

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ShiftedLaplacian : public BaseClass {

#undef MUELU_SHIFTEDLAPLACIAN_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef Tpetra::Vector<SC,LO,GO,NO>                  TVEC;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>             TMV;
    typedef Tpetra::Operator<SC,LO,GO,NO>                OP;
    typedef Belos::LinearProblem<SC,TMV,OP>              LinearProblem;
    typedef Belos::SolverManager<SC,TMV,OP>              SolverManager;
    typedef Belos::SolverFactory<SC,TMV,OP>              SolverFactory;

  public:

    //! Constructors
    ShiftedLaplacian()
      : Problem_("acoustic"), numPDEs_(1), Smoother_("schwarz"), Aggregation_("uncoupled"), Nullspace_("constant"), numLevels_(5), coarseGridSize_(100),
	omega_(2.0*M_PI), ashift1_((SC) 0.0), ashift2_((SC) -1.0), pshift1_((SC) 0.0), pshift2_((SC) -1.0), iters_(500), blksize_(1),
	tol_(1.0e-4), nsweeps_(5), ncycles_(1), cycles_(8), subiters_(10), option_(1), nproblems_(0), solverType_(1), restart_size_(100), recycle_size_(25),
	smoother_sweeps_(4), smoother_damping_((SC)1.0), krylov_type_(1), krylov_iterations_(5), krylov_preconditioner_(1),
	ilu_leveloffill_(5.0), ilu_abs_thresh_(0.0), ilu_rel_thresh_(1.0), ilu_diagpivotthresh_(0.1), ilu_drop_tol_(0.01), ilu_fill_tol_(0.01), ilu_relax_val_(1.0),
	ilu_rowperm_("LargeDiag"), ilu_colperm_("COLAMD"), ilu_drop_rule_("DROP_BASIC"), ilu_normtype_("INF_NORM"), ilu_milutype_("SILU"),
	schwarz_overlap_(0), schwarz_usereorder_(true), schwarz_combinemode_(Tpetra::ADD), schwarz_ordermethod_("rcm"),
	GridTransfersExist_(false), UseLaplacian_(true), VariableShift_(false),
	LaplaceOperatorSet_(false), ProblemMatrixSet_(false), PreconditioningMatrixSet_(false),
	StiffMatrixSet_(false), MassMatrixSet_(false), DampMatrixSet_(false),
	LevelShiftsSet_(false), isSymmetric_(true), useKrylov_(true)
    { }

    // Destructor
    virtual ~ShiftedLaplacian();

    // Parameters
    void setParameters(Teuchos::RCP< Teuchos::ParameterList > paramList);

    // Set matrices
    void setLaplacian(RCP<Matrix>& L);
    void setProblemMatrix(RCP<Matrix>& A);
    void setProblemMatrix(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >& TpetraA);
    void setPreconditioningMatrix(RCP<Matrix>& P);
    void setPreconditioningMatrix(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >& TpetraP);
    void setstiff(RCP<Matrix>& K);
    void setstiff(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >& TpetraK);
    void setmass(RCP<Matrix>& M);
    void setmass(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >& TpetraM);
    void setdamp(RCP<Matrix>& C);
    void setdamp(RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >& TpetraC);
    void setcoords(RCP<MultiVector>& Coords);
    void setNullSpace(RCP<MultiVector> NullSpace);
    void setProblemShifts(Scalar ashift1, Scalar ashift2);
    void setPreconditioningShifts(Scalar pshift1, Scalar pshift2);
    void setLevelShifts(std::vector<Scalar> levelshifts);

    // various initialization/setup functions
    void initialize();
    void setupFastRAP();
    void setupSlowRAP();
    void setupNormalRAP();
    void resetLinearProblem();

    // Solve phase
    int solve(const RCP<TMV> B, RCP<TMV>& X);
    void multigrid_apply(const RCP<MultiVector> B, RCP<MultiVector>& X);
    void multigrid_apply(const RCP<Tpetra::MultiVector<SC,LO,GO,NO> > B, RCP<Tpetra::MultiVector<SC,LO,GO,NO> >& X);
    int GetIterations();
    double GetResidual();

  private:

    // Problem options
    // Problem  -> acoustic, elastic, acoustic-elastic
    // numPDEs_ -> number of DOFs at each node

    std::string Problem_;
    int numPDEs_, numSetups_;

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
    SC         ashift1_, ashift2_, pshift1_, pshift2_;
    std::vector<SC> levelshifts_;

    // Krylov solver inputs
    // iters  -> max number of iterations
    // tol    -> residual tolerance
    // FMGRES -> if true, FGMRES is chosen as solver

    int    iters_, blksize_;
    double tol_;
    int    nsweeps_, ncycles_;
    int    cycles_, subiters_, option_, nproblems_, solverType_;
    int    restart_size_, recycle_size_;

    // Smoother parameters
    int    smoother_sweeps_;
    Scalar smoother_damping_;
    int    krylov_type_;
    int    krylov_iterations_;
    int    krylov_preconditioner_;
    double ilu_leveloffill_, ilu_abs_thresh_, ilu_rel_thresh_, ilu_diagpivotthresh_;
    double ilu_drop_tol_, ilu_fill_tol_, ilu_relax_val_;
    std::string ilu_rowperm_, ilu_colperm_, ilu_drop_rule_, ilu_normtype_, ilu_milutype_;
    int    schwarz_overlap_;
    bool   schwarz_usereorder_;
    Tpetra::CombineMode schwarz_combinemode_;
    std::string schwarz_ordermethod_;

    // flags for setup
    bool GridTransfersExist_;
    bool UseLaplacian_, VariableShift_;
    bool LaplaceOperatorSet_, ProblemMatrixSet_, PreconditioningMatrixSet_;
    bool StiffMatrixSet_, MassMatrixSet_, DampMatrixSet_, LevelShiftsSet_;
    bool isSymmetric_, useKrylov_;

    // Xpetra matrices
    // K_ -> stiffness matrix
    // C_ -> damping matrix
    // M_ -> mass matrix
    // L_ -> Laplacian
    // A_ -> Problem matrix
    // P_ -> Preconditioning matrix
    RCP<Matrix>                       K_, C_, M_, L_, A_, P_;
    RCP<MultiVector>                  Coords_, NullSpace_;

    // Multigrid Hierarchy and Factory Manager
    RCP<Hierarchy>                    Hierarchy_;
    RCP<FactoryManager>               Manager_;

    // Factories and prototypes
    RCP<TentativePFactory>            TentPfact_;
    RCP<PFactory>                     Pfact_;
    RCP<PgPFactory>                   PgPfact_;
    RCP<TransPFactory>                TransPfact_;
    RCP<GenericRFactory>              Rfact_;
    RCP<RAPFactory>                   Acfact_;
    RCP<RAPShiftFactory>              Acshift_;
    RCP<CoalesceDropFactory>          Dropfact_;
    RCP<CoupledAggregationFactory>    Aggfact_;
    RCP<UncoupledAggregationFactory>  UCaggfact_;
    RCP<SmootherPrototype>            smooProto_, coarsestSmooProto_;
    RCP<SmootherFactory>              smooFact_,  coarsestSmooFact_;
    Teuchos::ParameterList            coarsestSmooList_;
    std::string                       precType_;
    Teuchos::ParameterList            precList_;

    // Operator and Preconditioner
    RCP< MueLu::ShiftedLaplacianOperator<SC,LO,GO,NO> > MueLuOp_;
    RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >           TpetraA_;

    // Belos Linear Problem and Solver
    RCP<LinearProblem>                LinearProblem_;
    RCP<SolverManager>                SolverManager_;
    RCP<SolverFactory>                SolverFactory_;
    RCP<Teuchos::ParameterList>       BelosList_;

  };

}

#define MUELU_SHIFTEDLAPLACIAN_SHORT
#endif // MUELU_SHIFTEDLAPLACIAN_DECL_HPP
