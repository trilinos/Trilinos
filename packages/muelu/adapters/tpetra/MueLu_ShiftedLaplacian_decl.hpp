// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#if defined(HAVE_MUELU_IFPACK2)

#include <MueLu_BaseClass.hpp>
#include <MueLu_AmalgamationFactory_fwd.hpp>
#include <MueLu_CoalesceDropFactory_fwd.hpp>
#include <MueLu_CoarseMapFactory_fwd.hpp>
#include <MueLu_CoupledRBMFactory_fwd.hpp>
#include <MueLu_DirectSolver_fwd.hpp>
#include <MueLu_GenericRFactory_fwd.hpp>
#include <MueLu_Hierarchy_fwd.hpp>
#include <MueLu_Ifpack2Smoother_fwd.hpp>
#include <MueLu_PFactory_fwd.hpp>
#include <MueLu_PgPFactory_fwd.hpp>
#include <MueLu_RAPFactory_fwd.hpp>
#include <MueLu_RAPShiftFactory_fwd.hpp>
#include <MueLu_SaPFactory_fwd.hpp>
#include <MueLu_ShiftedLaplacian_fwd.hpp>
#include <MueLu_ShiftedLaplacianOperator.hpp>
#include <MueLu_SmootherFactory_fwd.hpp>
#include <MueLu_SmootherPrototype_fwd.hpp>
#include <MueLu_TentativePFactory_fwd.hpp>
#include <MueLu_TransPFactory_fwd.hpp>
#include <MueLu_UncoupledAggregationFactory_fwd.hpp>
#include <MueLu_Utilities_fwd.hpp>

// Belos
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#endif

namespace MueLu {

/*!
  @brief Shifted Laplacian Helmholtz solver

  This class provides a black box solver for indefinite Helmholtz problems.
  An AMG-Shifted Laplacian is used as a preconditioner for Krylov iterative
  solvers in Belos.

  @ingroup MueLuAdapters
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class ShiftedLaplacian : public BaseClass {
#undef MUELU_SHIFTEDLAPLACIAN_SHORT
#include "MueLu_UseShortNames.hpp"

  typedef Tpetra::Vector<SC, LO, GO, NO> TVEC;
  typedef Tpetra::MultiVector<SC, LO, GO, NO> TMV;
  typedef Tpetra::Operator<SC, LO, GO, NO> OP;
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  typedef Belos::LinearProblem<SC, TMV, OP> LinearProblem;
  typedef Belos::SolverManager<SC, TMV, OP> SolverManager;
  typedef Belos::SolverFactory<SC, TMV, OP> SolverFactory;
#endif

 public:
  /*
    FIXME  26-June-2015 JJH:  This contructor is setting numerous defaults.  However, they don't match the defaults
    FIXME  int the method setParameters().  There also isn't any parameter validation that I can see.
  */

  //! Constructors
  ShiftedLaplacian()
    : numPDEs_(1)
    , Smoother_("schwarz")
    , Aggregation_("uncoupled")
    , Nullspace_("constant")
    , numLevels_(5)
    , coarseGridSize_(100)
    , omega_(2.0 * M_PI)
    , iters_(500)
    , blksize_(1)
    , tol_(1.0e-4)
    , nsweeps_(5)
    , ncycles_(1)
    , cycles_(8)
    , subiters_(10)
    , option_(1)
    , nproblems_(0)
    , solverType_(1)
    , restart_size_(100)
    , recycle_size_(25)
    , smoother_sweeps_(4)
    , smoother_damping_((SC)1.0)
    , krylov_type_(1)
    , krylov_iterations_(5)
    , krylov_preconditioner_(1)
    , ilu_leveloffill_(5.0)
    , ilu_abs_thresh_(0.0)
    , ilu_rel_thresh_(1.0)
    , ilu_diagpivotthresh_(0.1)
    , ilu_drop_tol_(0.01)
    , ilu_fill_tol_(0.01)
    , ilu_relax_val_(1.0)
    , ilu_rowperm_("LargeDiag")
    , ilu_colperm_("COLAMD")
    , ilu_drop_rule_("DROP_BASIC")
    , ilu_normtype_("INF_NORM")
    , ilu_milutype_("SILU")
    , schwarz_overlap_(0)
    , schwarz_usereorder_(true)
    , schwarz_combinemode_(Tpetra::ADD)
    , schwarz_ordermethod_("rcm")
    , GridTransfersExist_(false)
    , isSymmetric_(true) {}

  // Destructor
  virtual ~ShiftedLaplacian();

  // Parameters
  void setParameters(Teuchos::RCP<Teuchos::ParameterList> paramList);

  // Set matrices
  void setProblemMatrix(RCP<Matrix>& A);
  void setProblemMatrix(RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >& TpetraA);
  void setPreconditioningMatrix(RCP<Matrix>& P);
  void setPreconditioningMatrix(RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >& TpetraP);
  void setstiff(RCP<Matrix>& K);
  void setstiff(RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >& TpetraK);
  void setmass(RCP<Matrix>& M);
  void setmass(RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> >& TpetraM);
  void setcoords(RCP<MultiVector>& Coords);
  void setNullSpace(RCP<MultiVector> NullSpace);
  void setLevelShifts(std::vector<Scalar> levelshifts);

  // initialize: set parameters and factories, construct
  // prolongation and restriction matrices
  void initialize();
  // setupFastRAP: setup hierarchy with
  // prolongators of the stiffness matrix
  // constant complex shifts
  void setupFastRAP();
  // setupSlowRAP: setup hierarchy with
  // prolongators of the stiffness matrix
  // variable complex shifts
  void setupSlowRAP();
  // setupNormalRAP: setup hierarchy with
  // prolongators of the preconditioning matrix
  void setupNormalRAP();
  // setupSolver: initialize Belos solver
  void setupSolver();
  // resetLinearProblem: for multiple frequencies;
  // reset the Belos operator if the frequency changes
  void resetLinearProblem();

  // Solve phase
  int solve(const RCP<TMV> B, RCP<TMV>& X);
  void multigrid_apply(const RCP<MultiVector> B,
                       RCP<MultiVector>& X);
  void multigrid_apply(const RCP<Tpetra::MultiVector<SC, LO, GO, NO> > B,
                       RCP<Tpetra::MultiVector<SC, LO, GO, NO> >& X);
  int GetIterations();
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType GetResidual();

  RCP<FactoryManager> Manager_;

 private:
  // Problem options
  // numPDEs_ -> number of DOFs at each node
  int numPDEs_;

  // Multigrid options
  // numLevels_      -> number of Multigrid levels
  // coarseGridSize_ -> size of coarsest grid (if current level has less DOFs, stop coarsening)
  std::string Smoother_, Aggregation_, Nullspace_;
  int numLevels_, coarseGridSize_;

  // Shifted Laplacian/Helmholtz parameters
  double omega_;
  std::vector<SC> levelshifts_;

  // Krylov solver inputs
  // iters  -> max number of iterations
  // tol    -> residual tolerance
  int iters_, blksize_;
  double tol_;
  int nsweeps_, ncycles_;
  int cycles_, subiters_, option_, nproblems_, solverType_;
  int restart_size_, recycle_size_;

  // Smoother parameters
  int smoother_sweeps_;
  Scalar smoother_damping_;
  int krylov_type_;
  int krylov_iterations_;
  int krylov_preconditioner_;
  double ilu_leveloffill_, ilu_abs_thresh_, ilu_rel_thresh_, ilu_diagpivotthresh_;
  double ilu_drop_tol_, ilu_fill_tol_, ilu_relax_val_;
  std::string ilu_rowperm_, ilu_colperm_, ilu_drop_rule_, ilu_normtype_, ilu_milutype_;
  int schwarz_overlap_;
  bool schwarz_usereorder_;
  Tpetra::CombineMode schwarz_combinemode_;
  std::string schwarz_ordermethod_;

  // flags for setup
  bool GridTransfersExist_;
  bool isSymmetric_;

  // Xpetra matrices
  // K_ -> stiffness matrix
  // M_ -> mass matrix
  // A_ -> Problem matrix
  // P_ -> Preconditioning matrix
  RCP<Matrix> K_, M_, A_, P_;
  RCP<MultiVector> Coords_, NullSpace_;

  // Multigrid Hierarchy
  RCP<Hierarchy> Hierarchy_;

  // Factories and prototypes
  RCP<TentativePFactory> TentPfact_;
  RCP<PFactory> Pfact_;
  RCP<PgPFactory> PgPfact_;
  RCP<TransPFactory> TransPfact_;
  RCP<GenericRFactory> Rfact_;
  RCP<RAPFactory> Acfact_;
  RCP<RAPShiftFactory> Acshift_;
  RCP<AmalgamationFactory> Amalgfact_;
  RCP<CoalesceDropFactory> Dropfact_;
  RCP<UncoupledAggregationFactory> UCaggfact_;
  RCP<CoarseMapFactory> CoarseMapfact_;
  RCP<SmootherPrototype> smooProto_, coarsestSmooProto_;
  RCP<SmootherFactory> smooFact_, coarsestSmooFact_;
  Teuchos::ParameterList coarsestSmooList_;
  std::string precType_;
  Teuchos::ParameterList precList_;

  // Operator and Preconditioner
  RCP<MueLu::ShiftedLaplacianOperator<SC, LO, GO, NO> > MueLuOp_;
  RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > TpetraA_;

#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  // Belos Linear Problem and Solver
  RCP<LinearProblem> LinearProblem_;
  RCP<SolverManager> SolverManager_;
  RCP<SolverFactory> SolverFactory_;
  RCP<Teuchos::ParameterList> BelosList_;
#endif
};

}  // namespace MueLu

#define MUELU_SHIFTEDLAPLACIAN_SHORT

#endif  // if defined(HAVE_MUELU_IFPACK2) and defined(HAVE_MUELU_TPETRA)

#endif  // MUELU_SHIFTEDLAPLACIAN_DECL_HPP
