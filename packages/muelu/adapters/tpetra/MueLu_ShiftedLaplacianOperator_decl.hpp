// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SHIFTEDLAPLACIANOPERATOR_DECL_HPP
#define MUELU_SHIFTEDLAPLACIANOPERATOR_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy_decl.hpp"
#include "MueLu_Utilities.hpp"

// TODO: Kokkos headers

namespace MueLu {

/*! @brief Wraps an existing MueLu::Hierarchy as a Tpetra::Operator, with an optional two-level correction.
    Intended to be used with MueLu::ShiftedLaplacian.
*/
template <class Scalar        = Tpetra::Operator<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::Operator<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class ShiftedLaplacianOperator
  : public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
  typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> OP;
  typedef MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node> MUtils;

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  ShiftedLaplacianOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H)
    : Hierarchy_(H)
    , option_(0) {}

  //! Auxiliary Constructor
  ShiftedLaplacianOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H,
                           const RCP<Matrix> A, int cycles, int iters, int option, double tol)
    : Hierarchy_(H)
    , A_(A)
    , cycles_(cycles)
    , iters_(iters)
    , option_(option)
    , tol_(tol) {
    // setup 2-level correction
    /*RCP< MueLu::Level > Level1 = H -> GetLevel(1);
    R_ = Level1 -> Get< RCP<Matrix> >("R");
    P_ = Level1 -> Get< RCP<Matrix> >("P");
    //RCP<Matrix> AP = Level1 -> Get< RCP<Matrix> >("AP graph");
    RCP<Matrix> AP;
    AP = MUtils::Multiply(*A_, false, *P_, false, AP);
    // Optimization storage option. If matrix is not changing later, allow this.
    bool doOptimizedStorage = true;
    // Reuse coarse matrix memory if available (multiple solve)
    //RCP<Matrix> Ac = Level1 -> Get< RCP<Matrix> >("RAP graph");
    RCP<Matrix> Ac;
    Ac = MUtils::Multiply(*R_, false, *AP, false, Ac, true, doOptimizedStorage);
    Ac_ = MUtils::Op2NonConstTpetraCrs(Ac);

    // Setup Belos for two-level correction
    BelosList_ = rcp( new Teuchos::ParameterList("GMRES") );
    BelosList_ -> set("Maximum Iterations", iters_ );
    BelosList_ -> set("Convergence Tolerance", tol_ );
    BelosLP_   = rcp( new Belos::LinearProblem<Scalar,MV,OP> );
    BelosLP_   -> setOperator ( Ac_ );
    BelosSM_   = rcp( new Belos::BlockGmresSolMgr<Scalar,MV,OP>(BelosLP_, BelosList_) );*/
  }

  //! Destructor.
  virtual ~ShiftedLaplacianOperator() {}

  //@}

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const;

  //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
  /*!
    \param[in] X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.

  */
  void apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::one()) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const;

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > R_, P_, A_;
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ac_;
  RCP<Teuchos::ParameterList> BelosList_;

  // RCP< Belos::LinearProblem<Scalar,MV,OP> > BelosLP_;
  // RCP< Belos::SolverManager<Scalar,MV,OP> > BelosSM_;

  // cycles -> number of V-cycles
  // iters  -> number of GMRES iterations per correction
  // option -> 0 if no correction is desired
  int cycles_, iters_, option_;
  double tol_;
};

}  // namespace MueLu

#endif  // MUELU_SHIFTEDLAPLACIANOPERATOR_DECL_HPP
