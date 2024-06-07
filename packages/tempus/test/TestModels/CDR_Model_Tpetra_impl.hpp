//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_CDR_MODEL_TPETRA_IMPL_HPP
#define TEMPUS_CDR_MODEL_TPETRA_IMPL_HPP

#include "CDR_Model_Functors.hpp"

// Thyra support
#include "Kokkos_Core_fwd.hpp"
#include "Teuchos_Assert.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_VectorStdOps.hpp"

// Tpetra support
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Tpetra_CrsGraph_def.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Tpetra_Import_def.hpp"
#include "Tpetra_Map_def.hpp"
#include "Tpetra_Vector_def.hpp"
#include "impl/Kokkos_HostThreadTeam.hpp"
#include <Teuchos_DefaultMpiComm.hpp>

namespace Tempus_Test {

// Constructor

template <typename SC, typename LO, typename GO, typename Node>
CDR_Model_Tpetra<SC, LO, GO, Node>::CDR_Model_Tpetra(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
    const GO numGlobalElements, const SC zMin, const SC zMax, const SC a,
    const SC k)
  : comm_(comm),
    numGlobalElements_(numGlobalElements),
    zMin_(zMin),
    zMax_(zMax),
    a_(a),
    k_(k),
    showGetInvalidArg_(false)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using ::Thyra::VectorBase;
  using MEB = ::Thyra::ModelEvaluatorBase;
  using ST  = Teuchos::ScalarTraits<SC>;

  TEUCHOS_ASSERT(nonnull(comm_));

  const auto num_nodes = numGlobalElements_ + 1;
  const auto myRank    = comm_->getRank();

  // owned space
  xOwnedMap_ = rcp(new tpetra_map(num_nodes, 0, comm_));
  xSpace_    = ::Thyra::createVectorSpace<SC, LO, GO, Node>(xOwnedMap_);

  // ghosted space
  if (comm_->getSize() == 1) {
    xGhostedMap_ = xOwnedMap_;
  }
  else {
    LO overlapNumMyElements;
    GO overlapGetMinGLobalIndex;
    overlapNumMyElements = xOwnedMap_->getLocalNumElements() + 2;
    if ((myRank == 0) || (myRank == (comm_->getSize() - 1))) {
      overlapNumMyElements--;
    }

    if (myRank == 0) {
      overlapGetMinGLobalIndex = xOwnedMap_->getMinGlobalIndex();
    }
    else {
      overlapGetMinGLobalIndex = xOwnedMap_->getMinGlobalIndex() - 1;
    }

    Teuchos::Array<GO> overlapMyGlobalNodes(overlapNumMyElements);

    GO getGlobalElement = overlapGetMinGLobalIndex;
    for (auto &globalElem : overlapMyGlobalNodes) {
      globalElem = overlapGetMinGLobalIndex;
      ++getGlobalElement;
    }

    const auto invalid =
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    xGhostedMap_ =
        Teuchos::rcp(new tpetra_map(invalid, overlapMyGlobalNodes, 0, comm_));
  }

  importer_ =
      Teuchos::rcp(new Tpetra::Import<LO, GO, Node>(xOwnedMap_, xGhostedMap_));

  // residual space
  fOwnedMap_ = xOwnedMap_;
  fSpace_    = xSpace_;

  x0_ = ::Thyra::createMember(xSpace_);
  V_S(x0_.ptr(), ST::zero());

  // Initialize the graph for W CrsMatrix object
  wGraph_ = createGraph();

  // Create the nodal coordinates
  nodeCoordinates_ = Teuchos::rcp(new tpetra_vec(xOwnedMap_));

  auto length      = zMax_ - zMin_;
  auto dx          = length / static_cast<SC>(numGlobalElements_ - 1);
  const auto minGI = xOwnedMap_->getMinGlobalIndex();
  {
    CoordFiller<tpetra_vec> coordFiller(*nodeCoordinates_, zMin_, dx, minGI);
    Kokkos::parallel_for("coords_fill", xOwnedMap_->getLocalNumElements(),
                         coordFiller);

    Kokkos::fence();
  }

  MEB::InArgsSetup<SC> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_t);
  inArgs.setSupports(MEB::IN_ARG_x);
  inArgs.setSupports(MEB::IN_ARG_beta);
  inArgs.setSupports(MEB::IN_ARG_x_dot);
  inArgs.setSupports(MEB::IN_ARG_alpha);
  prototypeInArgs_ = inArgs;

  MEB::OutArgsSetup<SC> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);

  prototypeOutArgs_ = outArgs;

  // Setup nominal values
  nominalValues_ = inArgs;
  nominalValues_.set_x(x0_);
  auto x_dot_init = Thyra::createMember(this->get_x_space());
  Thyra::put_scalar(SC(0.0), x_dot_init.ptr());
  nominalValues_.set_x_dot(x_dot_init);
}

// Initializers/Accessors

template <typename SC, typename LO, typename GO, typename Node>
Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, Node>>
CDR_Model_Tpetra<SC, LO, GO, Node>::createGraph()
{
  auto W_graph = Teuchos::rcp(new tpetra_graph(xOwnedMap_, 5));
  W_graph->resumeFill();

  auto overlapNumMyElements =
      static_cast<LO>(xGhostedMap_->getLocalNumElements());

  // Loop Over # of Finite Elements on Processor
  for (LO elem = 0; elem < overlapNumMyElements - 1; elem++) {
    // Loop over Nodes in Element
    for (LO i = 0; i < 2; i++) {
      auto row = xGhostedMap_->getGlobalElement(elem + i);

      // Loop over Trial Functions
      for (LO j = 0; j < 2; j++) {
        // If this row is owned by current processor, add the index
        if (xOwnedMap_->isNodeGlobalElement(row)) {
          auto colIndex = xGhostedMap_->getGlobalElement(elem + j);
          Teuchos::ArrayView<const GO> column(&colIndex, 1);
          W_graph->insertGlobalIndices(row, column);
        }
      }
    }
  }
  W_graph->fillComplete();

  return W_graph;
}

template <typename SC, typename LO, typename GO, typename Node>
void CDR_Model_Tpetra<SC, LO, GO, Node>::set_x0(
    const Teuchos::ArrayView<const SC> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(xSpace_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<SC> x0(x0_);
  x0.sv().values()().assign(x0_in);
}

template <typename SC, typename LO, typename GO, typename Node>
void CDR_Model_Tpetra<SC, LO, GO, Node>::setShowGetInvalidArgs(
    bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}

template <typename SC, typename LO, typename GO, typename Node>
void CDR_Model_Tpetra<SC, LO, GO, Node>::set_W_factory(
    const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<SC>>
        &W_factory)
{
  wFactory_ = W_factory;
}

// Public functions overridden from ModelEvaluator

template <typename SC, typename LO, typename GO, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
CDR_Model_Tpetra<SC, LO, GO, Node>::get_x_space() const
{
  return xSpace_;
}

template <typename SC, typename LO, typename GO, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>
CDR_Model_Tpetra<SC, LO, GO, Node>::get_f_space() const
{
  return fSpace_;
}

template <typename SC, typename LO, typename GO, typename Node>
Thyra::ModelEvaluatorBase::InArgs<SC>
CDR_Model_Tpetra<SC, LO, GO, Node>::getNominalValues() const
{
  return nominalValues_;
}

template <typename SC, typename LO, typename GO, typename Node>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>>
CDR_Model_Tpetra<SC, LO, GO, Node>::create_W() const
{
  auto W_factory = this->get_W_factory();

  TEUCHOS_TEST_FOR_EXCEPTION(
      is_null(W_factory), std::runtime_error,
      "W_factory in CDR_Model_Tpetra has a null W_factory!");

  auto matrix = this->create_W_op();

  return Thyra::linearOpWithSolve<SC>(*W_factory, matrix);
}

template <typename SC, typename LO, typename GO, typename Node>
Teuchos::RCP<Thyra::LinearOpBase<SC>>
CDR_Model_Tpetra<SC, LO, GO, Node>::create_W_op() const
{
  auto W_tpetra = Teuchos::rcp(new tpetra_matrix(wGraph_));

  return Thyra::tpetraLinearOp<SC, LO, GO, Node>(fSpace_, xSpace_, W_tpetra);
}

template <typename SC, typename LO, typename GO, typename Node>
Teuchos::RCP<::Thyra::PreconditionerBase<SC>>
CDR_Model_Tpetra<SC, LO, GO, Node>::create_W_prec() const
{
  auto W_op = create_W_op();
  auto prec = Teuchos::rcp(new Thyra::DefaultPreconditioner<SC>);

  prec->initializeRight(W_op);

  return prec;
}

template <typename SC, typename LO, typename GO, typename Node>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<SC>>
CDR_Model_Tpetra<SC, LO, GO, Node>::get_W_factory() const
{
  return wFactory_;
}

template <typename SC, typename LO, typename GO, typename Node>
Thyra::ModelEvaluatorBase::InArgs<SC>
CDR_Model_Tpetra<SC, LO, GO, Node>::createInArgs() const
{
  return prototypeInArgs_;
}

// Private functions overridden from ModelEvaluatorDefaultBase

template <typename SC, typename LO, typename GO, typename Node>
Thyra::ModelEvaluatorBase::OutArgs<SC>
CDR_Model_Tpetra<SC, LO, GO, Node>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}

template <typename SC, typename LO, typename GO, typename Node>
void CDR_Model_Tpetra<SC, LO, GO, Node>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));
  TEUCHOS_ASSERT(nonnull(inArgs.get_x_dot()));

  auto f_out      = outArgs.get_f();
  auto W_out      = outArgs.get_W_op();
  auto W_prec_out = outArgs.get_W_prec();

  if (nonnull(f_out) || nonnull(W_out) || nonnull(W_prec_out)) {
    // ****************
    // Get the underlying Tpetra objects
    // ****************

    RCP<tpetra_vec> f;
    if (nonnull(f_out)) {
      f = tpetra_extract::getTpetraVector(outArgs.get_f());
    }

    RCP<tpetra_matrix> J;
    if (nonnull(W_out)) {
      auto W_epetra = tpetra_extract::getTpetraOperator(W_out);
      J             = rcp_dynamic_cast<tpetra_matrix>(W_epetra);
      TEUCHOS_ASSERT(nonnull(J));
    }

    RCP<tpetra_matrix> M_inv;
    if (nonnull(W_prec_out)) {
      auto M_tpetra = tpetra_extract::getTpetraOperator(
          W_prec_out->getNonconstRightPrecOp());
      M_inv = rcp_dynamic_cast<tpetra_matrix>(M_tpetra);
      TEUCHOS_ASSERT(nonnull(M_inv));
      jDiag_ = Teuchos::rcp(new tpetra_vec(xOwnedMap_));
      jDiag_->putScalar(0.0);
    }

    // ****************
    // Create ghosted objects
    // ****************

    // Set the boundary condition directly.  Works for both x and xDot solves.
    if (comm_->getRank() == 0) {
      auto x      = Teuchos::rcp_const_cast<Thyra::VectorBase<SC>>(inArgs.get_x());
      auto xVec   = tpetra_extract::getTpetraVector(x);
      auto xView  = xVec->getLocalViewHost(Tpetra::Access::ReadWrite);
      xView(0, 0) = 1.0;
    }

    if (is_null(uPtr_)) {
      uPtr_ = Teuchos::rcp(new tpetra_vec(xGhostedMap_));
    }

    uPtr_->doImport(*(tpetra_extract::getConstTpetraVector(inArgs.get_x())),
                    *importer_, Tpetra::INSERT);

    if (is_null(uDotPtr_)) {
      uDotPtr_ = Teuchos::rcp(new tpetra_vec(xGhostedMap_));
    }

    uDotPtr_->doImport(
        *(tpetra_extract::getConstTpetraVector(inArgs.get_x_dot())), *importer_,
        Tpetra::INSERT);

    if (is_null(xPtr_)) {
      xPtr_ = Teuchos::rcp(new tpetra_vec(xGhostedMap_));
      xPtr_->doImport(*nodeCoordinates_, *importer_, Tpetra::INSERT);
    }

    auto overlapNumMyElements =
        static_cast<LO>(xGhostedMap_->getLocalNumElements());

    // Zero out the objects that will be filled
    if (nonnull(f)) {
      f->putScalar(0.0);
    }

    if (nonnull(J)) {
      J->setAllToScalar(0.0);
    }

    if (nonnull(M_inv)) {
      M_inv->setAllToScalar(0.0);
    }

    if (nonnull(f)) {
      DfDp2EvaluatorFunctor<tpetra_vec> fFunctor(*f, *xPtr_, *uPtr_, *uDotPtr_,
                                                 comm_->getRank(), a_, k_);
      Kokkos::parallel_for("DfDp2EvaluatorFunctor", overlapNumMyElements - 1,
                           fFunctor);
      Kokkos::fence();
    }

    const auto alpha = inArgs.get_alpha();
    const auto beta  = inArgs.get_beta();

    if (nonnull(J)) {
      JacobianEvaluatorFunctor<tpetra_vec, tpetra_matrix> jFunctor(
          *J, *xPtr_, *uPtr_, *uDotPtr_, comm_->getRank(), a_, k_, alpha, beta);

      J->resumeFill();
      Kokkos::parallel_for("JacobianEvaluatorFunctor", overlapNumMyElements - 1,
                           jFunctor);
      Kokkos::fence();
      J->fillComplete();
    }

    if (nonnull(M_inv)) {
      PreconditionerEvaluatorFunctor<tpetra_vec, tpetra_matrix> mFunctor(
          *M_inv, *xPtr_, *uPtr_, *uDotPtr_, comm_->getRank(), a_, k_, alpha,
          beta);

      M_inv->resumeFill();
      Kokkos::parallel_for("PreconditionerEvaluatorFunctor",
                           overlapNumMyElements - 1, mFunctor);
      Kokkos::fence();
      M_inv->fillComplete();
    }

    if (nonnull(M_inv)) {
      // Invert the Jacobian diagonal for the preconditioner
      auto &diag = *jDiag_;
      M_inv->getLocalDiagCopy(diag);
      diag.reciprocal(diag);

      M_inv->rightScale(diag);
      M_inv->rightScale(diag);
    }
  }
}

}  // namespace Tempus_Test

#endif  // TEMPUS_CDR_MODEL_TPETRA_IMPL_HPP
