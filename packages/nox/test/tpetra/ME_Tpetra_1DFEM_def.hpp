#ifndef NOX_TPETRA_ME_1DFEM_DEF_HPP
#define NOX_TPETRA_ME_1DFEM_DEF_HPP

// Thyra support
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Tpetra_Import_Util2.hpp"  //for sortCrsEntries

// Tpetra support
#include "Thyra_TpetraThyraWrappers.hpp"

// Kokkos support
#include "Kokkos_Core.hpp"

#include "1DFEM_Functors.hpp"

// Nonmember constuctors

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar, LO, GO, Node> >
evaluatorTpetra1DFEM(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const Tpetra::global_size_t numGlobalElements,
                     const Scalar zMin,
                     const Scalar zMax)
{
  return Teuchos::rcp(new EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>(comm,numGlobalElements,zMin,zMax));
}

// Constructor

template<class Scalar, class LO, class GO, class Node>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::
EvaluatorTpetra1DFEM(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const Tpetra::global_size_t numGlobalElements,
                     const Scalar zMin,
                     const Scalar zMax) :
  comm_(comm),
  numGlobalElements_(numGlobalElements),
  zMin_(zMin),
  zMax_(zMax),
  showGetInvalidArg_(false)
{
  typedef ::Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<scalar_type> ST;

  TEUCHOS_ASSERT(nonnull(comm_));

  const Tpetra::global_size_t numNodes = numGlobalElements_ + 1;

  // owned space
  GO indexBase = 0;
  xOwnedMap_ = Teuchos::rcp(new const tpetra_map(numNodes, indexBase, comm_));
  xSpace_ = ::Thyra::createVectorSpace<Scalar, LO, GO, Node>(xOwnedMap_);

  // ghosted space
  if (comm_->getSize() == 1) {
    xGhostedMap_ = xOwnedMap_;
  }
  else {
    std::size_t overlapNumMyNodes;
    GO overlapMinMyGID;
    overlapNumMyNodes = xOwnedMap_->getNodeNumElements() + 2;
    if ( (comm_->getRank() == 0) || (comm_->getRank() == (comm_->getSize() - 1)) )
      --overlapNumMyNodes;

    if (comm_->getRank() == 0)
      overlapMinMyGID = xOwnedMap_->getMinGlobalIndex();
    else
      overlapMinMyGID = xOwnedMap_->getMinGlobalIndex() - 1;

    Teuchos::Array<GO> overlapMyGlobalNodes(overlapNumMyNodes);
    GO gid = overlapMinMyGID;
    for (auto gidIter = overlapMyGlobalNodes.begin(); gidIter != overlapMyGlobalNodes.end(); ++gidIter) {
      *gidIter = gid;
      ++gid;
    }

    const Tpetra::global_size_t invalid = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    xGhostedMap_ = Teuchos::rcp(new tpetra_map(invalid, overlapMyGlobalNodes, indexBase, comm_));
  }

  importer_ = Teuchos::rcp(new Tpetra::Import<LO,GO,Node>(xOwnedMap_, xGhostedMap_));

  // residual space
  fOwnedMap_ = xOwnedMap_;
  fSpace_ = xSpace_;

  x0_ = ::Thyra::createMember(xSpace_);
  V_S(x0_.ptr(), ST::zero());

  // Initialize the graph for W CrsMatrix object
  W_graph_ = createGraph();

  // Create the nodal coorinates
  std::size_t numLocalNodes = xOwnedMap_->getNodeNumElements();
  GO minGID = xOwnedMap_->getMinGlobalIndex();
  Scalar dz = (zMax_ - zMin_)/static_cast<Scalar>(numGlobalElements_);
  nodeCoordinates_ = Teuchos::rcp(new tpetra_vec(xOwnedMap_));
  nodeCoordinates_->template modify<typename tpetra_vec::execution_space>();

  MeshFillFunctor<tpetra_vec> functor(*nodeCoordinates_, zMin_, dz, minGID);
  Kokkos::parallel_for("coords fill", numLocalNodes, functor);

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  nominalValues_.set_x(x0_);

  residTimer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Residual Evaluation");
  jacTimer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Jacobian Evaluation");
}

// Initializers/Accessors

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, Node> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::createGraph()
{
  typedef typename tpetra_graph::local_graph_type::size_type size_type;

  // Compute graph offset array
  int numProcs = comm_->getSize();
  int myRank = comm_->getRank();
  std::size_t numMyNodes = xOwnedMap_->getNodeNumElements();
  std::size_t numLocalEntries = 0;
  //Kokkos::View<std::size_t*> counts("row counts", numMyNodes);
  Kokkos::View<size_type*> counts("row counts", numMyNodes);
  {
    RowCountsFunctor<size_type, LO> functor(counts, numMyNodes, numProcs, myRank);
    Kokkos::parallel_reduce("row counts comp", numMyNodes, functor, numLocalEntries);
  }

  //Kokkos::View<std::size_t*> offsets("row offsets", numMyNodes+1);
  Kokkos::View<size_type*> offsets("row offsets", numMyNodes+1);
  {
    RowOffsetsFunctor<size_type, LO> functor(offsets, counts, numMyNodes);
    Kokkos::parallel_scan("row offsets comp", numMyNodes+1, functor);
  }

  // Create array of non-zero entry column indices
  Kokkos::View<LO*> indices("column indices", numLocalEntries);
  //typename local_graph_type::entries_type::non_const_type indices("column indices", numLocalEntries);
  {
    ColumnIndexCompFunctor<size_type, LO> functor(indices, offsets, counts, numMyNodes, numProcs, myRank);
    Kokkos::parallel_for("column indices comp", numMyNodes, functor);
  }

  //Sort the indices within each row.
  Tpetra::Import_Util::sortCrsEntries(offsets, indices);

  // Construct the graph
  Teuchos::RCP<tpetra_graph> W_graph =
    Teuchos::rcp(new tpetra_graph(xOwnedMap_, xGhostedMap_, offsets, indices));
  W_graph->fillComplete();
  return W_graph;
}

template<class Scalar, class LO, class GO, class Node>
void EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::
set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(xSpace_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


template<class Scalar, class LO, class GO, class Node>
void EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::
setShowGetInvalidArgs(bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}

template<class Scalar, class LO, class GO, class Node>
void EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::
set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory)
{
  W_factory_ = W_factory;
}

// Public functions overridden from ModelEvaulator


template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_x_space() const
{
  return xSpace_;
}


template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_f_space() const
{
  return fSpace_;
}


template<class Scalar, class LO, class GO, class Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::create_W_op() const
{
  Teuchos::RCP<tpetra_matrix> W_tpetra = Teuchos::rcp(new tpetra_matrix(W_graph_));

  return Thyra::tpetraLinearOp<Scalar, LO, GO, Node>(fSpace_, xSpace_, W_tpetra);
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP< ::Thyra::PreconditionerBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::create_W_prec() const
{
  Teuchos::RCP<tpetra_matrix> W_tpetra = Teuchos::rcp(new tpetra_matrix(W_graph_));

  Teuchos::RCP<thyra_op> W_op = Thyra::tpetraLinearOp<Scalar, LO, GO, Node>(fSpace_, xSpace_, W_tpetra);

  Teuchos::RCP<Thyra::DefaultPreconditioner<Scalar> > prec =
    Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>);

  prec->initializeRight(W_op);

  return prec;
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar, class LO, class GO, class Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar, class LO, class GO, class Node>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar, class LO, class GO, class Node>
void EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));

  const Teuchos::RCP<thyra_vec> f_out = outArgs.get_f();
  const Teuchos::RCP<thyra_op> W_out = outArgs.get_W_op();
  const Teuchos::RCP<thyra_prec> W_prec_out = outArgs.get_W_prec();

  const bool fill_f = nonnull(f_out);
  const bool fill_W = nonnull(W_out);
  const bool fill_W_prec = nonnull(W_prec_out);

  typedef Tpetra::Operator<Scalar,LO,GO,Node> tpetra_op;
  typedef ::Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node> tpetra_extract;

  if ( fill_f || fill_W || fill_W_prec ) {

    // Get the underlying tpetra objects
    Teuchos::RCP<tpetra_vec> f;
    if (fill_f) {
      f = tpetra_extract::getTpetraVector(f_out);
    }

    Teuchos::RCP<tpetra_matrix> J;
    if (fill_W) {
      Teuchos::RCP<tpetra_op> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
      J = Teuchos::rcp_dynamic_cast<tpetra_matrix>(W_tpetra);
      TEUCHOS_ASSERT(nonnull(J));
      J->resumeFill();
    }

    Teuchos::RCP<tpetra_matrix> M_inv;
    if (fill_W_prec) {
      Teuchos::RCP<tpetra_op> M_tpetra = tpetra_extract::getTpetraOperator(W_prec_out->getNonconstRightPrecOp());
      M_inv = Teuchos::rcp_dynamic_cast<tpetra_matrix>(M_tpetra);
      TEUCHOS_ASSERT(nonnull(M_inv));
      M_inv->resumeFill();

      if (is_null(J_diagonal_))
        J_diagonal_ = Teuchos::rcp(new tpetra_vec(xOwnedMap_));
    }

    //typedef Kokkos::HostSpace host_space;
    typedef typename tpetra_vec::execution_space execution_space;

    // Create ghosted objects
    if (is_null(uPtr_))
      uPtr_ = Teuchos::rcp(new tpetra_vec(xGhostedMap_));

    uPtr_->doImport(*(tpetra_extract::getConstTpetraVector(inArgs.get_x())), *importer_, Tpetra::REPLACE);

    if (is_null(xPtr_)) {
      xPtr_ = Teuchos::rcp(new tpetra_vec(xGhostedMap_));
      xPtr_->doImport(*nodeCoordinates_, *importer_, Tpetra::INSERT);
    }

    // Zero out the objects that will be filled
    if (fill_f) {
      f->putScalar(0.0);
    }
    if (fill_W) {
      J->setAllToScalar(0.0);
    }
    if (fill_W_prec) {
      M_inv->setAllToScalar(0.0);
      J_diagonal_->putScalar(0.0);
    }

    // Get local Views of data
    int myRank = comm_->getRank();
    std::size_t numMyElements = xGhostedMap_->getNodeNumElements()-1;

    xPtr_->template sync<execution_space>();
    uPtr_->template sync<execution_space>();

    // Residual fill
    if (fill_f) {
      Teuchos::TimeMonitor timer(*residTimer_);
      f->template sync<execution_space>();
      f->template modify<execution_space>();

      ResidualEvaluatorFunctor<tpetra_vec> functor(*f, *xPtr_, *uPtr_, myRank);
      Kokkos::parallel_for("residual evaluation", numMyElements, functor);
      Kokkos::fence();
    }

    // Jacobian fill
    if (fill_W) {
      Teuchos::TimeMonitor timer(*jacTimer_);
      JacobianEvaluatorFunctor<tpetra_vec, tpetra_matrix> functor(*J, *xPtr_, *uPtr_, myRank);
      Kokkos::parallel_for("jacobian evaluation", numMyElements, functor);
      Kokkos::fence();
    }

    // Preconditioner fill
    if (fill_W_prec) {
      PreconditionerEvaluatorFunctor<tpetra_vec, tpetra_matrix> functor(*M_inv, *xPtr_, *uPtr_, myRank);
      Kokkos::parallel_for("prec evaluation", numMyElements, functor);
    }

    if (fill_W) {
      J->fillComplete();
    }

    if (fill_W_prec) {
      // Invert the Jacobian diagonal for the preconditioner
      // For some reason the matrix must be fill complete before calling rightScale
      M_inv->fillComplete();
      tpetra_vec& diag = *J_diagonal_;
      M_inv->getLocalDiagCopy(diag);
      diag.reciprocal(diag);
      M_inv->rightScale(diag);
      M_inv->rightScale(diag);
    }

  }

}

#endif
