// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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
#include "NOX_TpetraTypedefs.hpp"

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
  Np_(5),
  Ng_(7),
  printDebug_(false),
  showGetInvalidArg_(false),
  pNames_(Np_),
  gNames_(Ng_)
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
    overlapNumMyNodes = xOwnedMap_->getLocalNumElements() + 2;
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
  std::size_t numLocalNodes = xOwnedMap_->getLocalNumElements();
  GO minGID = xOwnedMap_->getMinGlobalIndex();
  Scalar dz = (zMax_ - zMin_)/static_cast<Scalar>(numGlobalElements_);
  nodeCoordinates_ = Teuchos::rcp(new tpetra_vec(xOwnedMap_));

  MeshFillFunctor<tpetra_vec> functor(*nodeCoordinates_, zMin_, dz, minGID);
  Kokkos::parallel_for("coords fill", numLocalNodes, functor);

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  inArgs.set_Np_Ng(Np_,Ng_);
  prototypeInArgs_ = inArgs;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);
  outArgs.set_Np_Ng(Np_,Ng_);
  outArgs.setSupports(MEB::OUT_ARG_DfDp,2,MEB::DerivativeSupport(MEB::DERIV_MV_JACOBIAN_FORM));
  outArgs.setSupports(MEB::OUT_ARG_DgDx,4,MEB::DerivativeSupport(MEB::DERIV_MV_GRADIENT_FORM));
  outArgs.setSupports(MEB::OUT_ARG_DgDp,4,2,MEB::DerivativeSupport(MEB::DERIV_MV_JACOBIAN_FORM));

  outArgs.setSupports(MEB::OUT_ARG_DfDp,4,MEB::DerivativeSupport(MEB::DERIV_MV_JACOBIAN_FORM));
  outArgs.setSupports(MEB::OUT_ARG_DgDx,6,MEB::DerivativeSupport(MEB::DERIV_MV_GRADIENT_FORM));
  outArgs.setSupports(MEB::OUT_ARG_DgDp,4,4,MEB::DerivativeSupport(MEB::DERIV_MV_JACOBIAN_FORM));
  outArgs.setSupports(MEB::OUT_ARG_DgDp,6,2,MEB::DerivativeSupport(MEB::DERIV_MV_JACOBIAN_FORM));
  outArgs.setSupports(MEB::OUT_ARG_DgDp,6,4,MEB::DerivativeSupport(MEB::DERIV_MV_JACOBIAN_FORM));

  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  nominalValues_.set_x(x0_);

  residTimer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Residual Evaluation");
  jacTimer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Jacobian Evaluation");

  // Parameter and response support. There exists one parameter and one response.
  for (auto& p : pNames_)
    p = Teuchos::rcp(new Teuchos::Array<std::string>);
  pNames_[0]->push_back("Dummy p(0)");
  pNames_[1]->push_back("Dummy p(1)");
  pNames_[2]->push_back("k");
  pNames_[3]->push_back("Dummy p(3)");
  pNames_[4]->push_back("T_left");
  pMap_ = Teuchos::rcp(new const tpetra_map(1, 0, comm_, Tpetra::LocallyReplicated));
  pSpace_ = ::Thyra::createVectorSpace<Scalar, LO, GO, Node>(pMap_);
  p2_ = ::Thyra::createMember(pSpace_);
  V_S(p2_.ptr(),1.0);
  nominalValues_.set_p(2,p2_);
  p4_ = ::Thyra::createMember(pSpace_);
  V_S(p4_.ptr(),1.0);
  nominalValues_.set_p(4,p4_);

  for (auto& g : gNames_)
    g.clear();
  gNames_[0].push_back("Dummy g(0)");
  gNames_[1].push_back("Dummy g(1)");
  gNames_[2].push_back("Dummy g(2)");
  gNames_[3].push_back("Dummy g(3)");
  gNames_[4].push_back("Constraint: T_right=2");
  gNames_[5].push_back("Dummy g(5)");
  gNames_[6].push_back("Constraint: 2*T_left=T_right");
  gMap_ = Teuchos::rcp(new const tpetra_map(1, 0, comm_, Tpetra::LocallyReplicated));
  gSpace_ = ::Thyra::createVectorSpace<Scalar, LO, GO, Node>(gMap_);
  dgdpMap_ = Teuchos::rcp(new const tpetra_map(1, 0, comm_, Tpetra::LocallyReplicated));
  dgdpSpace_ = ::Thyra::createVectorSpace<Scalar, LO, GO, Node>(dgdpMap_);

  p_name_to_index_["Dummy p(0)"] = std::make_pair(0,0);
  p_name_to_index_["Dummy p(1)"] = std::make_pair(1,0);
  p_name_to_index_["k"] = std::make_pair(2,0);
  p_name_to_index_["Dummy p(3)"] = std::make_pair(3,0);
  p_name_to_index_["T_left"] = std::make_pair(4,0);

  g_name_to_index_["Dummy g(0)"] = std::make_pair(0,0);
  g_name_to_index_["Dummy g(1)"] = std::make_pair(1,0);
  g_name_to_index_["Dummy g(2)"] = std::make_pair(2,0);
  g_name_to_index_["Dummy g(3)"] = std::make_pair(3,0);
  g_name_to_index_["Constraint: T_right=2"] = std::make_pair(4,0);
  g_name_to_index_["Dummy g(5)"] = std::make_pair(5,0);
  g_name_to_index_["Constraint: 2*T_left=T_right"] = std::make_pair(6,0);
}

// Initializers/Accessors

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, Node> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::createGraph()
{
  using size_type = typename tpetra_graph::local_graph_device_type::size_type;

  // Compute graph offset array
  int numProcs = comm_->getSize();
  int myRank = comm_->getRank();
  std::size_t numMyNodes = xOwnedMap_->getLocalNumElements();
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
int EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::Np() const
{return Np_;}

template<class Scalar, class LO, class GO, class Node>
int EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::Ng() const
{return Ng_;}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_p_space(int /* l */) const
{
  // All parameters are locally replicated scalars of size 1
  return pSpace_;
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Teuchos::Array<std::string> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_p_names(int l) const
{
  return pNames_[l];
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_g_space(int /* j */) const
{
  // All parameters are locally replicated scalars of size 1
  return gSpace_;
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::ArrayView<const std::string>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_g_names(int j) const
{
  return gNames_[j];
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<::Thyra::LinearOpBase<Scalar>>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::create_DfDp_op(int l) const
{
  TEUCHOS_ASSERT( (l == 2) || (l == 4) );
  return ::Thyra::createMembers(xSpace_,1,"LOCA::DgDx");
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<::Thyra::LinearOpBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::create_DgDx_op(int j) const
{
  TEUCHOS_ASSERT( (j == 4) || (j == 6) );
  return ::Thyra::createMembers(xSpace_,1,"LOCA::DgDx");
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<::Thyra::LinearOpBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::create_DgDx_dot_op(int j) const
{
  TEUCHOS_ASSERT( (j == 4) || (j == 6) );
  return ::Thyra::createMembers(xSpace_,1,"LOCA::DgDx_dot");
}

template<class Scalar, class LO, class GO, class Node>
::Teuchos::RCP<::Thyra::LinearOpBase<Scalar> >
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::create_DgDp_op( int j, int l ) const
{
  TEUCHOS_ASSERT( (j == 4) || (j == 6) );
  TEUCHOS_ASSERT( (l == 2) || (l == 4) );
  // Instead of using a dense serial matrix, we use a locally
  // replicated multivector since it provides thyra wrappers.
  return ::Thyra::createMembers(gSpace_,1,"LOCA::DgDp");
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >
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

template<class Scalar, class LO, class GO, class Node>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::createOutArgs() const
{
  return prototypeOutArgs_;
}

template<class Scalar, class LO, class GO, class Node>
void EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::
evalModel(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
          const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));

  const Teuchos::RCP<thyra_vec> f_out = outArgs.get_f();
  const Teuchos::RCP<thyra_op> W_out = outArgs.get_W_op();
  const Teuchos::RCP<thyra_prec> W_prec_out = outArgs.get_W_prec();
  const Teuchos::RCP<thyra_vec> g4_out = outArgs.get_g(4);
  const Teuchos::RCP<thyra_vec> g6_out = outArgs.get_g(6);
  const Teuchos::RCP<thyra_mvec> Dg4Dx_out = outArgs.get_DgDx(4).getMultiVector();
  const Teuchos::RCP<thyra_mvec> Dg6Dx_out = outArgs.get_DgDx(6).getMultiVector();
  const Teuchos::RCP<thyra_mvec> DfDp2_out = outArgs.get_DfDp(2).getMultiVector();
  const Teuchos::RCP<thyra_mvec> DfDp4_out = outArgs.get_DfDp(4).getMultiVector();
  const Teuchos::RCP<thyra_mvec> Dg4Dp2_out = outArgs.get_DgDp(4,2).getMultiVector();
  const Teuchos::RCP<thyra_mvec> Dg4Dp4_out = outArgs.get_DgDp(4,4).getMultiVector();
  const Teuchos::RCP<thyra_mvec> Dg6Dp2_out = outArgs.get_DgDp(6,2).getMultiVector();
  const Teuchos::RCP<thyra_mvec> Dg6Dp4_out = outArgs.get_DgDp(6,4).getMultiVector();

  const bool fill_f = nonnull(f_out);
  const bool fill_W = nonnull(W_out);
  const bool fill_W_prec = nonnull(W_prec_out);
  const bool fill_g4 = nonnull(g4_out);
  const bool fill_g6 = nonnull(g6_out);
  const bool fill_Dg4Dx = nonnull(Dg4Dx_out);
  const bool fill_Dg6Dx = nonnull(Dg6Dx_out);
  const bool fill_DfDp2 = nonnull(DfDp2_out);
  const bool fill_DfDp4 = nonnull(DfDp4_out);
  const bool fill_Dg4Dp2 = nonnull(Dg4Dp2_out);
  const bool fill_Dg4Dp4 = nonnull(Dg4Dp4_out);
  const bool fill_Dg6Dp2 = nonnull(Dg6Dp2_out);
  const bool fill_Dg6Dp4 = nonnull(Dg6Dp4_out);

  if (printDebug_) {
    std::cout << "DEBUG: In evalModel: f=" << fill_f
              << ",W=" << fill_W
              << ",W_prec=" << fill_W_prec
              << ",g4=" << fill_g4
              << ",Dg4Dx=" << fill_Dg4Dx
              << ",Dg6Dx=" << fill_Dg6Dx
              << ",DfDp2=" << fill_DfDp2
              << ",DfDp4=" << fill_DfDp4
              << ",Dg4Dp2=" << fill_Dg4Dp2
              << ",Dg4Dp4=" << fill_Dg4Dp4
              << ",Dg6Dp2=" << fill_Dg6Dp2
              << ",Dg6Dp4=" << fill_Dg6Dp4
              << std::endl;
  }

  using tpetra_op = Tpetra::Operator<Scalar,LO,GO,Node>;
  using tpetra_extract = ::Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node>;

  // Create ghosted objects
  if (is_null(uPtr_))
    uPtr_ = Teuchos::rcp(new tpetra_vec(xGhostedMap_));

  uPtr_->doImport(*(tpetra_extract::getConstTpetraVector(inArgs.get_x())), *importer_, Tpetra::REPLACE);

  if (is_null(xPtr_)) {
    xPtr_ = Teuchos::rcp(new tpetra_vec(xGhostedMap_));
    xPtr_->doImport(*nodeCoordinates_, *importer_, Tpetra::INSERT);
  }

  // Sizes for functors
  int myRank = comm_->getRank();
  std::size_t numMyElements = xGhostedMap_->getLocalNumElements()-1;

  // Get parameters, default is from nominal values
  auto k_tpetra = tpetra_extract::getConstTpetraMultiVector(nominalValues_.get_p(2));
  if (nonnull(inArgs.get_p(2)))
    k_tpetra = tpetra_extract::getConstTpetraMultiVector(inArgs.get_p(2));

  Scalar k_val = (k_tpetra->getLocalViewHost(Tpetra::Access::ReadOnly))(0,0);

  auto p4_tpetra = tpetra_extract::getConstTpetraMultiVector(nominalValues_.get_p(4));
  if (nonnull(inArgs.get_p(4)))
    p4_tpetra = tpetra_extract::getConstTpetraMultiVector(inArgs.get_p(4));

  Scalar p4_val = (p4_tpetra->getLocalViewHost(Tpetra::Access::ReadOnly))(0,0);

  if (printDebug_) {
    if (nonnull(inArgs.get_p(2)))
      std::cout << "*** p2: k (NOT nominal)=" << k_val << std::endl;
    else
      std::cout << "*** p2: k (IS nominal)=" << k_val << std::endl;
    if (nonnull(inArgs.get_p(4)))
      std::cout << "*** p4: T_left (NOT nominal)=" << p4_val << std::endl;
    else
      std::cout << "*** p4: T_left (IS nominal)=" << p4_val << std::endl;
  }

  // Residual fill
  if (fill_f) {
    Teuchos::TimeMonitor timer(*residTimer_);
    Teuchos::RCP<tpetra_vec> f = tpetra_extract::getTpetraVector(f_out);
    f->putScalar(0.0);
    ResidualEvaluatorFunctor<tpetra_vec> functor(*f, *xPtr_, *uPtr_, myRank, k_val, p4_val);
    Kokkos::parallel_for("residual evaluation", numMyElements, functor);
    NOX::DeviceSpace().fence();
  }

  // Jacobian fill
  if (fill_W) {
    Teuchos::TimeMonitor timer(*jacTimer_);
    Teuchos::RCP<tpetra_op> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
    Teuchos::RCP<tpetra_matrix> J = Teuchos::rcp_dynamic_cast<tpetra_matrix>(W_tpetra);
    TEUCHOS_ASSERT(nonnull(J));
    J->resumeFill();
    J->setAllToScalar(0.0);
    JacobianEvaluatorFunctor<tpetra_vec, tpetra_matrix> functor(*J, *xPtr_, *uPtr_, myRank, k_val);
    Kokkos::parallel_for("jacobian evaluation", numMyElements, functor);
    NOX::DeviceSpace().fence();
    J->fillComplete();
  }

  // Preconditioner fill
  if (fill_W_prec) {
    Teuchos::RCP<tpetra_matrix> M_inv;
    Teuchos::RCP<tpetra_op> M_tpetra = tpetra_extract::getTpetraOperator(W_prec_out->getNonconstRightPrecOp());
    M_inv = Teuchos::rcp_dynamic_cast<tpetra_matrix>(M_tpetra);
    TEUCHOS_ASSERT(nonnull(M_inv));
    M_inv->resumeFill();
    if (is_null(J_diagonal_))
      J_diagonal_ = Teuchos::rcp(new tpetra_vec(xOwnedMap_));

    M_inv->setAllToScalar(0.0);
    J_diagonal_->putScalar(0.0);
    PreconditionerEvaluatorFunctor<tpetra_vec, tpetra_matrix> functor(*M_inv, *xPtr_, *uPtr_, myRank, k_val);
    Kokkos::parallel_for("prec evaluation", numMyElements, functor);
    NOX::DeviceSpace().fence();

    // Invert the Jacobian diagonal for the preconditioner
    // For some reason the matrix must be fill complete before calling rightScale
    M_inv->fillComplete();
    tpetra_vec& diag = *J_diagonal_;
    M_inv->getLocalDiagCopy(diag);
    diag.reciprocal(diag);
    M_inv->rightScale(diag);
    M_inv->rightScale(diag);
  }

  // Fill Responses. These are so small, we will do it on host. Don't
  // waste time with parallel dispatch.
  auto x = tpetra_extract::getConstTpetraVector(inArgs.get_x())->getLocalViewHost(Tpetra::Access::ReadOnly);

  if (fill_g4) {
    // g4 is locally replicated.
    auto g4_tpetra = tpetra_extract::getTpetraMultiVector(g4_out);

    // g4 = T(Zmax) - 2.0
    Scalar T_right = x(x.extent(0)-1,0);
    Teuchos::broadcast(*comm_,comm_->getSize()-1,&T_right);
    auto g4_host = g4_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
    g4_host(0,0) = T_right - 2.0;
    if (printDebug_)
      std::cout << "evalModel: g(4)= T_right - 2.0 =" << g4_host(0,0) << " T_right=" << T_right << std::endl;
  }

  if (fill_g6) {
    Scalar T_left = x(0,0);
    Teuchos::broadcast(*comm_,0,&T_left);

    Scalar T_right = x(x.extent(0)-1,0);
    Teuchos::broadcast(*comm_,comm_->getSize()-1,&T_right);

    // g6 = 2* T(Zmin) - T(Zmax)
    // g6 is locally replicated.
    auto g6_tpetra = tpetra_extract::getTpetraMultiVector(g6_out);
    auto g6_host = g6_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
    g6_host(0,0) = 2.0 * T_left - T_right;
    if (printDebug_)
      std::cout << "evalModel: g(6)= 2 * T_left - T_right =" << g6_host(0,0) << " T_left=" << T_left << " T_right="
                <<  T_right << std::endl;
  }

  if (fill_Dg4Dx) {
    auto Dg4Dx_tpetra = tpetra_extract::getTpetraMultiVector(Dg4Dx_out);
    Dg4Dx_tpetra->putScalar(0.0);

    // Right most value
    if (comm_->getRank() == (comm_->getSize()-1)) {
      auto Dg4Dx_host = Dg4Dx_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
      Dg4Dx_host(Dg4Dx_host.extent(0)-1,0) = 1.0;
    }
  }

  if (fill_Dg6Dx) {
    auto Dg6Dx_tpetra = tpetra_extract::getTpetraMultiVector(Dg6Dx_out);
    Dg6Dx_tpetra->putScalar(0.0);

    // Left most value
    if (comm_->getRank() == 0) {
      auto Dg6Dx_host = Dg6Dx_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
      Dg6Dx_host(0,0) = 2.0;
    }
    // Right most value
    if (comm_->getRank() == (comm_->getSize()-1)) {
      auto Dg6Dx_host = Dg6Dx_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
      Dg6Dx_host(Dg6Dx_host.extent(0)-1,0) = -1.0;
    }
  }

  if (fill_DfDp2) {
    auto DfDp2_tpetra = tpetra_extract::getTpetraMultiVector(DfDp2_out);
    DfDp2_tpetra->putScalar(0.0);

    DfDp2EvaluatorFunctor<NOX::TMultiVector> functor(*DfDp2_tpetra, *xPtr_, *uPtr_, myRank, k_val);
    Kokkos::parallel_for("DfDp2 evaluation", numMyElements, functor);
    NOX::DeviceSpace().fence();
  }

  if (fill_DfDp4) {
    auto DfDp4_tpetra = tpetra_extract::getTpetraMultiVector(DfDp4_out);
    DfDp4_tpetra->putScalar(0.0);

    // Dirichlet BC on left is the equation:
    // f(0) = T_left - p(4)
    if (comm_->getRank() == 0) {
      auto DfDp4_host = DfDp4_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
      DfDp4_host(0,0) = -1.0;
    }
  }

  if (fill_Dg4Dp2) {
    auto Dg4Dp2_tpetra = tpetra_extract::getTpetraMultiVector(Dg4Dp2_out);
    Dg4Dp2_tpetra->putScalar(0.0);
  }

  if (fill_Dg4Dp4) {
    auto Dg4Dp4_tpetra = tpetra_extract::getTpetraMultiVector(Dg4Dp4_out);
    Dg4Dp4_tpetra->putScalar(0.0);
  }

  if (fill_Dg6Dp2) {
    auto Dg6Dp2_tpetra = tpetra_extract::getTpetraMultiVector(Dg6Dp2_out);
    Dg6Dp2_tpetra->putScalar(0.0);
  }

  if (fill_Dg6Dp4) {
    auto Dg6Dp4_tpetra = tpetra_extract::getTpetraMultiVector(Dg6Dp4_out);
    Dg6Dp4_tpetra->putScalar(0.0);
  }

}

template<class Scalar, class LO, class GO, class Node>
std::pair<int,int>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_p_index(const std::string& p_name) const
{
  const auto search = p_name_to_index_.find(p_name);
  TEUCHOS_ASSERT(search != p_name_to_index_.end());
  return search->second;
}

template<class Scalar, class LO, class GO, class Node>
std::pair<int,int>
EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::get_g_index(const std::string& g_name) const
{
  const auto search = g_name_to_index_.find(g_name);
  TEUCHOS_ASSERT(search != g_name_to_index_.end());
  return search->second;
}

template<class Scalar, class LO, class GO, class Node>
void EvaluatorTpetra1DFEM<Scalar, LO, GO, Node>::
reportFinalPoint (const ::Thyra::ModelEvaluatorBase::InArgs<Scalar> &finalPoint, const bool wasSolved)
{}

#endif
