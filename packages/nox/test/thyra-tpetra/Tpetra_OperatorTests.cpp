// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "NOX_Thyra.H"
#include "NOX_Thyra_Vector.H"
#include "NOX_Thyra_MultiVector.H"

// Typedefs

typedef Tpetra::Vector<>::scalar_type Scalar;
typedef Tpetra::Vector<>::local_ordinal_type LO;
typedef Tpetra::Vector<>::global_ordinal_type GO;
typedef Tpetra::Vector<>::node_type Node;

typedef Teuchos::ScalarTraits<Scalar> ST;
typedef Tpetra::Map<LO,GO,Node> Map;
typedef Tpetra::Vector<Scalar,LO,GO,Node> TV;
typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TMV;
typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> CRSM;
typedef Thyra::VectorSpaceBase<Scalar> TVSB;
typedef Thyra::VectorBase<Scalar> TVB;
typedef Thyra::MultiVectorBase<Scalar> TMVB;
typedef Thyra::LinearOpBase<Scalar> TLOB;
typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node> TOVE;
typedef NOX::Thyra::Vector NTV;
typedef NOX::Thyra::MultiVector NTMV;
typedef typename TMV::mag_type mag_type;

#if defined HAVE_TPETRACORE_CUDA
#define NUM_LOCAL 10000
#else
#define NUM_LOCAL 100
#endif
const std::size_t numLocalElements = NUM_LOCAL;

// Model evaluator for returning precomputed Jacobian operator

class JacobianEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:
  JacobianEvaluator(const Teuchos::RCP<TLOB>& lop) :
    lop_(lop)
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(MEB::IN_ARG_x);
    inArgs_ = inArgs;

    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(MEB::OUT_ARG_f);
    outArgs.setSupports(MEB::OUT_ARG_W_op);
    outArgs_ = outArgs;

    Teuchos::RCP<TVB> x0 = Thyra::createMember(lop_->domain());
    Thyra::assign(x0.ptr(), ST::zero());
    nominalValues_ = inArgs;
    nominalValues_.set_x(x0);
  }

  Teuchos::RCP<const TVSB> get_x_space() const { return lop_->domain(); }
  Teuchos::RCP<const TVSB> get_f_space() const { return lop_->range(); }
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const { return nominalValues_; }
  Teuchos::RCP<TLOB> create_W_op() const { return lop_; }
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const { return inArgs_; }

private:
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const { return outArgs_; }
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const
  { if (nonnull(outArgs.get_W_op())) outArgs.get_W_op() = lop_; }

  Teuchos::RCP<TLOB> lop_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
};

// Tpetra matrix creation

Teuchos::RCP<CRSM> createTpetraOp(const Teuchos::RCP<const Map>& map)
{
  Teuchos::RCP<CRSM> mat = Teuchos::rcp(new CRSM(map, 2));

  Scalar posOne = ST::one();
  Scalar negOne = -1.0*posOne;
  for (std::size_t i = 0; i < map->getLocalNumElements(); ++i) {
    const GO row = map->getGlobalElement(static_cast<LO>(i));
    if (row == 0) {
      mat->insertGlobalValues(row,
                              Teuchos::tuple(row+1),
                              Teuchos::tuple(posOne));
    } else if (row == static_cast<GO>(map->getGlobalNumElements()-1)) {
      mat->insertGlobalValues(row,
                              Teuchos::tuple(row-1),
                              Teuchos::tuple(negOne));
    } else {
      mat->insertGlobalValues(row,
                              Teuchos::tuple(row-1, row+1),
                              Teuchos::tuple(negOne, posOne));
    }
  }

  mat->fillComplete();
  return mat;
}

// Routines for checking solution

void checkMultiVectors(const Teuchos::RCP<TMV>& a,
                       const Teuchos::RCP<TMV>& b,
                       const Teuchos::ArrayView<mag_type>& expectedNorms,
                       Teuchos::FancyOStream& out,
                       bool& success)
{
  TEUCHOS_ASSERT(a->getNumVectors() == b->getNumVectors());
  TEUCHOS_ASSERT(static_cast<size_t>(a->getNumVectors()) == static_cast<size_t>(expectedNorms.size()));
  Scalar one = ST::one();
  b->update(-1.0*one, *a, one);
  std::vector<ST::magnitudeType> norms(expectedNorms.size());
  b->norm2(Teuchos::arrayViewFromVector(norms));

  ST::magnitudeType tol = 1.0e-14;
  for (auto norm = norms.begin(); norm != norms.end(); ++norm) {
    TEUCHOS_TEST_EQUALITY(*norm < tol, true, out, success);
  }
  a->norm2(Teuchos::arrayViewFromVector(norms));
  auto normIter = norms.begin();
  auto expNormIter = expectedNorms.begin();
  for (; normIter != norms.end(); ++normIter, ++expNormIter) {
    TEUCHOS_TEST_FLOATING_EQUALITY(*normIter, *expNormIter, tol, out, success);
  }
}

// Unit tests

TEUCHOS_UNIT_TEST(Tpetra_OpTests, VecNoTrans)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create vector objects
  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Teuchos::RCP<TVB> z = Thyra::createMember(space);

  Scalar one = ST::one();
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);

  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);
  Teuchos::RCP<TV> z_tpetra = TOVE::getTpetraVector(z);

  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y));
  Teuchos::RCP<NTV> z_nox = Teuchos::rcp(new NTV(z));

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);
  TEUCHOS_TEST_EQUALITY(Thyra::opSupported(*op, Thyra::NOTRANS), true, out, success);

  // Apply op with Tpetra directly
  op_tpetra->apply(*z_tpetra, *x_tpetra, Teuchos::NO_TRANS);

  // Apply op through NOX Group interface
  Teuchos::RCP<JacobianEvaluator> model = Teuchos::rcp(new JacobianEvaluator(op));
  // Need to pass a lows factory to the group to call computeJacobian
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsf = builder.createLinearSolveStrategy("");
  Teuchos::RCP<NOX::Thyra::Group> group =
    Teuchos::rcp(new NOX::Thyra::Group(*z_nox, model, op, lowsf, Teuchos::null, Teuchos::null, Teuchos::null));
  group->computeJacobian();
  group->applyJacobian(*z_nox, *y_nox);

  // Check for correct results
  mag_type val = static_cast<mag_type>(ST::squareroot(static_cast<Scalar>(2.0)));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);

}

TEUCHOS_UNIT_TEST(Tpetra_OpTests, MultiVecNoTrans)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create vector objects
  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> z = Thyra::createMembers(space, numCols);

  Scalar one = ST::one();
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);

  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);
  Teuchos::RCP<TMV> z_tpetra = TOVE::getTpetraMultiVector(z);

  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y));
  Teuchos::RCP<NTMV> z_nox = Teuchos::rcp(new NTMV(z));

  auto scales = Teuchos::tuple(one, 2.0*one, 3.0*one, 4.0*one, 5.0*one);
  x_tpetra->scale(scales);
  y_tpetra->scale(scales);
  z_tpetra->scale(scales);

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);
  TEUCHOS_TEST_EQUALITY(Thyra::opSupported(*op, Thyra::NOTRANS), true, out, success);

  // Apply op with Tpetra directly
  op_tpetra->apply(*z_tpetra, *x_tpetra, Teuchos::NO_TRANS);

  // Apply op through NOX Group interface
  Teuchos::RCP<JacobianEvaluator> model = Teuchos::rcp(new JacobianEvaluator(op));
  // Need to pass a lows factory to the group to call computeJacobian
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsf = builder.createLinearSolveStrategy("");
  Teuchos::RCP<NOX::Thyra::Group> group =
    Teuchos::rcp(new NOX::Thyra::Group(dynamic_cast<NTV&>((*z_nox)[0]), model, op, lowsf, Teuchos::null, Teuchos::null, Teuchos::null));
  group->computeJacobian();
  group->applyJacobianMultiVector(*z_nox, *y_nox);

  // Check for correct results
  mag_type val = static_cast<mag_type>(ST::squareroot(static_cast<Scalar>(2.0)));
  auto ans = Teuchos::tuple(val, 2.0*val, 3.0*val, 4.0*val, 5.0*val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Tpetra_OpTests, VecTrans)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create vector objects
  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Teuchos::RCP<TVB> z = Thyra::createMember(space);

  Scalar one = ST::one();
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);

  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);
  Teuchos::RCP<TV> z_tpetra = TOVE::getTpetraVector(z);

  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y));
  Teuchos::RCP<NTV> z_nox = Teuchos::rcp(new NTV(z));

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);
  TEUCHOS_TEST_EQUALITY(Thyra::opSupported(*op, Thyra::TRANS), true, out, success);

  // Apply op with Tpetra directly
  op_tpetra->apply(*z_tpetra, *x_tpetra, Teuchos::TRANS);

  // Apply op through NOX Group interface
  Teuchos::RCP<JacobianEvaluator> model = Teuchos::rcp(new JacobianEvaluator(op));
  // Need to pass a lows factory to the group to call computeJacobian
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsf = builder.createLinearSolveStrategy("");
  Teuchos::RCP<NOX::Thyra::Group> group =
    Teuchos::rcp(new NOX::Thyra::Group(*z_nox, model, op, lowsf, Teuchos::null, Teuchos::null, Teuchos::null));
  group->computeJacobian();
  // Need this to get the group to call updateLOWS, or else shared_jacobian is unitialized
  (void) group->getJacobian();
  group->applyJacobianTranspose(*z_nox, *y_nox);

  // Check for correct results
  mag_type val = static_cast<mag_type>(ST::squareroot(static_cast<Scalar>(2.0)));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Tpetra_OpTests, MultiVecTrans)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  // Create vector objects
  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> z = Thyra::createMembers(space, numCols);

  Scalar one = ST::one();
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);

  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);
  Teuchos::RCP<TMV> z_tpetra = TOVE::getTpetraMultiVector(z);

  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y));
  Teuchos::RCP<NTMV> z_nox = Teuchos::rcp(new NTMV(z));

  auto scales = Teuchos::tuple(one, 2.0*one, 3.0*one, 4.0*one, 5.0*one);
  x_tpetra->scale(scales);
  y_tpetra->scale(scales);
  z_tpetra->scale(scales);

  // Create operator objects
  Teuchos::RCP<CRSM> op_tpetra = createTpetraOp(map);
  Teuchos::RCP<TLOB> op = Thyra::tpetraLinearOp<Scalar,LO,GO,Node>(space, space, op_tpetra);
  TEUCHOS_TEST_EQUALITY(Thyra::opSupported(*op, Thyra::TRANS), true, out, success);

  // Apply op with Tpetra directly
  op_tpetra->apply(*z_tpetra, *x_tpetra, Teuchos::TRANS);

  // Apply op through NOX Group interface
  Teuchos::RCP<JacobianEvaluator> model = Teuchos::rcp(new JacobianEvaluator(op));
  // Need to pass a lows factory to the group to call computeJacobian even though it isn't used
  Stratimikos::DefaultLinearSolverBuilder builder;
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::parameterList();
  p->set("Linear Solver Type", "Belos");
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsf = builder.createLinearSolveStrategy("");
  Teuchos::RCP<NOX::Thyra::Group> group =
    Teuchos::rcp(new NOX::Thyra::Group(dynamic_cast<NTV&>((*z_nox)[0]), model, op, lowsf, Teuchos::null, Teuchos::null, Teuchos::null));
  group->computeJacobian();
  // Need this to get the group to call updateLOWS, or else shared_jacobian is unitialized
  (void) group->getJacobian();
  group->applyJacobianTransposeMultiVector(*z_nox, *y_nox);

  // Check for correct results
  mag_type val = static_cast<mag_type>(ST::squareroot(static_cast<Scalar>(2.0)));
  auto ans = Teuchos::tuple(val, 2.0*val, 3.0*val, 4.0*val, 5.0*val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}
