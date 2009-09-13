
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_UnitTestHelpers.hpp"


namespace {


//
// Helper code and declarations
//


bool dumpAll = false;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-all", "no-dump-all", &dumpAll,
    "Dump lots of data" );
}


using Teuchos::as;
using Teuchos::null;
using Teuchos::is_null;
using Teuchos::updateSuccess;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Teuchos::ptrFromRef;
using Teuchos::fancyOStream;
using Teuchos::get_extra_data;
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using Thyra::MultiVectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::createMembers;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
using Thyra::DefaultBlockedLinearOp;
using Thyra::defaultBlockedLinearOp;
using Thyra::ProductVectorSpaceBase;
using Thyra::PhysicallyBlockedLinearOpBase;
using Thyra::LinearOpTester;
typedef Thyra::Ordinal Ordinal;


//
// Unit Tests
//


const int m = 5;
const int n = 3;
const int p = 2;


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedLinearOp, defaultConstruct,
  Scalar )
{

  typedef Scalar S;

  const RCP<PhysicallyBlockedLinearOpBase<S> > M = defaultBlockedLinearOp<S>();

  TEST_ASSERT(is_null(M->range()));
  TEST_ASSERT(is_null(M->domain()));
  TEST_ASSERT(is_null(M->productRange()));
  TEST_ASSERT(is_null(M->productDomain()));
 
  std::ostringstream describe_msg;

  describe_msg << "'";
  M->describe(*fancyOStream(rcpFromRef(describe_msg)), Teuchos::VERB_LOW);
  describe_msg << "'";

  std::ostringstream expected_msg;
  expected_msg
    << "' " << M->Describable::description() << "{"
    << "numRowBlocks="<<0
    << ",numColBlocks="<<0
    << "}\n'";

  TEST_EQUALITY_CONST( describe_msg.str(), expected_msg.str() );
  
}

THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedLinearOp,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedLinearOp, block1x1,
  Scalar )
{

  typedef Scalar S;  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const RCP<const VectorSpaceBase<S> > vsm =  defaultSpmdVectorSpace<S>(m);
  const RCP<const VectorSpaceBase<S> > vsn =  defaultSpmdVectorSpace<S>(n);

  const RCP<MultiVectorBase<S> > A = createMembers(vsm, n);
  Thyra::assign(A.ptr(), as<S>(1.0));

  const RCP<PhysicallyBlockedLinearOpBase<S> > M = defaultBlockedLinearOp<S>();

  TEST_EQUALITY_CONST(M->blockFillIsActive(), false);
  M->beginBlockFill(1, 1);
  TEST_EQUALITY_CONST(M->blockFillIsActive(), true);
  TEST_EQUALITY_CONST(M->acceptsBlock(0, 0), true);
  M->setBlock(0, 0, A);
  M->endBlockFill();
  TEST_EQUALITY_CONST(M->blockFillIsActive(), false);
 
  Thyra::LinearOpTester<S> linearOpTester;
  linearOpTester.set_all_error_tol(1e3*SMT::eps());
  linearOpTester.dump_all(dumpAll);

  updateSuccess(linearOpTester.check(*M, ptrFromRef(out)), success);

  const RCP<VectorBase<S> >  x = createMember<S>(vsn), y = createMember<S>(vsm);
  Thyra::assign(x.ptr(), as<Scalar>(2.0));
  Thyra::apply<S>( *M, Thyra::NOTRANS, *x, y.ptr() );
  TEST_FLOATING_EQUALITY( Thyra::sum(*y), as<Scalar>(m*n*1.0*2.0),
    as<ScalarMag>(SMT::eps() / (n*m)) );

  const RCP<const LinearOpBase<S> > M2 = Thyra::block1x1<S>(A);
  updateSuccess(linearOpTester.compare(*M, *M2, ptrFromRef(out)), success);

  const RCP<const LinearOpBase<S> > M3 = Thyra::nonconstBlock1x1<S>(A);
  updateSuccess(linearOpTester.compare(*M, *M3, ptrFromRef(out)), success);
  
}

THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedLinearOp,
  block1x1 )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultBlockedLinearOp, nestedBlock2x2,
  Scalar )
{

  typedef Scalar S;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  
  using Teuchos::describe;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::block1x2;
  using Thyra::block2x1;
  using Thyra::block2x2;
  using Thyra::adjoint;

  const RCP<const VectorSpaceBase<S> > vsm = defaultSpmdVectorSpace<S>(m);
  const RCP<const VectorSpaceBase<S> > vsn = defaultSpmdVectorSpace<S>(n);
  const RCP<const VectorSpaceBase<S> > vsp = defaultSpmdVectorSpace<S>(p);

  const RCP<MultiVectorBase<S> >
    A = createMembers(vsm, m, "A"),
    B = createMembers(vsn, m, "B"),
    Pm = createMembers(vsm, p, "Pm"),
    Pn = createMembers(vsn, p, "Pn"),
    Q = createMembers(vsp, p, "Q");
  
  assign<S>(A.ptr(), as<S>(1.0));
  assign<S>(B.ptr(), as<S>(2.0));
  assign<S>(Pm.ptr(), as<S>(3.0));
  assign<S>(Pn.ptr(), as<S>(4.0));
  assign<S>(Q.ptr(), as<S>(5.0));
  
  const RCP<const LinearOpBase<S> >
    M00 = block2x2<S>(A, adjoint<S>(B), B, null),
    M01 = block2x1<S>(Pm, Pn),
    M10 = block1x2<S>(adjoint<S>(Pm), adjoint<S>(Pn));

  const RCP<const LinearOpBase<S> >
    M = block2x2<S>( M00, M01, M10, Q, "M" );

  out << "M = " << describe(*M, Teuchos::VERB_MEDIUM);
  out << "M->range() = " << describe(*M->range(), Teuchos::VERB_MEDIUM);
  out << "M->domain() = " << describe(*M->range(), Teuchos::VERB_MEDIUM);

  const RCP<const PhysicallyBlockedLinearOpBase<S> >
    pbM = rcp_dynamic_cast<const PhysicallyBlockedLinearOpBase<S> >(M, true);

  TEST_NOTHROW(
    rcp_dynamic_cast<const ProductVectorSpaceBase<S> >(
      pbM->productRange()->getBlock(0), true)
    );
  TEST_NOTHROW(
    rcp_dynamic_cast<const DefaultSpmdVectorSpace<S> >(
      pbM->productRange()->getBlock(1), true)
    );
  TEST_NOTHROW(
    rcp_dynamic_cast<const ProductVectorSpaceBase<S> >(
      pbM->productRange()->getBlock(0), true)
    );
  TEST_NOTHROW(
    rcp_dynamic_cast<const DefaultSpmdVectorSpace<S> >(
      pbM->productDomain()->getBlock(1), true)
    );
 
  Thyra::LinearOpTester<S> linearOpTester;
  linearOpTester.set_all_error_tol(1e3*SMT::eps());
  linearOpTester.dump_all(dumpAll);

  updateSuccess(linearOpTester.check(*M, ptrFromRef(out)), success);

  const RCP<const LinearOpBase<S> >
    M2 = block2x2<S>( Q, M10, M01, M00, "M2" );

  out << "M2 = " << describe(*M2, Teuchos::VERB_MEDIUM);
  out << "M2->range() = " << describe(*M2->range(), Teuchos::VERB_MEDIUM);
  out << "M2->domain() = " << describe(*M2->range(), Teuchos::VERB_MEDIUM);

  const RCP<const PhysicallyBlockedLinearOpBase<S> >
    pbM2 = rcp_dynamic_cast<const PhysicallyBlockedLinearOpBase<S> >(M2, true);

  TEST_NOTHROW(
    rcp_dynamic_cast<const DefaultSpmdVectorSpace<S> >(
      pbM2->productRange()->getBlock(0), true)
    );
  TEST_NOTHROW(
    rcp_dynamic_cast<const ProductVectorSpaceBase<S> >(
      pbM2->productRange()->getBlock(1), true)
    );
  TEST_NOTHROW(
    rcp_dynamic_cast<const DefaultSpmdVectorSpace<S> >(
      pbM2->productDomain()->getBlock(0), true)
    );
  TEST_NOTHROW(
    rcp_dynamic_cast<const ProductVectorSpaceBase<S> >(
      pbM2->productRange()->getBlock(1), true)
    );
 
  updateSuccess(linearOpTester.check(*M2, ptrFromRef(out)), success);
  
}

THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultBlockedLinearOp,
  nestedBlock2x2 )


} // namespace
