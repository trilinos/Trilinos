
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_MultiVectorTester.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::inOutArg;
using Teuchos::Range1D;
using Thyra::VectorSpaceBase;
using Thyra::MultiVectorBase;
using Thyra::createMembers;
using Thyra::DefaultSpmdMultiVector;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
typedef Thyra::Ordinal Ordinal;


const int g_localDim = 6; // ToDo: Make variable!
const int g_numCols = 5;  // ToDo: Make variable!


template<class Scalar>
RCP<VectorSpaceBase<Scalar> > 
createSpmdVectorSpace(const Teuchos_Ordinal localDim)
{
  return defaultSpmdVectorSpace<Scalar>(
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm(),
    localDim, -1 );
}


//
// Unit Tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, defaultConstruct,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<const MultiVectorBase<Scalar> > mv = createMembers(*vs, g_numCols);
  Teuchos::Array<Scalar> mv_sums(g_numCols);
  Thyra::sums<Scalar>(*mv, mv_sums());
  out << "sums(mv) = " << mv_sums << "\n";
#ifdef THYRA_INITIALIZE_VECS_MULTIVECS_WITH_NANS
  for (int i = 0; i < g_numCols; ++i) {
    TEST_ASSERT(ST::isnaninf(mv_sums[i]));
  }
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, defaultTester,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  Thyra::MultiVectorTester<Scalar> mvTester;
  const bool mvTesterResult = mvTester.checkMultiVector(*vs, inOutArg(out));
  TEST_ASSERT(mvTesterResult);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  defaultTester )


// Make sure the currect public member access protections are in place
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, memberAccess,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, g_numCols);
  Thyra::assign<Scalar>(mv.ptr(), ST::zero());
  RCP<DefaultSpmdMultiVector<Scalar> > spmdMv = 
    rcp_dynamic_cast<DefaultSpmdMultiVector<Scalar> >(mv, true);
  TEST_ASSERT(nonnull(spmdMv->spmdSpace()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  memberAccess )


//
// Make sure that DefaultSpmdMultiVector::subview(...) return a
// DefaultSpmdMultiVector object!
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, constContigSubViewImpl,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, g_numCols);
  Thyra::assign<Scalar>(mv.ptr(), ST::zero());
  RCP<const MultiVectorBase<Scalar> > mvv = mv.getConst()->subView(Range1D(0, 1));
  RCP<const DefaultSpmdMultiVector<Scalar> > spmdMvv = 
    rcp_dynamic_cast<const DefaultSpmdMultiVector<Scalar> >(mvv, true);
  TEST_ASSERT(nonnull(spmdMvv->spmdSpace()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  constContigSubViewImpl )


//
// Make sure that DefaultSpmdMultiVector::subview(...) return a
// DefaultSpmdMultiVector object!
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, nonconstContigSubViewImpl,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, g_numCols);
  Thyra::assign<Scalar>(mv.ptr(), ST::zero());
  RCP<MultiVectorBase<Scalar> > mvv = mv->subView(Range1D(0, 1));
  RCP<DefaultSpmdMultiVector<Scalar> > spmdMvv = 
    rcp_dynamic_cast<DefaultSpmdMultiVector<Scalar> >(mvv, true);
  TEST_ASSERT(nonnull(spmdMvv->spmdSpace()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  nonconstContigSubViewImpl )


//
// Test that a dangling non-const column subviews don't write back their data.
// If the parent object is gone then why write back the data?  This is a
// performance optimization.
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, danglingSubViews,
  Scalar )
{
  using Teuchos::tuple;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(*vs, g_numCols);
  RCP<MultiVectorBase<Scalar> > mvView = mv->subView(tuple<int>(0, 1));
  mv = null;
#ifdef THYRA_DEBUG
  const int startingNumCopyBack = DefaultSpmdMultiVector<Scalar>::numSkipCopyBack;
#endif
  mvView = null; // Should not write back data since parent mv is gone now.
#ifdef THYRA_DEBUG
  TEST_EQUALITY(DefaultSpmdMultiVector<Scalar>::numSkipCopyBack, startingNumCopyBack+1);
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  danglingSubViews )


#ifdef THYRA_DEBUG


//
// Test that invalid column indexes throw the right exception messages.
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, invalidSubviews,
  Scalar )
{
  using Teuchos::tuple;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<const MultiVectorBase<Scalar> > mv = createMembers(*vs, 1);
  TEST_THROW(mv->subView(tuple<int>(0, 1)), std::invalid_argument);
  TEST_THROW(mv->subView(tuple<int>(-1)), std::out_of_range);
  TEST_THROW(mv->subView(tuple<int>(1)), std::out_of_range);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  invalidSubviews )


#endif // THYRA_DEBUG


} // namespace
