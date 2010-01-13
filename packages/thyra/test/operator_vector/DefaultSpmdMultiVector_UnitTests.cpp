
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


// Make sure the currect public member access protections are in place
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdMultiVector, memberAccess,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<const VectorSpaceBase<Scalar> > vs = createSpmdVectorSpace<Scalar>(g_localDim);
  RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, g_numCols);
  Thyra::assign<Scalar>(mv.ptr(), ST::zero());
  RCP<DefaultSpmdMultiVector<Scalar> > spmdMv = 
    rcp_dynamic_cast<DefaultSpmdMultiVector<Scalar> >(mv);
  TEST_ASSERT(nonnull(spmdMv->spmdSpace()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdMultiVector,
  memberAccess )


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


} // namespace
