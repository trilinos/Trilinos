
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
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
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using Thyra::MultiVectorBase;
using Thyra::createMember;
using Thyra::createMembers;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
using Thyra::ProductVectorSpaceBase;
using Thyra::DefaultProductVectorSpace;
using Thyra::productVectorSpace;
using Thyra::ConstDetachedVectorView;
using Thyra::DetachedVectorView;
using Thyra::ConstDetachedSpmdVectorView;
using Thyra::DetachedSpmdVectorView;
typedef Thyra::Ordinal Ordinal;


const int g_localDim = 4; // ToDo: Make variable!


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


//
// Test the default product vector space constructor
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, defaultConstruct,
  Scalar )
{

  RCP<const DefaultProductVectorSpace<Scalar> > vs = productVectorSpace<Scalar>();
  TEST_EQUALITY(vs->numBlocks(), -1);
  TEST_EQUALITY(vs->dim(), as<Ordinal>(-1));
  out << "vs = " << *vs;
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  defaultConstruct )


//
// Test dynmaic cast from a VectorSpaceBase object to a ProductVectorSpaceBase
// object
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, dynamicCast,
  Scalar )
{
  //using Thyra::productVectorSpace;
  RCP<const VectorSpaceBase<Scalar> > vs =  productVectorSpace<Scalar>();
  RCP<const ProductVectorSpaceBase<Scalar> > pvs = productVectorSpaceBase(vs);
  TEST_INEQUALITY(pvs, null);
  TEST_EQUALITY(pvs->numBlocks(), -1);
  TEST_EQUALITY(pvs->dim(), as<Ordinal>(-1));
  out << "vs = " << *vs;
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  dynamicCast )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, dynamicCastNonconst,
  Scalar )
{
  //using Thyra::nonconstProductVectorSpaceBase;
  RCP<VectorSpaceBase<Scalar> > vs =  productVectorSpace<Scalar>();
  RCP<ProductVectorSpaceBase<Scalar> > pvs = nonconstProductVectorSpaceBase(vs);
  TEST_INEQUALITY(pvs, null);
  TEST_EQUALITY(pvs->numBlocks(), -1);
  TEST_EQUALITY(pvs->dim(), as<Ordinal>(-1));
  out << "vs = " << *vs;
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  dynamicCastNonconst )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, dynamicCast_fail,
  Scalar )
{
  //using Thyra::nonconstProductVectorSpaceBase;
  RCP<const VectorSpaceBase<Scalar> > vs =  defaultSpmdVectorSpace<Scalar>(g_localDim);
  TEST_THROW(productVectorSpaceBase(vs), Teuchos::m_bad_cast);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  dynamicCast_fail )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, dynamicCastNonconst_fail,
  Scalar )
{
  //using Thyra::nonconstProductVectorSpaceBase;
  RCP<VectorSpaceBase<Scalar> > vs =  defaultSpmdVectorSpace<Scalar>(g_localDim);
  TEST_THROW(nonconstProductVectorSpaceBase(vs), Teuchos::m_bad_cast);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  dynamicCastNonconst_fail )


} // namespace
