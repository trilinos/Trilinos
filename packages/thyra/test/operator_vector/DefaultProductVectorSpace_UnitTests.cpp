
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_VectorStdOps.hpp"
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
using Thyra::V_S;
using Thyra::V_V;
using Thyra::assign;
using Thyra::sum;
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
// Test dynamic casting using the inline double helper functions. In this
// case, we don't need to employ a using declaration or explicit template
// arguments.
//

TEUCHOS_UNIT_TEST( DefaultProductVectorSpace, dynamicCast_double)
{

  const RCP<DefaultProductVectorSpace<double> > dpvs =  productVectorSpace<double>();
  const RCP<const DefaultProductVectorSpace<double> > cdpvs = dpvs;
  const RCP<VectorSpaceBase<double> > vs =  dpvs;
  const RCP<const VectorSpaceBase<double> > cvs =  vs;
    
  // VS -> PVSB
  RCP<ProductVectorSpaceBase<double> > pvs1 = nonconstProductVectorSpaceBase(vs);
  TEST_INEQUALITY(pvs1, null);
  out << "pvs1 = " << *pvs1;

  // DPVS -> PVSB
  RCP<ProductVectorSpaceBase<double> > pvs2 = nonconstProductVectorSpaceBase(dpvs);
  TEST_INEQUALITY(pvs2, null);
  out << "pvs2 = " << *pvs2;
  
  // VS -> const PVSB
  RCP<const ProductVectorSpaceBase<double> > pvs3 = productVectorSpaceBase(vs);
  TEST_INEQUALITY(pvs3, null);
  out << "pvs3 = " << *pvs3;

  // DPVS -> const PVSB
  RCP<const ProductVectorSpaceBase<double> > pvs4 = productVectorSpaceBase(dpvs);
  TEST_INEQUALITY(pvs4, null);
  out << "pvs4 = " << *pvs4;
  
  // const VS -> const PVSB
  RCP<const ProductVectorSpaceBase<double> > cpvs5 = productVectorSpaceBase(cvs);
  TEST_INEQUALITY(cpvs5, null);
  out << "cpvs5 = " << *cpvs5;

  // const DPVS -> const PVSB
  RCP<const ProductVectorSpaceBase<double> > cpvs6 = productVectorSpaceBase(cdpvs);
  TEST_INEQUALITY(cpvs6, null);
  out << "cpvs6 = " << *cpvs6;

}


//
// Test dynmaic cast from a VectorSpaceBase object to a ProductVectorSpaceBase
// object using the templated function.  In this case, we have to inject the
// names of the non-member conversion functions into the local namespace and
// we have to use the explicit template argument so that the templated
// function will be selected over the non-templated double version of the
// function.
//
  
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, dynamicCast,
  Scalar )
{

  using Thyra::nonconstProductVectorSpaceBase;
  using Thyra::productVectorSpaceBase;

  const RCP<DefaultProductVectorSpace<Scalar> > dpvs =  productVectorSpace<Scalar>();
  const RCP<const DefaultProductVectorSpace<Scalar> > cdpvs = dpvs;
  const RCP<VectorSpaceBase<Scalar> > vs =  dpvs;
  const RCP<const VectorSpaceBase<Scalar> > cvs =  vs;
  
  // VS -> PVSB (exact match so template arg is *not* needed)
  RCP<ProductVectorSpaceBase<Scalar> > 
    pvs1 = nonconstProductVectorSpaceBase(vs);
  TEST_INEQUALITY(pvs1, null);
  out << "pvs1 = " << *pvs1;

  // DPVS -> PVSB (base conversion so template arg is needed)
  RCP<ProductVectorSpaceBase<Scalar> >
    pvs2 = nonconstProductVectorSpaceBase<Scalar>(dpvs);
  TEST_INEQUALITY(pvs2, null);
  out << "pvs2 = " << *pvs2;
  
  // VS -> const PVSB (const conversion so template arg is needed) 
  RCP<const ProductVectorSpaceBase<Scalar> >
    pvs3 = productVectorSpaceBase<Scalar>(vs);
  TEST_INEQUALITY(pvs3, null);
  out << "pvs3 = " << *pvs3;

  // DPVS -> const PVSB (base & const conversion so template arg is needed)
  RCP<const ProductVectorSpaceBase<Scalar> >
    pvs4 = productVectorSpaceBase<Scalar>(dpvs);
  TEST_INEQUALITY(pvs4, null);
  out << "pvs4 = " << *pvs4;
  
  // const VS -> const PVSB (exact match so template arg is *not* needed)
  RCP<const ProductVectorSpaceBase<Scalar> >
    cpvs5 = productVectorSpaceBase(cvs);
  TEST_INEQUALITY(cpvs5, null);
  out << "cpvs5 = " << *cpvs5;

  // const DPVS -> const PVSB (base conversion so template arg is needed)
  RCP<const ProductVectorSpaceBase<Scalar> >
    cpvs6 = productVectorSpaceBase<Scalar>(cdpvs);
  TEST_INEQUALITY(cpvs6, null);
  out << "cpvs6 = " << *cpvs6;

}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  dynamicCast )


//
// Make sure that invalid dynamic casts fail and throw exceptions
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, dynamicCast_fail,
  Scalar )
{
  using Thyra::nonconstProductVectorSpaceBase;
  using Thyra::productVectorSpaceBase;
  RCP<VectorSpaceBase<Scalar> > vs = defaultSpmdVectorSpace<Scalar>(g_localDim);
  TEST_THROW(productVectorSpaceBase<Scalar>(vs), Teuchos::m_bad_cast);
  TEST_THROW(nonconstProductVectorSpaceBase<Scalar>(vs), Teuchos::m_bad_cast);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  dynamicCast_fail )


//
// Test that a product vector space with one embedded space is compatiable
// with that single space and visa versa.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, singleBlockCompatibility,
  Scalar )
{

  using Teuchos::describe;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const RCP<const VectorSpaceBase<Scalar> >
    vs = defaultSpmdVectorSpace<Scalar>(g_localDim),
    pvs = productVectorSpace<Scalar>(tuple(vs)());

  out << "vs=" << describe(*vs);
  out << "pvs=" << describe(*pvs);

  TEST_ASSERT(pvs->isCompatible(*vs));
  TEST_ASSERT(vs->isCompatible(*pvs));

  const ScalarMag dim_scalar = vs->dim();

  const RCP<VectorBase<Scalar> >
    v1 = createMember(vs),
    pv1 = createMember(pvs);

  out << "Test that you can copy from single vector to product vector (1) ...\n";
  const Scalar val1 = 2.0;
  V_S(v1.ptr(), val1);
  V_V(pv1.ptr(), *v1);
  TEST_FLOATING_EQUALITY( sum(*pv1), as<Scalar>(dim_scalar * val1),
    as<ScalarMag>(SMT::eps() / dim_scalar * 1e+2) );

  out << "Test that you can copy from product vector (1) to single vector ...\n";
  const Scalar val2 = 3.0;
  V_S(pv1.ptr(), val2);
  V_V(v1.ptr(), *pv1);
  TEST_FLOATING_EQUALITY( sum(*v1), as<Scalar>(dim_scalar * val2),
    as<ScalarMag>(SMT::eps() / dim_scalar * 1e+2) );

  const RCP<MultiVectorBase<Scalar> >
    mv1 = createMembers(vs, 1),
    pmv1 = createMembers(pvs, 1);

  out << "Test that you can copy from single multi-vector to product multi-vector (1) ...\n";
  assign(mv1.ptr(), val1);
  assign(pmv1.ptr(), *mv1);
  TEST_FLOATING_EQUALITY( sum(*pmv1->col(0)), as<Scalar>(dim_scalar * val1),
    as<ScalarMag>(SMT::eps() / dim_scalar * 1e+2) );

  out << "Test that you can copy from product multi-vector (1) to single multi-vector ...\n";
  assign(pmv1.ptr(), val2);
  assign(mv1.ptr(), *pmv1);
  TEST_FLOATING_EQUALITY( sum(*mv1->col(0)), as<Scalar>(dim_scalar * val2),
    as<ScalarMag>(SMT::eps() / dim_scalar * 1e+2) );

}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  singleBlockCompatibility )


} // namespace
