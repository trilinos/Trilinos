/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_UnitTestHelpers.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;


const int g_localDim = 4; // ToDo: Make variable!


template<class Scalar>
RCP<VectorSpaceBase<Scalar> > 
createSpmdVectorSpace(const Ordinal localDim)
{
  return defaultSpmdVectorSpace<Scalar>(
    Teuchos::DefaultComm<Ordinal>::getComm(),
    localDim, -1 );
}


template<class Scalar>
RCP<VectorSpaceBase<Scalar> > 
createProductVectorSpace(const Ordinal localDim, const int numBlocks)
{
  return productVectorSpace<Scalar>(
    createSpmdVectorSpace<Scalar>(localDim),
    numBlocks);
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
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
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
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
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
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  dynamicCast_fail )


//
// Make sure that DefaultProductVectorSpace::createMember() returns a
// DefaultProductVector object.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, DefaultProductVector,
  Scalar )
{
  using Teuchos::rcp_dynamic_cast;
  const RCP<const VectorSpaceBase<Scalar> >
    vs = defaultSpmdVectorSpace<Scalar>(g_localDim),
    pvs = productVectorSpace<Scalar>(tuple(vs)());
  const RCP<VectorBase<Scalar> > v1 = createMember(pvs);
  const RCP<DefaultProductVector<Scalar> > pv1 =
    rcp_dynamic_cast<DefaultProductVector<Scalar> >(v1, true);
  TEST_EQUALITY(pvs.get(), pv1->space().get());
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  DefaultProductVector )


//
// Make sure that DefaultProductVectorSpace::createMembers() returns a
// DefaultProductMultiVector object.
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, DefaultProductMultiVector,
  Scalar )
{
  using Teuchos::as;
  using Teuchos::rcp_dynamic_cast;
  const RCP<const VectorSpaceBase<Scalar> >
    vs = defaultSpmdVectorSpace<Scalar>(g_localDim),
    pvs = productVectorSpace<Scalar>(tuple(vs)());
  const int numCols = 2;
  const RCP<MultiVectorBase<Scalar> > mv1 = createMembers(pvs, numCols);
  const RCP<DefaultProductMultiVector<Scalar> > pmv1 =
    rcp_dynamic_cast<DefaultProductMultiVector<Scalar> >(mv1, true);
  TEST_EQUALITY(pvs.get(), pmv1->range().get());
  TEST_EQUALITY(as<int>(pmv1->domain()->dim()), numCols);
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  DefaultProductMultiVector )


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
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  singleBlockCompatibility )


//
// Test that a product vector space with more than one block is incompatible
// with a vector space of the same size always!
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVectorSpace, twoBlockIncompatible,
  Scalar )
{
  using Teuchos::describe;

  const RCP<const VectorSpaceBase<Scalar> >
    vs1 = defaultSpmdVectorSpace<Scalar>(g_localDim),
    vs2 = defaultSpmdVectorSpace<Scalar>(2*g_localDim),
    pvs = productVectorSpace<Scalar>(tuple(vs1,vs1)());
  
  out << "vs1=" << describe(*vs1);
  out << "vs2=" << describe(*vs2);
  out << "pvs=" << describe(*pvs);

  TEST_EQUALITY(pvs->dim(), vs2->dim());
  TEST_ASSERT(!pvs->isCompatible(*vs1));
  TEST_ASSERT(!vs1->isCompatible(*pvs));
  TEST_ASSERT(!pvs->isCompatible(*vs2));
  TEST_ASSERT(!vs2->isCompatible(*pvs));

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVectorSpace,
  twoBlockIncompatible )


//
// DefaultProductVector
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVector,
  castOrCreateNonconstProductVectorBase_cast, Scalar )
{

  const int numBlocks = 3;

  const RCP<const VectorSpaceBase<Scalar> > vs =
    createProductVectorSpace<Scalar>(g_localDim, numBlocks);

  const RCP<VectorBase<Scalar> > v = createMember(vs);

  const RCP<ProductVectorBase<Scalar> > prod_v =
    castOrCreateNonconstProductVectorBase(v);

  out << "prod_v = " << *prod_v;

  TEST_EQUALITY_CONST(prod_v->productSpace()->numBlocks(), numBlocks);

  TEST_EQUALITY(
    dynamic_cast<void*>(v.getRawPtr()),
    dynamic_cast<void*>(prod_v.getRawPtr()));

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVector,
  castOrCreateNonconstProductVectorBase_cast )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVector,
  castOrCreateNonconstProductVectorBase_create, Scalar )
{

  const RCP<const VectorSpaceBase<Scalar> > vs =
    createSpmdVectorSpace<Scalar>(g_localDim);

  const RCP<VectorBase<Scalar> > v = createMember(vs);

  const RCP<ProductVectorBase<Scalar> > prod_v =
    castOrCreateNonconstProductVectorBase(v);

  out << "prod_v = " << *prod_v;

  TEST_EQUALITY_CONST(prod_v->productSpace()->numBlocks(), 1);

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVector,
  castOrCreateNonconstProductVectorBase_create )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVector,
  castOrCreateProductVectorBase_cast, Scalar )
{

  const int numBlocks = 3;

  const RCP<const VectorSpaceBase<Scalar> > vs =
    createProductVectorSpace<Scalar>(g_localDim, numBlocks);

  const RCP<const VectorBase<Scalar> > v = createMember(vs);

  const RCP<const ProductVectorBase<Scalar> > prod_v =
    castOrCreateProductVectorBase(v);

  out << "prod_v = " << *prod_v;

  TEST_EQUALITY_CONST(prod_v->productSpace()->numBlocks(), numBlocks);

  TEST_EQUALITY(
    dynamic_cast<const void*>(v.getRawPtr()),
    dynamic_cast<const void*>(prod_v.getRawPtr()));

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVector,
  castOrCreateProductVectorBase_cast )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultProductVector,
  castOrCreateProductVectorBase_create, Scalar )
{

  const RCP<const VectorSpaceBase<Scalar> > vs =
    createSpmdVectorSpace<Scalar>(g_localDim);

  const RCP<const VectorBase<Scalar> > v = createMember(vs);

  const RCP<const ProductVectorBase<Scalar> > prod_v =
    castOrCreateProductVectorBase(v);

  out << "prod_v = " << *prod_v;

  TEST_EQUALITY_CONST(prod_v->productSpace()->numBlocks(), 1);

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultProductVector,
  castOrCreateProductVectorBase_create )




} // namespace Thyra
