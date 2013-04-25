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
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_VectorStdOpsTester.hpp"


//#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP

#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
#  include "RTOpPack_SPMD_apply_op_decl.hpp"
#  include "Thyra_SpmdVectorBase.hpp"
#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::get_extra_data;
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using Thyra::MultiVectorBase;
using Thyra::createMember;
using Thyra::createMembers;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, defaultConstruct,
  Scalar )
{

  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>());
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(-1));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(-1));
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, serialConstruct,
  Scalar )
{

  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(g_localDim));
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  serialConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, serialConstructZeroSize,
  Scalar )
{

  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(0));
  TEST_EQUALITY_CONST(vs->getComm(), null);
  TEST_EQUALITY_CONST(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(vs->localSubDim(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(vs->mapCode(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(vs->dim(), as<Ordinal>(0));

  ECHO(const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs));
  out << "v = " << *v;
  
  TEST_ASSERT(vs->isCompatible(*v->space()));
  TEST_EQUALITY_CONST(v->space()->dim(), as<Ordinal>(0));

  ECHO(const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, 0));
  out << "mv = " << *mv;
  
  TEST_ASSERT(vs->isCompatible(*mv->range()));
  TEST_EQUALITY_CONST(mv->range()->dim(), as<Ordinal>(0));
  TEST_EQUALITY_CONST(mv->domain()->dim(), as<Ordinal>(0));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  serialConstructZeroSize )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, parallelConstruct,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(comm, g_localDim, -1));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank()*g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize()*g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  parallelConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, parallelConstructGlobalDim,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(comm, g_localDim, g_localDim * comm->getSize()));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank()*g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  //TEST_EQUALITY_CONST(vs->isLocallyReplicated(), false);
  TEST_EQUALITY(vs->isLocallyReplicated(), (comm->getSize()==1));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize()*g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  parallelConstructGlobalDim )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, locallyReplicatedParallelConstruct,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(comm, g_localDim, g_localDim));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY_CONST(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY_CONST(vs->isLocallyReplicated(), true);
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  locallyReplicatedParallelConstruct )


/****
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, parallelConstructEmptyProc,
  Scalar )
{

  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  const int procRank = comm->getRank();
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(comm, procRank == 0 ? 0 : g_localDim, -1));
  TEST_EQUALITY(vs->getComm(), comm);

  if (procRank == 0) {
    TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(0));
  }
  else {
    TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  }
  TEST_EQUALITY(vs->dim(), as<Ordinal>((comm->getSize()-1)*g_localDim));

  if (vs->dim()) {

    ECHO(const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs));
    ECHO(V_S(v.ptr(), as<Scalar>(1.0)));
    out << "v = " << *v;
    
    // ToDo: Fix MultiVector to work with empty processors
    //ECHO(const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, 1));
    //ECHO(assign(mv.ptr(), as<Scalar>(1.0)));
    //out << "mv = " << *mv;

  }

  // AGS: Turn on vector space testor for emptyProc case

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;
  Scalar tol = 1.0e-12;
  bool showAllTests=true, dumpAll=true;

  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.warning_tol(ScalarMag(0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);

  Thyra::VectorStdOpsTester<Scalar> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(ScalarMag(0.1)*tol);
  vectorStdOpsTester.error_tol(tol);

  out << "\nTesting the VectorSpaceBase interface of vs ...\n";
  TEUCHOS_TEST_ASSERT(vectorSpaceTester.check(*vs, &out), out, success);

  out << "\nTesting standard vector ops for vs ...\n";
  TEUCHOS_TEST_ASSERT(vectorStdOpsTester.checkStdOps(*vs, &out), out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  parallelConstructEmptyProc )
****/


//#ifndef THYRA_HIDE_DEPRECATED_CODE


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, deprecatedSerialConstruct,
  Scalar )
{

  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(g_localDim));
  TEST_EQUALITY(vs->getComm(), null);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(0));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->mapCode(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(g_localDim));
  ECHO(const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs));
  TEST_ASSERT(vs->isCompatible(*v->space()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  deprecatedSerialConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, deprecatedParallelConstruct,
  Scalar )
{
  ECHO(const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm());
  ECHO(RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, g_localDim, -1));
  TEST_EQUALITY(vs->getComm(), comm);
  TEST_EQUALITY(vs->localOffset(), as<Ordinal>(comm->getRank()*g_localDim));
  TEST_EQUALITY(vs->localSubDim(), as<Ordinal>(g_localDim));
  TEST_EQUALITY(vs->dim(), as<Ordinal>(comm->getSize()*g_localDim));
  ECHO(const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs));
  TEST_ASSERT(vs->isCompatible(*v->space()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  deprecatedParallelConstruct )


//#endif // THYRA_HIDE_DEPRECATED_CODE


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, parallelFullExtract,
  Scalar )
{

  const RCP<const VectorSpaceBase<Scalar> > vs =
    createSpmdVectorSpace<Scalar>(g_localDim);

  const RCP<VectorBase<Scalar> > v = createMember(vs);
#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
    RTOpPack::show_spmd_apply_op_dump = true;
#endif
  {
    out << "\nSetting up v[i] = i, i=0...n-1 ...\n";
    const DetachedSpmdVectorView<Scalar> dv(v);
    const Ordinal localOffset = dv.spmdSpace()->localOffset();
    const Ordinal localSubDim = dv.spmdSpace()->localSubDim();
    for (Ordinal i = 0; i < localSubDim; ++i) {
      dv[i] = as<Scalar>(localOffset + i);
    }
  }

  out << "\nv = " << *v;

  {

    const ConstDetachedVectorView<Scalar> dv(v);

    TEST_EQUALITY(dv.subDim(), vs->dim());
    
    out << "\nTest that dv[i] == i, i=0...n-1 ... ";
    bool local_success = true;
    for (Ordinal i = 0; i < dv.subDim(); ++i) {
      TEST_ARRAY_ELE_EQUALITY( dv, i, as<Scalar>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;

  }

#ifdef THYRA_DEFAULT_SPMD_VECTOR_SPACE_UNIT_TESTS_DUMP
    RTOpPack::show_spmd_apply_op_dump = false;
#endif
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  parallelFullExtract)



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultSpmdVectorSpace, dangling_vs,
  Scalar )
{
  RCP<const Thyra::VectorSpaceBase<Scalar> > vs = 
    defaultSpmdVectorSpace<Scalar>(g_localDim);
  const int globalDim = vs->dim();
  RCP<VectorBase<Scalar> > x1 = createMember(vs);
  {
    // x1 owns a false RCP to vs
    TEST_EQUALITY_CONST( x1->space().has_ownership(), true );
    // RCP<> for x1 owns a true RCP to vs
    const std::string label = "VectorSpaceBase";
    RCP<const VectorSpaceBase<Scalar> > extra_data_x1 = 
      get_extra_data<RCP<const VectorSpaceBase<Scalar> >, VectorBase<Scalar> >(x1, label);
    TEST_EQUALITY_CONST( extra_data_x1.has_ownership(), true );
  }
  RCP<Thyra::VectorBase<Scalar> > x0 = x1->clone_v();
  {
    // x0 owns a false RCP to vs
    TEST_EQUALITY_CONST( x0->space().has_ownership(), true );
    // RCP<> for x0 owns a true RCP to a _DIFFERENT_ VectorSpaceBase
    // object because the one used to clone x1 is a false RCP, so the
    // VectorSpaceBase was cloned and that is the one that was set on the RCP.
    std::string label = "VectorSpaceBase";
    RCP<const VectorSpaceBase<Scalar> > extra_data_x0 = 
      get_extra_data<RCP<const
      VectorSpaceBase<Scalar> >, VectorBase<Scalar> >(x0, label );
    TEST_EQUALITY_CONST( extra_data_x0.has_ownership(), true );
    TEST_EQUALITY( extra_data_x0.ptr(), vs.ptr() );
  }
  vs = null; // vs still around because x1's RCP owns it
  x1 = null; // vs deleted
  {
    RCP<const VectorSpaceBase<Scalar> > vs_old = x0->space();
    TEST_EQUALITY_CONST( vs_old->dim(), globalDim );
  }
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultSpmdVectorSpace,
  dangling_vs)



} // namespace
