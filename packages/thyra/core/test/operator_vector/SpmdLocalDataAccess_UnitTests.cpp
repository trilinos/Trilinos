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
#include "Thyra_SpmdLocalDataAccess.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"


namespace {


int g_localDim = 4;
int g_numCols = 3;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "local-dim", &g_localDim, "Local dimension of each vector." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "num-cols", &g_numCols, "Number of columns in each multi-vector." );
}


} // namespace 


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::get_extra_data;


#define PRINT_VAR(varName) \
  out << #varName" = " << (varName) << "\n"


template<class Scalar>
RCP<const DefaultSpmdVectorSpace<Scalar> >
createVS(const Ordinal localDim)
{
  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();

  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, localDim, -1);

  return vs;
}


template<class Scalar>
RCP<const DefaultSpmdVectorSpace<Scalar> >
createZeroEleProcVS(const Ordinal localDim)
{
  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();

  const Ordinal thisLocalDim = comm->getRank()==0 ? 0 : localDim;
  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, thisLocalDim, -1);

  return vs;
}


template<class Scalar>
RCP<const DefaultSpmdVectorSpace<Scalar> >
createProcDimVS()
{
  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();

  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, comm->getRank()+1, -1);

  return vs;
}


//
// Unit Tests
//


//
// Test getLocalSubVectorView
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( getLocalSubVectorView, basic,
  Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createProcDimVS<Scalar>();
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(v.ptr(), val);
  const ScalarMag tol = 100.0*ScalarTraits<Scalar>::eps();
  TEST_FLOATING_EQUALITY(sum<Scalar>(*v), as<Scalar>(val * vs->dim()), tol);
  out << "*** Test that we get the view correctly ...\n";
  RTOpPack::ConstSubVectorView<Scalar> lsv = 
    getLocalSubVectorView<Scalar>(v);
  TEST_EQUALITY(lsv.globalOffset(), as<Ordinal>((procRank*(procRank+1))/2));
  TEST_EQUALITY(lsv.subDim(), procRank+1);
  TEST_EQUALITY_CONST(lsv.stride(), 1);
  for (int k = 0; k < lsv.subDim(); ++k) {
    TEST_EQUALITY(lsv[k], val);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( getLocalSubVectorView,
  basic)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( getLocalSubVectorView, empty_p0,
  Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(g_localDim);
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(v.ptr(), val);
  out << "*** Test that we get the view correctly including an empty view on p0 ...\n";
  RTOpPack::ConstSubVectorView<Scalar> lsv = 
    getLocalSubVectorView<Scalar>(v);
  if (procRank == 0) {
    TEST_EQUALITY_CONST(lsv.globalOffset(), 0);
    TEST_EQUALITY_CONST(lsv.subDim(), 0);
    TEST_EQUALITY_CONST(lsv.values(), null);
  }
  else {
    TEST_EQUALITY(lsv.globalOffset(), as<Ordinal>((procRank-1)*g_localDim));
    TEST_EQUALITY(lsv.subDim(), g_localDim);
  }
  TEST_EQUALITY_CONST(lsv.stride(), 1);
  for (int k = 0; k < lsv.subDim(); ++k) {
    TEST_EQUALITY(lsv[k], val);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( getLocalSubVectorView,
  empty_p0)


//
// Test getNonconstLocalSubVectorView
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( getNonconstLocalSubVectorView, basic,
  Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createProcDimVS<Scalar>();
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(v.ptr(), val);
  const ScalarMag tol = 100.0*ScalarTraits<Scalar>::eps();
  TEST_FLOATING_EQUALITY(sum<Scalar>(*v), as<Scalar>(val * vs->dim()), tol);
  {
    out << "*** Test that we get and change the nonconst view correctly ...\n";
    RTOpPack::SubVectorView<Scalar> lsv = 
      getNonconstLocalSubVectorView<Scalar>(v);
    TEST_EQUALITY(lsv.subDim(), procRank+1);
    TEST_EQUALITY_CONST(lsv.stride(), 1);
    for (int k = 0; k < lsv.subDim(); ++k) {
      lsv[k] = lsv.globalOffset() + k + 1;
    }
    const Ordinal n = vs->dim();
    TEST_FLOATING_EQUALITY(sum<Scalar>(*v), as<Scalar>((n*(n+1))/2.0), tol);
  }
  {
    out << "*** Test that we get the same values when we grab const view ...\n";
    RTOpPack::ConstSubVectorView<Scalar> lsv = 
      getLocalSubVectorView<Scalar>(v);
    TEST_EQUALITY(lsv.subDim(), procRank+1);
    TEST_EQUALITY_CONST(lsv.stride(), 1);
    for (int k = 0; k < lsv.subDim(); ++k) {
      TEST_EQUALITY(lsv[k], lsv.globalOffset() + k + 1);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( getNonconstLocalSubVectorView,
  basic)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( getNonconstLocalSubVectorView, empty_p0,
  Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(g_localDim);
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(v.ptr(), val);
  out << "*** Test that we get the view correctly including an empty view on p0 ...\n";
  RTOpPack::SubVectorView<Scalar> lsv = 
    getNonconstLocalSubVectorView<Scalar>(v);
  if (procRank == 0) {
    TEST_EQUALITY_CONST(lsv.subDim(), 0);
    TEST_EQUALITY_CONST(lsv.values(), null);
  }
  else {
    TEST_EQUALITY(lsv.subDim(), g_localDim);
  }
  TEST_EQUALITY_CONST(lsv.stride(), 1);
  for (int k = 0; k < lsv.subDim(); ++k) {
    TEST_EQUALITY(lsv[k], val);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( getNonconstLocalSubVectorView,
  empty_p0)





//
// Test getLocalSubMultiVectorView
//

/*

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( getLocalSubMultiVectorView, basic,
  Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs = createProcDimVS<Scalar>();
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<MultiVectorBase<Scalar> > v = createMembers<Scalar>(vs, g_numCols);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(mv.ptr(), val);
  const ScalarMag tol = 100.0*ScalarTraits<Scalar>::eps();
  TEST_FLOATING_EQUALITY(sum<Scalar>(*v), as<Scalar>(val * vs->dim()), tol);
  out << "*** Test that we get the view correctly ...\n";
  RTOpPack::ConstSubMultiVectorView<Scalar> lsmv = 
    getLocalSubMultiVectorView<Scalar>(v);
  TEST_EQUALITY(lsmv.subDim(), procRank+1);
  for (int k = 0; k < lsmv.subDim(); ++k) {
    TEST_EQUALITY(lsmv[k], val);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( getLocalSubMultiVectorView,
  basic)

*/


// ToDo:

// Test a basic construction on all processes returns the right view object

// Test that getting a non-const view the modifying it will modify the
// underlying V or MV object once the view is released.

// Test that getting data from a non SPMD object throws a defined excetpion on
// every process.

// Test that getting data from an unsized VS V or MV object returns an empty
// view on every process.

// ToDo: Write a generic testing class that runs all of these unit tests on
// any set of SPMD object to see that it does the right thing.


} // namespace Thyra
