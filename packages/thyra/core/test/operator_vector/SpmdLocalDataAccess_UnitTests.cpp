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
createLocallyReplicatedVS(const Ordinal localDim)
{
  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();

  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, localDim, localDim);
  // ToDo: Pass in argument to state that space is locally replicated!

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
createProcRankLocalDimVS()
{
  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();

  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, comm->getRank()+1, -1);

  return vs;
}


template<class Scalar>
RCP<const DefaultSpmdVectorSpace<Scalar> >
createZeroVS()
{
  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();

  RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm, 0, -1);

  return vs;
}


//
// Unit Tests
//


//
// Test getLocalSubVectorView
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess, 
  getLocalSubVectorView_procRankLocalDim, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createProcRankLocalDimVS<Scalar>();
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  TEST_EQUALITY(vs->isLocallyReplicated(), numProcs==1);
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
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  getLocalSubVectorView_procRankLocalDim)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  getLocalSubVectorView_empty_p0, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(g_localDim);
  TEST_ASSERT(!vs->isLocallyReplicated());
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
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  getLocalSubVectorView_empty_p0)


//
// Test getNonconstLocalSubVectorView
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  getNonconstLocalSubVectorView_procRankLocalDim, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createProcRankLocalDimVS<Scalar>();
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
    TEST_EQUALITY(lsv.globalOffset(), as<Ordinal>((procRank*(procRank+1))/2));
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
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess, 
  getNonconstLocalSubVectorView_procRankLocalDim)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  getNonconstLocalSubVectorView_empty_p0, Scalar )
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
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  getNonconstLocalSubVectorView_empty_p0)


//
// Test getLocalSubMultiVectorView
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  getLocalSubMultiVectorView_procRankLocalDim, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createProcRankLocalDimVS<Scalar>();
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, g_numCols);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(mv.ptr(), val);
  const ScalarMag tol = 100.0*ScalarTraits<Scalar>::eps();
  TEST_EQUALITY_CONST(mv->domain()->dim(), g_numCols);
  for (int j = 0; j < g_numCols; ++j) {
    TEST_FLOATING_EQUALITY(sum<Scalar>(*mv->col(0)), as<Scalar>(val * vs->dim()), tol);
  }
  out << "*** Test that we get the view correctly ...\n";
  RTOpPack::ConstSubMultiVectorView<Scalar> lsmv = 
    getLocalSubMultiVectorView<Scalar>(mv);
  TEST_EQUALITY(lsmv.globalOffset(), as<Ordinal>((procRank*(procRank+1))/2));
  TEST_EQUALITY(lsmv.subDim(), procRank+1);
  TEST_EQUALITY(lsmv.leadingDim(), lsmv.subDim());
  TEST_EQUALITY_CONST(lsmv.colOffset(), 0);
  TEST_EQUALITY(lsmv.numSubCols(), g_numCols);
  for (int i = 0; i < lsmv.subDim(); ++i) {
    for (int j = 0; j < lsmv.numSubCols(); ++j) {
      TEST_EQUALITY(lsmv(i,j), val);
    }
  }
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  getLocalSubMultiVectorView_procRankLocalDim)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  getLocalSubMultiVectorView_empty_p0, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(g_localDim);
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, g_numCols);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(mv.ptr(), val);
  out << "*** Test that we get the view correctly including an empty view on p0 ...\n";
  RTOpPack::ConstSubMultiVectorView<Scalar> lsmv = 
    getLocalSubMultiVectorView<Scalar>(mv);
  if (procRank == 0) {
    TEST_EQUALITY_CONST(lsmv.globalOffset(), 0);
    TEST_EQUALITY_CONST(lsmv.subDim(), 0);
    TEST_EQUALITY_CONST(lsmv.values(), null);
  }
  else {
    TEST_EQUALITY(lsmv.globalOffset(), as<Ordinal>((procRank-1)*g_localDim));
    TEST_EQUALITY(lsmv.subDim(), g_localDim);
  }
  TEST_EQUALITY(lsmv.leadingDim(), lsmv.subDim());
  TEST_EQUALITY_CONST(lsmv.colOffset(), 0);
  TEST_EQUALITY(lsmv.numSubCols(), g_numCols);
  for (int i = 0; i < lsmv.subDim(); ++i) {
    for (int j = 0; j < lsmv.numSubCols(); ++j) {
      TEST_EQUALITY(lsmv(i,j), val);
    }
  }
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  getLocalSubMultiVectorView_empty_p0)


//
// Test getNonconstLocalSubMultiVectorView
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  getNonconstLocalSubMultiVectorView_procRankLocalDim, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createProcRankLocalDimVS<Scalar>();
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, g_numCols);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(mv.ptr(), val);
  {
    out << "*** Test that we get and change the nonconst view correctly ...\n";
    RTOpPack::SubMultiVectorView<Scalar> lsmv = 
      getNonconstLocalSubMultiVectorView<Scalar>(mv);
    TEST_EQUALITY(lsmv.globalOffset(), as<Ordinal>((procRank*(procRank+1))/2));
    TEST_EQUALITY(lsmv.subDim(), procRank+1);
    TEST_EQUALITY(lsmv.leadingDim(), lsmv.subDim());
    TEST_EQUALITY_CONST(lsmv.colOffset(), 0);
    TEST_EQUALITY(lsmv.numSubCols(), g_numCols);
    for (int i = 0; i < lsmv.subDim(); ++i) {
      for (int j = 0; j < lsmv.numSubCols(); ++j) {
        lsmv(i,j) = lsmv.globalOffset() + i + 0.1 * j;
      }
    }
  }
  {
    out << "*** Test that we get the same values when we grab const view ...\n";
    RTOpPack::ConstSubMultiVectorView<Scalar> lsmv = 
      getLocalSubMultiVectorView<Scalar>(mv);
    for (int i = 0; i < lsmv.subDim(); ++i) {
      for (int j = 0; j < lsmv.numSubCols(); ++j) {
        TEST_EQUALITY(lsmv(i,j), as<Scalar>(lsmv.globalOffset() + i + 0.1 * j));
      }
    }
  }
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  getNonconstLocalSubMultiVectorView_procRankLocalDim)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  getNonconstLocalSubMultiVectorView_empty_p0, Scalar )
{
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroEleProcVS<Scalar>(g_localDim);
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const int procRank = comm->getRank();
  PRINT_VAR(procRank);
  const int numProcs = comm->getSize();
  PRINT_VAR(numProcs);
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, g_numCols);
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);
  assign<Scalar>(mv.ptr(), val);
  {
    out << "*** Test that we get and change the nonconst view correctly ...\n";
    RTOpPack::SubMultiVectorView<Scalar> lsmv = 
      getNonconstLocalSubMultiVectorView<Scalar>(mv);
    if (procRank == 0) {
      TEST_EQUALITY_CONST(lsmv.globalOffset(), 0);
      TEST_EQUALITY_CONST(lsmv.subDim(), 0);
      TEST_EQUALITY_CONST(lsmv.values(), null);
    }
    else {
      TEST_EQUALITY(lsmv.globalOffset(), as<Ordinal>((procRank-1)*g_localDim));
      TEST_EQUALITY(lsmv.subDim(), g_localDim);
    }
    TEST_EQUALITY_CONST(lsmv.colOffset(), 0);
    TEST_EQUALITY(lsmv.leadingDim(), lsmv.subDim());
    TEST_EQUALITY(lsmv.numSubCols(), g_numCols);
    for (int i = 0; i < lsmv.subDim(); ++i) {
      for (int j = 0; j < lsmv.numSubCols(); ++j) {
        lsmv(i,j) = lsmv.globalOffset() + i + 0.1 * j;
      }
    }
  }
  {
    out << "*** Test that we get the same values when we grab const view ...\n";
    RTOpPack::ConstSubMultiVectorView<Scalar> lsmv = 
      getLocalSubMultiVectorView<Scalar>(mv);
    for (int i = 0; i < lsmv.subDim(); ++i) {
      for (int j = 0; j < lsmv.numSubCols(); ++j) {
        TEST_EQUALITY(lsmv(i,j), as<Scalar>(lsmv.globalOffset() + i + 0.1 * j));
      }
    }
  }
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  getNonconstLocalSubMultiVectorView_empty_p0)


//
// Locally replicated objects
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  locallyReplicated, Scalar )
{
  out << "Create a locally replicated vector space ...\n";
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createLocallyReplicatedVS<Scalar>(g_localDim);
  TEST_EQUALITY(vs->dim(), g_localDim);
  TEST_ASSERT(vs->isLocallyReplicated());
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);

  out << "Test locally replicated Vector ...\n";
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  {
    ECHO(RTOpPack::SubVectorView<Scalar> lsv = 
      getNonconstLocalSubVectorView<Scalar>(v));
    TEST_EQUALITY_CONST(lsv.globalOffset(), 0);
    TEST_EQUALITY(lsv.subDim(), g_localDim);
    TEST_EQUALITY_CONST(lsv.stride(), 1);
    for (int k = 0; k < lsv.subDim(); ++k) {
      lsv[k] = k + 1;
    }
  }
  {
    ECHO(RTOpPack::ConstSubVectorView<Scalar> lsv = 
      getLocalSubVectorView<Scalar>(v));
    TEST_EQUALITY_CONST(lsv.globalOffset(), 0);
    TEST_EQUALITY(lsv.subDim(), g_localDim);
    TEST_EQUALITY_CONST(lsv.stride(), 1);
    for (int k = 0; k < lsv.subDim(); ++k) {
      TEST_EQUALITY(lsv[k], lsv.globalOffset() + k + 1);
    }
  }
  
  out << "Test locally replicated MultiVector ...\n";
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, g_numCols);
  {
    ECHO(RTOpPack::SubMultiVectorView<Scalar> lsmv = 
      getNonconstLocalSubMultiVectorView<Scalar>(mv));
    TEST_EQUALITY(lsmv.globalOffset(), 0);
    TEST_EQUALITY(lsmv.subDim(), g_localDim);
    TEST_EQUALITY(lsmv.leadingDim(), lsmv.subDim());
    TEST_EQUALITY_CONST(lsmv.colOffset(), 0);
    TEST_EQUALITY(lsmv.numSubCols(), g_numCols);
    for (int i = 0; i < lsmv.subDim(); ++i) {
      for (int j = 0; j < lsmv.numSubCols(); ++j) {
        lsmv(i,j) = i + 0.1 * j;
      }
    }
  }
  {
    ECHO(RTOpPack::ConstSubMultiVectorView<Scalar> lsmv = 
      getLocalSubMultiVectorView<Scalar>(mv));
    for (int i = 0; i < lsmv.subDim(); ++i) {
      for (int j = 0; j < lsmv.numSubCols(); ++j) {
        TEST_EQUALITY(lsmv(i,j), as<Scalar>(i + 0.1 * j));
      }
    }
  }
  
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  locallyReplicated)


//
// Objects for a zero-sized vector space
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SpmdLocalDataAccess,
  zeroVS, Scalar )
{
  out << "Create a locally replicated vector space ...\n";
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  const RCP<const DefaultSpmdVectorSpace<Scalar> > vs =
    createZeroVS<Scalar>();
  TEST_ASSERT(!vs->isLocallyReplicated());
  TEST_EQUALITY_CONST(vs->dim(), 0);
  const RCP<const Teuchos::Comm<Ordinal> > comm = vs->getComm();
  const Scalar val = as<Scalar>(1.5);
  PRINT_VAR(val);

  out << "Test locally replicated Vector ...\n";
  const RCP<VectorBase<Scalar> > v = createMember<Scalar>(vs);
  {
    ECHO(RTOpPack::SubVectorView<Scalar> lsv = 
      getNonconstLocalSubVectorView<Scalar>(v));
    TEST_EQUALITY_CONST(lsv.globalOffset(), 0);
    TEST_EQUALITY(lsv.subDim(), 0);
    TEST_EQUALITY_CONST(lsv.stride(), 1);
  }
  {
    ECHO(RTOpPack::ConstSubVectorView<Scalar> lsv = 
      getLocalSubVectorView<Scalar>(v));
    TEST_EQUALITY_CONST(lsv.globalOffset(), 0);
    TEST_EQUALITY(lsv.subDim(), 0);
    TEST_EQUALITY_CONST(lsv.stride(), 1);
  }
  assign<Scalar>(v.ptr(), val);
  TEST_EQUALITY(sum<Scalar>(*v), as<Scalar>(0.0));
  
  out << "Test locally replicated MultiVector ...\n";
  const RCP<MultiVectorBase<Scalar> > mv = createMembers<Scalar>(vs, g_numCols);
  {
    ECHO(RTOpPack::SubMultiVectorView<Scalar> lsmv = 
      getNonconstLocalSubMultiVectorView<Scalar>(mv));
    TEST_EQUALITY(lsmv.globalOffset(), 0);
    TEST_EQUALITY(lsmv.subDim(), 0);
    TEST_EQUALITY(lsmv.leadingDim(), lsmv.subDim());
    TEST_EQUALITY_CONST(lsmv.colOffset(), 0);
    TEST_EQUALITY(lsmv.numSubCols(), g_numCols);
  }
  {
    ECHO(RTOpPack::ConstSubMultiVectorView<Scalar> lsmv = 
      getLocalSubMultiVectorView<Scalar>(mv));
  }
  assign<Scalar>(mv.ptr(), val);
  for (int j = 0; j < g_numCols; ++j) {
    TEST_EQUALITY(sum<Scalar>(*mv->col(0)), as<Scalar>(0.0));
  }
  
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( SpmdLocalDataAccess,
  zeroVS)


} // namespace Thyra
