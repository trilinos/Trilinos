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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "Teuchos::as IS" AND ANY
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


#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_VectorStdOpsTester.hpp"
#include "Thyra_MultiVectorStdOpsTester.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_RowStatLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Thyra_EpetraVectorSpace.hpp"
#include "Thyra_EpetraMultiVector.hpp"
#include "Thyra_EpetraVector.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "Epetra_SerialComm.h"
#include "Teuchos_DefaultSerialComm.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"


#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayView;
using Teuchos::rcp_dynamic_cast;
using Teuchos::inOutArg;
using Teuchos::Comm;
using Teuchos::tuple;
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  using GO = long long;
#else
  using GO = int;
#endif

const int g_localDim = 4; // ToDo: Make variable!


RCP<const Epetra_Comm>
convertTeuchosCommToEpetraComm(const RCP<const Teuchos::Comm<int>> commT)
{
#ifdef HAVE_MPI
  auto mpiCommT = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(commT);
  if (!mpiCommT.is_null()) {
    return Teuchos::rcp( new Epetra_MpiComm((*mpiCommT->getRawMpiComm())()) );
  }
#endif
  auto serialCommT = Teuchos::rcp_dynamic_cast<const Teuchos::SerialComm<int>>(commT);
  if (!serialCommT.is_null()) {
    return Teuchos::rcp( new Epetra_SerialComm() );    
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error! Could not cast the Teuchos Comm to either MpiComm nor SerialComm.\n");

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

RCP<const Epetra_Map>
createEpetraMap(const int localDim)
{
  auto comm = convertTeuchosCommToEpetraComm(Teuchos::DefaultComm<int>::getComm());
  return Teuchos::rcp( new Epetra_Map(Teuchos::as<GO>(-1), localDim, 0, *comm) );
  // ToDo: Pass in the comm?
}


RCP<const VectorSpaceBase<double>>
createEpetraVectorSpace(const int localDim)
{
  return createVectorSpace(createEpetraMap(localDim));
}


RCP<Epetra_Operator>
createTriDiagonalEpetraOperator(const int numLocalRows)
{
  RCP<const Epetra_Map> map = createEpetraMap(numLocalRows);
  
  const int numMyElements = map->NumMyElements();
  const bool isLongLong = map->GlobalIndicesLongLong();
  GO numGlobalElements;
  if (isLongLong) {
    numGlobalElements = map->NumGlobalPoints64();
  } else {
#ifndef EPETRA_NO_32BIT_INDICES
    numGlobalElements = map->NumGlobalElements();
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Error! The indices in the Epetra_Map are not long long, "
                               "and yet 32bits indices are not enabled.\n");
#endif
  }

  GO* myGlobalElements;
  map->MyGlobalElementsPtr(myGlobalElements);

  // Create an OTeger vector numNz that is used to build the Petra Matrix.
  // numNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
  // on this processor

  Teuchos::ArrayRCP<int> numNz = Teuchos::arcp<int>(numMyElements);

  // We are building a tridiagonal matrix where each row has (-1 2 1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (int i=0; i < numMyElements; ++i) {
    if (myGlobalElements[i] == 0 || myGlobalElements[i] == numGlobalElements-1) {
      // boundary
      numNz[i] = 2;
    }
    else {
      numNz[i] = 3;
    }
  }

  // Create a Epetra_Matrix using the Map, with a static allocation dictated by numNz
  RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy, *map, numNz.getRawPtr(), true) );

  // We are done with NumNZ
  numNz = Teuchos::null;

  // Add  rows one-at-a-time
  // Off diagonal values will always be -1 on the left and +1 on the right
  const double two    =  2.0;
  const double posOne = +1.0;
  const double negOne = -1.0;
  Teuchos::Array<double> values(3);
  Teuchos::Array<GO> indices(3);
  int numEntries;
  for (int i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      numEntries = 2;
      indices[0] = myGlobalElements[i];
      indices[1] = myGlobalElements[i]+1;
      values[0] = two;
      values[1] = posOne;
    } else if (myGlobalElements[i] == (numGlobalElements-1)) {
      numEntries = 2;
      indices[0] = myGlobalElements[i]-1;
      indices[1] = myGlobalElements[i];
      values[0] = negOne;
      values[1] = two;
    } else {
      numEntries = 3;
      indices[0] = myGlobalElements[i]-1;
      indices[1] = myGlobalElements[i];
      indices[2] = myGlobalElements[i]+1;
      values[0] = negOne;
      values[1] = two;
      values[2] = posOne;
    }
    int err_code = A->InsertGlobalValues( myGlobalElements[i], numEntries, values().getConst().getRawPtr(), indices().getConst().getRawPtr());
    TEUCHOS_TEST_FOR_EXCEPTION(err_code!=0, std::runtime_error, "Error! Something went wrong while inserting global entries on row " << myGlobalElements[i] << ".\n");
  }

  // Finish up
  A->FillComplete();

  return A;
}


bool showAllTests = false;
bool dumpAll = false;
bool runLinearOpTester = true;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "show-all-tests", "no-show-all-tests", &showAllTests, "Show all tests or not" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-all", "no-dump-all", &dumpAll, "Dump all objects being tested" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "run-linear-op-tester", "no-run-linear-op-tester", &runLinearOpTester, "..." );
}


//
// Unit Tests
//


//
// convertEpetraToThyraComm
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, convertEpetraToThyraComm)
{
  RCP<const Epetra_Comm> epetraComm = convertTeuchosCommToEpetraComm(Teuchos::DefaultComm<int>::getComm());
  RCP<const Comm<Ordinal>> thyraComm = Thyra::convertEpetraToThyraComm(*epetraComm);
  TEST_ASSERT(nonnull(thyraComm));
}


//
// createVectorSpace
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createVectorSpace)
{
  const RCP<const Epetra_Map> epetraMap = createEpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<double>> vs = createVectorSpace(epetraMap);
  TEST_ASSERT(nonnull(vs));
  out << "vs = " << *vs;
  const RCP<const SpmdVectorSpaceBase<double>> vs_spmd =
    rcp_dynamic_cast<const SpmdVectorSpaceBase<double>>(vs, true);
  TEST_EQUALITY(vs_spmd->localSubDim(), g_localDim);
  TEST_EQUALITY(vs->dim(), Teuchos::as<Ordinal>(epetraMap->NumGlobalElements64()));
}


//
// createVector
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createVector)
{
  typedef EpetraOperatorVectorExtraction ConverterE;

  const RCP<const Epetra_Map> epetraMap = createEpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<double>> vs = createVectorSpace(epetraMap);

  const RCP<Epetra_Vector> epetraVector( new Epetra_Vector(*epetraMap) );

  // Test createVector with and without optional vector space
  {
    const RCP<VectorBase<double>> thyraVector = createVector(epetraVector, vs);
    TEST_EQUALITY(thyraVector->space(), vs);
    const RCP<Epetra_Vector> epetraVector2 = ConverterE::getEpetraVector(thyraVector);
    TEST_EQUALITY(epetraVector2, epetraVector);
  }

  {
    const RCP<VectorBase<double>> thyraVector = createVector(epetraVector);
    TEST_INEQUALITY(thyraVector->space(), vs);
    TEST_ASSERT(thyraVector->space()->isCompatible(*vs));
    const RCP<Epetra_Vector> epetraVector2 = ConverterE::getEpetraVector(thyraVector);
    TEST_EQUALITY(epetraVector2, epetraVector);
  }
}


//
// createConstVector
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createConstVector)
{
  typedef EpetraOperatorVectorExtraction ConverterE;

  const RCP<const Epetra_Map> epetraMap = createEpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<double>> vs = createVectorSpace(epetraMap);

  const RCP<const Epetra_Vector> epetraVector( new Epetra_Vector(*epetraMap) );

  // Test createConstVector with and without optional vector space
  {
    const RCP<const VectorBase<double>> thyraVector = createConstVector(epetraVector, vs);
    TEST_EQUALITY(thyraVector->space(), vs);
    const RCP<const Epetra_Vector> epetraVector2 = ConverterE::getConstEpetraVector(thyraVector);
    TEST_EQUALITY(epetraVector2, epetraVector);
  }

  {
    const RCP<const VectorBase<double>> thyraVector = createConstVector(epetraVector);
    TEST_INEQUALITY(thyraVector->space(), vs);
    TEST_ASSERT(thyraVector->space()->isCompatible(*vs));
    const RCP<const Epetra_Vector> epetraVector2 = ConverterE::getConstEpetraVector(thyraVector);
    TEST_EQUALITY(epetraVector2, epetraVector);
  }
}


//
// createMultiVector
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createMultiVector)
{
  typedef EpetraOperatorVectorExtraction ConverterE;

  const int numCols = 3;

  const RCP<const Epetra_Map> epetraMap = createEpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<double>> rangeVs = createVectorSpace(epetraMap);

  const RCP<const Epetra_Map> epetraLocRepMap( new Epetra_LocalMap(numCols,0,epetraMap->Comm()) );
  const RCP<const VectorSpaceBase<double>> domainVs = createVectorSpace(epetraLocRepMap);

  const RCP<Epetra_MultiVector> epetraMultiVector( new Epetra_MultiVector(*epetraMap, numCols) );

  // Test createMultiVector with and without optional vector spaces
  {
    const RCP<MultiVectorBase<double>> thyraMultiVector = createMultiVector(epetraMultiVector, rangeVs, domainVs);
    TEST_EQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_EQUALITY(thyraMultiVector->domain(), domainVs);
    const RCP<Epetra_MultiVector> epetraMultiVector2 = ConverterE::getEpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(epetraMultiVector2, epetraMultiVector);
  }

  {
    const RCP<MultiVectorBase<double>> thyraMultiVector = createMultiVector(epetraMultiVector);
    TEST_INEQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_INEQUALITY(thyraMultiVector->domain(), domainVs);
    TEST_ASSERT(thyraMultiVector->range()->isCompatible(*rangeVs));
    TEST_ASSERT(thyraMultiVector->domain()->isCompatible(*domainVs));
    const RCP<Epetra_MultiVector> epetraMultiVector2 = ConverterE::getEpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(epetraMultiVector2, epetraMultiVector);
  }
}


//
// createConstMultiVector
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createConstMultiVector)
{
  typedef EpetraOperatorVectorExtraction ConverterE;

  const int numCols = 3;

  const RCP<const Epetra_Map> epetraMap = createEpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<double>> rangeVs = createVectorSpace(epetraMap);

  const RCP<const Epetra_Map> epetraLocRepMap( new Epetra_LocalMap(numCols, 0, epetraMap->Comm()) );
  const RCP<const VectorSpaceBase<double>> domainVs = createVectorSpace(epetraLocRepMap);

  const RCP<const Epetra_MultiVector> epetraMultiVector( new Epetra_MultiVector(*epetraMap, numCols) );

  // Test createMultiVector with and without optional vector spaces
  {
    const RCP<const MultiVectorBase<double>> thyraMultiVector = createConstMultiVector(epetraMultiVector, rangeVs, domainVs);
    TEST_EQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_EQUALITY(thyraMultiVector->domain(), domainVs);
    const RCP<const Epetra_MultiVector> epetraMultiVector2 = ConverterE::getConstEpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(epetraMultiVector2, epetraMultiVector);
  }

  {
    const RCP<const MultiVectorBase<double>> thyraMultiVector = createConstMultiVector(epetraMultiVector);
    TEST_INEQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_INEQUALITY(thyraMultiVector->domain(), domainVs);
    TEST_ASSERT(thyraMultiVector->range()->isCompatible(*rangeVs));
    TEST_ASSERT(thyraMultiVector->domain()->isCompatible(*domainVs));
    const RCP<const Epetra_MultiVector> epetraMultiVector2 = ConverterE::getConstEpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(epetraMultiVector2, epetraMultiVector);
  }
}


//
// EpetraVectorSpace
//


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, TeptraVectorSpace)
{
  const RCP<const VectorSpaceBase<double>> vs = createVectorSpace(createEpetraMap(g_localDim));
  const RCP<VectorBase<double>> v = createMember(vs);
  TEST_ASSERT(nonnull(v));
  TEST_EQUALITY(v->space(), vs);
}


//
// vectorSpaceTester
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, vectorSpaceTester)
{
  const RCP<const VectorSpaceBase<double>> vs = createEpetraVectorSpace(g_localDim);
  Thyra::VectorSpaceTester<double> vectorSpaceTester;
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);
  TEST_ASSERT(vectorSpaceTester.check(*vs, &out));
}


//
// vectorStdOpsTester
//


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, vectorStdOpsTester)
{
  const RCP<const VectorSpaceBase<double>> vs = createVectorSpace(createEpetraMap(g_localDim));
  Thyra::VectorStdOpsTester<double> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(5.0e-13);
  vectorStdOpsTester.error_tol(5.0e-14);
  TEST_ASSERT(vectorStdOpsTester.checkStdOps(*vs, &out));
}


//
// multiVectorStdOpsTester
//


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, multiVectorStdOpsTester)
{
  const RCP<const VectorSpaceBase<double>> vs = createVectorSpace(createEpetraMap(g_localDim));
  Thyra::MultiVectorStdOpsTester<double> mvStdOpsTester;
  mvStdOpsTester.warning_tol(5.0e-13);
  mvStdOpsTester.error_tol(5.0e-14);
  TEST_ASSERT(mvStdOpsTester.checkStdOps(*vs, &out));
}


//
// getEpetraMultiVector
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, getEpetraMultiVector)
{
  typedef EpetraOperatorVectorExtraction ConverterE;

  const int numCols = 3;
  const RCP<const VectorSpaceBase<double>> vs = createEpetraVectorSpace(g_localDim);

  {
    const RCP<MultiVectorBase<double>> mv = createMembers(vs, numCols);
    const RCP<Epetra_MultiVector> emv = ConverterE::getEpetraMultiVector(mv);
    TEST_ASSERT(nonnull(emv));
    TEST_EQUALITY(Teuchos::as<Ordinal>(emv->Map().NumGlobalElements64()), vs->dim());
  }

  {
    const RCP<VectorBase<double>> v = createMember(vs);
    const RCP<Epetra_MultiVector> emv = ConverterE::getEpetraMultiVector(v);
    TEST_ASSERT(nonnull(emv));
    TEST_EQUALITY(Teuchos::as<Ordinal>(emv->Map().NumGlobalElements64()), vs->dim());
  }

#ifdef THYRA_DEBUG
  const RCP<VectorBase<double>> pv = Thyra::defaultProductVector<double>();
  TEST_THROW(ConverterE::getEpetraMultiVector(pv), std::logic_error);
#endif

}


//
// getConstEpetraMultiVector
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, getConstEpetraMultiVector)
{
  typedef EpetraOperatorVectorExtraction ConverterE;

  const int numCols = 3;
  const RCP<const VectorSpaceBase<double>> vs = createEpetraVectorSpace(g_localDim);

  {
    const RCP<const MultiVectorBase<double>> mv = createMembers(vs, numCols);
    const RCP<const Epetra_MultiVector> emv = ConverterE::getConstEpetraMultiVector(mv);
    TEST_ASSERT(nonnull(emv));
    TEST_EQUALITY(Teuchos::as<Ordinal>(emv->Map().NumGlobalPoints64()), vs->dim());
  }

  {
    const RCP<const VectorBase<double>> v = createMember(vs);
    const RCP<const Epetra_MultiVector> emv = ConverterE::getConstEpetraMultiVector(v);
    TEST_ASSERT(nonnull(emv));
    TEST_EQUALITY(Teuchos::as<Ordinal>(emv->Map().NumGlobalPoints64()), vs->dim());
  }

#ifdef THYRA_DEBUG
  const RCP<const VectorBase<double>> pv = Thyra::defaultProductVector<double>();
  TEST_THROW(ConverterE::getConstEpetraMultiVector(pv), std::logic_error);
#endif

}


//
// EpetraLinearOp
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, EpetraLinearOp)
{
  using Teuchos::rcpFromRef;
  typedef Teuchos::ScalarTraits<double> ST;

  const RCP<Epetra_Operator> epetraOp = createTriDiagonalEpetraOperator(g_localDim);
  out << "epetraOp = " << epetraOp->Label() << std::endl;
  TEST_ASSERT(nonnull(epetraOp));

  const RCP<const VectorSpaceBase<double>> rangeSpace  = createVectorSpace(rcpFromRef(epetraOp->OperatorRangeMap()));
  const RCP<const VectorSpaceBase<double>> domainSpace = createVectorSpace(rcpFromRef(epetraOp->OperatorDomainMap()));
  const RCP<const LinearOpBase<double>> thyraLinearOp = epetraLinearOp(rangeSpace, domainSpace, epetraOp);
  TEST_ASSERT(nonnull(thyraLinearOp));

  out << "\nCheck that operator returns the right thing ...\n";
  const RCP<VectorBase<double>> x = createMember(thyraLinearOp->domain());
  Thyra::V_S(x.ptr(), ST::one());
  const RCP<VectorBase<double>> y = createMember(thyraLinearOp->range());
  Thyra::apply<double>(*thyraLinearOp, Thyra::NOTRANS, *x, y.ptr());
  const double sum_y = sum(*y);
  TEST_FLOATING_EQUALITY( sum_y, Teuchos::as<double>(3+1+2*(y->space()->dim()-2)),
    100.0 * ST::eps() );

  out << "\nCheck the general LinearOp interface ...\n";
  Thyra::LinearOpTester<double> linearOpTester;
  linearOpTester.show_all_tests(showAllTests);
  linearOpTester.dump_all(dumpAll);
  if (runLinearOpTester) {
    TEST_ASSERT(linearOpTester.check(*thyraLinearOp, Teuchos::inOutArg(out)));
  }
}


//
// createLinearOp
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createLinearOp)
{
  using Teuchos::rcpFromRef;
  typedef EpetraOperatorVectorExtraction ConverterE;

  const RCP<Epetra_Operator> epetraOp = createTriDiagonalEpetraOperator(g_localDim);
  out << "epetraOp = " << epetraOp->Label() << std::endl;

  const RCP<const VectorSpaceBase<double>> rangeSpace  = createVectorSpace(rcpFromRef(epetraOp->OperatorRangeMap()));
  const RCP<const VectorSpaceBase<double>> domainSpace = createVectorSpace(rcpFromRef(epetraOp->OperatorDomainMap()));

  {
    const RCP<LinearOpBase<double>> thyraOp = createLinearOp(epetraOp, rangeSpace, domainSpace);
    TEST_EQUALITY(thyraOp->range(), rangeSpace);
    TEST_EQUALITY(thyraOp->domain(), domainSpace);
    const RCP<Epetra_Operator> epetraOp2 = ConverterE::getEpetraOperator(thyraOp);
    TEST_EQUALITY(epetraOp2, epetraOp);
  }

  {
    const RCP<LinearOpBase<double>> thyraOp = createLinearOp(epetraOp);
    TEST_INEQUALITY(thyraOp->range(), rangeSpace);
    TEST_INEQUALITY(thyraOp->domain(), domainSpace);
    TEST_ASSERT(thyraOp->range()->isCompatible(*rangeSpace));
    TEST_ASSERT(thyraOp->domain()->isCompatible(*domainSpace));
    const RCP<Epetra_Operator> epetraOp2 = ConverterE::getEpetraOperator(thyraOp);
    TEST_EQUALITY(epetraOp2, epetraOp);
  }
}


//
// createConstLinearOp
//

TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createConstLinearOp)
{
  using Teuchos::rcpFromRef;
  typedef EpetraOperatorVectorExtraction ConverterE;

  const RCP<const Epetra_Operator> epetraOp = createTriDiagonalEpetraOperator(g_localDim);
  out << "epetraOp = " << epetraOp->Label() << std::endl;

  const RCP<const VectorSpaceBase<double>> rangeSpace  = createVectorSpace(rcpFromRef(epetraOp->OperatorRangeMap()));
  const RCP<const VectorSpaceBase<double>> domainSpace = createVectorSpace(rcpFromRef(epetraOp->OperatorDomainMap()));

  {
    const RCP<const LinearOpBase<double>> thyraOp = createConstLinearOp(epetraOp, rangeSpace, domainSpace);
    TEST_EQUALITY(thyraOp->range(), rangeSpace);
    TEST_EQUALITY(thyraOp->domain(), domainSpace);
    const RCP<const Epetra_Operator> epetraOp2 = ConverterE::getConstEpetraOperator(thyraOp);
    TEST_EQUALITY(epetraOp2, epetraOp);
  }

  {
    const RCP<const LinearOpBase<double>> thyraOp = createConstLinearOp(epetraOp);
    TEST_INEQUALITY(thyraOp->range(), rangeSpace);
    TEST_INEQUALITY(thyraOp->domain(), domainSpace);
    TEST_ASSERT(thyraOp->range()->isCompatible(*rangeSpace));
    TEST_ASSERT(thyraOp->domain()->isCompatible(*domainSpace));
    const RCP<const Epetra_Operator> epetraOp2 = ConverterE::getConstEpetraOperator(thyraOp);
    TEST_EQUALITY(epetraOp2, epetraOp);
  }
}


//
// EpetraLinearOp_RowStatLinearOpBase
//


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, EpetraLinearOp_RowStatLinearOpBase)
{
  using Teuchos::rcpFromRef;
  using ST = Teuchos::ScalarTraits<double>;

  const RCP<Epetra_Operator> epetraOp = createTriDiagonalEpetraOperator(g_localDim);
  out << "epetraOp = " << epetraOp->Label() << std::endl;
  TEST_ASSERT(nonnull(epetraOp));

  const RCP<Epetra_RowMatrix> epetraRowMatrix =
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetraOp,true);

  const RCP<const VectorSpaceBase<double>> rangeSpace  = createVectorSpace(rcpFromRef(epetraOp->OperatorRangeMap()));
  const RCP<const VectorSpaceBase<double>> domainSpace = createVectorSpace(rcpFromRef(epetraOp->OperatorDomainMap()));
  const RCP<LinearOpBase<double>> thyraLinearOp = epetraLinearOp(rangeSpace, domainSpace, epetraOp);
  TEST_ASSERT(nonnull(thyraLinearOp));

  const Teuchos::RCP<Thyra::RowStatLinearOpBase<double>> rowStatOp =
    Teuchos::rcp_dynamic_cast<Thyra::RowStatLinearOpBase<double>>(thyraLinearOp, true);

  // Get the inverse row sums

  const RCP<VectorBase<double>> inv_row_sums = createMember(thyraLinearOp->range());
  const RCP<VectorBase<double>> row_sums = createMember(thyraLinearOp->range());

  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());
  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  out << "inv_row_sums = " << *inv_row_sums;
  out << "row_sums = " << *row_sums;

  TEST_FLOATING_EQUALITY(
    Thyra::sum<double>(*row_sums),
    Teuchos::as<double>(4.0 * thyraLinearOp->domain()->dim() - 2.0),
    Teuchos::as<double>(10.0 * ST::eps())
    );

  TEST_FLOATING_EQUALITY(
    Thyra::sum<double>(*inv_row_sums),
    Teuchos::as<double>( 1.0 / 4.0 * (thyraLinearOp->domain()->dim() - 2) + 2.0 / 3.0 ),
    Teuchos::as<double>(10.0 * ST::eps())
    );
}


//
// EpetraLinearOp_RowStatLinearOpBase
//


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, EpetraLinearOp_ScaledLinearOpBase)
{
  using Teuchos::rcpFromRef;

  const RCP<Epetra_Operator> epetraOp = createTriDiagonalEpetraOperator(g_localDim);
  out << "epetraOp = " << epetraOp->Label() << std::endl;
  TEST_ASSERT(nonnull(epetraOp));

  const RCP<Epetra_RowMatrix> epetraRowMatrix =
    Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetraOp,true);

  const RCP<const VectorSpaceBase<double>> rangeSpace  = createVectorSpace(rcpFromRef(epetraOp->OperatorRangeMap()));
  const RCP<const VectorSpaceBase<double>> domainSpace = createVectorSpace(rcpFromRef(epetraOp->OperatorDomainMap()));
  const RCP<LinearOpBase<double>> thyraLinearOp = epetraLinearOp(rangeSpace, domainSpace, epetraOp);
  TEST_ASSERT(nonnull(thyraLinearOp));

  const Teuchos::RCP<Thyra::RowStatLinearOpBase<double>> rowStatOp =
    Teuchos::rcp_dynamic_cast<Thyra::RowStatLinearOpBase<double>>(thyraLinearOp, true);

  // Get the inverse row sums

  const RCP<VectorBase<double>> inv_row_sums = createMember(thyraLinearOp->range());
  const RCP<VectorBase<double>> row_sums = createMember(thyraLinearOp->range());

  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());
  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  out << "inv_row_sums = " << *inv_row_sums;
  out << "row_sums = " << *row_sums;

  const Teuchos::RCP<Thyra::ScaledLinearOpBase<double>> scaledOp =
    Teuchos::rcp_dynamic_cast<Thyra::ScaledLinearOpBase<double>>(thyraLinearOp, true);

  TEUCHOS_ASSERT(scaledOp->supportsScaleLeft());

  scaledOp->scaleLeft(*inv_row_sums);

  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  out << "row_sums after left scaling by inv_row_sum = " << *row_sums;

  // scaled row sums should be one for each entry
  TEST_FLOATING_EQUALITY(
    Teuchos::as<double>(row_sums->space()->dim()),
    Thyra::sum<double>(*row_sums),
    Teuchos::as<double>(10.0 * Teuchos::ScalarTraits<double>::eps())
    );

  // Don't currently check the results of right scaling.  Epetra tests
  // already check this.  Once column sums are supported in epetra
  // adapters, this can be checked easily.
  TEUCHOS_ASSERT(scaledOp->supportsScaleRight());
  scaledOp->scaleRight(*inv_row_sums);
  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,row_sums.ptr());
  out << "row_sums after right scaling by inv_row_sum = " << *row_sums;
}


} // namespace Thyra
