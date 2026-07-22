// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_TpetraThyraWrappers.hpp"
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
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Tuple.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayView;
using Teuchos::rcp_dynamic_cast;
using Teuchos::inOutArg;
using Teuchos::Comm;
using Teuchos::tuple;


const int g_localDim = 4; // ToDo: Make variable!


typedef Tpetra::Map<> TpetraMap_t;


RCP<const TpetraMap_t>
createTpetraMap(const int localDim)
{
  typedef Teuchos::OrdinalTraits<Tpetra::global_size_t> OT;
  return Teuchos::rcp(new TpetraMap_t(OT::invalid(), localDim, 0,
      Teuchos::DefaultComm<int>::getComm()));
  // ToDo: Pass in the comm?
}

// ToDo: Above: Vary the LocalOrdinal and GlobalOrdinal types?


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
createTpetraVectorSpace(const int localDim)
{
  return Thyra::createVectorSpace<Scalar>(createTpetraMap(g_localDim));
}


template<class Scalar>
RCP<Tpetra::Operator<Scalar> >
createTriDiagonalTpetraOperator(const int numLocalRows)
{
  typedef Tpetra::global_size_t global_size_t;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  RCP<const Tpetra::Map<> > map = createTpetraMap(numLocalRows);

  const size_t numMyElements = map->getLocalNumElements();
  const global_size_t numGlobalElements = map->getGlobalNumElements();

  ArrayView<const GO> myGlobalElements = map->getLocalElementList();

  // Create an OTeger vector numNz that is used to build the Petra Matrix.
  // numNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
  // on this processor

  Teuchos::Array<size_t> numNz (numMyElements);

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (size_t i=0; i < numMyElements; ++i) {
    if (myGlobalElements[i] == 0 || static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
      // boundary
      numNz[i] = 2;
    }
    else {
      numNz[i] = 3;
    }
  }

  // Create a Tpetra::Matrix using the Map, with a static allocation dictated by numNz
  RCP< Tpetra::CrsMatrix<Scalar> > A =
    Teuchos::rcp( new Tpetra::CrsMatrix<Scalar>(map, numNz ()) );

  // We are done with NumNZ
  {
    Teuchos::Array<size_t> empty;
    swap (empty, numNz); // classic idiom for freeing a container
  }

  // Add  rows one-at-a-time
  // Off diagonal values will always be -1
  const Scalar two    = static_cast<Scalar>( 2.0);
  const Scalar posOne = static_cast<Scalar>(+1.0);
  const Scalar negOne = static_cast<Scalar>(-1.0);
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues( myGlobalElements[i],
        tuple<GO>(myGlobalElements[i], myGlobalElements[i]+1)(),
        tuple<Scalar> (two, posOne)()
        );
    }
    else if (static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
      A->insertGlobalValues( myGlobalElements[i],
        tuple<GO>(myGlobalElements[i]-1, myGlobalElements[i])(),
        tuple<Scalar> (negOne, two)()
        );
    }
    else {
      A->insertGlobalValues( myGlobalElements[i],
        tuple<GO>(myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1)(),
        tuple<Scalar> (negOne, two, posOne)()
        );
    }
  }

  // Finish up
  A->fillComplete();

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
// convertTpetraToThyraComm
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, convertTpetraToThyraComm,
  Scalar )
{
  RCP<const Comm<int> > tpetraComm = Teuchos::DefaultComm<int>::getComm();
  RCP<const Comm<Ordinal> > thyraComm = Thyra::convertTpetraToThyraComm(tpetraComm);
  TEST_ASSERT(nonnull(thyraComm));
}


//
// createVectorSpace
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, createVectorSpace,
  Scalar )
{
  const RCP<const TpetraMap_t> tpetraMap = createTpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<Scalar> > vs =
    Thyra::createVectorSpace<Scalar>(tpetraMap);
  TEST_ASSERT(nonnull(vs));
  out << "vs = " << *vs;
  const RCP<const SpmdVectorSpaceBase<Scalar> > vs_spmd =
    rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(vs, true);
  TEST_EQUALITY(vs_spmd->localSubDim(), g_localDim);
  TEST_EQUALITY(vs->dim(), as<Ordinal>(tpetraMap->getGlobalNumElements()));

  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;
  RCP<const TpetraMap_t> tpetraMap2 = ConverterT::getTpetraMap(vs);
  TEST_EQUALITY(tpetraMap2, tpetraMap);
}


//
// createVector
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, createVector,
  Scalar )
{

  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const RCP<const TpetraMap_t> tpetraMap = createTpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<Scalar> > vs =
    Thyra::createVectorSpace<Scalar>(tpetraMap);

  const RCP<Tpetra::Vector<Scalar> > tpetraVector =
    rcp(new Tpetra::Vector<Scalar>(tpetraMap));

  {
    const RCP<VectorBase<Scalar> > thyraVector = createVector(tpetraVector, vs);
    TEST_EQUALITY(thyraVector->space(), vs);
    const RCP<Tpetra::Vector<Scalar> > tpetraVector2 =
      ConverterT::getTpetraVector(thyraVector);
    TEST_EQUALITY(tpetraVector2, tpetraVector);
  }

  {
    const RCP<VectorBase<Scalar> > thyraVector = Thyra::createVector(tpetraVector);
    TEST_INEQUALITY(thyraVector->space(), vs);
    TEST_ASSERT(thyraVector->space()->isCompatible(*vs));
    const RCP<Tpetra::Vector<Scalar> > tpetraVector2 =
      ConverterT::getTpetraVector(thyraVector);
    TEST_EQUALITY(tpetraVector2, tpetraVector);
  }

}


//
// createConstVector
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, createConstVector,
  Scalar )
{

  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const RCP<const TpetraMap_t> tpetraMap = createTpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<Scalar> > vs =
    Thyra::createVectorSpace<Scalar>(tpetraMap);

  const RCP<const Tpetra::Vector<Scalar> > tpetraVector =
    rcp(new Tpetra::Vector<Scalar>(tpetraMap));

  {
    const RCP<const VectorBase<Scalar> > thyraVector =
      createConstVector(tpetraVector, vs);
    TEST_EQUALITY(thyraVector->space(), vs);
    const RCP<const Tpetra::Vector<Scalar> > tpetraVector2 =
      ConverterT::getConstTpetraVector(thyraVector);
    TEST_EQUALITY(tpetraVector2, tpetraVector);
  }

  {
    const RCP<const VectorBase<Scalar> > thyraVector =
      Thyra::createConstVector(tpetraVector);
    TEST_INEQUALITY(thyraVector->space(), vs);
    TEST_ASSERT(thyraVector->space()->isCompatible(*vs));
    const RCP<const Tpetra::Vector<Scalar> > tpetraVector2 =
      ConverterT::getConstTpetraVector(thyraVector);
    TEST_EQUALITY(tpetraVector2, tpetraVector);
  }

}


//
// createMultiVector
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, createMultiVector,
  Scalar )
{
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NODE;
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const int numCols = 3;

  const RCP<const TpetraMap_t> tpetraMap = createTpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<Scalar> > rangeVs =
    Thyra::createVectorSpace<Scalar>(tpetraMap);

  const RCP<const TpetraMap_t> tpetraLocRepMap =
    Tpetra::createLocalMapWithNode<LO,GO,NODE>(
      numCols, tpetraMap->getComm());
  const RCP<const VectorSpaceBase<Scalar> > domainVs =
    Thyra::createVectorSpace<Scalar>(tpetraLocRepMap);

  const RCP<Tpetra::MultiVector<Scalar> > tpetraMultiVector =
    rcp(new Tpetra::MultiVector<Scalar>(tpetraMap, numCols));

  {
    const RCP<MultiVectorBase<Scalar> > thyraMultiVector =
      createMultiVector(tpetraMultiVector, rangeVs, domainVs);
    TEST_EQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_EQUALITY(thyraMultiVector->domain(), domainVs);
    const RCP<Tpetra::MultiVector<Scalar> > tpetraMultiVector2 =
      ConverterT::getTpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(tpetraMultiVector2, tpetraMultiVector);
  }

  {
    const RCP<MultiVectorBase<Scalar> > thyraMultiVector =
      Thyra::createMultiVector(tpetraMultiVector);
    TEST_INEQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_INEQUALITY(thyraMultiVector->domain(), domainVs);
    TEST_ASSERT(thyraMultiVector->range()->isCompatible(*rangeVs));
    TEST_ASSERT(thyraMultiVector->domain()->isCompatible(*domainVs));
    const RCP<Tpetra::MultiVector<Scalar> > tpetraMultiVector2 =
      ConverterT::getTpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(tpetraMultiVector2, tpetraMultiVector);
  }

}


//
// createConstMultiVector
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, createConstMultiVector,
  Scalar )
{
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NODE;
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const int numCols = 3;

  const RCP<const TpetraMap_t> tpetraMap = createTpetraMap(g_localDim);
  const RCP<const VectorSpaceBase<Scalar> > rangeVs =
    Thyra::createVectorSpace<Scalar>(tpetraMap);

  const RCP<const TpetraMap_t> tpetraLocRepMap =
    Tpetra::createLocalMapWithNode<LO,GO,NODE>(
      numCols, tpetraMap->getComm());
  const RCP<const VectorSpaceBase<Scalar> > domainVs =
    Thyra::createVectorSpace<Scalar>(tpetraLocRepMap);

  const RCP<const Tpetra::MultiVector<Scalar> > tpetraMultiVector =
    rcp(new Tpetra::MultiVector<Scalar>(tpetraMap, numCols));

  {
    const RCP<const MultiVectorBase<Scalar> > thyraMultiVector =
      createConstMultiVector(tpetraMultiVector, rangeVs, domainVs);
    TEST_EQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_EQUALITY(thyraMultiVector->domain(), domainVs);
    const RCP<const Tpetra::MultiVector<Scalar> > tpetraMultiVector2 =
      ConverterT::getConstTpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(tpetraMultiVector2, tpetraMultiVector);
  }

  {
    const RCP<const MultiVectorBase<Scalar> > thyraMultiVector =
      Thyra::createConstMultiVector(tpetraMultiVector);
    TEST_INEQUALITY(thyraMultiVector->range(), rangeVs);
    TEST_INEQUALITY(thyraMultiVector->domain(), domainVs);
    TEST_ASSERT(thyraMultiVector->range()->isCompatible(*rangeVs));
    TEST_ASSERT(thyraMultiVector->domain()->isCompatible(*domainVs));
    const RCP<const Tpetra::MultiVector<Scalar> > tpetraMultiVector2 =
      ConverterT::getConstTpetraMultiVector(thyraMultiVector);
    TEST_EQUALITY(tpetraMultiVector2, tpetraMultiVector);
  }

}


//
// TeptraVectorSpace
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, TeptraVectorSpace,
  Scalar )
{
  const RCP<const VectorSpaceBase<Scalar> > vs =
    Thyra::createVectorSpace<Scalar>(createTpetraMap(g_localDim));
  const RCP<VectorBase<Scalar> > v = createMember(vs);
  TEST_ASSERT(nonnull(v));
  TEST_EQUALITY(v->space(), vs);
}


//
// vectorSpaceTester
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, vectorSpaceTester,
  Scalar )
{
  const RCP<const VectorSpaceBase<Scalar> > vs
    = createTpetraVectorSpace<Scalar>(g_localDim);
  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);
  TEST_ASSERT(vectorSpaceTester.check(*vs, &out));
}


//
// vectorStdOpsTester
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, vectorStdOpsTester,
  Scalar )
{
  const RCP<const VectorSpaceBase<Scalar> > vs =
    Thyra::createVectorSpace<Scalar>(createTpetraMap(g_localDim));
  Thyra::VectorStdOpsTester<Scalar> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(5.0e-13);
  vectorStdOpsTester.error_tol(5.0e-14);
  TEST_ASSERT(vectorStdOpsTester.checkStdOps(*vs, &out));
}


//
// multiVectorStdOpsTester
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, multiVectorStdOpsTester,
  Scalar )
{
  const RCP<const VectorSpaceBase<Scalar> > vs =
    Thyra::createVectorSpace<Scalar>(createTpetraMap(g_localDim));
  Thyra::MultiVectorStdOpsTester<Scalar> mvStdOpsTester;
  mvStdOpsTester.warning_tol(5.0e-13);
  mvStdOpsTester.error_tol(5.0e-14);
  TEST_ASSERT(mvStdOpsTester.checkStdOps(*vs, &out));
}


//
// getTpetraMultiVector
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, getTpetraMultiVector,
  Scalar )
{
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const int numCols = 3;
  const RCP<const VectorSpaceBase<Scalar> > vs
    = createTpetraVectorSpace<Scalar>(g_localDim);

  {
    const RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, numCols);
    const RCP<Tpetra::MultiVector<Scalar> > tmv =
      ConverterT::getTpetraMultiVector(mv);
    TEST_ASSERT(nonnull(tmv));
    TEST_EQUALITY(as<Ordinal>(tmv->getMap()->getGlobalNumElements()), vs->dim());
  }

  {
    const RCP<VectorBase<Scalar> > v = createMember(vs);
    const RCP<Tpetra::MultiVector<Scalar> > tmv =
      ConverterT::getTpetraMultiVector(v);
    TEST_ASSERT(nonnull(tmv));
    TEST_EQUALITY(as<Ordinal>(tmv->getMap()->getGlobalNumElements()), vs->dim());
  }

#ifdef THYRA_DEBUG
  const RCP<VectorBase<Scalar> > pv = Thyra::defaultProductVector<Scalar>();
  TEST_THROW(ConverterT::getTpetraMultiVector(pv), std::logic_error);
#endif

}


//
// getConstTpetraMultiVector
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, getConstTpetraMultiVector,
  Scalar )
{
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const int numCols = 3;
  const RCP<const VectorSpaceBase<Scalar> > vs
    = createTpetraVectorSpace<Scalar>(g_localDim);

  {
    const RCP<const MultiVectorBase<Scalar> > mv = createMembers(vs, numCols);
    const RCP<const Tpetra::MultiVector<Scalar> > tmv =
      ConverterT::getConstTpetraMultiVector(mv);
    TEST_ASSERT(nonnull(tmv));
    TEST_EQUALITY(as<Ordinal>(tmv->getMap()->getGlobalNumElements()), vs->dim());
  }

  {
    const RCP<const VectorBase<Scalar> > v = createMember(vs);
    const RCP<const Tpetra::MultiVector<Scalar> > tmv =
      ConverterT::getConstTpetraMultiVector(v);
    TEST_ASSERT(nonnull(tmv));
    TEST_EQUALITY(as<Ordinal>(tmv->getMap()->getGlobalNumElements()), vs->dim());
  }

#ifdef THYRA_DEBUG
  const RCP<const VectorBase<Scalar> > pv = Thyra::defaultProductVector<Scalar>();
  TEST_THROW(ConverterT::getConstTpetraMultiVector(pv), std::logic_error);
#endif

}


//
// TpetraLinearOp
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, TpetraLinearOp,
  Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::as;

  const RCP<Tpetra::Operator<Scalar> > tpetraOp =
    createTriDiagonalTpetraOperator<Scalar>(g_localDim);
  out << "tpetraOp = " << Teuchos::describe(*tpetraOp, Teuchos::VERB_HIGH) << std::endl;
  TEST_ASSERT(nonnull(tpetraOp));

  const RCP<const VectorSpaceBase<Scalar> > rangeSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getRangeMap());
  const RCP<const VectorSpaceBase<Scalar> > domainSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getDomainMap());
  const RCP<const LinearOpBase<Scalar> > thyraLinearOp =
    Thyra::tpetraLinearOp(rangeSpace, domainSpace, tpetraOp);
  TEST_ASSERT(nonnull(thyraLinearOp));

  out << "\nCheck that operator returns the right thing ...\n";
  const RCP<VectorBase<Scalar> > x = createMember(thyraLinearOp->domain());
  Thyra::V_S(x.ptr(), ST::one());
  const RCP<VectorBase<Scalar> > y = createMember(thyraLinearOp->range());
  Thyra::apply<Scalar>(*thyraLinearOp, Thyra::NOTRANS, *x, y.ptr());
  const Scalar sum_y = sum(*y);
  TEST_FLOATING_EQUALITY( sum_y, as<Scalar>(3+1+2*(y->space()->dim()-2)),
    100.0 * ST::eps() );

  out << "\nCheck the general LinearOp interface ...\n";
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.show_all_tests(showAllTests);
  linearOpTester.dump_all(dumpAll);
  if (runLinearOpTester) {
    TEST_ASSERT(linearOpTester.check(*thyraLinearOp, Teuchos::inOutArg(out)));
  }

}


//
// createLinearOp
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, createLinearOp,
  Scalar )
{

  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const RCP<Tpetra::Operator<Scalar> > tpetraOp =
    createTriDiagonalTpetraOperator<Scalar>(g_localDim);
  out << "tpetraOp = " << Teuchos::describe(*tpetraOp, Teuchos::VERB_HIGH) << std::endl;

  const RCP<const VectorSpaceBase<Scalar> > rangeSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getRangeMap());

  const RCP<const VectorSpaceBase<Scalar> > domainSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getDomainMap());

  {
    const RCP<LinearOpBase<Scalar> > thyraOp =
      createLinearOp(tpetraOp, rangeSpace, domainSpace);
    TEST_EQUALITY(thyraOp->range(), rangeSpace);
    TEST_EQUALITY(thyraOp->domain(), domainSpace);
    const RCP<Tpetra::Operator<Scalar> > tpetraOp2 =
      ConverterT::getTpetraOperator(thyraOp);
    TEST_EQUALITY(tpetraOp2, tpetraOp);
  }

  {
    const RCP<LinearOpBase<Scalar> > thyraOp =
      Thyra::createLinearOp(tpetraOp);
    TEST_INEQUALITY(thyraOp->range(), rangeSpace);
    TEST_INEQUALITY(thyraOp->domain(), domainSpace);
    TEST_ASSERT(thyraOp->range()->isCompatible(*rangeSpace));
    TEST_ASSERT(thyraOp->domain()->isCompatible(*domainSpace));
    const RCP<Tpetra::Operator<Scalar> > tpetraOp2 =
      ConverterT::getTpetraOperator(thyraOp);
    TEST_EQUALITY(tpetraOp2, tpetraOp);
  }

}


//
// createConstLinearOp
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, createConstLinearOp,
  Scalar )
{

  typedef Thyra::TpetraOperatorVectorExtraction<Scalar> ConverterT;

  const RCP<const Tpetra::Operator<Scalar> > tpetraOp =
    createTriDiagonalTpetraOperator<Scalar>(g_localDim);
  out << "tpetraOp = " << Teuchos::describe(*tpetraOp, Teuchos::VERB_HIGH) << std::endl;

  const RCP<const VectorSpaceBase<Scalar> > rangeSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getRangeMap());

  const RCP<const VectorSpaceBase<Scalar> > domainSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getDomainMap());

  {
    const RCP<const LinearOpBase<Scalar> > thyraOp =
      createConstLinearOp(tpetraOp, rangeSpace, domainSpace);
    TEST_EQUALITY(thyraOp->range(), rangeSpace);
    TEST_EQUALITY(thyraOp->domain(), domainSpace);
    const RCP<const Tpetra::Operator<Scalar> > tpetraOp2 =
      ConverterT::getConstTpetraOperator(thyraOp);
    TEST_EQUALITY(tpetraOp2, tpetraOp);
  }

  {
    const RCP<const LinearOpBase<Scalar> > thyraOp =
      Thyra::createConstLinearOp(tpetraOp);
    TEST_INEQUALITY(thyraOp->range(), rangeSpace);
    TEST_INEQUALITY(thyraOp->domain(), domainSpace);
    TEST_ASSERT(thyraOp->range()->isCompatible(*rangeSpace));
    TEST_ASSERT(thyraOp->domain()->isCompatible(*domainSpace));
    const RCP<const Tpetra::Operator<Scalar> > tpetraOp2 =
      ConverterT::getConstTpetraOperator(thyraOp);
    TEST_EQUALITY(tpetraOp2, tpetraOp);
  }

}


//
// Tpetra-implemented methods
//


Teuchos::RCP<Teuchos::Time> lookupAndAssertTimer(const std::string &label)
{
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter(label);
  TEUCHOS_TEST_FOR_EXCEPTION(timer == null,
    std::runtime_error,
    "lookupAndAssertTimer(): timer \"" << label << "\" was not present in Teuchos::TimeMonitor."
    " Unit test not valid.");
  return timer;
}


#define CHECK_TPETRA_FUNC_CALL_INCREMENT( timerStr, tpetraCode, thyraCode ) \
{ \
  out << "\nTesting that Thyra calls down to " << timerStr << "\n"; \
  ECHO(tpetraCode); \
  const RCP<const Time> timer = lookupAndAssertTimer(timerStr); \
  const int countBefore = timer->numCalls();  \
  ECHO(thyraCode); \
  const int countAfter = timer->numCalls(); \
  TEST_EQUALITY( countAfter, countBefore+1 ); \
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, UseTpetraImplementations,
  Scalar )
{
  using Teuchos::Time;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType Magnitude;
  typedef VectorSpaceBase<Scalar> VectorSpace;
  typedef MultiVectorBase<Scalar> MultiVec;
  //typedef Tpetra::Map<int> TpetraMap;
  typedef Tpetra::MultiVector<Scalar> TpetraMultiVec;
  typedef TpetraOperatorVectorExtraction<Scalar> TOVE;

  const int numCols = 3;

  const RCP<const VectorSpace> vs =
    createTpetraVectorSpace<Scalar>(g_localDim);
  const RCP<MultiVec>
    A = createMembers(vs, numCols),
    B = createMembers(vs, numCols);
  const RCP<TpetraMultiVec>
    tA = TOVE::getTpetraMultiVector(A),
    tB = TOVE::getTpetraMultiVector(B);
  Array<Scalar> C(numCols*numCols,ST::zero());

  Teuchos::Array<Magnitude> avMag(numCols);
  Teuchos::Array<Scalar> avScal(numCols);

  CHECK_TPETRA_FUNC_CALL_INCREMENT(
    "Tpetra::MultiVector::putScalar()",
    tA->putScalar(ST::zero()),
    Thyra::assign(A.ptr(), ST::zero())
    );

  CHECK_TPETRA_FUNC_CALL_INCREMENT(
    "Tpetra::MultiVector::dot()",
    tA->dot(*tB, avScal() ),      // norms calls scalarProd, which calls Tpetra::dot
    Thyra::norms( *A, avMag() )
    );

  CHECK_TPETRA_FUNC_CALL_INCREMENT(
    "Tpetra::MultiVector::dot()",
    tA->dot(*tB, avScal() ),
    A->range()->scalarProds(*A, *B, avScal() )
    );

  /*
  CHECK_TPETRA_FUNC_CALL_INCREMENT(
    "Tpetra::MultiVector::scale(alpha)",
    tB->scale( ST::zero() ),
    Thyra::scale( ST::zero(), B.ptr() )
    );
  */

  /*
  CHECK_TPETRA_FUNC_CALL_INCREMENT(
    "Tpetra::MultiVector::operator=()",
    (*tA) = *tB,
    Thyra::assign( A.ptr(), *B )
    );
  */

  /*
  {
    RCP<MultiVec> Ctmvb = Thyra::createMembersView(
        A->domain(),
        RTOpPack::SubMultiVectorView<Scalar>(
          0, numCols, 0, numCols,
          Teuchos::arcpFromArrayView(C()), numCols
        )
      );
    CHECK_TPETRA_FUNC_CALL_INCREMENT(
      "Tpetra::multiplyOntoHost()",
      Tpetra::multiplyOntoHost(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ST::one(),*tA,*tB, C()),
      A->apply(Thyra::CONJTRANS,*B,Ctmvb.ptr(),ST::one(),ST::zero())
      );
  }
  */

  /*
  {
    RCP<const MultiVec> Ctmvb = Thyra::createMembersView(
        A->domain(),
        RTOpPack::ConstSubMultiVectorView<Scalar>(
          0, numCols, 0, numCols,
          Teuchos::arcpFromArrayView(C()), numCols
        )
      );
    const RCP<const TpetraMultiVec>
      tC = TOVE::getConstTpetraMultiVector(Ctmvb);
    CHECK_TPETRA_FUNC_CALL_INCREMENT(
      "Tpetra::MultiVector::multiply()",
      tB->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ST::one(),*tA,*tC,ST::zero()),
      A->apply(Thyra::NOTRANS,*Ctmvb,B.ptr(),ST::one(),ST::zero())
    );
  }
  */

/*
  RCP<Time>
    timerUpdate = lookupAndAssertTimer("Tpetra::MultiVector::update(alpha,A,beta,B,gamma)"),
    timerUpdate2   = lookupAndAssertTimer("Tpetra::MultiVector::update(alpha,A,beta)"),

  // TODO: test update(two vector)
  // TODO: test update(three vector)
*/
}


#ifdef TPETRA_TEUCHOS_TIME_MONITOR
#  define TPETRA_TIMER_TESTS(SCALAR)  \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, UseTpetraImplementations, SCALAR )
#else
#  define TPETRA_TIMER_TESTS(SCALAR)
#endif


//
// TpetraLinearOp_EpetraRowMatrix
//


#ifdef HAVE_THYRA_TPETRA_EPETRA


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, TpetraLinearOp_EpetraRowMatrix,
  Scalar )
{

  using Teuchos::as;
  using Teuchos::outArg;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::Array;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  const RCP<Tpetra::Operator<Scalar> > tpetraOp =
    createTriDiagonalTpetraOperator<Scalar>(g_localDim);

  const RCP<LinearOpBase<Scalar> > thyraOp =
    Thyra::createLinearOp(tpetraOp);

  const RCP<Thyra::TpetraLinearOp<Scalar> > thyraTpetraOp =
    Teuchos::rcp_dynamic_cast<Thyra::TpetraLinearOp<Scalar> >(thyraOp);

  RCP<const Epetra_Operator> epetraOp;
  Thyra::EOpTransp epetraOpTransp;
  Thyra::EApplyEpetraOpAs epetraOpApplyAs;
  Thyra::EAdjointEpetraOp epetraOpAdjointSupport;

  thyraTpetraOp->getEpetraOpView( outArg(epetraOp), outArg(epetraOpTransp),
    outArg(epetraOpApplyAs), outArg(epetraOpAdjointSupport) );

  if (typeid(Scalar) == typeid(double)) {
    TEST_ASSERT(nonnull(epetraOp));
    const RCP<const Epetra_RowMatrix> epetraRowMatrix =
      rcp_dynamic_cast<const Epetra_RowMatrix>(epetraOp, true);
    int numRowEntries = -1;
    epetraRowMatrix->NumMyRowEntries(1, numRowEntries);
    TEST_EQUALITY_CONST(numRowEntries, 3);
    Array<double> row_values(numRowEntries);
    Array<int> row_indices(numRowEntries);
    epetraRowMatrix->ExtractMyRowCopy(1, numRowEntries, numRowEntries,
      row_values.getRawPtr(), row_indices.getRawPtr());
    TEST_EQUALITY_CONST(row_values[0], -1.0);
    TEST_EQUALITY_CONST(row_values[1], 2.0);
    TEST_EQUALITY_CONST(row_values[2], 1.0);
    // ToDo: Test column indices!
  }
  else {
    TEST_ASSERT(is_null(epetraOp));
  }

}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, TpetraLinearOp_EpetraRowMatrix,
  Scalar )
{
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, TpetraLinearOp_RowStatLinearOpBase,
  Scalar )
{
  using Teuchos::as;

  const RCP<Tpetra::Operator<Scalar> > tpetraOp =
    createTriDiagonalTpetraOperator<Scalar>(g_localDim);
  out << "tpetraOp = " << Teuchos::describe(*tpetraOp, Teuchos::VERB_HIGH) << std::endl;
  TEST_ASSERT(nonnull(tpetraOp));

  const RCP<Tpetra::CrsMatrix<Scalar> > tpetraCrsMatrix =
    Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar> >(tpetraOp,true);

  const RCP<const VectorSpaceBase<Scalar> > rangeSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getRangeMap());
  const RCP<const VectorSpaceBase<Scalar> > domainSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getDomainMap());
  const RCP<LinearOpBase<Scalar> > thyraLinearOp =
    Thyra::tpetraLinearOp(rangeSpace, domainSpace, tpetraOp);
  TEST_ASSERT(nonnull(thyraLinearOp));

  const Teuchos::RCP<Thyra::RowStatLinearOpBase<Scalar> > rowStatOp =
    Teuchos::rcp_dynamic_cast<Thyra::RowStatLinearOpBase<Scalar> >(thyraLinearOp, true);

  // Get the inverse row sums

  const RCP<VectorBase<Scalar> > inv_row_sums =
    createMember<Scalar>(thyraLinearOp->range());
  const RCP<VectorBase<Scalar> > row_sums =
    createMember<Scalar>(thyraLinearOp->range());

  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());
  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  out << "inv_row_sums = " << *inv_row_sums;
  out << "row_sums = " << *row_sums;

  TEST_FLOATING_EQUALITY(
    Thyra::sum<Scalar>(*row_sums),
    Teuchos::as<Scalar>(4.0 * thyraLinearOp->domain()->dim() - 2.0),
    Teuchos::as<Scalar>(10.0 * Teuchos::ScalarTraits<Scalar>::eps())
    );

  TEST_FLOATING_EQUALITY(
    Thyra::sum<Scalar>(*inv_row_sums),
    Teuchos::as<Scalar>( 1.0 / 4.0 * (thyraLinearOp->domain()->dim() - 2) + 2.0 / 3.0 ),
    Teuchos::as<Scalar>(10.0 * Teuchos::ScalarTraits<Scalar>::eps())
    );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraThyraWrappers, TpetraLinearOp_ScaledLinearOpBase,
  Scalar )
{
  const RCP<Tpetra::Operator<Scalar> > tpetraOp =
    createTriDiagonalTpetraOperator<Scalar>(g_localDim);
  out << "tpetraOp = " << Teuchos::describe(*tpetraOp, Teuchos::VERB_HIGH) << std::endl;
  TEST_ASSERT(nonnull(tpetraOp));

  const RCP<Tpetra::CrsMatrix<Scalar> > tpetraCrsMatrix =
    Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar> >(tpetraOp,true);

  const RCP<const VectorSpaceBase<Scalar> > rangeSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getRangeMap());
  const RCP<const VectorSpaceBase<Scalar> > domainSpace =
    Thyra::createVectorSpace<Scalar>(tpetraOp->getDomainMap());
  const RCP<LinearOpBase<Scalar> > thyraLinearOp =
    Thyra::tpetraLinearOp(rangeSpace, domainSpace, tpetraOp);
  TEST_ASSERT(nonnull(thyraLinearOp));

  const Teuchos::RCP<Thyra::RowStatLinearOpBase<Scalar> > rowStatOp =
    Teuchos::rcp_dynamic_cast<Thyra::RowStatLinearOpBase<Scalar> >(thyraLinearOp, true);

  // Get the inverse row sums

  const RCP<VectorBase<Scalar> > inv_row_sums =
    createMember<Scalar>(thyraLinearOp->range());
  const RCP<VectorBase<Scalar> > row_sums =
    createMember<Scalar>(thyraLinearOp->range());

  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM,
    inv_row_sums.ptr());
  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  out << "inv_row_sums = " << *inv_row_sums;
  out << "row_sums = " << *row_sums;

  const Teuchos::RCP<Thyra::ScaledLinearOpBase<Scalar> > scaledOp =
    Teuchos::rcp_dynamic_cast<Thyra::ScaledLinearOpBase<Scalar> >(thyraLinearOp, true);

  TEUCHOS_ASSERT(scaledOp->supportsScaleLeft());

  scaledOp->scaleLeft(*inv_row_sums);

  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,
    row_sums.ptr());

  out << "row_sums after left scaling by inv_row_sum = " << *row_sums;

  // scaled row sums should be one for each entry
  TEST_FLOATING_EQUALITY(
    Scalar(row_sums->space()->dim()),
    Thyra::sum<Scalar>(*row_sums),
    as<Scalar>(10.0 * Teuchos::ScalarTraits<Scalar>::eps())
    );

  // Don't currently check the results of right scaling.  Tpetra tests
  // already check this.  Once column sums are supported in tpetra
  // adapters, this can be checked easily.
  TEUCHOS_ASSERT(scaledOp->supportsScaleRight());
  scaledOp->scaleRight(*inv_row_sums);
  rowStatOp->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM,row_sums.ptr());
  out << "row_sums after right scaling by inv_row_sum = " << *row_sums;
}

#endif // HAVE_THYRA_TPETRA_EPETRA


//
// Unit test instantiations
//

#define THYRA_TPETRA_THYRA_WRAPPERS_INSTANT(SCALAR) \
 \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers,  \
    convertTpetraToThyraComm, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    createVectorSpace, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    createVector, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    createConstVector, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers,  \
    createMultiVector, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    createConstMultiVector, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    TeptraVectorSpace, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    vectorSpaceTester, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    vectorStdOpsTester, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    multiVectorStdOpsTester, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    getTpetraMultiVector, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    getConstTpetraMultiVector, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    TpetraLinearOp, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    createLinearOp, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    createConstLinearOp, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    TpetraLinearOp_EpetraRowMatrix, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    TpetraLinearOp_RowStatLinearOpBase, SCALAR ) \
   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraThyraWrappers, \
    TpetraLinearOp_ScaledLinearOpBase, SCALAR ) \


// We can currently only explicitly instantiate with double support because
// Tpetra only supports explicit instantaition with double.  As for implicit
// instantation, g++ 3.4.6 on my Linux machine was taking more than 30 minutes
// to compile this file when all of the types double, float, complex<double>,
// and complex<float> where enabled.  Therefore, we will only test double for
// now until explicit instantation with other types are supported by Tpetra.

THYRA_TPETRA_THYRA_WRAPPERS_INSTANT(double)


} // namespace Thyra
