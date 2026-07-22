// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tNeumannSeries.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_NeumannSeriesPreconditionerFactory.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

static RCP<Teuchos::ParameterList> buildLibPL(int numTerms, std::string scalingType) {
  RCP<Teuchos::ParameterList> pl    = rcp(new Teuchos::ParameterList());
  RCP<Teuchos::ParameterList> nList = rcp(new Teuchos::ParameterList());

  nList->set("Type", "Neumann Series");
  nList->set("Number of Terms", numTerms);
  nList->set("Scaling Type", scalingType);

  pl->set("Neumann", *nList);

  return pl;
}

RCP<Thyra::LinearOpBase<double> > buildExampleOp(int type, const Epetra_Comm& comm) {
  Epetra_Map map(3, 0, comm);
  RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, map, 3));

  double values[3];
  int indices[3] = {0, 1, 2};

  switch (type) {
    case 2:
      values[0] = 5.0;
      values[1] = 0.0;
      values[2] = 0.0;
      mat->InsertGlobalValues(0, 3, values, indices);

      values[0] = 2.0;
      values[1] = 6.0;
      values[2] = 0.0;
      mat->InsertGlobalValues(1, 3, values, indices);

      values[0] = 0.0;
      values[1] = 3.0;
      values[2] = 7.0;
      mat->InsertGlobalValues(2, 3, values, indices);
      break;
    case 1:
      values[0] = 1.0;
      values[1] = 0.0;
      values[2] = 0.0;
      mat->InsertGlobalValues(0, 3, values, indices);

      values[0] = 2.0;
      values[1] = 1.0;
      values[2] = 0.0;
      mat->InsertGlobalValues(1, 3, values, indices);

      values[0] = 0.0;
      values[1] = 3.0;
      values[2] = 1.0;
      mat->InsertGlobalValues(2, 3, values, indices);
    default: break;
  }

  mat->FillComplete();

  return Thyra::nonconstEpetraLinearOp(mat);
}

void tNeumannSeries::initializeTest() { tolerance_ = 1.0e-15; }

int tNeumannSeries::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                            int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tNeumannSeries";

  status = test_simpleOp(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"test_simpleOp\" ... PASSED", "   \"test_simpleOp\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_scaledOp(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"test_scaledOp\" ... PASSED", "   \"test_scaledOp\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tNeumannSeries...PASSED", "tNeumannSeries...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tNeumannSeries...FAILED");
  }

  return failcount;
}

bool tNeumannSeries::test_simpleOp(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  // perform actual test
  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  std::stringstream ss;
  Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");

  {
    // build library parameter list
    Teuchos::RCP<Teuchos::ParameterList> pl = buildLibPL(1, "None");
    if (verbosity >= 10) {
      os << "   tNeumannSeries::test_simpleOp :"
         << " printing library parameter list" << std::endl;
      pl->print(os);
    }

    RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(*pl);
    RCP<Teko::InverseFactory> neumann = invLib->getInverseFactory("Neumann");
    RCP<Teko::InverseFactory> direct  = invLib->getInverseFactory("Amesos");

    Teko::LinearOp op = buildExampleOp(1, *GetComm());

    Teko::LinearOp neuInv = Teko::buildInverse(*neumann, op);
    Teko::LinearOp dirInv = Teko::buildInverse(*direct, op);

    const bool result = tester.compare(*neuInv, *dirInv, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(not result,
                std::endl
                    << "   tNeumannSeries::test_simpleOp "
                    << ": Comparing underresolved factory generated operator to correct operator");
    if (result || verbosity >= 10) os << ss.str();
  }

  {
    // build library parameter list
    Teuchos::RCP<Teuchos::ParameterList> pl = buildLibPL(2, "None");
    if (verbosity >= 10) {
      os << "   tNeumannSeries::test_simpleOp :"
         << " printing library parameter list" << std::endl;
      pl->print(os);
    }

    RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(*pl);
    RCP<Teko::InverseFactory> neumann = invLib->getInverseFactory("Neumann");
    RCP<Teko::InverseFactory> direct  = invLib->getInverseFactory("Amesos");

    Teko::LinearOp op = buildExampleOp(1, *GetComm());

    Teko::LinearOp neuInv = Teko::buildInverse(*neumann, op);
    Teko::LinearOp dirInv = Teko::buildInverse(*direct, op);

    const bool result = tester.compare(*neuInv, *dirInv, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tNeumannSeries::test_simpleOp "
                            << ": Comparing factory generated operator to correct operator");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tNeumannSeries::test_scaledOp(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  // perform actual test
  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  std::stringstream ss;
  Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");

  {
    // build library parameter list
    Teuchos::RCP<Teuchos::ParameterList> pl = buildLibPL(2, "Diagonal");
    if (verbosity >= 10) {
      os << "   tNeumannSeries::test_scaledOp :"
         << " printing library parameter list" << std::endl;
      pl->print(os);
    }

    RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(*pl);
    RCP<Teko::InverseFactory> neumann = invLib->getInverseFactory("Neumann");
    RCP<Teko::InverseFactory> direct  = invLib->getInverseFactory("Amesos");

    Teko::LinearOp op = buildExampleOp(2, *GetComm());

    Teko::LinearOp neuInv = Teko::buildInverse(*neumann, op);
    Teko::LinearOp dirInv = Teko::buildInverse(*direct, op);

    const bool result = tester.compare(*neuInv, *dirInv, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tNeumannSeries::test_scaledOp "
                            << ": Comparing factory generated operator to correct operator");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
