// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"

#define SS_ECHO(ops)      \
  {                       \
    std::stringstream ss; \
    ss << ops;            \
    ECHO(ss.str())        \
  };

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

// Build PreconditionerFactory testing rig
///////////////////////////////////////////////////////////

class TestFactory : public Teko::PreconditionerFactory {
 public:
  Teko::LinearOp buildPreconditionerOperator(Teko::LinearOp& lo,
                                             Teko::PreconditionerState& state) const;

  mutable double timestep_;
  mutable Teko::LinearOp pcdOp_;
  mutable std::string string_;
};

Teko::LinearOp TestFactory::buildPreconditionerOperator(Teko::LinearOp& lo,
                                                        Teko::PreconditionerState& state) const {
  Teko::RequestMesg pcdMesg("PCD Op"), tsMesg("timestep"), strMesg("name");

  pcdOp_    = callbackHandler_->request<Teko::LinearOp>(pcdMesg);
  timestep_ = callbackHandler_->request<double>(tsMesg);
  string_   = callbackHandler_->request<std::string>(strMesg);

  return Teuchos::null;
}

// test rigs callbacks
///////////////////////////////////////////////////////////

class PCDCallback : public Teko::RequestCallback<Teko::LinearOp> {
 public:
  Teko::LinearOp request(const Teko::RequestMesg& rd) {
    TEUCHOS_ASSERT(handlesRequest(rd));
    return Teuchos::null;
  }

  bool handlesRequest(const Teko::RequestMesg& rd) {
    if (rd.getName() == "PCD Op") return true;
    return false;
  }

  void preRequest(const Teko::RequestMesg& rm) {}
};

class TSCallback : public Teko::RequestCallback<double> {
 public:
  double request(const Teko::RequestMesg& rd) {
    TEUCHOS_ASSERT(handlesRequest(rd));
    return 0.1;
  }

  bool handlesRequest(const Teko::RequestMesg& rd) {
    if (rd.getName() == "timestep") return true;
    return false;
  }

  void preRequest(const Teko::RequestMesg& rm) {}
};

class StringCallback : public Teko::RequestCallback<std::string> {
 public:
  std::string request(const Teko::RequestMesg& rd) {
    TEUCHOS_ASSERT(handlesRequest(rd));
    return "the string";
  }

  bool handlesRequest(const Teko::RequestMesg& rd) {
    if (rd.getName() == "name") return true;
    return false;
  }

  void preRequest(const Teko::RequestMesg& rm) {}
};

///////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST(tRequestInterface, test_request_interface) {
  RCP<TestFactory> precFact            = rcp(new TestFactory);
  RCP<Teko::PreconditionerState> state = precFact->buildPreconditionerState();

  RCP<Teko::RequestHandler> rh = rcp(new Teko::RequestHandler);

  rh->addRequestCallback(rcp(new PCDCallback));
  rh->addRequestCallback(rcp(new TSCallback));
  rh->addRequestCallback(rcp(new StringCallback));

  {
    Teko::RequestMesg pcdMesg("PCD Op"), tsMesg("timestep"), strMesg("name");
    rh->preRequest<Teko::LinearOp>(pcdMesg);
    rh->preRequest<double>(tsMesg);
    rh->preRequest<std::string>(strMesg);
  }
  precFact->setRequestHandler(rh);

  Teko::LinearOp A;
  Teko::LinearOp result = precFact->buildPreconditionerOperator(A, *state);

  TEST_ASSERT(precFact->timestep_ == 0.1);
  TEST_ASSERT(precFact->pcdOp_ == Teuchos::null);
  TEST_ASSERT(precFact->string_ == "the string");

  try {
    Teko::RequestMesg intMesg("Int Op");
    rh->preRequest<int>(intMesg);
    out << "Found <int> with name \"Int Op\" in preRequest" << std::endl;
    TEST_ASSERT(false);
  } catch (std::exception& e) {
    out << "expected exception = " << e.what() << std::endl;
  }

  try {
    Teko::RequestMesg intMesg("Int Op");
    int size = rh->request<int>(intMesg);
    out << "Found <int> with name \"Int Op\" value=" << size << std::endl;
    TEST_ASSERT(false);
  } catch (std::exception& e) {
    out << "expected exception = " << e.what() << std::endl;
  }

  try {
    Teko::RequestMesg intMesg("PCD Op");
    int lo = rh->request<int>(intMesg);
    out << "Found <int> with name \"PCD Op\" value=" << lo << std::endl;
    TEST_ASSERT(false);
  } catch (std::exception& e) {
    out << "expected exception = " << e.what() << std::endl;
  }
}

// Test widget for the parameter list based call back:
//    used in preconditioner_request_interface unit test
///////////////////////////////////////////////////////////////////
class PLCallback : public Teko::RequestCallback<Teuchos::RCP<Teuchos::ParameterList> > {
 public:
  Teuchos::RCP<Teuchos::ParameterList> request(const Teko::RequestMesg& rm);

  bool handlesRequest(const Teko::RequestMesg& rm);

  void preRequest(const Teko::RequestMesg& rm);
};

Teuchos::RCP<Teuchos::ParameterList> PLCallback::request(const Teko::RequestMesg& rm) {
  Teuchos::RCP<const Teuchos::ParameterList> inputPL = rm.getParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(inputPL == Teuchos::null, std::runtime_error,
                             "Parameter list not included in request message");

  Teuchos::RCP<Teuchos::ParameterList> outputPL = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::ParameterList::ConstIterator itr;
  for (itr = inputPL->begin(); itr != inputPL->end(); ++itr) {
    if (itr->first == "cat")
      outputPL->set<int>("fact: level-of-fill", 7);
    else
      outputPL->setEntry(itr->first, itr->second);
  }

  return outputPL;
}

void PLCallback::preRequest(const Teko::RequestMesg& rm) {}

bool PLCallback::handlesRequest(const Teko::RequestMesg& rm) {
  if (rm.getName() == "Parameter List")
    return true;
  else
    return false;
}

TEUCHOS_UNIT_TEST(tRequestInterface, preconditioner_request_interface) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::ParameterList pl;
  {  // Ifpack2-Test requires a handler
    Teuchos::ParameterList& ifpack2List = pl.sublist("Ifpack2-Test");
    ifpack2List.set<std::string>("Type", "Ifpack2");
    ifpack2List.sublist("Ifpack2 Settings");
    ifpack2List.sublist("Required Parameters").set<std::string>("cat", "dog");
  }
  {  // Ifpack2-Test2 does not require a handler
    Teuchos::ParameterList& ifpack2List = pl.sublist("Ifpack2-Test2");
    ifpack2List.set<std::string>("Type", "Ifpack2");
    ifpack2List.sublist("Ifpack2 Settings").set<int>("fact: level-of-fill", 3);
  }

  {
    RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(pl);
    TEST_NOTHROW(invLib->getInverseFactory("Ifpack2-Test2"));
    TEST_THROW(invLib->getInverseFactory("Ifpack2-Test"), std::runtime_error);
  }

  {
    RCP<Teko::RequestHandler> reqHandler = rcp(new Teko::RequestHandler);
    reqHandler->addRequestCallback(rcp(new PLCallback));

    RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(pl);
    invLib->setRequestHandler(reqHandler);

    RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory("Ifpack2-Test");

    RCP<Teko::PreconditionerInverseFactory> pInvFact =
        Teuchos::rcp_dynamic_cast<Teko::PreconditionerInverseFactory>(invFact, true);
    Teuchos::RCP<const Teuchos::ParameterList> pl2 = pInvFact->getPrecFactory()->getParameterList();

    TEST_ASSERT(pl2->sublist("Ifpack2 Settings").isParameter("fact: level-of-fill"));
    TEST_EQUALITY(pl2->sublist("Ifpack2 Settings").get<int>("fact: level-of-fill"), 7);
  }
}
