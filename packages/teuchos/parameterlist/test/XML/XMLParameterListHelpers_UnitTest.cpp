// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <sstream>
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Teuchos {


TEUCHOS_UNIT_TEST( XMLParameterListHelpers, toFromFile )
{

  ParameterList A;
  A.sublist("SublistA").set("param_a", "a");
  A.sublist("SublistA").sublist("SublistB").set("param_b", "b");
  out << "\nA:\n"; A.print(out);
  writeParameterListToXmlFile(A, "A.xml");

  ParameterList B;
  updateParametersFromXmlFile("A.xml", inoutArg(B));
  out << "\nB:\n"; B.print(out);
  TEST_ASSERT( A == B );


  TEST_EQUALITY(
    getParameter<std::string>(B.sublist("SublistA", true), "param_a"),
    "a");
  TEST_EQUALITY(
    getParameter<std::string>(B.sublist("SublistA", true).sublist("SublistB", true), "param_b"),
    "b");

}

TEUCHOS_UNIT_TEST( XMLParameterListHelpers, OverwriteTest )
{
  ParameterList A;  
  A.set("conflicting param","a");
  A.sublist("SublistA").set("param_a", "a");
  A.sublist("SublistA").sublist("SublistB").set("param_b", "b");
  out << "\nA:\n"; A.print(out);
  std::stringstream Astream;  
  writeParameterListToXmlOStream(A,Astream);

  // Overwrite
  ParameterList B;
  B.set("conflicting param","b");
  updateParametersFromXmlString(Astream.str(), inoutArg(B),true);
  out << "\nB:\n"; B.print(out);
  TEST_ASSERT( A == B );
  TEST_EQUALITY(getParameter<std::string>(B.sublist("SublistA", true), "param_a"),"a");
  TEST_EQUALITY(getParameter<std::string>(B.sublist("SublistA", true).sublist("SublistB", true), "param_b"),"b");


  // No Overwrite
  ParameterList C;
  C.set("conflicting param","c");
  updateParametersFromXmlString(Astream.str(), inoutArg(B),false);
  out << "\nC:\n"; C.print(out);
  TEST_ASSERT( C.get("conflicting param","x") == "c" )
  TEST_EQUALITY(getParameter<std::string>(B.sublist("SublistA", true), "param_a"),"a");
  TEST_EQUALITY(getParameter<std::string>(B.sublist("SublistA", true).sublist("SublistB", true), "param_b"),"b");

}



} // namespace Teuchos



