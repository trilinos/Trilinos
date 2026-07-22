// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RawParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos {


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, xmlUpdateAndBroadcast ) {
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  // Test the broadcast functionality to avoid unscalable I/O collisions
  std::string inputFile="input.xml";
  ParameterList A;
  ParameterList B;
  updateParametersFromXmlFile(inputFile, outArg(A));
  updateParametersFromXmlFileAndBroadcast(inputFile, outArg(B), *comm);
  out << "B = " << B;
  TEST_ASSERT( B.begin() != B.end() ); // Avoid false positive from empty lists

  // See if any process returned a failed (i.e. a non-zero local_failed)
  const int local_failed = !(A == B);
  int global_failed = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, local_failed, outArg(global_failed));
  TEST_EQUALITY_CONST(global_failed, 0);
}



TEUCHOS_UNIT_TEST( Teuchos_ParameterList, xmlUpdateAndBroadcastNoOverWrite ) {
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  // Test the broadcast functionality to avoid unscalable I/O collisions
  std::string inputFile="input.xml";
  ParameterList A;
  A.set("Transient",true);
  A.set("Phalanx Graph Visualization Detail",2);

  ParameterList B;
  B.set("Transient",true);
  B.set("Phalanx Graph Visualization Detail",2);


  updateParametersFromXmlFile(inputFile, outArg(A));
  updateParametersFromXmlFileAndBroadcast(inputFile, outArg(B), *comm);
  out << "B = " << B;
  TEST_ASSERT( B.begin() != B.end() ); // Avoid false positive from empty lists

  // See if any process returned a failed (i.e. a non-zero local_failed)
  const int local_failed = !(A == B);
  int global_failed = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, local_failed, outArg(global_failed));
  TEST_EQUALITY_CONST(global_failed, 0);

  // Check the assigned values.
  TEST_EQUALITY( B.get("Transient",false), true)
  TEST_EQUALITY_CONST( B.get("Phalanx Graph Visualization Detail",1), 2)

}

TEUCHOS_UNIT_TEST( Teuchos_ParameterList, rawUpdateAndBroadcast ) {
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  // Test the broadcast functionality to avoid unscalable I/O collisions
  std::string inputFile="input.xml";
  ParameterList A;
  ParameterList B;

  // Set root to something != 0, if we can
  int root = 0;
  if(comm->getSize() > 1) root = 1;

  updateParametersFromXmlFile(inputFile, outArg(A));
  updateParametersAndBroadcast(outArg(A),outArg(B),*comm,root);
  out << "B = " << B;
  TEST_ASSERT( B.begin() != B.end() ); // Avoid false positive from empty lists

  // See if any process returned a failed (i.e. a non-zero local_failed)
  const int local_failed = !(A == B);
  int global_failed = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, local_failed, outArg(global_failed));
  TEST_EQUALITY_CONST(global_failed, 0);
}


TEUCHOS_UNIT_TEST( Teuchos_ParameterList, rawUpdateAndBroadcastNoOverWrite ) {
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  // Test the broadcast functionality to avoid unscalable I/O collisions
  std::string inputFile="input.xml";
  ParameterList A;
  A.set("Transient",true);
  A.set("Phalanx Graph Visualization Detail",2);

  ParameterList B;
  B.set("Transient",true);
  B.set("Phalanx Graph Visualization Detail",2);

  // Set root to something != 0, if we can
  int root = 0;
  if(comm->getSize() > 1) root = 1;

  updateParametersFromXmlFile(inputFile, outArg(A));
  updateParametersAndBroadcast(outArg(A),outArg(B),*comm,root);
  out << "B = " << B;
  TEST_ASSERT( B.begin() != B.end() ); // Avoid false positive from empty lists

  // See if any process returned a failed (i.e. a non-zero local_failed)
  const int local_failed = !(A == B);
  int global_failed = -1;
  reduceAll(*comm, Teuchos::REDUCE_SUM, local_failed, outArg(global_failed));
  TEST_EQUALITY_CONST(global_failed, 0);

  // Check the assigned values.
  TEST_EQUALITY( B.get("Transient",false), true)
  TEST_EQUALITY_CONST( B.get("Phalanx Graph Visualization Detail",1), 2)

}

TEUCHOS_UNIT_TEST( Teuchos_ParameterList, xmlUpdateAndBroadcastSpecialChars ) {
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  // Test that the special characters '&' and '<' are correctly translated from '&amp;' and '&lt;'
  std::string inputFile="inputSpecialChars.xml";
  ParameterList A;
  updateParametersFromXmlFileAndBroadcast(inputFile, outArg(A), *comm);
  out << "A = " << A;
  TEST_ASSERT( A.begin() != A.end() ); // Avoid false positive from empty lists

  TEST_EQUALITY( A.get<std::string>("sigma"), "if (x >= 0.0 && y >= 0.0 && x <= 0.5 && y <= 0.5)" );
}

} // namespace Teuchos



