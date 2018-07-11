// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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

} // namespace Teuchos



