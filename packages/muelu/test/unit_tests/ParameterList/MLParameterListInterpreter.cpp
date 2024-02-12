// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Galeri_XpetraUtils.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_MLParameterListInterpreter.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MLParameterListInterpreter, SetParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  if (!TYPE_EQUAL(SC, double)) {
    out << "Skipping for SC != double" << std::endl;
    return;
  }
  out << "version: " << MueLu::Version() << std::endl;

  // TODO: this test can be done at compilation time
#if !defined(HAVE_MUELU_EPETRA) or !defined(HAVE_MUELU_IFPACK) or !defined(HAVE_MUELU_AMESOS)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Epetra, Ifpack, Amesos");
#endif

#if !defined(HAVE_MUELU_IFPACK2) or !defined(HAVE_MUELU_AMESOS2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Tpetra, Ifpack2, Amesos2");
#endif

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(99);
  Teuchos::ParameterList galeriParameters;
  galeriParameters.set<GO>("nx", 99);
  RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", A->getRowMap(), galeriParameters);

  ArrayRCP<std::string> fileList = TestHelpers::GetFileList(std::string("ParameterList/MLParameterListInterpreter/"), std::string(".xml"));

  for (int i = 0; i < fileList.size(); i++) {
    out << "Processing file: " << fileList[i] << std::endl;
    Teuchos::ParameterList myList;
    myList.set("xml parameter file", "ParameterList/MLParameterListInterpreter/" + fileList[i]);

    Teuchos::ArrayRCP<typename MultiVector::scalar_type> xcoord = coordinates->getDataNonConst(0);
    myList.set("x-coordinates", xcoord.get());

    MLParameterListInterpreter mueluFactory(myList, A->getRowMap()->getComm());

    RCP<Hierarchy> H = mueluFactory.CreateHierarchy();
    H->GetLevel(0)->template Set<RCP<Matrix> >("A", A);

    mueluFactory.SetupHierarchy(*H);

    // TODO: check no unused parameters
    // TODO: check results of Iterate()
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MLParameterListInterpreter, SetParameterList, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

// TODO: some tests of the parameter list parser can be done without building the Hierarchy.
