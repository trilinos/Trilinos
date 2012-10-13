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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_TestHelpers.hpp"

#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; // TODO: remove

  TEUCHOS_UNIT_TEST(ParameterListInterpreter, SetParameterList)
  {

    //TODO: this test can be done at compilation time
#if !defined(HAVE_MUELU_IFPACK)  or !defined(HAVE_MUELU_AMESOS)
    if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) {
      out << "Test skipped (dependencies not available)" << std::endl;
      return;
    }
#endif

#if !defined(HAVE_MUELU_IFPACK2) or !defined(HAVE_MUELU_AMESOS2)
    if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
      out << "Test skipped (dependencies not available)" << std::endl;
      return;
    }
#endif

    RCP<Matrix> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99);

    ArrayRCP<std::string> fileList = TestHelpers::GetFileList(std::string("ParameterList/ParameterListInterpreter/"), std::string(".xml"));

    for(int i=0; i< fileList.size(); i++) {
      out << "Processing file: " << fileList[i] << std::endl;
      ParameterListInterpreter mueluFactory("ParameterList/ParameterListInterpreter/" + fileList[i]);
      
      RCP<Hierarchy> H = mueluFactory.CreateHierarchy();
      H->GetLevel(0)->Set<RCP<Matrix> >("A", A);
     
      mueluFactory.SetupHierarchy(*H);
      
      //TODO: check no unused parameters
      //TODO: check results of Iterate()
    }
  }
  
} // namespace MueLuTests

// TODO: some tests of the parameter list parser can be done without building the Hierarchy. 
