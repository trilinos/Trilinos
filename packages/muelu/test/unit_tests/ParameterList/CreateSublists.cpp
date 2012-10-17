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

#include <MueLu_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#ifdef HAVE_MUELU_ML
#include <ml_epetra_utils.h>
#endif

#include "MueLu_TestHelpers.hpp"

#include "MueLu_MLParameterListInterpreter.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(MueLu_CreateSubLists, SetParameterList)
  {
    std::string dir("ParameterList/CreateSublists/");

    ArrayRCP<std::string> fileList = TestHelpers::GetFileList(dir, std::string(".xml"));

    for(int i=0; i< fileList.size(); i++) {
      out << "Processing file: " << fileList[i] << std::endl;

      Teuchos::RCP<const Teuchos::ParameterList> inputList = Teuchos::getParametersFromXmlFile(dir + fileList[i]);
      Teuchos::ParameterList outputList;
      
      MueLu::CreateSublists(*inputList, outputList);
      
      // Test against reference output (replace '.xml' by '.output' to get the filename)
      Teuchos::RCP<Teuchos::ParameterList> refOutputList = Teuchos::getParametersFromXmlFile(dir + fileList[i].substr(0, fileList[i].find_last_of(".")) + ".output");
      TEST_EQUALITY(outputList, *refOutputList);
    }
  }
 
#ifdef HAVE_MUELU_ML
  TEUCHOS_UNIT_TEST(ML_CreateSublists, SetParameterList)
  {
    std::string dir("ParameterList/CreateSublists/");

    ArrayRCP<std::string> fileList = TestHelpers::GetFileList(dir, std::string(".xml"));

    for(int i=0; i< fileList.size(); i++) {
      out << "Processing file: " << fileList[i] << std::endl;

      Teuchos::RCP<const Teuchos::ParameterList> inputList = Teuchos::getParametersFromXmlFile(dir + fileList[i]);
      Teuchos::ParameterList outputList;
      
      ML_CreateSublists(*inputList, outputList);

      // Test against reference output (replace '.xml' by '.output' to get the filename)
      Teuchos::RCP<Teuchos::ParameterList> refOutputList = Teuchos::getParametersFromXmlFile(dir + fileList[i].substr(0, fileList[i].find_last_of(".")) + ".output");
      TEST_EQUALITY(outputList, *refOutputList);
    }
  }
#endif

} // namespace MueLuTests
