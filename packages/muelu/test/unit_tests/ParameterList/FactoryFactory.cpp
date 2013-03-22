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

#include "MueLu_FactoryFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; // TODO: remove

  // This is not a real unit test, because output of BuildFactory is not verified. But anyway, it still useful.
  TEUCHOS_UNIT_TEST(FactoryFactory, BuildFactory)
  {
    ArrayRCP<std::string> fileList =  TestHelpers::GetFileList(std::string("ParameterList/FactoryFactory/"), std::string(".xml"));

    for(int i=0; i< fileList.size(); i++) {
      out << "Processing file: " << fileList[i] << std::endl;
      Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("ParameterList/FactoryFactory/" + fileList[i]);

      for (Teuchos::ParameterList::ConstIterator param = paramList->begin(); param != paramList->end(); ++param) {
        Teuchos::OSTab tab(out);
        const std::string             & paramName  = paramList->name(param);
        const Teuchos::ParameterEntry & paramValue = paramList->entry(param);

        out << "Building object '" << paramName << "'" << std::endl;

        const FactoryMap factoryMapIn;
        FactoryFactory().BuildFactory(paramValue, factoryMapIn);
      }
    }
  }

} // namespace MueLuTests


