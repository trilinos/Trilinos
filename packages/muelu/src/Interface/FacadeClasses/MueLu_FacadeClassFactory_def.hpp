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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DEF_HPP_

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>


#include "MueLu_Exceptions.hpp"

#include "MueLu_Facade_Simple_decl.hpp"

#include "MueLu_FacadeClassFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FacadeClassFactory() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Teuchos::ParameterList> FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const ParameterList& paramList) {

    FacadeSimple<Scalar,LocalOrdinal,GlobalOrdinal,Node> test;

    Teuchos::RCP<Teuchos::ParameterList> facadeParams = test.SetParameterList(paramList);

    return facadeParams;
    // obtain ParameterList with default input parameters for this facade class
    /*std::string defaultString = FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::defaultParams_;
    Teuchos::RCP<ParameterList> defaultList = Teuchos::getParametersFromXmlString(defaultString);

    // validate user input parameters (and set defaults if necessary)
    Teuchos::ParameterList inputParameters = paramList;
    inputParameters.validateParametersAndSetDefaults(*defaultList);

    TEUCHOS_TEST_FOR_EXCEPTION(inputParameters.get<std::string>("MueLu preconditioner") == "undefined", MueLu::Exceptions::RuntimeError, "FacadeClassFactory: undefined MueLu preconditioner. Set the \"MueLu preconditioner\" parameter correctly in your input file.");

    // create copy of template string which is updated with in-place string replacements
    std::string finalString = FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stringTemplate_;

    // loop over all input parameters
    for(Teuchos::ParameterList::ConstIterator it = inputParameters.begin(); it != inputParameters.end(); it++) {
      // form replacement string
      std::string par_name = inputParameters.name(it);
      std::stringstream ss;
      ss << "XXX" << par_name << "YYY";

      // update final string with parameters
      Teuchos::ParameterEntry par_entry = inputParameters.entry(it);
      ReplaceString(finalString,
              ss.str(), Teuchos::toString(par_entry.getAny()));
    }

    // logical code for more complicated distinctions
    if(inputParameters.get<bool>("Block 1: transfer smoothing") == true) {
      ReplaceString(finalString, "XXXBlock 1: prolongatorYYY", "myPFact1");
      ReplaceString(finalString, "XXXBlock 1: restrictor YYY", "myRFact1");
    } else {
      ReplaceString(finalString, "XXXBlock 1: prolongatorYYY", "myTentativePFact1");
      ReplaceString(finalString, "XXXBlock 1: restrictor YYY", "myTransPFact1");
    }
    if(inputParameters.get<bool>("Block 2: transfer smoothing") == true) {
      ReplaceString(finalString, "XXXBlock 2: prolongatorYYY", "myPFact2");
      ReplaceString(finalString, "XXXBlock 2: restrictor YYY", "myRFact2");
    } else {
      ReplaceString(finalString, "XXXBlock 2: prolongatorYYY", "myTentativePFact2");
      ReplaceString(finalString, "XXXBlock 2: restrictor YYY", "myTransPFact2");
    }
    // end logical code

    Teuchos::RCP<ParameterList> ret = Teuchos::getParametersFromXmlString(finalString);*/
    //return ret;
  }

} // end namespace MueLu

#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DEF_HPP_ */
