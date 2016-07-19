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

#include "MueLu_FacadeClassBase.hpp"
#include "MueLu_Facade_Simple_decl.hpp"
#include "MueLu_Facade_BGS2x2_decl.hpp"

#include "MueLu_FacadeClassFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FacadeClassFactory() {
    facadeClasses_["Simple"] = Teuchos::rcp(new FacadeSimple<Scalar,LocalOrdinal,GlobalOrdinal,Node>());
    facadeClasses_["BGS2x2"] = Teuchos::rcp(new FacadeBGS2x2<Scalar,LocalOrdinal,GlobalOrdinal,Node>());
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Teuchos::ParameterList> FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const ParameterList& paramList) {

    TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("MueLu preconditioner") == false, MueLu::Exceptions::RuntimeError, "FacadeClassFactory: undefined MueLu preconditioner. Set the \"MueLu preconditioner\" parameter correctly in your input file.");
    TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("MueLu preconditioner") == "undefined", MueLu::Exceptions::RuntimeError, "FacadeClassFactory: undefined MueLu preconditioner. Set the \"MueLu preconditioner\" parameter correctly in your input file.");

    std::string precMueLu = paramList.get<std::string>("MueLu preconditioner");

    // could not find requested facade class
    if(facadeClasses_.find(precMueLu) == facadeClasses_.end()) {
      GetOStream(Errors) << "FacadeClassFactory: Could not find facade class \"" << precMueLu << "\"!" << std::endl;
      GetOStream(Errors) << "The available facade classes are:" << std::endl;
      for(typename std::map<std::string, Teuchos::RCP<FacadeClassBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >::const_iterator it =facadeClasses_.begin(); it != facadeClasses_.end(); it++){
        GetOStream(Errors) << "   " << it->first << std::endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "FacadeClassFactory: Could not find facade class \"" << precMueLu << "\".");
    }

    return facadeClasses_[precMueLu]->SetParameterList(paramList);
  }

} // end namespace MueLu

#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DEF_HPP_ */
