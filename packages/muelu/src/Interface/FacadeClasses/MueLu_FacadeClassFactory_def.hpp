// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DEF_HPP_

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "MueLu_Exceptions.hpp"

#include "MueLu_FacadeClassBase.hpp"
#include "MueLu_FacadeSimple_def.hpp"
#include "MueLu_FacadeBGS2x2_def.hpp"

#include "MueLu_FacadeClassFactory_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FacadeClassFactory() {
  facadeClasses_["Simple"] = Teuchos::rcp(new FacadeSimple<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
  facadeClasses_["BGS2x2"] = Teuchos::rcp(new FacadeBGS2x2<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Teuchos::ParameterList> FacadeClassFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const ParameterList& paramList) {
  TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("MueLu preconditioner") == false, MueLu::Exceptions::RuntimeError, "FacadeClassFactory: undefined MueLu preconditioner. Set the \"MueLu preconditioner\" parameter correctly in your input file.");
  TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("MueLu preconditioner") == "undefined", MueLu::Exceptions::RuntimeError, "FacadeClassFactory: undefined MueLu preconditioner. Set the \"MueLu preconditioner\" parameter correctly in your input file.");

  std::string precMueLu = paramList.get<std::string>("MueLu preconditioner");

  // could not find requested facade class
  if (facadeClasses_.find(precMueLu) == facadeClasses_.end()) {
    GetOStream(Errors) << "FacadeClassFactory: Could not find facade class \"" << precMueLu << "\"!" << std::endl;
    GetOStream(Errors) << "The available facade classes are:" << std::endl;
    for (typename std::map<std::string, Teuchos::RCP<FacadeClassBase> >::const_iterator it = facadeClasses_.begin(); it != facadeClasses_.end(); it++) {
      GetOStream(Errors) << "   " << it->first << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "FacadeClassFactory: Could not find facade class \"" << precMueLu << "\".");
  }

  return facadeClasses_[precMueLu]->SetParameterList(paramList);
}

}  // end namespace MueLu

#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DEF_HPP_ */
