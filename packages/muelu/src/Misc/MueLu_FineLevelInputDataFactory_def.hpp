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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DEF_HPP_

#include "Xpetra_Matrix.hpp"

#include "MueLu_FineLevelInputDataFactory_decl.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> FineLevelInputDataFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< std::string >("Fine level variable", std::string("A"), "Variable name on finest level.");
    validParamList->set< std::string >("Variable", std::string("A"), "Variable name on all coarse levels (except the finest level).");

    validParamList->set< std::string >("Fine level factory", std::string("NoFactory"), "Generating factory of the fine level variable");
    validParamList->set< std::string >("Coarse level factory", std::string("NoFactory"), "Generating factory for data on all coarse levels (except the finest)");

    validParamList->set< RCP<const FactoryBase> >("FinestLevelDataFactory", Teuchos::null, "Generating factory for level data on finest level");
    validParamList->set< RCP<const FactoryBase> >("CoarseLevelDataFactory", Teuchos::null, "Generating factory for coarse level data");

    //return validParamList;
    return Teuchos::null;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FineLevelInputDataFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {

    const ParameterList & pL = GetParameterList();

    std::string variableName = "";
    if(pL.isParameter("Variable"))
      variableName = pL.get<std::string>("Variable");

    std::string factoryName  = "NoFactory";

    if (currentLevel.GetLevelID() == 0) {
      if(pL.isParameter("Fine level factory"))
        factoryName = pL.get<std::string>("Fine level factory");
    } else {
      if(pL.isParameter("Coarse level factory"))
        factoryName = pL.get<std::string>("Coarse level factory");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(variableName == "", MueLu::Exceptions::RuntimeError, "FineLevelInputDataFactory: no variable name provided. Please set \'Variable\' parameter in your input deck.");

    RCP<const FactoryBase> fact = GetFactory(factoryName);
    if (fact == Teuchos::null) { fact = currentLevel.GetFactoryManager()->GetFactory(factoryName); }
    currentLevel.DeclareInput(variableName, fact.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FineLevelInputDataFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "InputUserData", currentLevel);

    const ParameterList& pL = GetParameterList();

    std::string variableName       = "";
    if (pL.isParameter("Variable"))
      variableName = pL.get<std::string>("Variable");
    std::string factoryName        = "NoFactory";

    if (currentLevel.GetLevelID() == 0) {
      if(pL.isParameter("Fine level factory"))
        factoryName        = pL.get<std::string>("Fine level factory");
    } else {
      if(pL.isParameter("Coarse level factory"))
        factoryName        = pL.get<std::string>("Coarse level factory");
    }
    RCP<const FactoryBase> fact = GetFactory(factoryName);
    if (fact == Teuchos::null) { fact = currentLevel.GetFactoryManager()->GetFactory(factoryName); }

    GetOStream(Debug) << "Use " << variableName << " from " << factoryName << "(" << fact.get() << ")" << std::endl;

    // TODO make it work for other variable types
    RCP<Matrix> data = currentLevel.Get< RCP<Matrix> >(variableName, fact.get());
    Set(currentLevel, variableName, data);
      //Teuchos::TypeNameTraits<T>::name();

    //GetOStream(Runtime0) << "Lumping dropped entries" << std::endl;
  }

} //namespace MueLu

#endif /* PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DEF_HPP_ */
