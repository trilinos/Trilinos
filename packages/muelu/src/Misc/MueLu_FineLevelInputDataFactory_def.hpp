// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DEF_HPP_

#include "Xpetra_Matrix.hpp"

#include "MueLu_FineLevelInputDataFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> FineLevelInputDataFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  // Variable name (e.g. A or P or Coordinates)
  validParamList->set<std::string>("Variable", std::string("A"), "Variable name on all coarse levels (except the finest level).");

  // Names of generating factories (on finest level and coarse levels)
  validParamList->set<RCP<const FactoryBase> >("Fine level factory", Teuchos::null, "Generating factory of the fine level variable");
  validParamList->set<RCP<const FactoryBase> >("Coarse level factory", Teuchos::null, "Generating factory for data on all coarse levels (except the finest)");

  // Type of variable (see source code for a complete list of all available types)
  validParamList->set<std::string>("Variable type", std::string("Matrix"), "Type of variable");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FineLevelInputDataFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const ParameterList& pL = GetParameterList();

  std::string variableName = "";
  if (pL.isParameter("Variable"))
    variableName = pL.get<std::string>("Variable");

  std::string factoryName = "NoFactory";
  if (currentLevel.GetLevelID() == 0) {
    factoryName = "Fine level factory";
  } else {
    factoryName = "Coarse level factory";
  }

  TEUCHOS_TEST_FOR_EXCEPTION(variableName == "", MueLu::Exceptions::RuntimeError, "FineLevelInputDataFactory: no variable name provided. Please set \'Variable\' parameter in your input deck.");

  // data must be specified in factory! (not in factory manager)
  RCP<const FactoryBase> fact = GetFactory(factoryName);
  currentLevel.DeclareInput(variableName, fact.get(), this);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FineLevelInputDataFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "InputUserData", currentLevel);

  const ParameterList& pL = GetParameterList();

  std::string variableName = "";
  if (pL.isParameter("Variable"))
    variableName = pL.get<std::string>("Variable");

  std::string variableType = "";
  if (pL.isParameter("Variable type"))
    variableType = pL.get<std::string>("Variable type");

  std::string factoryName = "NoFactory";
  if (currentLevel.GetLevelID() == 0) {
    factoryName = "Fine level factory";
  } else {
    factoryName = "Coarse level factory";
  }
  RCP<const FactoryBase> fact = GetFactory(factoryName);

  GetOStream(Debug) << "Use " << variableName << " of type " << variableType << " from " << factoryName << "(" << fact.get() << ")" << std::endl;

  // check data type
  // std::string strType = currentLevel.GetTypeName(variableName, fact.get());
  if (variableType == "int") {
    int data = currentLevel.Get<int>(variableName, fact.get());
    Set(currentLevel, variableName, data);
  } else if (variableType == "double") {
    double data = currentLevel.Get<double>(variableName, fact.get());
    Set(currentLevel, variableName, data);
  } else if (variableType == "string") {
    std::string data = currentLevel.Get<std::string>(variableName, fact.get());
    Set(currentLevel, variableName, data);
  } else {
    size_t npos = std::string::npos;

    if (variableType.find("Aggregates") != npos) {
      RCP<Aggregates> data = currentLevel.Get<RCP<Aggregates> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("SmootherBase") != npos) {
      RCP<SmootherBase> data = currentLevel.Get<RCP<SmootherBase> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("SmootherPrototype") != npos) {
      RCP<SmootherPrototype> data = currentLevel.Get<RCP<SmootherPrototype> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("Export") != npos) {
      RCP<Export> data = currentLevel.Get<RCP<Export> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("Import") != npos) {
      RCP<Import> data = currentLevel.Get<RCP<Import> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("Map") != npos) {
      RCP<Map> data = currentLevel.Get<RCP<Map> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("Matrix") != npos) {
      RCP<Matrix> data = currentLevel.Get<RCP<Matrix> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("MultiVector") != npos) {
      RCP<MultiVector> data = currentLevel.Get<RCP<MultiVector> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else if (variableType.find("Operator") != npos) {
      RCP<Operator> data = currentLevel.Get<RCP<Operator> >(variableName, fact.get());
      Set(currentLevel, variableName, data);
    } else {
      // TAW: is this working with empty procs?
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "FineLevelInputDataFactory: cannot detect type of variable " << variableName << " generated by " << fact.get() << ". User provided type " << variableType);
    }
  }
}

}  // namespace MueLu

#endif /* PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DEF_HPP_ */
