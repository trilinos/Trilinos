// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LOWPRECISIONFACTORY_DEF_HPP
#define MUELU_LOWPRECISIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>

#include "MueLu_LowPrecisionFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> LowPrecisionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<std::string>("matrix key", "A", "");
  validParamList->set<RCP<const FactoryBase> >("R", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LowPrecisionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const ParameterList& pL = GetParameterList();
  std::string matrixKey   = pL.get<std::string>("matrix key");
  Input(currentLevel, matrixKey);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LowPrecisionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  using Teuchos::ParameterList;

  const ParameterList& pL = GetParameterList();
  std::string matrixKey   = pL.get<std::string>("matrix key");

  FactoryMonitor m(*this, "Converting " + matrixKey + " to half precision", currentLevel);

  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, matrixKey);

  GetOStream(Warnings) << "Matrix not converted to half precision. This only works for Tpetra and when both Scalar and HalfScalar have been instantiated." << std::endl;
  Set(currentLevel, matrixKey, A);
}

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> LowPrecisionFactory<double, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<std::string>("matrix key", "A", "");
  validParamList->set<RCP<const FactoryBase> >("R", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LowPrecisionFactory<double, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const ParameterList& pL = GetParameterList();
  std::string matrixKey   = pL.get<std::string>("matrix key");
  Input(currentLevel, matrixKey);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LowPrecisionFactory<double, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  using Teuchos::ParameterList;
  using HalfScalar = typename Teuchos::ScalarTraits<Scalar>::halfPrecision;

  const ParameterList& pL = GetParameterList();
  std::string matrixKey   = pL.get<std::string>("matrix key");

  FactoryMonitor m(*this, "Converting " + matrixKey + " to half precision", currentLevel);

  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, matrixKey);

  if ((A->getRowMap()->lib() == Xpetra::UseTpetra) && std::is_same<Scalar, double>::value) {
    auto tpA        = rcp_dynamic_cast<TpetraCrsMatrix>(rcp_dynamic_cast<CrsMatrixWrap>(A)->getCrsMatrix(), true)->getTpetra_CrsMatrix();
    auto tpLowA     = tpA->template convert<HalfScalar>();
    auto tpLowOpA   = rcp(new Tpetra::CrsMatrixMultiplyOp<Scalar, HalfScalar, LocalOrdinal, GlobalOrdinal, Node>(tpLowA));
    auto xpTpLowOpA = rcp(new TpetraOperator(tpLowOpA));
    auto xpLowOpA   = rcp_dynamic_cast<Operator>(xpTpLowOpA);
    Set(currentLevel, matrixKey, xpLowOpA);
    return;
  }

  GetOStream(Warnings) << "Matrix not converted to half precision. This only works for Tpetra and when both Scalar and HalfScalar have been instantiated." << std::endl;
  Set(currentLevel, matrixKey, A);
}
#endif

#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> LowPrecisionFactory<std::complex<double>, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<std::string>("matrix key", "A", "");
  validParamList->set<RCP<const FactoryBase> >("R", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory of the matrix A to be converted to lower precision");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LowPrecisionFactory<std::complex<double>, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const ParameterList& pL = GetParameterList();
  std::string matrixKey   = pL.get<std::string>("matrix key");
  Input(currentLevel, matrixKey);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LowPrecisionFactory<std::complex<double>, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  using Teuchos::ParameterList;
  using HalfScalar = typename Teuchos::ScalarTraits<Scalar>::halfPrecision;

  const ParameterList& pL = GetParameterList();
  std::string matrixKey   = pL.get<std::string>("matrix key");

  FactoryMonitor m(*this, "Converting " + matrixKey + " to half precision", currentLevel);

  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, matrixKey);

  if ((A->getRowMap()->lib() == Xpetra::UseTpetra) && std::is_same<Scalar, std::complex<double> >::value) {
    auto tpA        = rcp_dynamic_cast<TpetraCrsMatrix>(rcp_dynamic_cast<CrsMatrixWrap>(A)->getCrsMatrix(), true)->getTpetra_CrsMatrix();
    auto tpLowA     = tpA->template convert<HalfScalar>();
    auto tpLowOpA   = rcp(new Tpetra::CrsMatrixMultiplyOp<Scalar, HalfScalar, LocalOrdinal, GlobalOrdinal, Node>(tpLowA));
    auto xpTpLowOpA = rcp(new TpetraOperator(tpLowOpA));
    auto xpLowOpA   = rcp_dynamic_cast<Operator>(xpTpLowOpA);
    Set(currentLevel, matrixKey, xpLowOpA);
    return;
  }

  GetOStream(Warnings) << "Matrix not converted to half precision. This only works for Tpetra and when both Scalar and HalfScalar have been instantiated." << std::endl;
  Set(currentLevel, matrixKey, A);
}
#endif

}  // namespace MueLu

#endif  // MUELU_LOWPRECISIONFACTORY_DEF_HPP
