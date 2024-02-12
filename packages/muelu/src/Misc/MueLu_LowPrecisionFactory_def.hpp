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
