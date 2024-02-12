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
#ifndef MUELU_USERPFACTORY_DEF_HPP
#define MUELU_USERPFACTORY_DEF_HPP

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_UserPFactory_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> UserPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory of the nullspace");
  validParamList->set<std::string>("matrixFileName", "", "File with matrix data");
  validParamList->set<std::string>("mapFileName", "", "File with coarse map data");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UserPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
  Input(fineLevel, "A");
  Input(fineLevel, "Nullspace");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UserPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UserPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  RCP<Matrix> A                  = Get<RCP<Matrix> >(fineLevel, "A");
  RCP<MultiVector> fineNullspace = Get<RCP<MultiVector> >(fineLevel, "Nullspace");

  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() != 1, Exceptions::RuntimeError, "Block size > 1 has not been implemented");

  const Teuchos::ParameterList& pL = GetParameterList();

  std::string mapFile      = pL.get<std::string>("mapFileName");
  RCP<const Map> rowMap    = A->getRowMap();
  RCP<const Map> coarseMap = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMap(mapFile, rowMap->lib(), rowMap->getComm());
  Set(coarseLevel, "CoarseMap", coarseMap);

  std::string matrixFile = pL.get<std::string>("matrixFileName");
  RCP<Matrix> P          = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(matrixFile, rowMap, coarseMap, coarseMap, rowMap);
#if 1
  Set(coarseLevel, "P", P);
#else
  // Expand column map by 1
  RCP<Matrix> P1 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *P, false);
  P              = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(matrixFile, rowMap, P1->getColMap(), coarseMap, rowMap);
  Set(coarseLevel, "P", P);
#endif

  RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(coarseMap, fineNullspace->getNumVectors());
  P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(), Teuchos::ScalarTraits<SC>::zero());
  Set(coarseLevel, "Nullspace", coarseNullspace);

  // Coordinates transfer
  size_t n = Teuchos::as<size_t>(sqrt(coarseMap->getGlobalNumElements()));
  TEUCHOS_TEST_FOR_EXCEPTION(n * n != coarseMap->getGlobalNumElements(), Exceptions::RuntimeError, "Unfortunately, this is not the case, don't know what to do");

  RCP<MultiVector> coarseCoords = MultiVectorFactory::Build(coarseMap, 2);
  ArrayRCP<Scalar> x = coarseCoords->getDataNonConst(0), y = coarseCoords->getDataNonConst(1);
  for (size_t LID = 0; LID < coarseMap->getLocalNumElements(); ++LID) {
    GlobalOrdinal GID = coarseMap->getGlobalElement(LID) - coarseMap->getIndexBase();
    GlobalOrdinal i = GID % n, j = GID / n;
    x[LID] = i;
    y[LID] = j;
  }
  Set(coarseLevel, "Coordinates", coarseCoords);

  if (IsPrint(Statistics1)) {
    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*P, "P", params);
  }
}

}  // namespace MueLu

#define MUELU_USERPFACTORY_SHORT
#endif  // MUELU_USERPFACTORY_DEF_HPP
