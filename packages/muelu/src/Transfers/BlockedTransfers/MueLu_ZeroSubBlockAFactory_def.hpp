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

#ifndef MUELU_ZEROSUBBLOCKAFACTORY_DEF_HPP_
#define MUELU_ZEROSUBBLOCKAFACTORY_DEF_HPP_

#include "MueLu_ZeroSubBlockAFactory_decl.hpp"

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ZeroSubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", MueLu::NoFactory::getRCP(), "Generating factory for A.");

  validParamList->set<int>("block row", 0, "Block row of subblock in block matrix A");
  validParamList->set<int>("block col", 0, "Block column of subblock in block matrix A");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ZeroSubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ZeroSubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  const ParameterList& pL = GetParameterList();
  const size_t row        = Teuchos::as<size_t>(pL.get<int>("block row"));
  const size_t col        = Teuchos::as<size_t>(pL.get<int>("block col"));

  // Get the blocked input matrix
  RCP<Matrix> Ain         = currentLevel.Get<RCP<Matrix>>("A", this->GetFactory("A").get());
  RCP<BlockedCrsMatrix> A = rcp_dynamic_cast<BlockedCrsMatrix>(Ain);

  // Perform some basic sanity checks
  TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), Exceptions::BadCast, "Input matrix A is not a BlockedCrsMatrix.");
  TEUCHOS_TEST_FOR_EXCEPTION(row > A->Rows(), Exceptions::RuntimeError, "row [" << row << "] > A.Rows() [" << A->Rows() << "].");
  TEUCHOS_TEST_FOR_EXCEPTION(col > A->Cols(), Exceptions::RuntimeError, "col [" << col << "] > A.Cols() [" << A->Cols() << "].");

  // strided maps for range and domain map of sub matrix
  RCP<const StridedMap> stridedRangeMap  = Teuchos::null;
  RCP<const StridedMap> stridedDomainMap = Teuchos::null;

  // extract map information from MapExtractor
  RCP<const MapExtractor> rangeMapExtractor  = A->getRangeMapExtractor();
  RCP<const MapExtractor> domainMapExtractor = A->getDomainMapExtractor();

  // In case that both user-specified and internal striding information from the submaps
  // does not contain valid striding information, try to extract it from the global maps
  // in the map extractor.
  if (stridedRangeMap.is_null()) {
    TEUCHOS_TEST_FOR_EXCEPTION(rangeMapExtractor->getMap(row).is_null(), Exceptions::BadCast,
                               "Range map extractor contains non-strided maps in block row " << row << ". This should not be.");
    stridedRangeMap = rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getMap(row));
  }

  if (stridedDomainMap.is_null()) {
    TEUCHOS_TEST_FOR_EXCEPTION(domainMapExtractor->getMap(row).is_null(), Exceptions::BadCast,
                               "Domain map extractor contains non-strided maps in block row " << row << ". This should not be.");
    stridedDomainMap = rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getMap(col));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(stridedRangeMap.is_null(), Exceptions::BadCast, "rangeMap " << row << " is not a strided map.");
  TEUCHOS_TEST_FOR_EXCEPTION(stridedDomainMap.is_null(), Exceptions::BadCast, "domainMap " << col << " is not a strided map.");

  RCP<Matrix> Op = MatrixFactory::Build(stridedRangeMap, stridedDomainMap, static_cast<size_t>(0));
  TEUCHOS_ASSERT(!Op.is_null());

  Op->fillComplete(stridedDomainMap, stridedRangeMap);
  TEUCHOS_ASSERT(Op->isFillComplete());

  GetOStream(Statistics1) << "A(" << row << "," << col << ") is a single block and has strided maps:"
                          << "\n  range  map fixed block size = " << stridedRangeMap->getFixedBlockSize() << ", strided block id = " << stridedRangeMap->getStridedBlockId()
                          << "\n  domain map fixed block size = " << stridedDomainMap->getFixedBlockSize() << ", strided block id = " << stridedDomainMap->getStridedBlockId() << std::endl;
  GetOStream(Statistics2) << "A(" << row << "," << col << ") has " << Op->getGlobalNumRows() << "x" << Op->getGlobalNumCols() << " rows and columns." << std::endl;

  if (Op->IsView("stridedMaps") == true)
    Op->RemoveView("stridedMaps");
  Op->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

  currentLevel.Set("A", Op, this);
}

}  // namespace MueLu

#endif /* MUELU_ZEROSUBBLOCKAFACTORY_DEF_HPP_ */
