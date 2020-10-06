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
#include <Xpetra_CrsMatrixWrap.hpp>
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

    validParamList->set<int>("block row", 0, "Block row of subblock matrix A");
    validParamList->set<int>("block col", 0, "Block column of subblock matrix A");

    validParamList->set<int>("row map from block col", -1, "Block row of subblock matrix A");
    validParamList->set<int>("column map from block row", -1, "Block row of subblock matrix A");

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
    const size_t row = Teuchos::as<size_t>(pL.get<int>("block row"));
    const size_t col = Teuchos::as<size_t>(pL.get<int>("block col"));

    const size_t rowForColMap = Teuchos::as<size_t>(pL.get<int>("row map from block col"));
    const size_t colForRowMap = Teuchos::as<size_t>(pL.get<int>("column map from block row"));
    TEUCHOS_TEST_FOR_EXCEPTION(rowForColMap==-1, Exceptions::InvalidArgument,
        "Block row to grab the column map is not specified. Specify a valid row via the parameter \"row map from block col\".");
    TEUCHOS_TEST_FOR_EXCEPTION(colForRowMap==-1, Exceptions::InvalidArgument,
        "Block column to grab the row map is not specified. Specify a valid column via the parameter \"column map from block row\".");

    // Get the blocked input matrix
    RCP<Matrix> Ain = currentLevel.Get<RCP<Matrix>>("A", this->GetFactory("A").get());
    RCP<BlockedCrsMatrix> A = rcp_dynamic_cast<BlockedCrsMatrix>(Ain);

    // Perform some basic sanity checks
    TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), Exceptions::BadCast, "Input matrix A is not a BlockedCrsMatrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(row > A->Rows(), Exceptions::RuntimeError, "row [" << row << "] > A.Rows() [" << A->Rows() << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(col > A->Cols(), Exceptions::RuntimeError, "col [" << col << "] > A.Cols() [" << A->Cols() << "].");

    // Create the required sub-matrix as an emtpy matrix w/ proper maps
    RCP<const Map> rowMap = A->getMatrix(row, colForRowMap)->getRowMap();
    RCP<const Map> colMap = A->getMatrix(rowForColMap, col)->getColMap();

    RCP<const Map> rangeMap = A->getMatrix(row, colForRowMap)->getRangeMap();
    RCP<const Map> domainMap = A->getMatrix(rowForColMap, col)->getDomainMap();
    RCP<Matrix> Op = MatrixFactory::Build(rowMap, colMap, static_cast<size_t>(0));

    TEUCHOS_ASSERT(!Op.is_null());

    Op->fillComplete(domainMap, rangeMap);
    TEUCHOS_ASSERT(Op->isFillComplete());

    // strided maps for range and domain map of sub matrix
    RCP<const StridedMap> stridedRangeMap = Teuchos::null;
    RCP<const StridedMap> stridedDomainMap = Teuchos::null;

    // extract map information from MapExtractor
    RCP<const MapExtractor> rangeMapExtractor = A->getRangeMapExtractor();
    RCP<const MapExtractor> domainMapExtractor = A->getDomainMapExtractor();

    stridedRangeMap = rcp_dynamic_cast<const StridedMap>(rangeMap);
    stridedDomainMap = rcp_dynamic_cast<const StridedMap>(domainMap);

    // In case that both user-specified and internal striding information from the submaps
    // does not contain valid striding information, try to extract it from the global maps
    // in the map extractor.
    if (stridedRangeMap.is_null()) {
      
      if (rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getMap(colForRowMap)) != Teuchos::null) {
        stridedRangeMap = rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getMap(colForRowMap));
      }
      else
      {     
        RCP<const Map> fullRangeMap = rangeMapExtractor->getFullMap();
        RCP<const StridedMap> stridedFullRangeMap = rcp_dynamic_cast<const StridedMap>(fullRangeMap);
        TEUCHOS_TEST_FOR_EXCEPTION(stridedFullRangeMap.is_null(), Exceptions::BadCast, "Full rangeMap is not a strided map.");

        std::vector<size_t> stridedData = stridedFullRangeMap->getStridingData();
        if (stridedData.size() == 1 && row > 0) {
          // We have block matrices. use striding block information 0
          stridedRangeMap = StridedMapFactory::Build(rangeMap, stridedData, 0, stridedFullRangeMap->getOffset());

        } else {
          // We have strided matrices. use striding information of the corresponding block
          stridedRangeMap = StridedMapFactory::Build(rangeMap, stridedData, row, stridedFullRangeMap->getOffset());
        }
      }
    }

    if (stridedDomainMap.is_null()) {
      if (rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getMap(rowForColMap)) != Teuchos::null) {
        stridedDomainMap = rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getMap(rowForColMap));
      }
      else
      {        
        RCP<const Map> fullDomainMap = domainMapExtractor->getFullMap();
        RCP<const StridedMap> stridedFullDomainMap = rcp_dynamic_cast<const StridedMap>(fullDomainMap);
        TEUCHOS_TEST_FOR_EXCEPTION(stridedFullDomainMap.is_null(), Exceptions::BadCast, "Full domainMap is not a strided map");

        std::vector<size_t> stridedData = stridedFullDomainMap->getStridingData();
        if (stridedData.size() == 1 && col > 0) {
          // We have block matrices. use striding block information 0
          stridedDomainMap = StridedMapFactory::Build(domainMap, stridedData, 0, stridedFullDomainMap->getOffset());

        } else {
          // We have strided matrices. use striding information of the corresponding block
          stridedDomainMap = StridedMapFactory::Build(domainMap, stridedData, col, stridedFullDomainMap->getOffset());
        }
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(stridedRangeMap.is_null(), Exceptions::BadCast, "rangeMap " << row << " is not a strided map.");
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDomainMap.is_null(), Exceptions::BadCast, "domainMap " << col << " is not a strided map.");

    GetOStream(Statistics1) << "A(" << row << "," << col << ") is a single block and has strided maps:"
        << "\n  range  map fixed block size = " << stridedRangeMap ->getFixedBlockSize() << ", strided block id = " << stridedRangeMap ->getStridedBlockId()
        << "\n  domain map fixed block size = " << stridedDomainMap->getFixedBlockSize() << ", strided block id = " << stridedDomainMap->getStridedBlockId() << std::endl;
    GetOStream(Statistics2) << "A(" << row << "," << col << ") has " << Op->getGlobalNumRows() << "x" << Op->getGlobalNumCols() << " rows and columns." << std::endl;

    // TODO do we really need that? we moved the code to getMatrix...
    if (Op->IsView("stridedMaps") == true)
      Op->RemoveView("stridedMaps");
    Op->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);


    currentLevel.Set("A", Op, this);
  }

} // namespace MueLu

#endif /* MUELU_ZEROSUBBLOCKAFACTORY_DEF_HPP_ */
