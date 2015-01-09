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
/*
 * MueLu_SubBlockAFactory_def.hpp
 *
 *  Created on: 02.01.2012
 *      Author: tobias
 */

#ifndef MUELU_SUBBLOCKAFACTORY_DEF_HPP_
#define MUELU_SUBBLOCKAFACTORY_DEF_HPP_


#include "MueLu_SubBlockAFactory_decl.hpp"

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",          MueLu::NoFactory::getRCP(), "Generating factory for A.");
    validParamList->set< int >                   ("block row",                           0, "Block row of subblock matrix A");
    validParamList->set< int >                   ("block col",                           0, "Block column of subblock matrix A");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    const ParameterList& pL = GetParameterList();
    size_t row = Teuchos::as<size_t>(pL.get<int>("block row"));
    size_t col = Teuchos::as<size_t>(pL.get<int>("block col"));

    RCP<Matrix>           Ain = Get<RCP<Matrix> >(currentLevel, "A");
    RCP<BlockedCrsMatrix> A   = rcp_dynamic_cast<BlockedCrsMatrix>(Ain);

    TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(),     Exceptions::BadCast,      "Input matrix A is not a BlockedCrsMatrix.");
    TEUCHOS_TEST_FOR_EXCEPTION(row > A->Rows(), Exceptions::RuntimeError, "row [" << row << "] > A.Rows() [" << A->Rows() << "].");
    TEUCHOS_TEST_FOR_EXCEPTION(col > A->Cols(), Exceptions::RuntimeError, "col [" << col << "] > A.Cols() [" << A->Cols() << "].");

    RCP<CrsMatrixWrap> Op = Teuchos::rcp(new CrsMatrixWrap(A->getMatrix(row, col)));

    //////////////// EXPERIMENTAL
    // extract striding information from RangeMapExtractor

    RCP<const MapExtractor> rangeMapExtractor  = A->getRangeMapExtractor();
    RCP<const MapExtractor> domainMapExtractor = A->getDomainMapExtractor();

    RCP<const Map> rangeMap  = rangeMapExtractor ->getMap(row);
    RCP<const Map> domainMap = domainMapExtractor->getMap(col);

    RCP<const StridedMap> srangeMap  = rcp_dynamic_cast<const StridedMap>(rangeMap);
    RCP<const StridedMap> sdomainMap = rcp_dynamic_cast<const StridedMap>(domainMap);

    if (srangeMap.is_null()) {
      RCP<const Map>         fullRangeMap = rangeMapExtractor->getFullMap();
      RCP<const StridedMap> sFullRangeMap = rcp_dynamic_cast<const StridedMap>(fullRangeMap);
      TEUCHOS_TEST_FOR_EXCEPTION(sFullRangeMap.is_null(), Exceptions::BadCast, "Full rangeMap is not a strided map.");

      std::vector<size_t> stridedData = sFullRangeMap->getStridingData();
      if (stridedData.size() == 1 && row > 0) {
        // We have block matrices. use striding block information 0
        srangeMap = StridedMapFactory::Build(rangeMap, stridedData,   0, sFullRangeMap->getOffset());

      } else {
        // We have strided matrices. use striding information of the corresponding block
        srangeMap = StridedMapFactory::Build(rangeMap, stridedData, row, sFullRangeMap->getOffset());
      }
    }

    if (sdomainMap.is_null()) {
      RCP<const Map>         fullDomainMap = domainMapExtractor->getFullMap();
      RCP<const StridedMap> sFullDomainMap = rcp_dynamic_cast<const StridedMap>(fullDomainMap);
      TEUCHOS_TEST_FOR_EXCEPTION(sFullDomainMap.is_null(), Exceptions::BadCast, "Full domainMap is not a strided map");

      std::vector<size_t> stridedData = sFullDomainMap->getStridingData();
      if (stridedData.size() == 1 && col > 0) {
        // We have block matrices. use striding block information 0
        sdomainMap = StridedMapFactory::Build(domainMap, stridedData,   0, sFullDomainMap->getOffset());

      } else {
        // We have strided matrices. use striding information of the corresponding block
        sdomainMap = StridedMapFactory::Build(domainMap, stridedData, col, sFullDomainMap->getOffset());
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(srangeMap.is_null(),  Exceptions::BadCast, "rangeMap "  << row << " is not a strided map.");
    TEUCHOS_TEST_FOR_EXCEPTION(sdomainMap.is_null(), Exceptions::BadCast, "domainMap " << col << " is not a strided map.");

    GetOStream(Statistics1) << "A(" << row << "," << col << ") has strided maps:"
        << "\n  range  map fixed block size = " << srangeMap ->getFixedBlockSize() << ", strided block id = " << srangeMap ->getStridedBlockId()
        << "\n  domain map fixed block size = " << sdomainMap->getFixedBlockSize() << ", strided block id = " << sdomainMap->getStridedBlockId() << std::endl;

    if (Op->IsView("stridedMaps") == true)
      Op->RemoveView("stridedMaps");
    Op->CreateView("stridedMaps", srangeMap, sdomainMap);

    TEUCHOS_TEST_FOR_EXCEPTION(Op->IsView("stridedMaps") == false, Exceptions::RuntimeError, "Failed to set \"stridedMaps\" view.");

    //////////////// EXPERIMENTAL

    currentLevel.Set("A", rcp_dynamic_cast<Matrix>(Op), this);
  }

} // namespace MueLu

#endif /* MUELU_SUBBLOCKAFACTORY_DEF_HPP_ */
