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

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SubBlockAFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SubBlockAFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A", MueLu::NoFactory::getRCP(), "Generating factory for A.");
    validParamList->set< int >("block row", 0, "block row of subblock matrix A");
    validParamList->set< int >("block col", 0, "block column of subblock matrix A");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SubBlockAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OMatrix; //TODO
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixClass; //TODO
    typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixWrapClass; //TODO
    typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> BlockedCrsOMatrix; //TODO
    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;

    const ParameterList & pL = GetParameterList();
    size_t row = Teuchos::as<size_t>(pL.get<int>("block row"));
    size_t col = Teuchos::as<size_t>(pL.get<int>("block col"));

    RCP<OMatrix> Ain = Teuchos::null;
    Ain = Get< RCP<OMatrix> >(currentLevel, "A");

    RCP<BlockedCrsOMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsOMatrix>(Ain);

    TEUCHOS_TEST_FOR_EXCEPTION(bA==Teuchos::null, Exceptions::BadCast, "MueLu::SubBlockAFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(row>bA->Rows(), Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: A.Rows() > rows_! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(col>bA->Cols(), Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: A.Cols() > cols_! error.");
    TEUCHOS_TEST_FOR_EXCEPTION(row<0, Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: row_<0 error.");
    TEUCHOS_TEST_FOR_EXCEPTION(col<0, Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: col_<0 error.");

    Teuchos::RCP<CrsMatrixClass> A = bA->getMatrix(row, col);

    Teuchos::RCP<CrsMatrixWrapClass> Op = Teuchos::rcp(new CrsMatrixWrapClass(A));

    //////////////// EXPERIMENTAL
    // extract striding information from RangeMapExtractor

    Teuchos::RCP<const MapExtractorClass> rgMapExtractor = bA->getRangeMapExtractor();
    Teuchos::RCP<const MapExtractorClass> doMapExtractor = bA->getDomainMapExtractor();

    Teuchos::RCP<const Map> rgMap = rgMapExtractor->getMap(row);
    Teuchos::RCP<const Map> doMap = doMapExtractor->getMap(col);

    Teuchos::RCP<const StridedMap> srgMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rgMap);
    Teuchos::RCP<const StridedMap> sdoMap = Teuchos::rcp_dynamic_cast<const StridedMap>(doMap);

    if(srgMap == Teuchos::null) {
      Teuchos::RCP<const Map> fullRgMap = rgMapExtractor->getFullMap();
      Teuchos::RCP<const StridedMap> sFullRgMap = Teuchos::rcp_dynamic_cast<const StridedMap>(fullRgMap);
      TEUCHOS_TEST_FOR_EXCEPTION(sFullRgMap==Teuchos::null, Exceptions::BadCast, "MueLu::SubBlockAFactory::Build: full rangeMap is not a strided map");
      std::vector<size_t> stridedData = sFullRgMap->getStridingData();
      if(stridedData.size() == 1 && row > 0) // we have block matrices. use striding block information 0
        srgMap = StridedMapFactory::Build(rgMap, stridedData, 0, sFullRgMap->getOffset());
      else // we have strided matrices. use striding information of the corresponding block
        srgMap = StridedMapFactory::Build(rgMap, stridedData, row, sFullRgMap->getOffset());
    }

    if(sdoMap == Teuchos::null) {
      Teuchos::RCP<const Map> fullDoMap = doMapExtractor->getFullMap();
      Teuchos::RCP<const StridedMap> sFullDoMap = Teuchos::rcp_dynamic_cast<const StridedMap>(fullDoMap);
      TEUCHOS_TEST_FOR_EXCEPTION(sFullDoMap==Teuchos::null, Exceptions::BadCast, "MueLu::SubBlockAFactory::Build: full domainMap is not a strided map");
      std::vector<size_t> stridedData2 = sFullDoMap->getStridingData();
      if(stridedData2.size() == 1 && col > 0) // we have block matrices. use striding block information 0
        sdoMap = StridedMapFactory::Build(doMap, stridedData2, 0, sFullDoMap->getOffset());
      else // we have strided matrices. use striding information of the corresponding block
        sdoMap = StridedMapFactory::Build(doMap, stridedData2, col, sFullDoMap->getOffset());
    }

    TEUCHOS_TEST_FOR_EXCEPTION(srgMap==Teuchos::null, Exceptions::BadCast, "MueLu::SubBlockAFactory::Build: rangeMap " << row << " is not a strided map");
    TEUCHOS_TEST_FOR_EXCEPTION(sdoMap==Teuchos::null, Exceptions::BadCast, "MueLu::SubBlockAFactory::Build: domainMap " << col << " is not a strided map");

    GetOStream(Statistics1) << "A(" << row << "," << col << ") has strided maps: range map fixed block size=" << srgMap->getFixedBlockSize() << " strided block id = " << srgMap->getStridedBlockId() << ", domain map fixed block size=" << sdoMap->getFixedBlockSize() << ", strided block id=" << sdoMap->getStridedBlockId() << std::endl;

    if(Op->IsView("stridedMaps") == true) Op->RemoveView("stridedMaps");
    Op->CreateView("stridedMaps", srgMap, sdoMap);
    TEUCHOS_TEST_FOR_EXCEPTION(Op->IsView("stridedMaps")==false, Exceptions::RuntimeError, "MueLu::SubBlockAFactory::Build: failed to set stridedMaps");

    //////////////// EXPERIMENTAL

    currentLevel.Set("A", Teuchos::rcp_dynamic_cast<OMatrix>(Op), this);
  }

} // namespace MueLu

#endif /* MUELU_SUBBLOCKAFACTORY_DEF_HPP_ */
