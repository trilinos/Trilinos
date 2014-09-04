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
#ifndef MUELU_AMALGAMATIONFACTORY_DEF_HPP
#define MUELU_AMALGAMATIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_AmalgamationFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A"); // sub-block from blocked A
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const
  {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    LO fullblocksize    = 1;    // block dim for fixed size blocks
    GO offset           = 0;   // global offset of dof gids
    LO blockid          = -1;  // block id in strided map
    LO nStridedOffset   = 0;   // DOF offset for strided block id "blockid" (default = 0)
    LO stridedblocksize = fullblocksize; // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)
    GO indexBase        = A->getRowMap()->getIndexBase();  // index base for maps

    // 1) check for blocking/striding information

    if (A->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
      Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // NOTE: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
      RCP<const StridedMap> stridedRowMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
      TEUCHOS_TEST_FOR_EXCEPTION(stridedRowMap == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
      fullblocksize = stridedRowMap->getFixedBlockSize();
      offset        = stridedRowMap->getOffset();
      blockid       = stridedRowMap->getStridedBlockId();

      if (blockid > -1) {
        std::vector<size_t> stridingInfo = stridedRowMap->getStridingData();
        for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
          nStridedOffset += stridingInfo[j];
        stridedblocksize = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);

      } else {
        stridedblocksize = fullblocksize;
      }
      oldView = A->SwitchToView(oldView);
      GetOStream(Runtime1) << "AmalagamationFactory::Build():" << " found fullblocksize=" << fullblocksize << " and stridedblocksize=" << stridedblocksize << " from strided maps. offset=" << offset << std::endl;

    } else {
      GetOStream(Warnings0) << "AmalagamationFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;
    }
    // TODO: maybe no striding information on coarser levels -> misuse nullspace vector?

    // 2) prepare maps for amalgamated graph of A and
    //    setup unamalgamation information

    RCP<const Map> const &rowMap = A->getRowMap();
    Teuchos::ArrayView<const GlobalOrdinal> globalColDofs = A->getColMap()->getNodeElementList();
    LocalOrdinal nColEle = Teuchos::as<LocalOrdinal>(globalColDofs.size());
    RCP<std::vector<GlobalOrdinal> > gNodeIds; // contains global node ids on current proc
    gNodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);
    for (LocalOrdinal i = 0; i < nColEle; i++) {
      GlobalOrdinal gDofId = globalColDofs[i];
      GlobalOrdinal gNodeId = DOFGid2NodeId(gDofId, fullblocksize, offset, indexBase);
      if (rowMap->isNodeGlobalElement(gDofId))
        gNodeIds->push_back(gNodeId);
    }
    //make the gNodeIds unique
    std::sort( gNodeIds->begin(), gNodeIds->end() );
    gNodeIds->erase( std::unique( gNodeIds->begin(), gNodeIds->end() ), gNodeIds->end() );

    // store (un)amalgamation information on current level
    RCP<AmalgamationInfo> amalgamationData = rcp(new AmalgamationInfo(gNodeIds, A->getColMap(), fullblocksize, offset, blockid, nStridedOffset, stridedblocksize) );

    Set(currentLevel, "UnAmalgamationInfo", amalgamationData);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  const GlobalOrdinal AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DOFGid2NodeId(GlobalOrdinal gid, LocalOrdinal blockSize, const GlobalOrdinal offset, const GlobalOrdinal indexBase) {
    // here, the assumption is, that the node map has the same indexBase as the dof map
    GlobalOrdinal globalblockid = ((GlobalOrdinal) gid - offset - indexBase) / blockSize + indexBase;
    return globalblockid;
  }

} //namespace MueLu

#endif /* MUELU_SUBBLOCKUNAMALGAMATIONFACTORY_DEF_HPP */

