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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
#define MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP

#include "MueLu_MultiVectorTransferFactory_decl.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< std::string >           ("Vector name",   "undefined", "Name of the vector that will be transfered on the coarse grid (level key)"); // TODO: how to set a validator without default value?
    validParamList->set< RCP<const FactoryBase> >("R",           Teuchos::null, "Factory of the transfer operator (restriction)");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MultiVectorTransferFactory(std::string const & vectorName) {
    SetParameter("Vector name", ParameterEntry(vectorName));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    if (fineLevel.GetLevelID() == 0) fineLevel.DeclareInput(vectorName, MueLu::NoFactory::get(), this);
    else                             fineLevel.DeclareInput(vectorName, this                   , this);

    Input(coarseLevel, "R");

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    RCP<MultiVector> vector = fineLevel.Get< RCP<MultiVector> >(vectorName, MueLu::NoFactory::get()); //FIXME

    RCP<Matrix> transferOp = Get<RCP<Matrix> >(coarseLevel, "R");

    RCP<MultiVector> result = MultiVectorFactory::Build(transferOp->getRangeMap(), vector->getNumVectors());
    GetOStream(Runtime0, 0) << "Transferring multivector \"" << vectorName << "\"" << std::endl;

    transferOp->apply(*vector, *result);

    if (vectorName == "Coordinates") {
      GetOStream(Runtime0, 0) << "Averaging coordinates" << std::endl;
      size_t numLocalRows = transferOp->getNodeNumRows();
      Array<size_t> numEntriesPerRow(numLocalRows);
      for (size_t i=0; i<numLocalRows; ++i) {
        numEntriesPerRow.push_back(transferOp->getNumEntriesInLocalRow(i));
        if (numEntriesPerRow[i] == 0) numEntriesPerRow[i] = 1;
      }

      assert(numLocalRows == result->getLocalLength());
      for (size_t i=0; i<result->getNumVectors(); ++i) {
        ArrayRCP<Scalar> vals = result->getDataNonConst(i);
        for (size_t j=0; j<numLocalRows; ++j) {
          vals[j] /= numEntriesPerRow[j];
        }
      }
    }

    coarseLevel.Set<RCP<MultiVector> >(vectorName, result, this);
    coarseLevel.Set<RCP<MultiVector> >(vectorName, result); //FIXME

  } //Build

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayRCP<Scalar> MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::expandCoordinates(ArrayRCP<SC> coord, LocalOrdinal blksize) {
    if (blksize == 1)
      return coord;

    ArrayRCP<SC> expandCoord(coord.size()*blksize); //TODO: how to avoid automatic initialization of the vector? using arcp()?

    for(int i=0; i<coord.size(); i++) {
      for(int j=0; j< blksize; j++) {
        expandCoord[i*blksize + j] = coord[i];
      }
    }
    return expandCoord;

  } // expandCoord

} // namespace MueLu

#endif // MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
