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
  MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MultiVectorTransferFactory(
    std::string const & vectorName,
    std::string const & restrictionName)
    : vectorName_(vectorName),
      restrictionName_(restrictionName)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    if (fineLevel.GetLevelID() == 0) fineLevel.DeclareInput(vectorName_, MueLu::NoFactory::get(), this);
    else                             fineLevel.DeclareInput(vectorName_, this                   , this);
    coarseLevel.DeclareInput(restrictionName_, GetFactory(restrictionName_).get(), this);
    /*
    //FIXME ThreeLevels unit test dies
    fineLevel.DeclareInput(vectorName_, restrictionFact_.get(), this);
    coarseLevel.DeclareInput(restrictionName_, restrictionFact_.get(), this);
    */
  } // DeclareInput

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level &coarseLevel) const {

    FactoryMonitor m(*this, "Build", coarseLevel);

    RCP<MultiVector> vector = fineLevel.Get<RCP<MultiVector> >(vectorName_, MueLu::NoFactory::get()); //FIXME

    RCP<Matrix> transferOp = coarseLevel.Get<RCP<Matrix> >(restrictionName_, GetFactory(restrictionName_).get());

    RCP<MultiVector> result = MultiVectorFactory::Build(transferOp->getRangeMap(), vector->getNumVectors());
    GetOStream(Runtime0, 0) << "Transferring multivector \"" << vectorName_ << "\"" << std::endl;

    transferOp->apply(*vector, *result);

    if (vectorName_ == "Coordinates") {
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

    coarseLevel.Set<RCP<MultiVector> >(vectorName_, result, this);
    coarseLevel.Set<RCP<MultiVector> >(vectorName_, result); //FIXME

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
