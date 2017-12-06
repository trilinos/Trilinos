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
#ifndef MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
#define MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP

#include "MueLu_MultiVectorTransferFactory_decl.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< std::string >           ("Vector name",      "undefined", "Name of the vector that will be transferred on the coarse grid (level key)"); // TODO: how to set a validator without default value?
    validParamList->set< RCP<const FactoryBase> >("Vector factory", Teuchos::null, "Factory of the vector");
    validParamList->set< RCP<const FactoryBase> >("R",              Teuchos::null, "Factory of the transfer operator (restriction)");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVectorTransferFactory(std::string const & vectorName) {
    SetParameter("Vector name", ParameterEntry(vectorName));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    fineLevel.DeclareInput(vectorName, GetFactory("Vector factory").get(), this);
    Input(coarseLevel, "R");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    const ParameterList & pL = GetParameterList();
    std::string vectorName   = pL.get<std::string>("Vector name");

    RCP<MultiVector> fineVector = fineLevel.Get< RCP<MultiVector> >(vectorName, GetFactory("Vector factory").get());
    RCP<Matrix>      transferOp = Get<RCP<Matrix> >(coarseLevel, "R");

    RCP<MultiVector> coarseVector = MultiVectorFactory::Build(transferOp->getRangeMap(), fineVector->getNumVectors());
    GetOStream(Runtime0) << "Transferring multivector \"" << vectorName << "\"" << std::endl;

    RCP<MultiVector> onesVector = MultiVectorFactory::Build(transferOp->getDomainMap(), 1);
    onesVector->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    RCP<MultiVector> rowSumVector = MultiVectorFactory::Build(transferOp->getRangeMap(), 1);
    transferOp->apply(*onesVector, *rowSumVector);
    transferOp->apply(*fineVector, *coarseVector);

    if (vectorName == "Coordinates")
      TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"Use CoordinatesTransferFactory to transfer coordinates instead of MultiVectorTransferFactory.");

    Set<RCP<MultiVector> >(coarseLevel, vectorName, coarseVector);

  } // Build

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayRCP<Scalar> MultiVectorTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::expandCoordinates(ArrayRCP<SC> coordinates, LocalOrdinal blksize) {
    if (blksize == 1)
      return coordinates;

    ArrayRCP<SC> expandCoord(coordinates.size()*blksize); //TODO: how to avoid automatic initialization of the vector? using arcp()?

    for(int i=0; i<coordinates.size(); i++) {
      for(int j=0; j< blksize; j++) {
        expandCoord[i*blksize + j] = coordinates[i];
      }
    }
    return expandCoord;

  } // expandCoordinates

} // namespace MueLu

#endif // MUELU_MULTIVECTORTRANSFER_FACTORY_DEF_HPP
