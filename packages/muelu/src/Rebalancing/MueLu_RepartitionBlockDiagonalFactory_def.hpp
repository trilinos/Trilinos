// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2013 Sandia Corporation
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
#ifndef MUELU_REPARTITIONBLOCKDIAGONALFACTORY_DEF_HPP_
#define MUELU_REPARTITIONBLOCKDIAGONALFACTORY_DEF_HPP_

#include "MueLu_RepartitionBlockDiagonalFactory_decl.hpp"

#include <Teuchos_Utils.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

 template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
 RCP<const ParameterList> RepartitionBlockDiagonalFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >   ("A",             Teuchos::null, "Factory of the matrix A");

    return validParamList;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RepartitionBlockDiagonalFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level & currentLevel) const {
    Input(currentLevel, "A");
  } //DeclareInput()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RepartitionBlockDiagonalFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);
    typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>  BlockCrs;

    RCP<Matrix> originalA = Get< RCP<Matrix> >(currentLevel, "A");    
    RCP<BlockCrs> A = Teuchos::rcp_dynamic_cast<BlockCrs>(originalA);

    //    RCP<BlockCrs> A = Get< RCP<BlockCrs> >(currentLevel, "A");    
    TEUCHOS_TEST_FOR_EXCEPTION(A==Teuchos::null, Exceptions::BadCast, "MueLu::RepartitionBlockDiagonalFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");

    // Build the block diagonal
    RCP<BlockCrs> DiagonalMatrix = Teuchos::rcp(new BlockCrs(A->getBlockedRangeMap(),A->getBlockedDomainMap(),0,Xpetra::StaticProfile));    
    for(size_t i=0; i< A->Rows(); i++)
      DiagonalMatrix->setMatrix(i,i,A->getMatrix(i,i));

    Set(currentLevel,"A",Teuchos::rcp_dynamic_cast<Matrix>(DiagonalMatrix));

    
  } //Build()

} // end namespace MueLu

#endif /* MUELU_REPARTITIONBLOCKDIAGONALFACTORY_DEF_HPP_ */
