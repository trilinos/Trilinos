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
#ifndef MUELU_TRANSPFACTORY_DEF_HPP
#define MUELU_TRANSPFACTORY_DEF_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

#include <Xpetra_Operator.hpp>

#include "MueLu_TransPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TransPFactory(RCP<const FactoryBase> PFact)
    : PFact_(PFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TransPFactory() {}
 
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    coarseLevel.DeclareInput("P", PFact_.get(), this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    return BuildR(fineLevel,coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildR(Level & fineLevel, Level & coarseLevel) const {
    FactoryMonitor m(*this, "Transpose P", coarseLevel);

    Teuchos::OSTab tab(this->getOStream());
    Teuchos::ParameterList matrixList;
    RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PFact_.get());

    //doesn't work -- bug in EpetraExt?
    //RCP<Operator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getRangeMap(),matrixList);
    //      RCP<CrsOperator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getDomainMap(),matrixList);
    //RCP<Operator> R = Utils::TwoMatrixMultiply(P,true,I,false); //doesn't work -- bug in EpetraExt?
    //      RCP<Operator> R = Utils::TwoMatrixMultiply(I,false,P,true);

    RCP<Operator> R = Utils2::Transpose(P,true);

    coarseLevel.Set("R", R, this);
    
    ///////////////////////// EXPERIMENTAL
    if(P->IsView("stridedMaps")) R->CreateView("stridedMaps", P, true);  
    ///////////////////////// EXPERIMENTAL

  } //BuildR

} //namespace MueLu

#endif // MUELU_TRANSPFACTORY_DEF_HPP
