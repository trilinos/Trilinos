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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_


#include "MueLu_Monitor.hpp"

#include "MueLu_UnsmooshFactory_decl.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::UnsmooshFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("DofStatus",          Teuchos::null, "Generating factory for dofStatus array (usually the VariableDofLaplacdianFactory)");

    validParamList->set< int  >                  ("maxDofPerNode", 1,     "Maximum number of DOFs per node");
    validParamList->set< bool >                  ("fineIsPadded" , false, "true if finest level input matrix is padded");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    //const ParameterList& pL = GetParameterList();
    Input(currentLevel, "A");
    Input(currentLevel, "DofStatus");
    //Input(currentLevel, "Coordinates");

    /*if (currentLevel.IsAvailable("DofPresent", NoFactory::get())) {
      currentLevel.DeclareInput("DofPresent", NoFactory::get(), this);
    }*/
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);
    typedef Teuchos::ScalarTraits<SC> STS;

    const ParameterList  & pL = GetParameterList();

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    Teuchos::RCP< const Teuchos::Comm< int > > comm = A->getRowMap()->getComm();
    Xpetra::UnderlyingLib lib = A->getRowMap()->lib();

    Teuchos::Array<char> dofStatus = Get<Teuchos::Array<char> >(currentLevel, "DofStatus");

    int maxDofPerNode = pL.get<int> ("maxDofPerNode");
    bool fineIsPadded = pL.get<bool>("fineIsPadded");

    //Set(currentLevel,"A",lapMat);
  }


} /* MueLu */


#endif /* PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_ */
