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
#ifndef MUELU_PATTERNFACTORY_DEF_HPP
#define MUELU_PATTERNFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_PatternFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory for the matrix");
    validParamList->set< RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory for the matrix providing nonzero graph");
    validParamList->set<int>                     ("k",             0, "Polynomial degree: the resulting pattern is A^k*P [default = 0]");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(coarseLevel, "P");

    const ParameterList& pL = GetParameterList();
    if (pL.get<int>("k") > 0)
      Input(fineLevel, "A");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Ppattern", coarseLevel);

    RCP<Matrix> P = Get< RCP<Matrix> >(coarseLevel, "P");

    const ParameterList& pL = GetParameterList();
    int k = pL.get<int>("k");

    if (k > 0) {
      RCP<Matrix> A = Get< RCP<Matrix> >(fineLevel, "A");
      RCP<Matrix> AP;

      bool doFillComplete  = true;
      bool optimizeStorage = true;
      bool allowMLMultiply = false;

      for (int i = 0; i < k; i++) {
        AP = Utils::Multiply(*A, false, *P, false, doFillComplete, optimizeStorage, allowMLMultiply);
        P.swap(AP);
      }
    }

    Set(coarseLevel, "Ppattern", P->getCrsGraph());
  }


} // namespace MueLu

#endif // MUELU_PATTERNFACTORY_DEF_HPP
