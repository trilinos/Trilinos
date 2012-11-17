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
#ifndef MUELU_NULLSPACEPRESMOOTHFACTORY_DEF_HPP
#define MUELU_NULLSPACEPRESMOOTHFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_NullspacePresmoothFactory_decl.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NullspacePresmoothFactory()
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~NullspacePresmoothFactory() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "Nullspace");

    if (currentLevel.GetLevelID() == 0)
      Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void NullspacePresmoothFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Nullspace presmooth factory", currentLevel);

    RCP<MultiVector> nullspace;

    RCP<MultiVector> newB;
    if (currentLevel.GetLevelID() == 0) {
      RCP<Matrix>      A = Get< RCP<Matrix> >     (currentLevel, "A");
      RCP<MultiVector> B = Get< RCP<MultiVector> >(currentLevel, "Nullspace");
      newB = MultiVectorFactory::Build(B->getMap(), B->getNumVectors());

      Teuchos::ArrayRCP<SC> D = Utils::GetMatrixDiagonal(*A);

      SC damping = 4./3;
      damping /= Utils::PowerMethod(*A, true, (LO) 10, 1e-4);

      A->apply(*B, *newB, Teuchos::NO_TRANS);

      size_t numVec  = newB->getNumVectors();
      LO numElements = newB->getLocalLength();
      for (size_t j = 0; j < numVec; j++) {
        Teuchos::ArrayRCP<const SC> Bj    = B->getData(j);
        Teuchos::ArrayRCP<SC>       newBj = newB->getDataNonConst(j);

        for (LO i = 0; i < numElements; i++) {
          if (newBj[i]/D[i])
            std::cout << "Update [" << i << "]: " << -damping*newBj[i]/D[i] << std::endl;
          newBj[i] = Bj[i] - damping*newBj[i]/D[i];
        }
      }
    } else {
      newB = Get< RCP<MultiVector> >(currentLevel, "Nullspace");
    }

    // provide "Nullspace" variable on current level
    Set(currentLevel, "Nullspace", newB);

  } // Build

} //namespace MueLu

#endif // MUELU_NULLSPACEPRESMOOTHFACTORY_DEF_HPP
