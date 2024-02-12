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
#ifndef MUELU_SCALEDNULLSPACEFACTORY_DEF_HPP
#define MUELU_SCALEDNULLSPACEFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ScaledNullspaceFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ScaledNullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<std::string>("Fine level nullspace", "Nullspace", "Variable name which is used to store null space multi vector on the finest level (default=\"Nullspace\").");

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the fine level matrix (only needed if default null space is generated)");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory of the fine level null space");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ScaledNullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  // Scaled Nullspace always needs A & a Nullspace from somewhere
  Input(currentLevel, "A");
  Input(currentLevel, "Nullspace");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ScaledNullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<MultiVector> nullspace, tentativeNullspace;

  // TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.GetLevelID() != 0, Exceptions::RuntimeError, "MueLu::ScaledNullspaceFactory::Build(): ScaledNullspaceFactory can be used for finest level (LevelID == 0) only.");

  if (currentLevel.GetLevelID() == 0) {
    if (currentLevel.IsAvailable("Nullspace", NoFactory::get())) {
      // When a fine nullspace have already been defined by user using Set("Nullspace", ...) or
      // Set("Nullspace1", ...), we use it.
      tentativeNullspace = currentLevel.Get<RCP<MultiVector> >("Nullspace", NoFactory::get());
    } else {
      // A User "Nullspace" (nspName) is not available (use factory)
      tentativeNullspace = currentLevel.Get<RCP<MultiVector> >("Nullspace", GetFactory("Nullspace").get());

    }  // end if "Nullspace" not available
  } else {
    tentativeNullspace = currentLevel.Get<RCP<MultiVector> >("Nullspace", GetFactory("Nullspace").get());
  }

  // Scale!
  RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");

  // determine numPDEs
  LocalOrdinal numPDEs = 1;
  if (A->IsView("stridedMaps") == true) {
    Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap()) == Teuchos::null, Exceptions::BadCast, "MueLu::ScaledNullspaceFactory::Build: cast to strided row map failed.");
    numPDEs = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
    oldView = A->SwitchToView(oldView);
  }

  GetOStream(Runtime1) << "ScaledNullspaceFactory: Generating scaled nullspace, blocksize = " << tentativeNullspace->getNumVectors() << std::endl;
  nullspace                      = MultiVectorFactory::Build(tentativeNullspace->getMap(), tentativeNullspace->getNumVectors());
  *nullspace                     = *tentativeNullspace;  // Copy the tentative nullspace
  RCP<MultiVector> blockDiagonal = MultiVectorFactory::Build(A->getDomainMap(), numPDEs);

  Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::extractBlockDiagonal(*A, *blockDiagonal);
  Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::inverseScaleBlockDiagonal(*blockDiagonal, true, *nullspace);

  // provide "Scaled Nullspace" variable on current level (used by TentativePFactory)
  Set(currentLevel, "Scaled Nullspace", nullspace);

}  // Build

}  // namespace MueLu

#endif  // MUELU_SCALEDNULLSPACEFACTORY_DEF_HPP
