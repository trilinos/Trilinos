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
#ifndef MUELU_PLANEDETECTIONFACTORY_DEF_HPP
#define MUELU_PLANEDETECTIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
//#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_PlaneDetectionFactory_decl.hpp"

#include "MueLu_FactoryManager.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> PlaneDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

// #define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
//     SET_VALID_ENTRY("linedetection: orientation");
//     SET_VALID_ENTRY("linedetection: num layers");
// #undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",               Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",     Teuchos::null, "Generating factory for coorindates");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void PlaneDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");

    // The factory needs the information about the number of z-layers. While this information is
    // provided by the user for the finest level, the factory itself is responsible to provide the
    // corresponding information on the coarser levels. Since a factory cannot be dependent on itself
    // we use the NoFactory class as generator class, but remove the UserData keep flag, such that
    // "NumZLayers" is part of the request/release mechanism.
    // Please note, that this prevents us from having several (independent) CoarsePFactory instances!
    // TODO: allow factory to dependent on self-generated data for TwoLevelFactories -> introduce ExpertRequest/Release in Level
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void PlaneDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Line detection (Ray style)", currentLevel);

    RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");
    const LO dofsPerNode = A->GetFixedBlockSize();
    const LO numNodes    = A->getNodeNumRows() / dofsPerNode;

    std::cout << "Nodes in the problem: " << numNodes << std::endl;

  } // Build

} //namespace MueLu

#endif // MUELU_PLANEDETECTIONFACTORY_DEF_HPP
