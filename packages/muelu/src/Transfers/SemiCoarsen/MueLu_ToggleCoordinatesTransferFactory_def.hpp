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
#ifndef MUELU_TOGGLECOORDINATESTRANSFER_FACTORY_DEF_HPP
#define MUELU_TOGGLECOORDINATESTRANSFER_FACTORY_DEF_HPP

#include "MueLu_ToggleCoordinatesTransferFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ToggleCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("Chosen P", Teuchos::null, "Name of TogglePFactory this ToggleCoordinatesTransferFactory is connected to. Parameter provides information which execution path (prolongator) has been chosen.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ToggleCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* fineLevel */, Level& coarseLevel) const {
  const ParameterList& pL = GetParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(!pL.isParameter("Chosen P"), Exceptions::RuntimeError, "MueLu::ToggleCoordinatesTransferFactory::DeclareInput: You have to set the 'Chosen P' parameter to a factory name of type TogglePFactory. The ToggleCoordinatesTransferFactory must be used together with a TogglePFactory!");
  Input(coarseLevel, "Chosen P");
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = coordFacts_.begin(); it != coordFacts_.end(); ++it) {
    coarseLevel.DeclareInput("Coordinates", (*it).get(), this);  // request/release coarse coordinates
    (*it)->CallDeclareInput(coarseLevel);                        // request dependencies
  }
  hasDeclaredInput_ = true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ToggleCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& /* fineLevel */, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Coordinate transfer toggle", coarseLevel);

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> xdMV;

  TEUCHOS_TEST_FOR_EXCEPTION(coordFacts_.size() != 2, Exceptions::RuntimeError, "MueLu::TogglePFactory::Build: ToggleCoordinatesTransferFactory needs two different transfer operator strategies for toggling.");

  int chosenP = Get<int>(coarseLevel, "Chosen P");
  GetOStream(Runtime1) << "Transfer Coordinates" << chosenP << " to coarse level" << std::endl;
  RCP<xdMV> coarseCoords = coarseLevel.Get<RCP<xdMV> >("Coordinates", (coordFacts_[chosenP]).get());
  Set(coarseLevel, "Coordinates", coarseCoords);

  // loop through all coord facts and check whether the coarse coordinates are available.
  // This is the coarse coordinate transfer factory which belongs to the execution path
  // chosen by the TogglePFactory
  /*RCP<xdMV> coarseCoords = Teuchos::null;
  for(size_t t=0; t<coordFacts_.size(); ++t) {
    bool bIsAv = coarseLevel.IsAvailable("Coordinates",(coordFacts_[t]).get());
    std::cout << "Coordinates generated by " << (coordFacts_[t]).get() << " available? " << bIsAv << std::endl;
    if ( coarseLevel.IsAvailable("Coordinates",(coordFacts_[t]).get()) ) {
      GetOStream(Runtime1) << "Choose factory " << t << " (" << (coordFacts_[t]).get() << ")" << std::endl;
      coarseCoords = coarseLevel.Get< RCP<xdMV> >("Coordinates",(coordFacts_[t]).get());
      Set(coarseLevel, "Coordinates", coarseCoords);
    }
  }*/

  // Release dependencies of all coordinate transfer factories
  for (size_t t = 0; t < coordFacts_.size(); ++t) {
    coarseLevel.Release(*(coordFacts_[t]));
  }

  // TODO: exception if coarseCoords == Teuchos::null
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ToggleCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddCoordTransferFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::ToggleCoordinatesTransferFactory::AddCoordTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. Make sure you provide the factory which generates the coarse coordinates. Usually this is a coordinate transfer factory."
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::ToggleCoordinatesTransferFactory::AddCoordTransferFactory: Factory is being added after we have already declared input");
  coordFacts_.push_back(factory);
}
}  // namespace MueLu

#endif  // MUELU_TOGGLECOORDINATESTRANSFER_FACTORY_DEF_HPP
