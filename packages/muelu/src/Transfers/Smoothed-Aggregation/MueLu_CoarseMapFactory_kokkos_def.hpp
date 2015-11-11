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
#ifndef MUELU_COARSEMAPFACTORY_KOKKOS_DEF_HPP_
#define MUELU_COARSEMAPFACTORY_KOKKOS_DEF_HPP_

#include "MueLu_CoarseMapFactory_kokkos_decl.hpp"

#include <Teuchos_Array.hpp>

#include <Xpetra_MultiVector.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> CoarseMapFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Generating factory for aggregates.");
    validParamList->set<RCP<const FactoryBase> >("Nullspace",  Teuchos::null, "Generating factory for null space.");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CoarseMapFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "Aggregates");
    Input(currentLevel, "Nullspace");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CoarseMapFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    auto aggregates = Get<RCP<Aggregates_kokkos> >(currentLevel, "Aggregates");
    auto nullspace  = Get<RCP<MultiVector> >      (currentLevel, "Nullspace");

    auto map = aggregates->GetMap();

    const GO     numAggs = aggregates->GetNumAggregates();
    const size_t NSDim   = nullspace->getNumVectors();

    GlobalOrdinal numCoarseDofs = numAggs * NSDim;

    std::vector<size_t> stridingInfo(1);
    stridingInfo[0] = NSDim;

    auto INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    RCP<const Map> coarseMap = StridedMapFactory::Build(map->lib(), INVALID, numCoarseDofs, map->getIndexBase(), stridingInfo, map->getComm());

    Set(currentLevel, "CoarseMap", coarseMap);
  }

} //namespace MueLu

#endif /* MUELU_COARSEMAPFACTORY_KOKKOS_DEF_HPP_ */
