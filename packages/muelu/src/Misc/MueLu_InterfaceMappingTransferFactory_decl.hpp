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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DECL_HPP_
#define MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DECL_HPP_

#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class InterfaceMappingTransferFactory
  @brief Transfer mapping data for interface aggregation to the coarse level

  The fine level knows about the mapping on the coarse level.
  This factory just pushes that pre-computed information to the coarse level.

  ## Input/output ##

  ### User parameters ###
  Parameter | type | default | master.xml | validated | requested | description
  ----------|------|---------|:----------:|:---------:|:---------:|------------
  CoarseDualNodeID2PrimalNodeID | Factory | null |   | * | * | Generating factory of the coarse dual-to-primal node mapping

  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see InterfaceAggregationFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see InterfaceAggregationFactory::DeclareInput).
*/
template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class InterfaceMappingTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_INTERFACEMAPPINGTRANSFERFACTORY_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"
 public:
  //! Constructor.
  InterfaceMappingTransferFactory() = default;

  //! Destructor.
  ~InterfaceMappingTransferFactory() {}

  RCP<const ParameterList> GetValidParameterList() const override;
  void DeclareInput(Level &fineLevel, Level &coarseLevel) const override;
  void Build(Level &fineLevel, Level &coarseLevel) const override;
};

}  // namespace MueLu
#define MUELU_INTERFACEMAPPINGTRANSFERFACTORY_SHORT
#endif /* MUELU_INTERFACEMAPPINGTRANSFERFACTORY_DECL_HPP_ */
