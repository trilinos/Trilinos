// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
