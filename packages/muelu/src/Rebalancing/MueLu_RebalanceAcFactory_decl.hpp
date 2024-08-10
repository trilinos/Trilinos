// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCEACFACTORY_DECL_HPP
#define MUELU_REBALANCEACFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"

namespace MueLu {
/*!
  @class RebalanceAcFactory
  @brief Factory for building coarse matrices.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RebalanceAcFactory : public TwoLevelFactoryBase {
#undef MUELU_REBALANCEACFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  RebalanceAcFactory();

  virtual ~RebalanceAcFactory();

  RCP<const ParameterList> GetValidParameterList() const;
  //@}

  //! @name Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level &fineLevel, Level &coarseLevel) const;
  //@}

  //@{
  /*! @brief Add rebalancing factory in the end of list of rebalancing factories in RebalanceAcFactory.

  Rebalancing factories are derived from SingleLevelFactoryBase and rebalance the underlaying object
  (e.g. map, vector,...) to fit to the rebalanced maps.
  */
  void AddRebalanceFactory(const RCP<const FactoryBase> &factory);

  //! Returns number of transfer factories.
  size_t NumRebalanceFactories() const { return rebalanceFacts_.size(); }

  //@}

 private:
  //! list of user-defined rebalancing Factories
  std::vector<RCP<const FactoryBase> > rebalanceFacts_;

};  // class RebalanceAcFactory

}  // namespace MueLu

#define MUELU_REBALANCEACFACTORY_SHORT
#endif  // MUELU_REBALANCEACFACTORY_DECL_HPP
