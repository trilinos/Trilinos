// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_USERAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_USERAGGREGATIONFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_UserAggregationFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class UserAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_USERAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  UserAggregationFactory(){};

  //! Destructor.
  virtual ~UserAggregationFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Set/get methods.
  //@{

  // Options shared by all aggregation algorithms

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level &currentLevel) const;

  //@}

 private:
};  // class UserAggregationFactory

}  // namespace MueLu

#define MUELU_USERAGGREGATIONFACTORY_SHORT
#endif /* MUELU_USERAGGREGATIONFACTORY_DECL_HPP_ */
