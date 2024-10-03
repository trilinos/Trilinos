// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_GENERICRFACTORY_DECL_HPP
#define MUELU_GENERICRFACTORY_DECL_HPP

/*
 * MueLu_GenericRFactory.hpp
 *
 *  Created on: 20.09.2011
 *      Author: tobias
 */

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_GenericRFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_PFactory_fwd.hpp"

namespace MueLu {

/*!
  @class GenericRFactory class.
  @brief Factory for building restriction operators using a prolongator factory
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class GenericRFactory : public TwoLevelFactoryBase {
#undef MUELU_GENERICRFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  GenericRFactory() {}

  //! Destructor.
  virtual ~GenericRFactory() {}
  //@}

  //! Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level &fineLevel, Level &coarseLevel) const;

  //@}

};  // class GenericRFactory

}  // namespace MueLu

#define MUELU_GENERICRFACTORY_SHORT
#endif  // MUELU_GENERICRFACTORY_DECL_HPP
