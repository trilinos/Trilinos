// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_EDGEPROLONGATORPATTERNFACTORY_DECL_HPP
#define MUELU_EDGEPROLONGATORPATTERNFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class EdgeProlongatorFactory class.
  @brief Factory for building nonzero patterns for energy minimization.
  @ingroup MueLuTransferClasses
  */

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class EdgeProlongatorPatternFactory : public TwoLevelFactoryBase {
#undef MUELU_EDGEPROLONGATORPATTERNFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! @brief Constructor.
  EdgeProlongatorPatternFactory() {}

  //! Destructor.
  virtual ~EdgeProlongatorPatternFactory() = default;

  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! @name Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds nonzero pattern (graph) and returns it in <tt>coarseLevel</tt>.
    */
  void Build(Level& fineLevel, Level& coarseLevel) const;

  //@}

};  // class EdgeProlongatorFactory

}  // namespace MueLu

#define MUELU_EDGEPROLONGATORPATTERNFACTORY_SHORT
#endif  // MUELU_EDGEPROLONGATORPATTERNFACTORY_DECL_HPP
