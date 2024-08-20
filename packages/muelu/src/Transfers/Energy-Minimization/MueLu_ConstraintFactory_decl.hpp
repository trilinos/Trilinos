// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CONSTRAINTFACTORY_DECL_HPP
#define MUELU_CONSTRAINTFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Constraint_fwd.hpp"
#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class ConstraintFactory class.
  @brief Factory for building the constraint operator

  Factory for creating the constraint operator.

  @ingroup MueLuTransferClasses
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class ConstraintFactory : public TwoLevelFactoryBase {
#undef MUELU_CONSTRAINTFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.
    User can supply a factory for generating the nonzero pattern. The nullspace vectors (both fine and coarse) will
    be taken from the corresponding level factories
    */
  ConstraintFactory() {}

  //! Destructor.
  virtual ~ConstraintFactory() {}

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

    Builds Constraint and returns it in <tt>coarseLevel</tt>.
    */
  void Build(Level& fineLevel, Level& coarseLevel) const;

  //@}
};  // class ConstraintFactory

}  // namespace MueLu

#define MUELU_CONSTRAINTFACTORY_SHORT
#endif  // MUELU_CONSTRAINTFACTORY_DECL_HPP
