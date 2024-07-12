// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_STRUCTUREDLINEDETECTIONFACTORY_DECL_HPP
#define MUELU_STRUCTUREDLINEDETECTIONFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_StructuredLineDetectionFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class StructuredLineDetectionFactory class.
  @brief Factory building line detection information on structured meshes
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class StructuredLineDetectionFactory : public SingleLevelFactoryBase {
#undef MUELU_STRUCTUREDLINEDETECTIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  StructuredLineDetectionFactory() {}

  //! Destructor.
  virtual ~StructuredLineDetectionFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds line detection information and stores it in currentLevel
    */
  void Build(Level& currentLevel) const;

  //@}

 private:
};  // class StructuredLineDetectionFactory

}  // namespace MueLu

#define MUELU_STRUCTUREDLINEDETECTIONFACTORY_SHORT
#endif  // MUELU_STRUCTUREDLINEDETECTIONFACTORY_DECL_HPP
