// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_NULLSPACEPRESMOOTHFACTORY_DECL_HPP
#define MUELU_NULLSPACEPRESMOOTHFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class NullspacePresmoothFactory : public SingleLevelFactoryBase {
#undef MUELU_NULLSPACEPRESMOOTHFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  NullspacePresmoothFactory() {}

  //! Destructor
  virtual ~NullspacePresmoothFactory() {}

  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! @name Input
  //@{
  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level &currentLevel) const;

  //@}

 private:
};  // class NullspacePresmoothFactory

}  // namespace MueLu

#define MUELU_NULLSPACEPRESMOOTHFACTORY_SHORT
#endif  // MUELU_NULLSPACEPRESMOOTHFACTORY_DECL_HPP
