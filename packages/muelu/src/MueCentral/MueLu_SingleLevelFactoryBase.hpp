// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SINGLELEVELFACTORY_HPP
#define MUELU_SINGLELEVELFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_Factory.hpp"
#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class SingleLevelFactoryBase
  @brief Base class for factories that use one level (currentLevel).

  @ingroup MueLuBaseClasses
*/
class SingleLevelFactoryBase : public Factory {
 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  SingleLevelFactoryBase();

  //! Destructor.
  virtual ~SingleLevelFactoryBase();

  //@}

  //! @name Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */
  virtual void DeclareInput(Level& currentLevel) const = 0;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  virtual void Build(Level& currentLevel) const = 0;

  //!
  virtual void CallBuild(Level& requestedLevel) const;

  //!
  virtual void CallDeclareInput(Level& requestedLevel) const;
  //@}

};  // class SingleLevelFactoryBase

}  // namespace MueLu

#define MUELU_SINGLELEVELFACTORY_SHORT
#endif  // ifndef MUELU_SINGLELEVELFACTORY_HPP

// TODO: code factorization between SingleLevelFactoryBase and TwoLevelFactoryBase
