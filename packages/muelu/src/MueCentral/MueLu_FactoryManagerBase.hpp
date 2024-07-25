// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_FACTORYMANAGERBASE_HPP
#define MUELU_FACTORYMANAGERBASE_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {

/*!
   @class FactoryManagerBase
   @brief Class that provides default factories within Needs class.
   @ingroup MueLuBaseClasses
*/
class FactoryManagerBase : public BaseClass {
 public:
  //@{ Constructors/Destructors.
  FactoryManagerBase()
    : bIgnoreUserData_(false) {}

  //! Destructor.
  virtual ~FactoryManagerBase() {}

  //@}

  //@{ Get/Set functions.

  //! Get
  // Return ref because user also give ref to the Hierarchy.
  const virtual RCP<const FactoryBase> GetFactory(const std::string& varName) const = 0;
  //@}

  //! Check
  // Return true if Factory associated with varName is registered
  virtual bool hasFactory(const std::string& varName) const = 0;

  // Free temporarily hold data at the end of Hierarchy::Setup()
  // This method is const because the clean concerns only mutable data.
  virtual void Clean() const {}  // TODO: should be used inside of MueLu::Hierarchy

#ifdef HAVE_MUELU_DEBUG
  virtual void ResetDebugData() const = 0;
#endif

  //! get IgnoreUserData flag
  bool IgnoreUserData() const { return bIgnoreUserData_; }

  //! set IgnoreUserData flag
  void SetIgnoreUserData(bool bIgnoreUserData = false) { bIgnoreUserData_ = bIgnoreUserData; }

 private:
  //! boolean flag that controls behaviour of Level::GetFactory
  //! if bIgnoreUserData == true,  the Level::GetFactory function always asks the Factory manager for a valid factory given a variable name
  //! if bIgnoreUserData == false, the Level::GetFactory prefers user-provided data for a variable name if available. Otherwise the factory manager is asked for a valid factory
  //! default: bIgnoreUserData = false;
  bool bIgnoreUserData_;

};  // class FactoryManagerBase

}  // namespace MueLu

#define MUELU_FACTORYMANAGERBASE_SHORT
#endif  // ifndef MUELU_FACTORYMANAGERBASE_HPP
