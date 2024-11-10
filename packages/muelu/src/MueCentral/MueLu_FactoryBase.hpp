// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_FACTORYBASE_HPP
#define MUELU_FACTORYBASE_HPP

#include "MueLu_config.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class FactoryBase
  @brief Base class for factories (e.g., R, P, and A_coarse).
  @ingroup MueLuBaseClasses
*/
class FactoryBase : public virtual BaseClass {
 public:
  //@{ Constructors/Destructors.

  //! Constructor.
  FactoryBase()
    : id_(FactoryBase::GenerateUniqueId()) {}

  //! Destructor.
  virtual ~FactoryBase() {}
  //@}

  //@{
  //! @name Build methods.

  virtual void CallBuild(Level& requestedLevel) const = 0;

  virtual void CallDeclareInput(Level& requestedLevel) const = 0;
  //@}

  //@{
  //! @name Access factory properties

  /// return unique factory id
  int GetID() const { return id_; };

    //@}

#ifdef HAVE_MUELU_DEBUG
  virtual void ResetDebugData() const = 0;
#endif

 private:
  static int GenerateUniqueId();

  const int id_;

};  // class FactoryBase

}  // namespace MueLu

#define MUELU_FACTORYBASE_SHORT
#endif  // ifndef MUELU_FACTORYBASE_HPP
