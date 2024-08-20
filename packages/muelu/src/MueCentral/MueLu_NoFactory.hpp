// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_NOFACTORY_HPP
#define MUELU_NOFACTORY_HPP

#include <algorithm>            // for swap
#include "Teuchos_RCPDecl.hpp"  // for RCP
#include "Teuchos_RCP.hpp"      // for RCP::RCP<T>, RCP::operator=, etc
#include "MueLu_config.hpp"     // for HAVE_MUELU_DEBUG
#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class NoFactory class.
  @brief NoFactory that is used for data stored in level class for that no generating factory is available/necessary.

  This should be used as the "generating" factory for user-data.  Uses Singleton pattern.
*/
class NoFactory : public FactoryBase {
  //! Constructor.
  NoFactory();

 public:
  //! Destructor.
  virtual ~NoFactory();

  //! Implementation of FactoryBase interface
  //@{

  //!
  void CallBuild(Level& requestedLevel) const;

  //!
  void CallDeclareInput(Level& /* requestedLevel */) const;

  //@}

  //! Static Get() functions
  //@{

  //!
  static const RCP<const NoFactory> getRCP();

  //!
  static const NoFactory* get();

  //@}
#ifdef HAVE_MUELU_DEBUG
  void ResetDebugData() const {}
#endif

 private:
  static RCP<const NoFactory> noFactory_;  // static NoFactory instance for user defined "factories"

};  // class NoFactory

}  // namespace MueLu

#endif  // MUELU_NOFACTORY_HPP
