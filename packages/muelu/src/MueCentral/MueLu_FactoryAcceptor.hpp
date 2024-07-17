// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_FACTORYACCEPTOR_HPP
#define MUELU_FACTORYACCEPTOR_HPP

#include <string>

#include "Teuchos_RCP.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"

namespace MueLu {

class FactoryAcceptor {
 public:
  virtual ~FactoryAcceptor() {}

  //@{
  //! Configuration

  //! SetFactory is for expert users only. To change configuration of the preconditioner, use a factory manager.
  virtual void SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) = 0;

  virtual const RCP<const FactoryBase> GetFactory(const std::string& varName) const = 0;

  // SetParameterList(...);

  // GetParameterList(...);

  //@}

};  // class FactoryAcceptor

}  // namespace MueLu

#define MUELU_FACTORYACCEPTOR_SHORT
#endif  // ifndef MUELU_FACTORYACCEPTOR_HPP
