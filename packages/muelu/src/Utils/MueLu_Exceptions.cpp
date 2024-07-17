// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_Exceptions.hpp>
#include <MueLu_Exceptions.hpp>

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {
namespace Exceptions {

BadCast::BadCast(const std::string& what_arg)
  : Teuchos::ExceptionBase(what_arg) {}

NotImplemented::NotImplemented(const std::string& what_arg)
  : Teuchos::ExceptionBase(what_arg) {}

RuntimeError::RuntimeError(const std::string& what_arg)
  : Teuchos::ExceptionBase(what_arg) {}

Overflow::Overflow(const std::string& what_arg)
  : Teuchos::ExceptionBase(what_arg) {}

Incompatible::Incompatible(const std::string& what_arg)
  : Teuchos::ExceptionBase(what_arg) {}

DependencyError::DependencyError(const std::string& what_arg)
  : Teuchos::ExceptionBase(what_arg) {}

InvalidArgument::InvalidArgument(const std::string& what_arg)
  : Teuchos::ExceptionBase(what_arg) {}

}  // namespace Exceptions
}  // namespace MueLu
