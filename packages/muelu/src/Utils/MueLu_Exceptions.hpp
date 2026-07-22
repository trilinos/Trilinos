// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_EXCEPTIONS_HPP
#define MUELU_EXCEPTIONS_HPP

#include <Teuchos_Exceptions.hpp>

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {
namespace Exceptions {

//! Exception indicating invalid cast attempted
class BadCast : public Teuchos::ExceptionBase {
 public:
  BadCast(const std::string& what_arg);
  ~BadCast();
};

//! Exception throws when you call an unimplemented method of MueLu
/** Mainly use for development in progress. **/
class NotImplemented : public Teuchos::ExceptionBase {
 public:
  NotImplemented(const std::string& what_arg);
  ~NotImplemented();
};

//! Exception throws to report errors in the internal logical of the program.
class RuntimeError : public Teuchos::ExceptionBase {
 public:
  RuntimeError(const std::string& what_arg);
  ~RuntimeError();
};

//! Exception throws to report overflows.
class Overflow : public Teuchos::ExceptionBase {
 public:
  Overflow(const std::string& what_arg);
  ~Overflow();
};

//! Exception throws to report incompatible objects (like maps).
class Incompatible : public Teuchos::ExceptionBase {
 public:
  Incompatible(const std::string& what_arg);
  ~Incompatible();
};

//! Exception throws to report data dependency problems between factories.
class DependencyError : public Teuchos::ExceptionBase {
 public:
  DependencyError(const std::string& what_arg);
  ~DependencyError();
};

//! Exception throws to report invalid user entry
class InvalidArgument : public Teuchos::ExceptionBase {
 public:
  InvalidArgument(const std::string& what_arg);
  ~InvalidArgument();
};

}  // namespace Exceptions
}  // namespace MueLu

#define MUELU_TPETRA_ETI_EXCEPTION(cl, obj, go) TEUCHOS_TEST_FOR_EXCEPTION(1, ::MueLu::Exceptions::BadCast, "Problem in " #cl "! Cannot create new object " #obj " with GO=" #go ". MueLu has been compiled with Tpetra enabled bug GO!=" #go ". Please add TPETRA_INST_INT_INT to your configuration.");

#endif  // ifndef MUELU_EXCEPTIONS_HPP
