// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EXCEPTIONS_HPP
#define XPETRA_EXCEPTIONS_HPP

#include <exception>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_Assert.hpp>

#include "Xpetra_ConfigDefs.hpp"

// Dynamically cast 'obj' to 'type newObj'. newObj is declared inside of the macro.
// If the dynamic cast failed, throw an exception of type Xpetra::Exception::Bad_Cast (using the message exceptionMsg).
#define XPETRA_COMMA ,
#define XPETRA_DYNAMIC_CAST(type, obj, newObj, exceptionMsg)                                                                             \
  type* newObj##_pt = dynamic_cast<type*>(&obj);                                                                                         \
  TEUCHOS_TEST_FOR_EXCEPTION(newObj##_pt == NULL, Xpetra::Exceptions::BadCast, "Cannot cast '" #obj "' to a " #type ". " #exceptionMsg); \
  type& newObj = *newObj##_pt;

// Dynamically cast the RCP 'obj' to 'RCP<type> newObj'. newObj is declared inside of the macro.
// If the dynamic cast failed, throw an exception of type Xpetra::Exception::Bad_Cast (using the message exceptionMsg).
#define XPETRA_RCP_DYNAMIC_CAST(type, obj, newObj, exceptionMsg)  \
  const RCP<type>& newObj = Teuchos::rcp_dynamic_cast<type>(obj); \
  TEUCHOS_TEST_FOR_EXCEPTION(newObj == Teuchos::null, Xpetra::Exceptions::BadCast, "Cannot cast '" #obj "' to a " #type ". " #exceptionMsg);

#ifdef HAVE_XPETRA_EPETRA
#define XPETRA_FACTORY_ERROR_IF_EPETRA(lib) \
  if ((lib) == ::Xpetra::UseEpetra)         \
    TEUCHOS_TEST_FOR_EXCEPTION(1, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and LocalOrdinal=GlobalOrdinal=int on a serial node");
#else
#define XPETRA_FACTORY_ERROR_IF_EPETRA(lib)
#endif

#define XPETRA_FACTORY_END TEUCHOS_TEST_FOR_EXCEPTION(1, ::Xpetra::Exceptions::BadCast, "Unknown map->lib() type. Did you compile with Epetra and Tpetra support?");

#define XPETRA_TPETRA_ETI_EXCEPTION(cl, obj, go, node) TEUCHOS_TEST_FOR_EXCEPTION(1, ::Xpetra::Exceptions::BadCast, "Problem in " #cl "! Cannot create new object " #obj " with GO=" #go " and NO=" #node ". Xpetra has been compiled with Tpetra enabled but GO=" #go " disabled and/or NO=" #node " disabled. Please check your Trilinos configuration.");

namespace Xpetra {
namespace Exceptions {

//! Exception indicating invalid cast attempted
/** For instance, this exception is throw when you try to mix Epetra and Tpetra objects: for example, a Xpetra::EpetraCrsMatrix cannot be built from a Xpetra::TpetraMap. **/
class BadCast : public Teuchos::ExceptionBase {
 public:
  BadCast(const std::string& what_arg)
    : Teuchos::ExceptionBase(what_arg) {}
};

//! Exception throws when you call an unimplemented method of Xpetra
/** Mainly use for development in progress. **/
class NotImplemented : public Teuchos::ExceptionBase {
 public:
  NotImplemented(const std::string& what_arg)
    : Teuchos::ExceptionBase(what_arg) {}
};

//! Exception throws to report errors in the internal logical of the program.
class RuntimeError : public Teuchos::ExceptionBase {
 public:
  RuntimeError(const std::string& what_arg)
    : Teuchos::ExceptionBase(what_arg) {}
};

//! Exception throws to report incompatible objects (like maps).
class Incompatible : public Teuchos::ExceptionBase {
 public:
  Incompatible(const std::string& what_arg)
    : Teuchos::ExceptionBase(what_arg) {}
};

}  // namespace Exceptions
}  // namespace Xpetra

#endif
