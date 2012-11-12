// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_EXCEPTIONS_HPP
#define XPETRA_EXCEPTIONS_HPP

#include <exception>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_Assert.hpp>

#include "Xpetra_ConfigDefs.hpp"

// Dynamically cast 'obj' to 'type newObj'. newObj is declared inside of the macro.
// If the dynamic cast failed, throw an exception of type Xpetra::Exception::Bad_Cast (using the message exceptionMsg).
#define XPETRA_DYNAMIC_CAST(type, obj, newObj, exceptionMsg)            \
  type * newObj ## _pt = dynamic_cast<type *>(&obj);                    \
  TEUCHOS_TEST_FOR_EXCEPTION(newObj ## _pt == NULL, Xpetra::Exceptions::BadCast, "Cannot cast '" #obj "' to a " #type ". " #exceptionMsg); \
  type & newObj = *newObj ## _pt;

// Dynamically cast the RCP 'obj' to 'RCP<type> newObj'. newObj is declared inside of the macro.
// If the dynamic cast failed, throw an exception of type Xpetra::Exception::Bad_Cast (using the message exceptionMsg).
#define XPETRA_RCP_DYNAMIC_CAST(type, obj, newObj, exceptionMsg)        \
  const RCP<type > & newObj = Teuchos::rcp_dynamic_cast<type >(obj);    \
  TEUCHOS_TEST_FOR_EXCEPTION(newObj == Teuchos::null, Xpetra::Exceptions::BadCast, "Cannot cast '" #obj "' to a " #type ". " #exceptionMsg);

#ifdef HAVE_XPETRA_EPETRA
#define XPETRA_FACTORY_ERROR_IF_EPETRA(lib)                            \
  if ((lib) == ::Xpetra::UseEpetra)                                     \
    TEUCHOS_TEST_FOR_EXCEPTION(1, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");
#else
#define XPETRA_FACTORY_ERROR_IF_EPETRA(lib)
#endif

#define XPETRA_FACTORY_END TEUCHOS_TEST_FOR_EXCEPTION(1, ::Xpetra::Exceptions::BadCast, "Unknown map->lib() type. Did you compile with Epetra and Tpetra support?");

namespace Xpetra {
  namespace Exceptions {

    //! Exception indicating invalid cast attempted
    /** For instance, this exception is throw when you try to mix Epetra and Tpetra objects: for example, a Xpetra::EpetraCrsMatrix cannot be built from a Xpetra::TpetraMap. **/
    class BadCast : public Teuchos::ExceptionBase
    {
    public:
      BadCast(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

    //! Exception throws when you call an unimplemented method of Xpetra
    /** Mainly use for development in progress. **/
    class NotImplemented : public Teuchos::ExceptionBase
    {
    public:
      NotImplemented(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

    //! Exception throws to report errors in the internal logical of the program.
    class RuntimeError : public Teuchos::ExceptionBase
    {
    public:
      RuntimeError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

  }
}

#endif
