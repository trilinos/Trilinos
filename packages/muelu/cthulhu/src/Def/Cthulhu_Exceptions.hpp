#ifndef CTHULHU_EXCEPTIONS_HPP
#define CTHULHU_EXCEPTIONS_HPP

#include <exception>
#include "Teuchos_Exceptions.hpp"

namespace Cthulhu {
  namespace Exceptions {
    
    //! Exception indicating invalid cast attempted
    /** For instance, this exception is throw when you try to mix Epetra and Tpetra objects: a Cthulhu::EpetraCrsMatrix cannot be built from a Cthulhu::TpetraMap. **/
    class BadCast : public Teuchos::ExceptionBase
    {
    public:
      BadCast(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

    //! Exception throws when you call an unimplemented method of Cthulhu
    /** Mainly use for developement work. **/
    class NotImplemented : public Teuchos::ExceptionBase
    {
    public:
      NotImplemented(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };
    
  }
}

#endif
