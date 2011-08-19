#ifndef MUELU_EXCEPTIONS_HPP
#define MUELU_EXCEPTIONS_HPP

#include <Teuchos_Exceptions.hpp>

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {
  namespace Exceptions {
    
    //! Exception indicating invalid cast attempted
    class BadCast : public Teuchos::ExceptionBase
    {
    public:
      BadCast(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

    //! Exception throws when you call an unimplemented method of MueLu
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

    //! Exception throws to report overflows.
    class Overflow : public Teuchos::ExceptionBase
    {
    public:
      Overflow(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

    //! Exception throws to report incompatible objects (like maps).
    class Incompatible : public Teuchos::ExceptionBase
    {
    public:
      Incompatible(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

  }
}

#endif //ifndef MUELU_EXCEPTIONS_HPP
