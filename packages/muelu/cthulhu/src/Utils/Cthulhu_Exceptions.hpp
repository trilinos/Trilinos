#ifndef CTHULHU_EXCEPTIONS_HPP
#define CTHULHU_EXCEPTIONS_HPP

#include <exception>
#include <Teuchos_Exceptions.hpp>

#include "Cthulhu_ConfigDefs.hpp"

// Dynamically cast 'obj' to 'type newObj'. newObj is declared inside of the macro.
// If the dynamic cast failed, throw an exception of type Cthulhu::Exception::Bad_Cast (using the message exceptionMsg).
#define CTHULHU_DYNAMIC_CAST(type, obj, newObj, exceptionMsg)           \
  type * newObj ## _pt = dynamic_cast<type *>(&obj);                    \
  TEST_FOR_EXCEPTION(newObj ## _pt == NULL, Cthulhu::Exceptions::BadCast, "Cannot cast '" #obj "' to a " #type ". " #exceptionMsg); \
  type & newObj = *newObj ## _pt;

// Dynamically cast the RCP 'obj' to 'RCP<type> newObj'. newObj is declared inside of the macro.
// If the dynamic cast failed, throw an exception of type Cthulhu::Exception::Bad_Cast (using the message exceptionMsg).
#define CTHULHU_RCP_DYNAMIC_CAST(type, obj, newObj, exceptionMsg)       \
  const RCP<type > & newObj = Teuchos::rcp_dynamic_cast<type >(obj);    \
  TEST_FOR_EXCEPTION(newObj == Teuchos::null, Cthulhu::Exceptions::BadCast, "Cannot cast '" #obj "' to a " #type ". " #exceptionMsg); 

#ifdef HAVE_CTHULHU_EPETRA
#define CTHULHU_FACTORY_ERROR_IF_EPETRA(lib)                            \
  if ((lib) == UseEpetra)                                               \
        TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");
#else
#define CTHULHU_FACTORY_ERROR_IF_EPETRA(lib)
#endif

#define CTHULHU_FACTORY_END TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::BadCast, "Unknown map->lib() type. Did you compile with Epetra and Tpetra support?");

namespace Cthulhu {
  namespace Exceptions {
    
    //! Exception indicating invalid cast attempted
    /** For instance, this exception is throw when you try to mix Epetra and Tpetra objects: for example, a Cthulhu::EpetraCrsMatrix cannot be built from a Cthulhu::TpetraMap. **/
    class BadCast : public Teuchos::ExceptionBase
    {
    public:
      BadCast(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
    };

    //! Exception throws when you call an unimplemented method of Cthulhu
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
