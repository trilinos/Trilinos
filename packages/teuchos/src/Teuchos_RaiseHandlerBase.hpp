#ifndef TEUCHOS_RAISEHANDLERBASE_H
#define TEUCHOS_RAISEHANDLERBASE_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos
{
  using std::string;

  /** 
   * Teuchos::RaiseHandlerBase is a base class for objects that specify what is to
   * happen upon a call to Teuchos::Error::raise(). When Teuchos::Error::raise() is called,
   * the handleRaise() method of a raise handler subclass will be called.
   * The behavior can be customized by supplying a custom raise hander.
   * By default, the Teuchos::DefaultRaiseHandler is used and the handleRaise()
   * method throws an exception.
   */
  class RaiseHandlerBase
    {
    public:
      /** empty ctor */
      RaiseHandlerBase(){;}

      /** virtual dtor */
      virtual ~RaiseHandlerBase(){;}

      /** the handleRaise() method is called inside Error::raise() */
      virtual void handleRaise(const char* msg) = 0 ;

      /** Return a RefCountPtr containing self. This is easy to implement: {return rcp(this);} */
      virtual RefCountPtr<RaiseHandlerBase> getRcp() = 0 ;

    private:

    };
}
#endif
