#ifndef TEUCHOS_ERROR_H
#define TEUCHOS_ERROR_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RaiseHandlerBase.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos
{
  using std::exception;
  using std::string;
  
  /**
   * Class Error provides customizable error handling for Trilinos. 
   * All errors should be reported by a call to Error::raise(). 
   */

  class Error
  {
  public:

    /** raise() is to be called when an error is detected. */
    static void raise(const string& msg);

    /** report a violation of a range [low, high) */
    static void boundsError(const string& msg, 
                            int badValue, int low, int high);

    /** install a custom raise handler */
    static void setRaiseHandler(RaiseHandlerBase* h) {handler() = h->getRcp();}

    /** pass along a caught exception with an additional message */
    static void trace(const exception& e, const string& msg);
    

  private:
    static RefCountPtr<RaiseHandlerBase>& handler();

    /** Returns the default raise handler. Default action is to
     * throw an exception. */
    static RefCountPtr<RaiseHandlerBase>& defaultHandler() ;
  };
}

#endif
