#ifndef TEUCHOS_DEFAULTRAISEHANDLER_H
#define TEUCHOS_DEFAULTRAISEHANDLER_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RaiseHandlerBase.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include <stdexcept>

namespace Teuchos
{
  using std::string;

  /** 
   * The default raise handler throws a std::runtime_error() exception.
   */
  class DefaultRaiseHandler : public RaiseHandlerBase
    {
    public:
      /** empty ctor */
      DefaultRaiseHandler(){;}

      /** TUVD */
      virtual ~DefaultRaiseHandler(){;}

      /** When Error::raise() is called, throw a std::runtime_error
       * exception. */
      virtual void handleRaise(const char* msg);

      /** Return a RefCountPtr containing self. */
      virtual RefCountPtr<RaiseHandlerBase> getRcp() {return rcp(this, true);}

    private:
    };
}

#endif
