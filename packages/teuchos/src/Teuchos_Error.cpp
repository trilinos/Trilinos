#include "Teuchos_Error.hpp"
#include "Teuchos_DefaultRaiseHandler.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_Out.hpp"

using namespace Teuchos;
using std::exception;
using std::string;

/* initialize the static raise handler object to the default handler. This
 * can be changed later with a call to setRaiseHandler() */


RefCountPtr<RaiseHandlerBase>& Error::defaultHandler() 
{
  static RefCountPtr<RaiseHandlerBase> h = rcp(new DefaultRaiseHandler());
  return h;
}

RefCountPtr<RaiseHandlerBase>& Error::handler() 
{
  static RefCountPtr<RaiseHandlerBase> h = rcp(new DefaultRaiseHandler());
  return h;
}

void Error::raise(const std::string& msg)
{
	handler()->handleRaise(msg.c_str());
}

void Error::boundsError(const std::string& msg, int badValue, int low, int high)
{
  std::string s = msg + ": bounds error: index=" + Teuchos::toString(badValue)
    + ", range=[" + Teuchos::toString(low) + ", " + Teuchos::toString(high)
    + ")";
  raise(s);
}

void Error::trace(const std::exception& e, const std::string& where)
{
	std::string msg = std::string(e.what()) + " at " + where;

	Out::println(msg);

	throw runtime_error(msg);
}
