#include "Teko_RequestHandler.hpp"

namespace Teko {

RequestHandler::RequestHandler() 
{ }

void RequestHandler::addRequestCallback(const Teuchos::RCP<RequestCallbackBase> & callback)
{
   callbacks_.push_back(callback);
}

} // end namespace Teko
