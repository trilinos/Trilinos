#include "Teuchos_DefaultRaiseHandler.hpp"

using namespace Teuchos;

void DefaultRaiseHandler::handleRaise(const char* msg)
{
	throw std::runtime_error(msg);
}


