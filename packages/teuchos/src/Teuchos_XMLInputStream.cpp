#include "Teuchos_XMLInputStream.hpp"
#include "Teuchos_Error.hpp"

using namespace Teuchos;


unsigned int XMLInputStream::curPos() const 
{
	Error::raise("XMLInputStream::curPos() should never be called. It is"
							 "there only for compatibility with Xerces");
	return 0;
}
