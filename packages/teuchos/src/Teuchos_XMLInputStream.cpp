#include "Teuchos_XMLInputStream.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;


unsigned int XMLInputStream::curPos() const 
{
	TEST_FOR_EXCEPTION(true, std::logic_error,
                     "XMLInputStream::curPos() should never be called. It is"
                     "there only for compatibility with Xerces");
	return 0;
}
