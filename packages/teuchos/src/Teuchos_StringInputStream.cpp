#include "Teuchos_StringInputStream.hpp"

using namespace Teuchos;


unsigned int StringInputStream::readBytes(unsigned char* const toFill, 
																					const unsigned int maxToRead)
{
	if (pos_ == text_.length()) return 0;
	
	int toRead = text_.length() - pos_;
	if ((int) maxToRead < toRead) toRead = maxToRead;

	strncpy((char*) toFill, text_.c_str(), toRead);

	pos_ += toRead;
	
	return (size_t) toRead;
}

