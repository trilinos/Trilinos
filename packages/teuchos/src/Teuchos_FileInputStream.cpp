#include "Teuchos_FileInputStream.hpp"
#include "Teuchos_Error.hpp"

using namespace Teuchos;

FileInputStream::FileInputStream(const string& filename)
	: XMLInputStream(), file_(fopen(filename.c_str(), "rb"))
{
	if (file_ == NULL)
		{
			Error::raise("FileInputStream ctor failed to open file: " 
									 + filename);
		}
}

unsigned int FileInputStream::readBytes(unsigned char* const toFill, 
																				const unsigned int maxToRead)
{
	if (feof(file_)) return 0;
	int n = ::fread((void*) toFill, sizeof(char), maxToRead, file_);
	if (n <= 0 || (n<(int) maxToRead && !feof(file_)))
		Error::raise("FileInputStream::readBytes error");
	
	return (size_t) n;
}

