#include "Teuchos_FileInputStream.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;

FileInputStream::FileInputStream(const string& filename)
	: XMLInputStream(), file_(fopen(filename.c_str(), "rb"))
{
  TEST_FOR_EXCEPTION(file_ == NULL,
                     runtime_error,
                     "FileInputStream ctor failed to open file: " 
                     << filename);
}

unsigned int FileInputStream::readBytes(unsigned char* const toFill, 
																				const unsigned int maxToRead)
{
	if (feof(file_)) return 0;
	int n = ::fread((void*) toFill, sizeof(char), maxToRead, file_);
  
	TEST_FOR_EXCEPTION(n <= 0 || (n<(int) maxToRead && !feof(file_)),
                     runtime_error,
                     "FileInputStream::readBytes error");
	
	return (size_t) n;
}

