#include "Teuchos_Out.hpp"
#include "Teuchos_DefaultWriter.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_MPISession.hpp"

using namespace Teuchos;



/* initialize the static raise handler object to the default handler. This
 * can be changed later with a call to setRaiseHandler() */

RefCountPtr<WriterBase> Out::writer_ = rcp(new DefaultWriter(), true);



void Out::print(const std::string& msg)
{
	writer_->print(msg);
}

void Out::println(const std::string& msg)
{
	writer_->println(msg);
}

void Out::rootPrintln(const std::string& msg)
{
	int rank=0;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (rank==0)
		{
			writer_->println(msg);
		}
}

void Out::rootPrintf(const char* format, ...)
{
	int rank=0;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (rank==0)
		{
			va_list args;
			va_start(args, format);
			Out::vprintf(format, args);
		}
}

void Out::printf(const char* format, ...)
{
	va_list args;
	va_start(args, format);

	int bufSize = 100;

	while (bufSize > 0)
		{
			char* str = new char[bufSize+1];
			int rtn = hack_vsnprintf(str, bufSize, format, args);
			if (rtn > 0)
				{
					writer_->print(str);
					va_end(args);
					delete [] str;
					return;
				} 
			else if (rtn==0)
				{
					va_end(args);
					delete [] str;
					return;
				}
			else
				{
					bufSize *= 2;
				}
		}

  TEST_FOR_EXCEPTION(true, std::length_error, 
                     "buffer overflow in Out::printf()");
}

void Out::vprintf(const char* format, va_list args)
{
	int bufSize = 100;

	while (bufSize > 0)
		{
			char* str = new char[bufSize+1];
			int rtn = hack_vsnprintf(str, bufSize, format, args);
			if (rtn > 0)
				{
					writer_->print(str);
					va_end(args);
					delete [] str;
					return;
				} 
			else if (rtn==0)
				{
					va_end(args);
					delete [] str;
					return;
				}
			else
				{
					bufSize *= 2;
				}
		}

	TEST_FOR_EXCEPTION(true, std::length_error, 
                     "buffer overflow in Out::vprintf()");
}

int Out::hack_vsnprintf(char* str, size_t size, const char* format, va_list args)
{
#ifndef TFLOP
	return vsnprintf(str, size, format, args);
#else 
	static FILE* devnull = fopen("/dev/null", "w");

	if (devnull==0)
		{
			if (size < 0)
				{
					str[0] = '\0';
				}
			return 0;
		}
	else
		{
			/* count the characters */
			int n = vfprintf(devnull, format, args);
			if (n <= size-1)
				{
					vsprintf(str, format, args);
					return n;
				}
			else
				{
					return -1;
				}
		}
#endif
}

void Out::setWriter(const RefCountPtr<WriterBase>& writer )
{
	writer_ = writer;
}
