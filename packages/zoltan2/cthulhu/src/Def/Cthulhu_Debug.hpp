#ifndef CTHULHU_DEBUG_HPP
#define CTHULHU_DEBUG_HPP

// #include "Cthulhu_ConfigDefs.hpp" -> TODO: add this header line. But at the moment, Cthulhu_ConfigDefs use Cthulhu_Debug

// #ifdef HAVE_MUEMAT_DEBUG

// TODO: Description
#include <string>

void cthulhu_debug_me(const std::string & file, const std::string & funcName);
void cthulhu_debug_me_print();

//#define DEBUG_ME
#ifdef DEBUG_ME
#define CTHULHU_DEBUG_ME       { cthulhu_debug_me(__FILE__,__FUNCTION__); };
#define CTHULHU_DEBUG_ME_PRINT { cthulhu_debug_me_print();                };
#else
#define CTHULHU_DEBUG_ME      
#define CTHULHU_DEBUG_ME_PRINT 
#endif

// #else

// #error TEST

// #endif

//#define CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
// #define HAVE_CTHULHU_EPETRA
// #define HAVE_CTHULHU_TPETRA

//#define CTHULHU_NOT_IMPLEMENTED



#endif
