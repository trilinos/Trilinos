#ifndef CTHULHU_DEBUG_HPP
#define CTHULHU_DEBUG_HPP

// #ifdef HAVE_MUEMAT_DEBUG

// TODO: Description
#include <string>

void cthulhu_debug_me(const std::string & file, const std::string & funcName);
void cthulhu_debug_me_print();

#define CTHULHU_DEBUG_ME       { cthulhu_debug_me(__FILE__,__FUNCTION__); };
#define CTHULHU_DEBUG_ME_PRINT { cthulhu_debug_me_print();                };

// #else

// #error TEST

// #endif

#endif
