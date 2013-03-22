#ifndef _fei_ErrMacros_hpp_
#define _fei_ErrMacros_hpp_
/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_iostream.hpp"

//
//This file simply holds the macros used to check error returns
//and print appropriate output using cerr.
//

#ifdef CHK_ERR
#undef CHK_ERR
#endif

#ifndef fei_file
#define fei_file "unknown_fei_file"
#endif

#define CHK_ERR(a) { int fei_ErrorCode = a; \
                    if (fei_ErrorCode != 0) { \
                      fei::console_out() << " FEI ERROR, " << fei_file << ", line " \
                           << __LINE__ << " " << fei_ErrorCode << FEI_ENDL; \
                      return(fei_ErrorCode); \
                   } }

#ifdef ERReturn
#undef ERReturn
#endif

#define ERReturn(a) { fei::console_out() << " FEI ERROR, " << fei_file << ", line " \
                           << __LINE__ << FEI_ENDL; \
			   return(-1); }

#ifdef voidERReturn
#undef voidERReturn
#endif

#define voidERReturn { fei::console_out() << " FEI ERROR, " << fei_file \
			   << ", line " << __LINE__ << FEI_ENDL; \
			   return; }

#endif

