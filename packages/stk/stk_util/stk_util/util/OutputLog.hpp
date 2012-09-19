/*--------------------------------------------------------------------*/
/*    Copyright 2003 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_UTIL_OutputLog_h
#define STK_UTIL_UTIL_OutputLog_h

#include <iosfwd>

namespace sierra {

std::ostream &out();                ///< Normal output stream
std::ostream &dout();               ///< Diagnostic output stream
std::ostream &pout();               ///< Per-processor output stream (See RuntimeDeferredx)
std::ostream &tout();               ///< Regression test textual output stream

std::ostream &dwout();              ///< Diagnostic writer stream

} // namespace sierra

#endif // STK_UTIL_UTIL_OutputLog_h


