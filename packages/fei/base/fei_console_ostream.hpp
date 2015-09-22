/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_console_ostream_hpp_
#define _fei_console_ostream_hpp_

#include "fei_macros.hpp"
#include "fei_iostream.hpp"

namespace fei {

/** Set the output stream that fei writes output to.
  fei almost never writes console output, except when errors occur.
  This ostream is the ostream that is returned by the 'console_out()'
  function below.
*/
void set_console_ostream(std::ostream& os);


/** Obtain an output-stream to write 'screen' output to.
  By default this output-stream is std::cerr, but can be set
  to any std::ostream using the above 'set_console_ostream' function.
*/
std::ostream& console_out();

}//namespace fei

#endif
//endif for ifndef _fei_console_ostream_hpp_


