/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"
#include "fei_iostream.hpp"

namespace fei {

std::ostream* console_ostream_ptr(std::ostream* osptr=NULL)
{
  static std::ostream* fei_ostream_ptr = &std::cerr;
  if (osptr) fei_ostream_ptr = osptr;
  return fei_ostream_ptr;
}

void set_console_ostream(std::ostream& os)
{
  console_ostream_ptr(&os);
}

std::ostream& console_out()
{
  return *console_ostream_ptr();
}

}//namespace fei

