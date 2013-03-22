/**   ------------------------------------------------------------
 *    Copyright 2000-2008 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <stk_util/util/Null_Streambuf.hpp>

/*--------------------------------------------------------------------*/

null_streambuf::null_streambuf() : std::streambuf()
{
  setp( buf , buf + sizeof(buf) );
}

null_streambuf::~null_streambuf() {}

/*--------------------------------------------------------------------*/
/* Overflow */

int null_streambuf::overflow( int c )
{
  setp( buf , buf + sizeof(buf) );

  return c ;
}

/*--------------------------------------------------------------------*/

int null_streambuf::sync()
{
  return 0 ;
}

std::streambuf * null_streambuf::setbuf( char * s , std::streamsize n )
{
  return this ;
}




/*--------------------------------------------------------------------*/
