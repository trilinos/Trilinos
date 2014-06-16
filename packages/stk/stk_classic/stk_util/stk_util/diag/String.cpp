/**   ------------------------------------------------------------
 *    Copyright 2002-2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <stk_util/diag/String.hpp>

#include <cstddef>
#include <cctype>
#include <stdexcept>
#include <iostream>

//----------------------------------------------------------------------

namespace sierra {

namespace {

static const char arraylower_t[] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F,
                                    0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F,
                                    0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2A, 0x2B, 0x2C, 0x2D, 0x2E, 0x2F,
                                    0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x3B, 0x3C, 0x3D, 0x3E, 0x3F, 
                                    0x40,  'a',  'b',  'c',  'd',  'e',  'f',  'g',  'h',  'i',  'j',  'k',  'l',  'm',  'n',  'o',
                                     'p',  'q',  'r',  's',  't',  'u',  'v',  'w',  'x',  'y',  'z', 0x5B, 0x5C, 0x5D, 0x5E, 0x5F,
                                    0x60,  'a',  'b',  'c',  'd',  'e',  'f',  'g',  'h',  'i',  'j',  'k',  'l',  'm',  'n',  'o',
                                     'p',  'q',  'r',  's',  't',  'u',  'v',  'w',  'x',  'y',  'z', 0x7B, 0x7C, 0x7D, 0x7E, 0x7F
};

inline int arraylower(int c) {
  return arraylower_t[c];
}

} // namespace <unnamed>

int to_label( int c )
{
  return isspace(c) ? '_' : c;
}

size_t char_simple_traits::length( const char * c1 )
{
  const char * c = c1 ;
  if ( c ) while ( *c ) ++c ;
  return c - c1 ;
}

int char_simple_traits::compare( const char * c1 , const char * c2 )
{
  for ( ; *c1 && arraylower(*c1) == arraylower(*c2); c1++, c2++)
    ;

  return arraylower(*c1) - arraylower(*c2);
}

size_t char_label_traits::length( const char * c1 )
{
  const char * c = c1 ;
  if ( c ) while ( *c ) ++c ;
  return c - c1 ;
}

void char_label_traits::convert( char * c , size_t n )
{
  for ( char * const e = c + n ; c != e ; ++c )
    *c = to_label(*c);
}

int char_label_traits::compare( const char * c1 , const char * c2 )
{
  for ( ; *c1 && arraylower(to_label(*c1)) == arraylower(to_label(*c2)); c1++, c2++)
    ;

  return arraylower(to_label(*c1)) - arraylower(to_label(*c2));
}

} // namespace sierra

//----------------------------------------------------------------------

namespace sierra {

namespace implementation {

/** @union StringData
 * @par Objectives:
 * @li  Don't allocate for short strings, strings <= max_len
 * @li  buf_len == sizeof(String)
 * @li  buf_len % sizeof(unsigned long) == 0
 *
 * @par Limitations:
 * @li  sizeof(StringData::Large) + 2 <= buf_len < ( 127 == 0x7f )
 *
 * @par Memory layout for short strings that are <= max_len
 * @li  buf[ 0 .. max_len - 1 ] = buffer for characters
 * @li  buf[ max_len ]          = null
 * @li  buf[ off_len ]          = length of string
 * @li  data                    = must not be used
 *
 * @par Memory layout for long strings that are > max_len
 * @li  data.ptr                  = pointer to allocated memory
 * @li  data.len                  = length of string
 * @li  data.siz                  = allocated size
 * @li  buf[ max_len ]           != null
 * @li  buf[ 0 .. max_len - 1 ]   = must not be used
 */

// Manage memory but do not invalidate current memory
// until operation is complete.

char * StringData::mem( const char * cs , size_t n )
{
  enum { inc = 4 * sizeof(long) };
  static std::string::allocator_type a ;

  const bool is_allocated = small[ max_len ] == 1;
  const bool to_allocated = max_len < n ;

  size_t new_alloc = 0 ;

  if ( to_allocated && ( ! is_allocated || large.siz <= n ) ) {
    const size_t n_total = n + 1 ;
    new_alloc = n_total % inc ? n_total + inc - n_total % inc : n_total ;
    if ( new_alloc == 0 || new_alloc - 1 < n ) {
      throw std::runtime_error("FAILED MEMORY ALLOCATION SIZING in sierra::String");
    }
  }

  char * dst = NULL ;
  char * del_ptr = NULL ;
  size_t del_size = 0 ;

  if ( is_allocated && ( new_alloc || ! to_allocated ) ) {
    // Deallocate currently allocated memory after copying input,
    // input might be a subset of currently allocated memory.
    del_ptr  = large.ptr ;
    del_size = large.siz ;
  }

  if ( to_allocated ) {
    // Needs to be allocated to hold input

    if ( new_alloc ) {
      // New allocation or reallocation to increase size

      { //----------------------------------------
	// Verify memory layout
	static bool first_pass = true ;

	if ( first_pass ) {
	  first_pass = false ;

	  if ( buf_len % sizeof(long) ||
	       sizeof(StringData) != buf_len ||
	       small + max_len < (char*)(&(large)) + sizeof(Large) ) {
	    throw std::logic_error("StringData memory layout error");
	  }
	}
      } //----------------------------------------

      try {
	large.siz = new_alloc ;
	large.ptr = (char *) a.allocate( new_alloc );
//	std::cout << "String allocated at " << (void *)large.ptr << " for " << new_alloc << std::endl;
      }
      catch (...) {
	throw std::runtime_error("FAILED MEMORY ALLOCATION in sierra::String");
      }
    }

    small[max_len] = 1 ;
    large.len      = n ;
    dst = large.ptr ;
  }
  else {
    small[max_len] = 0 ;
    small[off_len] = n ;
    dst = small ;
  }

  {
    const char * const cs_e = cs + n ;
    char * d = dst ;
    while ( cs != cs_e ) *d++ = *cs++ ;
    *d = 0 ;
  }

  if ( del_ptr != NULL ) {
    try {
//      std::cout << "String deallocated at " << (void *)del_ptr << " for " << del_size << std::endl;
      a.deallocate( del_ptr , del_size );
    }
    catch (...) {
      throw std::runtime_error("FAILED MEMORY DEALLOCATION in sierra::String");
    }
  }

  return dst ;
}

StringData::~StringData()
{ mem(NULL, 0); }

StringData::StringData()
{ small[ max_len ] = small[ off_len ] = small[ 0 ] = 0 ; }

size_t StringData::len() const
{ return small[ max_len ] ? large.len : small[ off_len ] ; }

char * StringData::c_str()
{ return small[ max_len ] ? large.ptr : small ; }

const char * StringData::c_str() const
{ return small[ max_len ] ? large.ptr : small ; }

} // namespace internal

std::ostream &
operator<<( std::ostream & os, const sierra::String & s)
{ return os << s.c_str(); }

std::istream &
operator>>( std::istream & is, sierra::String &	s )
{ std::string tmp; is >> tmp; s.assign(tmp); return is; }

std::ostream &
operator<<( std::ostream & os, const sierra::Identifier &s)
{ return os << s.c_str(); }

std::istream &
operator>>( std::istream & is, sierra::Identifier &s )
{ std::string tmp; is >> tmp; s.assign(tmp); return is; }

} // namespace sierra
