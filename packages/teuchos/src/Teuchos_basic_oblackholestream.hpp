// //////////////////////////////////////////////////////////////////////
// Teuchos_basic_oblackholeostream.h

// @HEADER 
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_BASIC_O_BLACK_HOLE_STREAM_H
#define TEUCHOS_BASIC_O_BLACK_HOLE_STREAM_H

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

///
/** <tt>basic_ostream<></tt> subclass that does nothing but discard output.
 *
 * Use the class anytime you must pass an <tt>basic_ostream<></tt> object
 * but don't want the output for any reason.
 */
//@{
template<typename _CharT, typename _Traits>
class basic_oblackholestream
	: virtual public std::basic_ostream<_CharT, _Traits>
{
public:

	// Types (inherited from basic_ios (27.4.4)):
	typedef _CharT                     		char_type;
	typedef typename _Traits::int_type 		int_type;
	typedef typename _Traits::pos_type 		pos_type;
	typedef typename _Traits::off_type 		off_type;
	typedef _Traits                    		traits_type;
	
	// Non-standard Types:
	typedef std::basic_streambuf<_CharT, _Traits>      __streambuf_type;
	typedef std::basic_ios<_CharT, _Traits>            __ios_type;
	typedef std::basic_ostream<_CharT, _Traits>        __ostream_type;

	// 27.6.2.2 Constructor/destructor:
	explicit basic_oblackholestream() : __ostream_type(NULL) {}
//	explicit basic_oblackholestream(__streambuf_type* __sb) : __ostream_type(NULL){}
	virtual ~basic_oblackholestream() { }

	// 27.6.2.5 Formatted output:
	// 27.6.2.5.3  basic_ostream::operator<<
	__ostream_type& operator<<(__ostream_type& (*__pf)(__ostream_type&)) { return *this; }
	__ostream_type& operator<<(__ios_type& (*__pf)(__ios_type&)) { return *this; }
	__ostream_type& operator<<(std::ios_base& (*__pf) (std::ios_base&)) { return *this; }

	// 27.6.2.5.2 Arithmetic Inserters
	__ostream_type& operator<<(long __n) { return *this; }
	__ostream_type& operator<<(unsigned long __n) { return *this; }
    __ostream_type& operator<<(bool __n) { return *this; }
	__ostream_type& operator<<(short __n) { return *this; }
	__ostream_type& operator<<(unsigned short __n) { return *this; }
	__ostream_type& operator<<(int __n) { return *this; }
	__ostream_type& operator<<(unsigned int __n) { return *this; }
	__ostream_type& operator<<(double __f) { return *this; }
	__ostream_type& operator<<(float __f) { return *this; }
	__ostream_type& operator<<(long double __f) { return *this; }
	__ostream_type& operator<<(const void* __p) { return *this; }
	__ostream_type& operator<<(__streambuf_type* __sb) { return *this; }

	// Unformatted output:
	__ostream_type& put(char_type __c) { return *this; }
	__ostream_type& write(const char_type* __s, std::streamsize __n) { return *this; }
	__ostream_type& flush() { return *this; }

	// Seeks:
	pos_type tellp() { return 0; }
	__ostream_type& seekp(pos_type) { return *this; }
	__ostream_type& seekp(off_type, std::ios_base::seekdir) { return *this; }

}; // end class basic_oblackholestream

} // end namespace Teuchos

#endif // TEUCHOS_BASIC_O_BLACK_HOLE_STREAM_H
