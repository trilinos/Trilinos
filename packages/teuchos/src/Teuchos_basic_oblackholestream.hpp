// //////////////////////////////////////////////////////////////////////
// Teuchos_basic_oblackholeostream.h

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
	__ostream_type& operator<<(__ostream_type& (*__pf)(__ostream_type&)) { *this; }
	__ostream_type& operator<<(__ios_type& (*__pf)(__ios_type&)) { *this; }
	__ostream_type& operator<<(std::ios_base& (*__pf) (std::ios_base&)) { *this; }

	// 27.6.2.5.2 Arithmetic Inserters
	__ostream_type& operator<<(long __n) { *this; }
	__ostream_type& operator<<(unsigned long __n) { *this; }
    __ostream_type& operator<<(bool __n) { *this; }
	__ostream_type& operator<<(short __n) { *this; }
	__ostream_type& operator<<(unsigned short __n) { *this; }
	__ostream_type& operator<<(int __n) { *this; }
	__ostream_type& operator<<(unsigned int __n) { *this; }
	__ostream_type& operator<<(double __f) { *this; }
	__ostream_type& operator<<(float __f) { *this; }
	__ostream_type& operator<<(long double __f) { *this; }
	__ostream_type& operator<<(const void* __p) { *this; }
	__ostream_type& operator<<(__streambuf_type* __sb) { *this; }

	// Unformatted output:
	__ostream_type& put(char_type __c) { *this; }
	__ostream_type& write(const char_type* __s, std::streamsize __n) { *this; }
	__ostream_type& flush() { *this; }

	// Seeks:
	pos_type tellp() { return 0; }
	__ostream_type& seekp(pos_type) { *this; }
	__ostream_type& seekp(off_type, std::ios_base::seekdir) { *this; }

}; // end class basic_oblackholestream

} // end namespace Teuchos

#endif // TEUCHOS_BASIC_O_BLACK_HOLE_STREAM_H
