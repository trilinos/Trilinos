// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_string_case_compare_hpp
#define stk_util_string_case_compare_hpp

#include <strings.h>
#include <string>
#include <functional>

namespace stk {

/** \addtogroup util_module
 *  \{
 */

//----------------------------------------------------------------------

/** \brief  Case-insensitive equality compare */
inline
bool equal_case( const char * lhs , const char * rhs )
{ return strcasecmp( lhs , rhs ) == 0 ; }

/** \brief  Case-insensitive inequality compare */
inline
bool not_equal_case( const char * lhs , const char * rhs )
{ return strcasecmp( lhs , rhs ) != 0 ; }

/** \brief  Case-insensitive less-than compare */
inline
bool less_case( const char * lhs , const char * rhs )
{ return strcasecmp( lhs , rhs ) < 0 ; }

/** \brief  Case-insensitive less-than-or-equal-to compare */
inline
bool less_equal_case( const char * lhs , const char * rhs )
{ return strcasecmp( lhs , rhs ) <= 0 ; }

/** \brief  Case-insensitive greater-than compare */
inline
bool greater_case( const char * lhs , const char * rhs )
{ return strcasecmp( lhs , rhs ) > 0 ; }

/** \brief  Case-insensitive greater-than-or-equal-to compare */
inline
bool greater_equal_case( const char * lhs , const char * rhs )
{ return strcasecmp( lhs , rhs ) >= 0 ; }

//----------------------------------------------------------------------

/** \brief  Case-insensitive equality compare */
inline
bool equal_case( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) == 0 ; }

/** \brief  Case-insensitive inequality compare */
inline
bool not_equal_case( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) != 0 ; }

/** \brief  Case-insensitive less-than compare */
inline
bool less_case( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) < 0 ; }

/** \brief  Case-insensitive less-than-or-equal-to compare */
inline
bool less_equal_case( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) <= 0 ; }

/** \brief  Case-insensitive greater-than compare */
inline
bool greater_case( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) > 0 ; }

/** \brief  Case-insensitive greater-than-or-equal-to compare */
inline
bool greater_equal_case( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) >= 0 ; }

//----------------------------------------------------------------------

/** \brief  Case-insensitive equality compare binary function object. */
struct EqualCase {
  /** \brief  Case-insensitive equality compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return equal_case( lhs , rhs ); }
};

/** \brief  Case-insensitive inequality compare binary function object. */
struct NotEqualCase {
  /** \brief  Case-insensitive inequality compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return not_equal_case( lhs , rhs ); }
};

/** \brief  Case-insensitive less-than compare binary function object. */
struct LessCase {
  /** \brief  Case-insensitive less-than compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return less_case( lhs , rhs ); }
};

/** \brief  Case-insensitive less-than-or-equal-to compare binary function object. */
struct LessEqualCase {
  /** \brief  Case-insensitive less-than-or-equal-to compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return less_equal_case( lhs , rhs ); }
};

/** \brief  Case-insensitive greater-than compare binary function object. */
struct GreaterCase {
  /** \brief  Case-insensitive greater-than compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return greater_case( lhs , rhs ); }
};

/** \brief  Case-insensitive greater-than-or-equal-to compare binary function object. */
struct GreaterEqualCase {
  /** \brief  Case-insensitive greater-than-or-equal-to compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return greater_equal_case( lhs , rhs ); }
};

//----------------------------------------------------------------------

/** \} */

} // namespace stk

#endif /* stk_util_string_case_compare_hpp */

