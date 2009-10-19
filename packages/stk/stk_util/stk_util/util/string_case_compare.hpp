
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
struct EqualCase : public std::binary_function<std::string,std::string,bool> {
  /** \brief  Case-insensitive equality compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return equal_case( lhs , rhs ); }
};

/** \brief  Case-insensitive inequality compare binary function object. */
struct NotEqualCase : public std::binary_function<std::string,std::string,bool> {
  /** \brief  Case-insensitive inequality compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return not_equal_case( lhs , rhs ); }
};

/** \brief  Case-insensitive less-than compare binary function object. */
struct LessCase : public std::binary_function<std::string,std::string,bool> {
  /** \brief  Case-insensitive less-than compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return less_case( lhs , rhs ); }
};

/** \brief  Case-insensitive less-than-or-equal-to compare binary function object. */
struct LessEqualCase : public std::binary_function<std::string,std::string,bool> {
  /** \brief  Case-insensitive less-than-or-equal-to compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return less_equal_case( lhs , rhs ); }
};

/** \brief  Case-insensitive greater-than compare binary function object. */
struct GreaterCase : public std::binary_function<std::string,std::string,bool> {
  /** \brief  Case-insensitive greater-than compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return greater_case( lhs , rhs ); }
};

/** \brief  Case-insensitive greater-than-or-equal-to compare binary function object. */
struct GreaterEqualCase : public std::binary_function<std::string,std::string,bool> {
  /** \brief  Case-insensitive greater-than-or-equal-to compare binary function object. */
  bool operator()( const std::string & lhs , const std::string & rhs ) const
    { return greater_equal_case( lhs , rhs ); }
};

//----------------------------------------------------------------------

/** \} */

} // namespace stk

#endif /* stk_util_string_case_compare_hpp */

