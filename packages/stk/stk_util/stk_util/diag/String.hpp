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

#ifndef STK_UTIL_DIAG_STRING_H
#define STK_UTIL_DIAG_STRING_H

#include <string.h>
#include <iostream>
#include <type_traits>

#ifdef USE_CISTRING

#include <stk_util/util/cistring.hpp>

namespace sierra {

typedef cistring String;
typedef cistring Identifier;
typedef cistring ParamId;

} // namespace sierra

#else 

#include <iosfwd>
#include <string>

namespace sierra {

template <typename T>
using If_String_View = std::enable_if_t<std::is_convertible_v<const T &, std::string_view> &&
                                        !std::is_convertible_v<const T &, const char *> >;

struct char_simple_traits ;

struct char_label_traits ;

/** Define precedence between two types */
template<typename T1, typename T2> struct Precedence ;

template<typename T> struct Precedence<T, T> { typedef T Type ; };

template<>
struct Precedence<char_simple_traits, char_label_traits> {
  typedef char_label_traits Type ;
};

template<>
struct Precedence<char_label_traits, char_simple_traits> {
  typedef char_label_traits Type ;
};

/**
 *  Template base class for don't allocate short strings class.
 */
template<class CharTraits> class StringBase ;

/**
 *  Ordinary characters.
 */
typedef StringBase< char_simple_traits > String ;
typedef StringBase< char_label_traits >  Identifier ;
typedef StringBase< char_label_traits >  ParamId ;

}

//----------------------------------------------------------------------
// Binary Operators:

namespace sierra {

template<class CT1, class CT2>
bool operator== ( const StringBase<CT1> &, const StringBase<CT2> & );

template<class CT1, class CT2>
bool operator!= ( const StringBase<CT1> &, const StringBase<CT2> & );

template<class CT1, class CT2>
bool operator<  ( const StringBase<CT1> &, const StringBase<CT2> & );

template<class CT1, class CT2>
bool operator>  ( const StringBase<CT1> &, const StringBase<CT2> & );

template<class CT1, class CT2>
bool operator<= ( const StringBase<CT1> &, const StringBase<CT2> & );

template<class CT1, class CT2>
bool operator>= ( const StringBase<CT1> &, const StringBase<CT2> & );


template<class CT1>
bool operator== ( const StringBase<CT1> &, std::string_view );

template<class CT1>
bool operator!= ( const StringBase<CT1> &, std::string_view );

template<class CT1>
bool operator<  ( const StringBase<CT1> &, std::string_view );

template<class CT1>
bool operator>  ( const StringBase<CT1> &, std::string_view );

template<class CT1>
bool operator<= ( const StringBase<CT1> &, std::string_view );

template<class CT1>
bool operator>= ( const StringBase<CT1> &, std::string_view );

template<class CT1>
bool operator== (std::string_view, const StringBase<CT1> &);

template<class CT1>
bool operator!= (std::string_view, const StringBase<CT1> &);

template<class CT1>
bool operator<  (std::string_view, const StringBase<CT1> &);

template<class CT1>
bool operator>  (std::string_view, const StringBase<CT1> &);

template<class CT1>
bool operator<= (std::string_view, const StringBase<CT1> &);

template<class CT1>
bool operator>= (std::string_view, const StringBase<CT1> &);



std::ostream &
operator<<( std::ostream & os, const sierra::String &s);

std::istream &
operator>>( std::istream & is, sierra::String &s );

std::ostream &
operator<<( std::ostream & os, const sierra::Identifier &s);

std::istream &
operator>>( std::istream & is, sierra::Identifier &s );

}

namespace sierra {

template<class CT>
class StringBase {
public:
  typedef const char * const_iterator;
  typedef char * iterator;

  typedef CT     traits_type ;
  typedef char   value_type ;
  typedef size_t size_type ;

  StringBase() {}
  ~StringBase() = default;
  StringBase(StringBase &&) = default;
  StringBase(const StringBase &) = default;
  StringBase &operator=(const StringBase &) = default;
  StringBase &operator=(StringBase &&) = default;

  //Allow construction from things implicitly convertable to std::string_view
  template <typename StringViewLike, typename = If_String_View<StringViewLike>> 
  explicit StringBase(const StringViewLike &sv) : data_(std::string_view(sv)) 
  { 
    traits_type::convert(begin(), length()); 
  }

  //Now we allow implicit conversion from std::string to StrinBase anywhere in code
  //This can silently create copy and potential heap allocations
  //This code path should be deprecated!!!
  StringBase(const std::string & s) : data_(s) 
  { 
    traits_type::convert(begin(), length()); 
  }

  StringBase(std::string && cs) : data_(std::move(cs)) 
  {
    traits_type::convert(begin(), length()); 
  }

  StringBase( const_iterator );
  template <class It>
  StringBase( It, It );
  StringBase( const char *, size_type );

  operator std::string_view() const noexcept {return data_;}

  template<class CT2>
  StringBase( const StringBase<CT2> & );

  template <class CT2>
  StringBase(StringBase<CT2> &&);

  StringBase<CT> & operator= ( const char * );

  template <typename StringViewLike, typename = If_String_View<StringViewLike>> 
  StringBase<CT> & operator= ( const StringViewLike & cs)
  { 
    return assign( cs.data(), cs.length() ); 
  }
  StringBase<CT> &operator=(std::string &&);

  template <class CT2>
  StringBase<CT> &operator=(StringBase<CT2> &&);

  StringBase<CT> & operator+= ( const char * );

  template <typename StringViewLike, typename = If_String_View<StringViewLike>> 
  StringBase<CT> & operator+= ( const StringViewLike & cs)
  { 
    return append( cs.data(), cs.length() ); 
  }

  // Capacity:
  constexpr size_type size() const noexcept { return data_.size();}
  constexpr size_type length() const noexcept { return data_.length();}
  constexpr bool empty() const noexcept { return data_.empty();}
  char * data() noexcept { return data_.data();}
  const char * data() const noexcept { return data_.data();}

  const_iterator begin() const;
  iterator begin();
  const_iterator end() const;
  iterator end();

  // Modifiers:
  StringBase<CT> & assign( const char * cs)
  {
    if (cs == nullptr) {
      data_ = std::string();
      return *this;
    }
    return assign(cs, strlen(cs));
  }

  StringBase<CT> & assign( const char * cs, const size_type n)
  {
    //cs must be non-null to prevent undefined behavior
    data_.assign(cs, n);
    traits_type::convert(begin(), length());
    return *this ;

  }

  template <typename StringViewLike, typename = If_String_View<StringViewLike>> 
  StringBase<CT> &assign(const StringViewLike & cs)
  {
    data_ = std::string_view(cs);
    traits_type::convert(begin(), length());
    return *this;
  }
  StringBase<CT> &assign(std::string && cs)
  {
    data_ = std::move(cs);
    traits_type::convert(begin(), length());
    return *this;
  }

  template <class CT2>
  StringBase<CT> &assign(StringBase<CT2> && cs)
  {
    data_ = std::move(cs).s_str();
    traits_type::convert(begin(), length());
    return *this;
  }

  StringBase<CT> & append( const char * );
  StringBase<CT> & append( const char *, const typename StringBase<CT>::size_type );

  template <typename StringViewLike, typename = If_String_View<StringViewLike>> 
  StringBase<CT> & append( const StringViewLike & cs)
  { 
    return append( cs.data(), cs.length() ); 
  }

  void swap( StringBase<CT> & );

  // string operations
  const char* c_str() const;
  const std::string &s_str() const &;
  std::string &&s_str() &&;

  int compare(const char *cs) const 
  { 
    const auto lhs_size = length();
    if (cs == nullptr)
    {
      if(lhs_size) return -1;
      return 0;
    }
    const auto rhs_size = strlen(cs);
    const auto result = CT::compare(c_str(), cs, std::min(lhs_size, rhs_size));
    if (result != 0) return result;
    if (lhs_size < rhs_size) return -1;
    if (lhs_size > rhs_size) return 1;

    return 0;
  }

  int compare( std::string_view str ) const
  {
    const auto lhs_size = length();
    const auto rhs_size = str.length();
    const auto result = CT::compare(c_str(), str.data(), std::min(lhs_size, rhs_size));
    if (result != 0) return result;
    if (lhs_size < rhs_size) return -1;
    if (lhs_size > rhs_size) return 1;

    return 0;
  }

private:
  std::string data_;
};

/** @class char_simple_traits
 *  Minimalist subset of character traits.
 */
struct char_simple_traits {
public:
  /** Convert 'n' characters, a no-op */
  static void convert( char *, size_t )
  {}

  /** Compare null-terminated strings */
  // len is min of two string lengths
  static int compare( const char * c1, const char * c2, size_t len ) noexcept;
};

/** @class char_label_traits
 *  All upper case, spaces and control characters converted to '_'.
 *  Does not support eof, stream types, or state types.
 */
struct char_label_traits {
public:
  /** Convert 'n' characters */
  static void convert( char * c, size_t n );

  /** Compare null-terminated strings as per conversion */
  // len is min of two string lengths
  static int compare( const char * c1, const char * c2, size_t len ) noexcept;
};


template<class CT>
typename StringBase<CT>::iterator
StringBase<CT>::begin()
{ return data_.data(); }

template<class CT>
typename StringBase<CT>::const_iterator
StringBase<CT>::begin() const
{ return data_.c_str(); }

template<class CT>
typename StringBase<CT>::iterator
StringBase<CT>::end()
{
  return begin() + size();
}

template<class CT>
typename StringBase<CT>::const_iterator
StringBase<CT>::end() const
{
  return begin() + size();
}

template<class CT>
const char* StringBase<CT>::c_str() const
{ return data_.c_str(); }

template <class CT>
const std::string &StringBase<CT>::s_str() const &
{
  return data_;
}

template <class CT>
std::string &&StringBase<CT>::s_str() &&
{
  return std::move(data_);
}

//----------------------------------------------------------------------

template <class CT>
template <class CT2>
StringBase<CT>::StringBase(const StringBase<CT2> &cs) : data_(cs.s_str())
{
  traits_type::convert(begin(), length());
}

template <class CT>
template <class CT2>
StringBase<CT>::StringBase(StringBase<CT2> &&cs) : data_(std::move(cs).s_str())
{
  traits_type::convert(begin(), length());
}

template <class CT>
StringBase<CT>::StringBase(const char *cs, typename StringBase<CT>::size_type n) : data_(cs, n)
{
  traits_type::convert(begin(), length());
}

template <class CT>
template <class It>
StringBase<CT>::StringBase(It l_begin, It l_end) : data_(l_begin, l_end)
{
  traits_type::convert(begin(), length());
}

template <class CT>
StringBase<CT>::StringBase(const char *cs)
{
  if (cs != nullptr) {
    data_ = std::string(cs);
    traits_type::convert(begin(), length());
  }
}

template <class CT>
template <class CT2>
StringBase<CT> &StringBase<CT>::operator=(StringBase<CT2> &&cs)
{
  return assign(std::move(cs));
}

template<class CT>
StringBase<CT>&
StringBase<CT>::operator= ( const char * cs ) { 
  if (cs == nullptr) {
    data_ = std::string();
    return *this;
  }
  return assign( cs, strlen(cs) ); 
}

template <class CT>
StringBase<CT> &StringBase<CT>::operator=(std::string &&cs)
{
  data_ = std::move(cs);
  traits_type::convert(begin(), length());

  return *this;
}

template<class CT>
StringBase<CT> &
StringBase<CT>::append( const char * cs, const typename StringBase<CT>::size_type n )
{
  auto orig_len = length();
  data_.append(cs, n);
  traits_type::convert(begin() + orig_len, length() - orig_len);
  return *this ;
}

template <>
inline StringBase<char_label_traits> &StringBase<char_label_traits>::append(
    const char *cs, const typename StringBase<char_label_traits>::size_type n)
{
  if (n > 0) {
    auto orig_len = length();
    data_.reserve(orig_len + n + 1);
    if (orig_len != 0) data_.append("_");
    data_.append(cs, n);
    traits_type::convert(begin() + orig_len, length() - orig_len);
  }
  return *this;
}

template<class CT>
StringBase<CT> & StringBase<CT>::append( const char * cs ) {
  return ( cs == nullptr ) ? *this : append( cs, strlen(cs) );
}
template<class CT>
StringBase<CT>&
StringBase<CT>::operator+= ( const char * cs ) {
  return (cs == nullptr) ? *this : append( cs, strlen(cs) ); 
}


template<class CT, class CT2>
bool operator== ( const StringBase<CT> & lhs,
		   const StringBase<CT2> & rhs )
{ return lhs.compare(rhs) == 0 ; }

template<class CT, class CT2>
bool operator!= ( const StringBase<CT> & lhs,
		   const StringBase<CT2> & rhs )
{ return lhs.compare(rhs) != 0 ; }

template<class CT, class CT2>
bool operator< ( const StringBase<CT> & lhs,
		  const StringBase<CT2> & rhs )
{ return lhs.compare(rhs) < 0 ; }

template<class CT, class CT2>
bool operator<= ( const StringBase<CT> & lhs,
		   const StringBase<CT2> & rhs )
{ return lhs.compare(rhs) <= 0 ; }

template<class CT, class CT2>
bool operator> ( const StringBase<CT> & lhs,
		  const StringBase<CT2> & rhs )
{ return lhs.compare(rhs) > 0 ; }

template<class CT, class CT2>
bool operator>= ( const StringBase<CT> & lhs,
		   const StringBase<CT2> & rhs )
{ return lhs.compare(rhs) >= 0 ; }


template<class CT>
bool operator== ( const StringBase<CT> & lhs,
		   std::string_view rhs )
{ return lhs.compare(rhs) == 0 ; }

template<class CT>
bool operator!= ( const StringBase<CT> & lhs,
		   std::string_view rhs )
{ return lhs.compare(rhs) != 0 ; }

template<class CT>
bool operator< ( const StringBase<CT> & lhs,
		  std::string_view rhs )
{ return lhs.compare(rhs) < 0 ; }

template<class CT>
bool operator<= ( const StringBase<CT> & lhs,
		   std::string_view rhs )
{ return lhs.compare(rhs) <= 0 ; }

template<class CT>
bool operator> ( const StringBase<CT> & lhs,
		  std::string_view rhs )
{ return lhs.compare(rhs) > 0 ; }

template<class CT>
bool operator>= ( const StringBase<CT> & lhs,
		   std::string_view rhs )
{ return lhs.compare(rhs) >= 0 ; }

template<class CT>
bool operator== ( std::string_view lhs,
		   const StringBase<CT> & rhs)
{ return rhs.compare(lhs) == 0 ; }

template<class CT>
bool operator!= ( std::string_view lhs,
		   const StringBase<CT> & rhs)
{ return rhs.compare(lhs) != 0 ; }

template<class CT>
bool operator< ( std::string_view lhs,
		  const StringBase<CT> & rhs)
{ return rhs.compare(lhs) > 0 ; }

template<class CT>
bool operator<= ( std::string_view lhs,
		   const StringBase<CT> & rhs)
{ return rhs.compare(lhs) >= 0 ; }

template<class CT>
bool operator> ( std::string_view lhs,
		  const StringBase<CT> & rhs)
{ return rhs.compare(lhs) < 0 ; }

template<class CT>
bool operator>= ( std::string_view lhs,
		   const StringBase<CT> & rhs)
{ return rhs.compare(lhs) <= 0 ; }

template<class CT, class CT2>
StringBase<CT>
operator+( const StringBase<CT> &cs1, const StringBase<CT2> &cs2) {
  StringBase<CT> t(cs1);
  t.append(cs2.c_str(), cs2.size());
  return t;
}

template<class CT>
StringBase<CT>
operator+( const StringBase<CT> &cs1, std::string_view cs2) {
  StringBase<CT> t(cs1);
  t.append(cs2.data(), cs2.length());
  return t;
}

template<class CT>
StringBase<CT>
operator+(std::string_view cs1, const StringBase<CT> & cs2 ) {
  StringBase<CT> t(cs1);
  t.append(cs2.c_str(), cs2.length());
  return t;
}

///
/// @}
///

} // namespace sierra

#endif // USE_CISTRING

#endif // STK_UTIL_DIAG_STRING_H
