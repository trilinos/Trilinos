// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef STK_UTIL_UTIL_TypeUtil_hpp
#define STK_UTIL_UTIL_TypeUtil_hpp

namespace sierra {

//-----------------------------------

template<typename T> struct IsFundamentalType ;

template<> struct IsFundamentalType<         char>  { enum { value = true }; };
template<> struct IsFundamentalType<unsigned char>  { enum { value = true }; };
template<> struct IsFundamentalType<signed   char>  { enum { value = true }; };

template<> struct IsFundamentalType<short> { enum { value = true }; };
template<> struct IsFundamentalType<int>   { enum { value = true }; };
template<> struct IsFundamentalType<long>  { enum { value = true }; };

template<> struct IsFundamentalType<unsigned short> { enum { value = true }; };
template<> struct IsFundamentalType<unsigned int>   { enum { value = true }; };
template<> struct IsFundamentalType<unsigned long>  { enum { value = true }; };

template<> struct IsFundamentalType<float>  { enum { value = true }; };
template<> struct IsFundamentalType<double> { enum { value = true }; };

template<typename T> struct IsFundamentalType { enum { value = false }; };

//-----------------------------------

template<typename T>
class TypeTraits
{
public: //private:
  template <class U> struct Traits
  {
    enum {is_pointer = false};
    enum {is_reference = false};
    typedef U Type;
  };

  template <class U> struct Traits<U*>
  {
    enum {is_pointer = true};
    enum {is_reference = false};
    typedef U Type;
  };

  template <class U> struct Traits<U&>
  {
    enum {is_pointer = false};
    enum {is_reference = true};
    typedef U Type;
  };

  template <class U> struct ConstTraits
  {
    enum {is_const = false};
    typedef U Type;
  };

  template <class U> struct ConstTraits<const U>
  {
    enum {is_const = true};
    typedef U Type;
  };

public:
  enum {isPointer = Traits<T>::is_pointer};
  enum {isReference = Traits<T>::is_reference};
  enum {isConst = ConstTraits<T>::is_const};
  typedef typename Traits<T>::Type BaseType;
};

} // namespace sierra

#endif // STK_UTIL_UTIL_TypeUtil_hpp
