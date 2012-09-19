/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
/**
 * @file
 * @author H. Carter Edwards
 * @date   January 2006
 *
 * Templates for static / compile-time checking
 */

#ifndef STK_UTIL_UTIL_TypeUtil_hpp
#define STK_UTIL_UTIL_TypeUtil_hpp

/**
 * @file
 *
 */

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
