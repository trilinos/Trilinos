#ifndef KOKKOS_CONTAINERS_PAIR_HPP
#define KOKKOS_CONTAINERS_PAIR_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Functional.hpp>
#include <utility>

namespace Kokkos {

template <class T1, class T2>
struct pair
{
  typedef T1 first_type;
  typedef T2 second_type;

  first_type  first;
  second_type second;

  KOKKOS_FORCEINLINE_FUNCTION
  pair()
    : first(), second()
  {}


  KOKKOS_FORCEINLINE_FUNCTION
  pair(const first_type & f, const second_type & s)
    : first(f), second(s)
  {}

  template <class U, class V>
  KOKKOS_FORCEINLINE_FUNCTION
  pair( const pair<U,V> &p)
    : first(p.first), second(p.second)
  {}

  template <class U, class V>
  KOKKOS_FORCEINLINE_FUNCTION
  pair<T1, T2> & operator=(const pair<U,V> &p)
  {
    first = p.first;
    second = p.second;
    return *this;
  }

  // from std::pair<U,V>
  template <class U, class V>
  pair( const std::pair<U,V> &p)
    : first(p.first), second(p.second)
  {}

  // NOTE: this is not a conversion operator since
  // a conversion operator made the relational operators
  // have ambiguous definitions
  std::pair<T1,T2> to_std_pair() const
  { return std::make_pair(first,second); }
};

// Relational operators
template <class T1, class T2>
KOKKOS_FORCEINLINE_FUNCTION
bool operator== (const pair<T1,T2>& lhs, const pair<T1,T2>& rhs)
{ return lhs.first==rhs.first && lhs.second==rhs.second; }

template <class T1, class T2>
KOKKOS_FORCEINLINE_FUNCTION
bool operator!= (const pair<T1,T2>& lhs, const pair<T1,T2>& rhs)
{ return !(lhs==rhs); }

template <class T1, class T2>
KOKKOS_FORCEINLINE_FUNCTION
bool operator<  (const pair<T1,T2>& lhs, const pair<T1,T2>& rhs)
{ return lhs.first<rhs.first || (!(rhs.first<lhs.first) && lhs.second<rhs.second); }

template <class T1, class T2>
KOKKOS_FORCEINLINE_FUNCTION
bool operator<= (const pair<T1,T2>& lhs, const pair<T1,T2>& rhs)
{ return !(rhs<lhs); }

template <class T1, class T2>
KOKKOS_FORCEINLINE_FUNCTION
bool operator>  (const pair<T1,T2>& lhs, const pair<T1,T2>& rhs)
{ return rhs<lhs; }

template <class T1, class T2>
KOKKOS_FORCEINLINE_FUNCTION
bool operator>= (const pair<T1,T2>& lhs, const pair<T1,T2>& rhs)
{ return !(lhs<rhs); }

// make_pair
template <class T1,class T2>
KOKKOS_FORCEINLINE_FUNCTION
pair<T1,T2> make_pair (T1 x, T2 y)
{ return ( pair<T1,T2>(x,y) ); }


template <class T1,class T2>
KOKKOS_FORCEINLINE_FUNCTION
pair<T1 &,T2 &> tie (T1 & x, T2 & y)
{ return ( pair<T1 &,T2 &>(x,y) ); }

template <class T1>
struct pair<T1,void>
{
  typedef T1 first_type;
  typedef void second_type;

  first_type  first;
  enum { second = 0 };

  KOKKOS_FORCEINLINE_FUNCTION
  pair()
    : first()
  {}

  KOKKOS_FORCEINLINE_FUNCTION
  pair(const first_type & f)
    : first(f)
  {}

  KOKKOS_FORCEINLINE_FUNCTION
  pair(const first_type & f, int)
    : first(f)
  {}

  template <class U>
  KOKKOS_FORCEINLINE_FUNCTION
  pair( const pair<U,void> &p)
    : first(p.first)
  {}

  template <class U>
  KOKKOS_FORCEINLINE_FUNCTION
  pair<T1, void> & operator=(const pair<U,void> &p)
  {
    first = p.first;
    return *this;
  }
};


// Relational operators
template <class T1>
KOKKOS_FORCEINLINE_FUNCTION
bool operator== (const pair<T1,void>& lhs, const pair<T1,void>& rhs)
{ return lhs.first==rhs.first; }

template <class T1>
KOKKOS_FORCEINLINE_FUNCTION
bool operator!= (const pair<T1,void>& lhs, const pair<T1,void>& rhs)
{ return !(lhs==rhs); }

template <class T1>
KOKKOS_FORCEINLINE_FUNCTION
bool operator<  (const pair<T1,void>& lhs, const pair<T1,void>& rhs)
{ return lhs.first<rhs.first; }

template <class T1>
KOKKOS_FORCEINLINE_FUNCTION
bool operator<= (const pair<T1,void>& lhs, const pair<T1,void>& rhs)
{ return !(rhs<lhs); }

template <class T1>
KOKKOS_FORCEINLINE_FUNCTION
bool operator>  (const pair<T1,void>& lhs, const pair<T1,void>& rhs)
{ return rhs<lhs; }

template <class T1>
KOKKOS_FORCEINLINE_FUNCTION
bool operator>= (const pair<T1,void>& lhs, const pair<T1,void>& rhs)
{ return !(lhs<rhs); }


template <class T1, class T2>
struct hash< pair<T1,T2> >
{
  typedef pair<T1,T2> argument_type;
  typedef pair<T1,T2> first_argument_type;
  typedef uint32_t second_argument_type;
  typedef uint32_t result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()( const pair<T1,T2> & p, uint32_t seed = 0u) const
  {
    typedef hash<T1> hash1;
    typedef hash<T2> hash2;
    return hash1(p.first, hash2(p.second,seed));
  }
};


} // namespace Kokkos


#endif //KOKKOS_CONTAINERS_PAIR_HPP
