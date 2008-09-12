/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   November 2006
 */

#ifndef util_SimpleArrayOps_hpp
#define util_SimpleArrayOps_hpp

namespace phdmesh {

// Basic operations for compile-time fixed length arrays
//
//  Example:  Sum<5>(x,y)   results in  x[i] += y[i]  for i=0..4

template< unsigned n, unsigned i=0> struct Copy ;
template< unsigned n, unsigned i=0> struct Sum ;
template< unsigned n, unsigned i=0> struct Prod ;
template< unsigned n, unsigned i=0> struct Max ;
template< unsigned n, unsigned i=0> struct Min ;
template< unsigned n, unsigned i=0> struct BitOr ;
template< unsigned n, unsigned i=0> struct BitAnd ;
template< unsigned n, unsigned i=0> struct InnerProduct ;
template< unsigned n, unsigned i=0> struct Compare ;

//----------------------------------------------------------------------

template<unsigned i>
struct Copy<0,i> {
  enum { N = 0 };
  Copy() {}
  template<typename T> inline Copy( T * const , const T * const ) {}
  template<typename T> inline Copy( T * const , const T ) {}
};

template<unsigned i>
struct Sum<0,i> {
  enum { N = 0 };
  Sum() {}
  template<typename T> inline Sum( T * const , const T * const ) {}
  template<typename T> inline Sum( T * const , const T , const T * const ) {}
};

template<unsigned i>
struct Prod<0,i> {
  enum { N = 0 };
  Prod() {}
  template<typename T> inline Prod( T * const , const T * const ) {}
};

template<unsigned i>
struct Max<0,i> {
  enum { N = 0 };
  Max() {}
  template<typename T> inline Max( T * const , const T * const ) {}
};

template<unsigned i>
struct Min<0,i> {
  enum { N = 0 };
  Min() {}
  template<typename T> inline Min( T * const , const T * const ) {}
};

template<unsigned i>
struct BitOr<0,i> {
  enum { N = 0 };
  BitOr() {}
  template<typename T> inline BitOr( T * const , const T * const ) {}
};

template<unsigned i>
struct BitAnd<0,i> {
  enum { N = 0 };
  BitAnd() {}
  template<typename T> inline BitAnd( T * const , const T * const ) {}
};

template<unsigned i>
struct InnerProduct<0,i> {
  enum { N = 0 };
  InnerProduct() {}
  template<typename T>
  inline InnerProduct( T & value , const T * const x , const T * const y ) {}
};

template<unsigned i>
struct Compare<0,i> {
  enum { N = 0 };
  Compare() {}

  template<typename T>
  inline static bool equal( const T * const , const T * const )
    { return true ; }

  template<typename T>
  inline static bool not_equal( const T * const , const T * const )
    { return false ; }

  template<typename T>
  inline static bool less( const T * const , const T * const )
    { return false ; }

  template<typename T>
  inline static bool less_equal( const T * const , const T * const )
    { return true ; }

  template<typename T>
  inline static bool greater( const T * const , const T * const )
    { return false ; }

  template<typename T>
  inline static bool greater_equal( const T * const , const T * const )
    { return true ; }
};


template<unsigned n,unsigned i>
struct Copy {
  enum { N = n , I = i };
  Copy() {}

  template<typename T>
  inline Copy( T * const dst , const T * const src )
    { dst[I] = src[I] ; Copy<N-1,I+1>(dst,src); }

  template<typename T>
  inline Copy( T * const dst , const T src )
    { dst[I] = src ; Copy<N-1,I+1>(dst,src); }
};

template<unsigned n,unsigned i>
struct Sum {
  enum { N = n , I = i };
  Sum() {}

  template<typename T>
  inline Sum( T * const dst , const T * const src )
    { dst[I] += src[I] ; Sum<N-1,I+1>(dst,src); }

  template<typename T>
  inline Sum( T * const dst , const T a , const T * const src )
    { dst[I] += a * src[I] ; Sum<N-1,I+1>(dst,a,src); }
};

template<unsigned n,unsigned i>
struct Prod {
  enum { N = n , I = i };
  Prod() {}
  template<typename T>
  inline Prod( T * const dst , const T * const src )
    { dst[I] *= src[I] ; Prod<N-1,I+1>(dst,src); }
};

template<unsigned n,unsigned i>
struct BitOr {
  enum { N = n , I = i };
  BitOr() {}
  template<typename T>
  inline BitOr( T * const dst , const T * const src )
    { dst[I] |= src[I] ; BitOr<N-1,I+1>(dst,src); }
};

template<unsigned n,unsigned i>
struct BitAnd {
  enum { N = n , I = i };
  BitAnd() {}
  template<typename T>
  inline BitAnd( T * const dst , const T * const src )
    { dst[I] &= src[I] ; BitAnd<N-1,I+1>(dst,src); }
};

template<unsigned n,unsigned i>
struct Max {
  enum { N = n , I = i };
  Max() {}
  template<typename T>
  inline Max( T * const dst , const T * const src )
    { if ( dst[I] < src[I] ) { dst[I] = src[I] ; } Max<N-1,I+1>(dst,src); }
};

template<unsigned n,unsigned i>
struct Min {
  enum { N = n , I = i };
  Min() {}
  template<typename T>
  inline Min( T * const dst , const T * const src )
    { if ( src[I] < dst[I] ) { dst[I] = src[I] ; } Min<N-1,I+1>(dst,src); }
};

template<unsigned n,unsigned i>
struct InnerProduct {
  enum { N = n , I = i };
  InnerProduct() {}
  template<typename T>
  inline InnerProduct( T & value , const T * const x , const T * const y )
    { value += x[I] * y[I] ; InnerProduct<N-1,I+1>( value , x , y ); }
};


template<unsigned n,unsigned i>
struct Compare {
  enum { N = n , I = i };
  Compare() {}

  template<typename T>
  inline static bool equal( const T * const dst , const T * const src )
    { return dst[I] == src[I] && Compare<N-1,I+1>::equal(dst,src); }

  template<typename T>
  inline static bool not_equal( const T * const dst , const T * const src )
    { return dst[I] != src[I] || Compare<N-1,I+1>::not_equal(dst,src); }

  template<typename T>
  inline static bool less( const T * const dst , const T * const src )
    {
      return dst[I] != src[I] ? dst[I] < src[I]
                              : Compare<N-1,I+1>::less(dst,src);
    }

  template<typename T>
  inline static bool less_equal( const T * const dst , const T * const src )
    {
      return dst[I] != src[I] ? dst[I] < src[I]
                              : Compare<N-1,I+1>::less_equal(dst,src);
    }

  template<typename T>
  inline static bool greater( const T * const dst , const T * const src )
    {
      return dst[I] != src[I] ? dst[I] > src[I]
                              : Compare<N-1,I+1>::greater(dst,src);
    }

  template<typename T>
  inline static bool greater_equal( const T * const dst , const T * const src )
    {
      return dst[I] != src[I] ? dst[I] > src[I]
                              : Compare<N-1,I+1>::greater_equal(dst,src);
    }
};

//-----------------------------------

} // namespace phdmesh

#endif

