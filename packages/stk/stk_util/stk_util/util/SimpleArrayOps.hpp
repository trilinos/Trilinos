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

#ifndef stk_util_util_SimpleArrayOps_hpp
#define stk_util_util_SimpleArrayOps_hpp

namespace stk {

  /// \todo REFACTOR: Add a sub-namespace?

// Basic operations for compile-time fixed length arrays
//
//  Example:  Sum<5>(x,y)   results in  x[i] += y[i]  for i=0..4

template< unsigned n, int i=0> struct Copy ;
template< unsigned n, int i=0> struct Sum ;
template< unsigned n, int i=0> struct Prod ;
template< unsigned n, int i=0> struct Max ;
template< unsigned n, int i=0> struct Min ;
template< unsigned n, int i=0> struct BitOr ;
template< unsigned n, int i=0> struct BitAnd ;
template< unsigned n, int i=0> struct InnerProduct ;
template< unsigned n, int i=0> struct Compare ;

//----------------------------------------------------------------------

template<int i>
struct Copy<0,i> {
  enum { N = 0 };
  Copy() {}
  template<typename T> inline Copy( T * const , const T * const ) {}
  template<typename T> inline Copy( T * const , const T ) {}
};

template<int i>
struct Sum<0,i> {
  enum { N = 0 };
  Sum() {}
  template<typename T> inline Sum( T * const , const T * const ) {}
  template<typename T> inline Sum( T * const , const T , const T * const ) {}
};

template<int i>
struct Prod<0,i> {
  enum { N = 0 };
  Prod() {}
  template<typename T> inline Prod( T * const , const T * const ) {}
};

template<int i>
struct Max<0,i> {
  enum { N = 0 };
  Max() {}
  template<typename T> inline Max( T * const , const T * const ) {}
};

template<int i>
struct Min<0,i> {
  enum { N = 0 };
  Min() {}
  template<typename T> inline Min( T * const , const T * const ) {}
};

template<int i>
struct BitOr<0,i> {
  enum { N = 0 };
  BitOr() {}
  template<typename T> inline BitOr( T * const , const T * const ) {}
};

template<int i>
struct BitAnd<0,i> {
  enum { N = 0 };
  BitAnd() {}
  template<typename T> inline BitAnd( T * const , const T * const ) {}
};

template<int i>
struct InnerProduct<0,i> {
  enum { N = 0 };
  InnerProduct() {}
  template<typename T>
  inline InnerProduct( T & /*value*/ , const T * const /*x*/ , const T * const /*y*/ ) {}
};

template<int i>
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


template<unsigned n>
struct Copy<n,-1> {
  enum { N = n };
  Copy() {}

  template<typename T>
  inline Copy( T * const dst , const T * const src )
    { for ( unsigned i = 0 ; i < N ; ++i ) { dst[i] = src[i]; } }

  template<typename T>
  inline Copy( T * const dst , const T src )
    { for ( unsigned i = 0 ; i < N ; ++i ) { dst[i] = src ; } }
};

template<unsigned n>
struct Sum<n,-1> {
  enum { N = n };
  Sum() {}

  template<typename T>
  inline Sum( T * const dst , const T * const src )
    { for ( unsigned i = 0 ; i < N ; ++i ) { dst[i] += src[i]; } }

  template<typename T>
  inline Sum( T * const dst , const T a , const T * const src )
    { for ( unsigned i = 0 ; i < N ; ++i ) { dst[i] += a * src[i]; } }
};

template<unsigned n>
struct Prod<n,-1> {
  enum { N = n };
  Prod() {}
  template<typename T>
  inline Prod( T * const dst , const T * const src )
    { for ( unsigned i = 0 ; i < N ; ++i ) { dst[i] *= src[i]; } }
};

template<unsigned n>
struct BitOr<n,-1> {
  enum { N = n };
  BitOr() {}
  template<typename T>
  inline BitOr( T * const dst , const T * const src )
    { for ( unsigned i = 0 ; i < N ; ++i ) { dst[i] |= src[i]; } }
};

template<unsigned n>
struct BitAnd<n,-1> {
  enum { N = n };
  BitAnd() {}
  template<typename T>
  inline BitAnd( T * const dst , const T * const src )
    { for ( unsigned i = 0 ; i < N ; ++i ) { dst[i] &= src[i]; } }
};

template<unsigned n>
struct Max<n,-1> {
  enum { N = n };
  Max() {}
  template<typename T>
  inline Max( T * const dst , const T * const src )
    {
      for ( unsigned i = 0 ; i < N ; ++i )
        { if ( dst[i] < src[i] ) { dst[i] = src[i] ; } }
    }
};

template<unsigned n>
struct Min<n,-1> {
  enum { N = n };
  Min() {}
  template<typename T>
  inline Min( T * const dst , const T * const src )
    {
      for ( unsigned i = 0 ; i < N ; ++i )
        { if ( src[i] < dst[i] ) { dst[i] = src[i] ; } }
    }
};

template<unsigned n>
struct InnerProduct<n,-1> {
  enum { N = n };
  InnerProduct() {}
  template<typename T>
  inline InnerProduct( T & value , const T * const x , const T * const y )
    { for ( unsigned i = 0 ; i < N ; ++i ) { value += x[i] * y[i] ; } }
};


template<unsigned n>
struct Compare<n,-1> {
  enum { N = n };
  Compare() {}

  template<typename T>
  inline static bool equal( const T * const dst , const T * const src )
    {
      bool result = true ;
      for ( unsigned i = 0 ; result && i < N ; ++i )
        { result = dst[i] == src[i] ; }
      return result ;
    }

  template<typename T>
  inline static bool not_equal( const T * const dst , const T * const src )
    { return ! equal( dst , src ); }

  template<typename T>
  inline static bool less( const T * const dst , const T * const src )
    {
      unsigned i = 0 ;
      for ( ; i < N && dst[i] == src[i] ; ++i );
      return i < N && dst[i] < src[i] ;
    }

  template<typename T>
  inline static bool less_equal( const T * const dst , const T * const src )
    { return ! less( src , dst ); }

  template<typename T>
  inline static bool greater( const T * const dst , const T * const src )
    { return less( src , dst ); }

  template<typename T>
  inline static bool greater_equal( const T * const dst , const T * const src )
    { return ! less( dst , src ); }
};


template<unsigned n,int i>
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

template<unsigned n,int i>
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

template<unsigned n,int i>
struct Prod {
  enum { N = n , I = i };
  Prod() {}
  template<typename T>
  inline Prod( T * const dst , const T * const src )
    { dst[I] *= src[I] ; Prod<N-1,I+1>(dst,src); }
};

template<unsigned n,int i>
struct BitOr {
  enum { N = n , I = i };
  BitOr() {}
  template<typename T>
  inline BitOr( T * const dst , const T * const src )
    { dst[I] |= src[I] ; BitOr<N-1,I+1>(dst,src); }
};

template<unsigned n,int i>
struct BitAnd {
  enum { N = n , I = i };
  BitAnd() {}
  template<typename T>
  inline BitAnd( T * const dst , const T * const src )
    { dst[I] &= src[I] ; BitAnd<N-1,I+1>(dst,src); }
};

template<unsigned n,int i>
struct Max {
  enum { N = n , I = i };
  Max() {}
  template<typename T>
  inline Max( T * const dst , const T * const src )
    { if ( dst[I] < src[I] ) { dst[I] = src[I] ; } Max<N-1,I+1>(dst,src); }
};

template<unsigned n,int i>
struct Min {
  enum { N = n , I = i };
  Min() {}
  template<typename T>
  inline Min( T * const dst , const T * const src )
    { if ( src[I] < dst[I] ) { dst[I] = src[I] ; } Min<N-1,I+1>(dst,src); }
};

template<unsigned n,int i>
struct InnerProduct {
  enum { N = n , I = i };
  InnerProduct() {}
  template<typename T>
  inline InnerProduct( T & value , const T * const x , const T * const y )
    { value += x[I] * y[I] ; InnerProduct<N-1,I+1>( value , x , y ); }
};


template<unsigned n,int i>
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

} // namespace stk

#endif

