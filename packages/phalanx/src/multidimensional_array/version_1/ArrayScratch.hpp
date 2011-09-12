/*
// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

*/

#ifndef util_ArrayScratch_hpp
#define util_ArrayScratch_hpp

#include <Array.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 >
class Array< Dimension , DataType[N1] > {
private:
  enum { Size = N1 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 ) const
    { return const_cast<data_type*>( m_data )[ m_dim(i1) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
};

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 , unsigned N2 >
class Array< Dimension , DataType[N1][N2] > {
private:
  enum { Size = N1 * N2 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 , unsigned i2 ) const
    { return const_cast<data_type*>( m_data )[ m_dim(i1,i2) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1,N2) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 , unsigned N2 ,
                              unsigned N3 >
class Array< Dimension , DataType[N1][N2][N3] > {
private:
  enum { Size = N1 * N2 * N3 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 ) const
    { return const_cast<data_type*>( m_data )[ m_dim(i1,i2,i3) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1,N2,N3) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 , unsigned N2 ,
                              unsigned N3 , unsigned N4 >
class Array< Dimension , DataType[N1][N2][N3][N4] > {
private:
  enum { Size = N1 * N2 * N3 * N4 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ) const
    { return const_cast<data_type*>( m_data )[ m_dim(i1,i2,i3,i4) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1,N2,N3,N4) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 , unsigned N2 ,
                              unsigned N3 , unsigned N4 ,
                              unsigned N5 >
class Array< Dimension , DataType[N1][N2][N3][N4][N5] > {
private:
  enum { Size = N1 * N2 * N3 * N4 * N5 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 ) const
    { return const_cast<data_type*>( m_data )[ m_dim(i1,i2,i3,i4,i5) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1,N2,N3,N4,N5) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 , unsigned N2 ,
                              unsigned N3 , unsigned N4 ,
                              unsigned N5 , unsigned N6 >
class Array< Dimension , DataType[N1][N2][N3][N4][N5][N6] > {
private:
  enum { Size = N1 * N2 * N3 * N4 * N5 * N6 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ) const
    { return const_cast<data_type*>( m_data )[ m_dim(i1,i2,i3,i4,i5,i6) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1,N2,N3,N4,N5,N6) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 , unsigned N2 ,
                              unsigned N3 , unsigned N4 ,
                              unsigned N5 , unsigned N6 ,
                              unsigned N7 >
class Array< Dimension , DataType[N1][N2][N3][N4][N5][N6][N7] > {
private:
  enum { Size = N1 * N2 * N3 * N4 * N5 * N6 * N7 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ,
                          unsigned i7 ) const
    { return const_cast<data_type*>( m_data )[ m_dim(i1,i2,i3,i4,i5,i6,i7) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1,N2,N3,N4,N5,N6,N7) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

template< class Dimension ,
          typename DataType , unsigned N1 , unsigned N2 ,
                              unsigned N3 , unsigned N4 ,
                              unsigned N5 , unsigned N6 ,
                              unsigned N7 , unsigned N8 >
class Array< Dimension , DataType[N1][N2][N3][N4][N5][N6][N7][N8] > {
private:
  enum { Size = N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8 };
  Dimension m_dim ;
  DataType  m_data[ Size ];
public:
  enum { DEEP_COPY = true };

  typedef DataType data_type ;
  typedef Dimension dimension_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return const_cast<data_type*>( m_data ); }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ,
                          unsigned i7 , unsigned i8 ) const
    { return const_cast<data_type*>( m_data )
                [ m_dim(i1,i2,i3,i4,i5,i6,i7,i8) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(N1,N2,N3,N4,N5,N6,N7,N8) {}

  Array( const Array & rhs )
   { for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; } }
 
  Array & operator = ( const Array & rhs )
   {
     for ( unsigned i = 0 ; i < Size ; ++i ) { m_data[i] = rhs.m_data[i] ; }
     return *this ;
   }
  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

}

#endif

