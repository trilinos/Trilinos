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

#ifndef util_Dimension_hpp
#define util_Dimension_hpp

#include <iosfwd>
#include <string>

namespace phdmesh {

//----------------------------------------------------------------------

/** \brief Maximum rank of DimNatural or DimFortran multi-index map. */
enum { DimMaxRank = 8 };

template< class Tag1        , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class DimNatural ;

template< class Tag1        , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class DimFortran ;

//----------------------------------------------------------------------
// Formatted output operators for the dimension types.

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
inline
std::ostream & operator <<
  ( std::ostream & ,
    const DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & );

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
inline
std::ostream & operator <<
  ( std::ostream & ,
    const DimNatural<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & );

//----------------------------------------------------------------------
/** \class DimTag
 *  \brief Virtual base class for DimNatural or DimFortran ordinate tags.
 *
 *  In addition to the three virtual methods a derived class 'DerivedDimTag'
 *  must provide a static method with the following signature:
 *  -  static const DerivedDimTag & descriptor();
 */
struct DimTag {

  virtual ~DimTag();

  /** The name of this tag, typically the name of the derived class. */
  virtual const char * name() const = 0 ;

  /** Generate a descriptive text label for a dimension and index.
   *  Default method simply converts the integer to a string.
   *  Throw an exception for invalid arguments.
   */
  virtual std::string to_string( unsigned size , unsigned index ) const ;

  /** Generate an index from a dimension and text label.
   *  Default method simply converts the string to an integer.
   *  Throw an exception for invalid arguments.
   */
  virtual unsigned to_index( unsigned size , const std::string & ) const ;

  //  static const DerivedDimTagType & descriptor();
};

//----------------------------------------------------------------------
/** \struct DimHelper
 *  \brief Internal helper class for copying and comparing.
 */
template< unsigned N , unsigned NGood > struct DimHelper ;


template<>
struct DimHelper<0,0> {
  static void zero( unsigned * dst ) {}
  static void copy( unsigned * dst , const unsigned * src ) {}
  static bool equal( const unsigned * lhs , const unsigned * rhs ) 
  {return true;}
};

template<>
struct DimHelper<1,1> {
  static void good() {}

  static void zero( unsigned * dst )
    { *dst = 0 ; }

  static void copy( unsigned * dst , const unsigned * src )
    { *dst = *src ; }

  static bool equal( const unsigned * lhs , const unsigned * rhs )
    { return *lhs == *rhs ; }
};

template< unsigned N >
struct DimHelper<N,N> {
private:
  enum { M = N - 1 };
public:
  static void good() {}

  static void zero( unsigned * dst )
    { *dst = 0 ; DimHelper<M,M>::zero( dst + 1 ); }

  static void copy( unsigned * dst , const unsigned * src )
    { *dst = *src ; DimHelper<M,M>::copy( dst + 1 , src + 1 ); }

  static bool equal( const unsigned * lhs , const unsigned * rhs )
    { return *lhs == *rhs && DimHelper<M,M>::equal( lhs + 1 , rhs + 1 ); }
};


//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** \class DimTagList
 *  \brief Internal helper class for list of ordinate tags.
 */
template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct DimTagList ;

template<>
struct DimTagList<void,void,void,void,void,void,void,void> {

  typedef DimNatural<void> ReverseNatural ;
  typedef DimFortran<void> ReverseFortran ;

  typedef DimNatural<void> TruncateNatural ;
  typedef DimFortran<void> TruncateFortran ;

  enum { Rank = 0 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = NULL ;
      tags[1] = NULL ;
      tags[2] = NULL ;
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 >
struct DimTagList<Tag1,void,void,void,void,void,void,void> {

  typedef DimNatural<Tag1> ReverseNatural ;
  typedef DimFortran<Tag1> ReverseFortran ;

  typedef DimNatural<void> TruncateNatural ;
  typedef DimFortran<void> TruncateFortran ;

  enum { Rank = 1 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = NULL ;
      tags[2] = NULL ;
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 , class Tag2 >
struct DimTagList<Tag1,Tag2,void,void,void,void,void,void> {

  typedef DimNatural<Tag2,Tag1> ReverseNatural ;
  typedef DimFortran<Tag2,Tag1> ReverseFortran ;

  typedef DimNatural<Tag2> TruncateNatural ;
  typedef DimFortran<Tag1> TruncateFortran ;

  enum { Rank = 2 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = & Tag2::descriptor();
      tags[2] = NULL ;
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 , class Tag2 , class Tag3 >
struct DimTagList<Tag1,Tag2,Tag3,void,void,void,void,void> {

  typedef DimNatural<Tag3,Tag2,Tag1> ReverseNatural ;
  typedef DimFortran<Tag3,Tag2,Tag1> ReverseFortran ;

  typedef DimNatural<     Tag2,Tag3> TruncateNatural ;
  typedef DimFortran<Tag1,Tag2>      TruncateFortran ;

  enum { Rank = 3 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = & Tag2::descriptor();
      tags[2] = & Tag3::descriptor();
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct DimTagList<Tag1,Tag2,Tag3,Tag4,void,void,void,void> {

  typedef DimNatural<Tag4,Tag3,Tag2,Tag1> ReverseNatural ;
  typedef DimFortran<Tag4,Tag3,Tag2,Tag1> ReverseFortran ;

  typedef DimNatural<     Tag2,Tag3,Tag4> TruncateNatural ;
  typedef DimFortran<Tag1,Tag2,Tag3>      TruncateFortran ;

  enum { Rank = 4 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = & Tag2::descriptor();
      tags[2] = & Tag3::descriptor();
      tags[3] = & Tag4::descriptor();
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
struct DimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> {

  typedef DimNatural<Tag5,Tag4,Tag3,Tag2,Tag1> ReverseNatural ;
  typedef DimFortran<Tag5,Tag4,Tag3,Tag2,Tag1> ReverseFortran ;

  typedef DimNatural<     Tag2,Tag3,Tag4,Tag5> TruncateNatural ;
  typedef DimFortran<Tag1,Tag2,Tag3,Tag4>      TruncateFortran ;

  enum { Rank = 5 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = & Tag2::descriptor();
      tags[2] = & Tag3::descriptor();
      tags[3] = & Tag4::descriptor();
      tags[4] = & Tag5::descriptor();
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct DimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> {

  typedef DimNatural<Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> ReverseNatural ;
  typedef DimFortran<Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> ReverseFortran ;

  typedef DimNatural<     Tag2,Tag3,Tag4,Tag5,Tag6> TruncateNatural ;
  typedef DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5>      TruncateFortran ;

  enum { Rank = 6 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = & Tag2::descriptor();
      tags[2] = & Tag3::descriptor();
      tags[3] = & Tag4::descriptor();
      tags[4] = & Tag5::descriptor();
      tags[5] = & Tag6::descriptor();
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct DimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> {

  typedef DimNatural<Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> ReverseNatural ;
  typedef DimFortran<Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> ReverseFortran ;

  typedef DimNatural<     Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> TruncateNatural ;
  typedef DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6>      TruncateFortran ;

  enum { Rank = 7 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = & Tag2::descriptor();
      tags[2] = & Tag3::descriptor();
      tags[3] = & Tag4::descriptor();
      tags[4] = & Tag5::descriptor();
      tags[5] = & Tag6::descriptor();
      tags[6] = & Tag7::descriptor();
      tags[7] = NULL ;
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct DimTagList {

  typedef DimNatural<Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> ReverseNatural ;
  typedef DimFortran<Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> ReverseFortran ;

  typedef DimNatural<     Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> TruncateNatural ;
  typedef DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>      TruncateFortran ;

  enum { Rank = 8 };

  const DimTag * tags[ DimMaxRank ];

  DimTagList()
    {
      tags[0] = & Tag1::descriptor();
      tags[1] = & Tag2::descriptor();
      tags[2] = & Tag3::descriptor();
      tags[3] = & Tag4::descriptor();
      tags[4] = & Tag5::descriptor();
      tags[5] = & Tag6::descriptor();
      tags[6] = & Tag7::descriptor();
      tags[7] = & Tag8::descriptor();
    }

private:
  DimTagList( const DimTagList & );
  DimTagList & operator = ( const DimTagList & );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** \class DimBase
 *  \brief Multi-index map with no compile-time knowledge of rank or tags.
 *
 *  This class allows the strongly typed DimNatural and DimFortran rank
 *  and size information to be stored and retrieved in weakly typed data.
 */
class DimBase {
private:
  typedef DimHelper< DimMaxRank , DimMaxRank > Helper ;
  unsigned m_rank ;
  unsigned m_stride[ DimMaxRank ];
public:

  /** \brief Rank of the multi-index */
  unsigned rank() const { return m_rank ; }

  /** \brief Size of the multi-index */
  unsigned size() const { return m_rank ? m_stride[ m_rank - 1 ] : 0 ; }

  /** \brief Natural multi-index map for indices[0 .. rank()-1] */
  unsigned natural_map( const unsigned * indices ) const ;

  /** \brief Fortran multi-index map for indices[0 .. rank()-1] */
  unsigned fortran_map( const unsigned * indices ) const ;

  /** \brief Natural multi-index inverse map for indices[0 .. rank()-1] */
  void natural_inv( unsigned offset , unsigned * indices ) const ;

  /** \brief Fortran multi-index inverse map for indices[0 .. rank()-1] */
  void fortran_inv( unsigned offset , unsigned * indices ) const ;

  /** \brief Natural multi-index map validity check for indices[0 .. rank()-1]
   */
  bool natural_valid( const unsigned * indices ) const ;

  /** \brief Fortran multi-index map validity check for indices[0 .. rank()-1]
   */
  bool fortran_valid( const unsigned * indices ) const ;

  /** \brief Natural multi-index map sizes[0 .. rank()-1] */
  void natural_size( unsigned * sizes ) const ;

  /** \brief Fortran multi-index map sizes[0 .. rank()-1] */
  void fortran_size( unsigned * sizes ) const ;

  //--------------------------------

  DimBase() : m_rank(0) { Helper::zero( m_stride ); }

  DimBase( const DimBase & rhs )
    : m_rank( rhs.m_rank )
    { Helper::copy( m_stride , rhs.m_stride ); }

  DimBase & operator = ( const DimBase & rhs )
    {
      m_rank = rhs.m_rank ;
      Helper::copy( m_stride , rhs.m_stride );
      return *this ;
    }

  bool operator == ( const DimBase & rhs ) const
    { return Helper::equal( m_stride , rhs.m_stride ); }

  bool operator != ( const DimBase & rhs ) const
    { return ! Helper::equal( m_stride , rhs.m_stride ); }

  //--------------------------------
  // 'DimFortran' and 'DimNatural' copy constructors:

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  explicit
  DimBase( const DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & rhs )
    {
      enum { N = DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::Rank };
      enum { M = DimMaxRank - N };
      m_rank = N ;
      DimHelper<N,N>::copy( m_stride , rhs.m_stride );
      DimHelper<M,M>::zero( m_stride + N );
    }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  explicit
  DimBase( const DimNatural<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & rhs )
    {
      enum { N = DimNatural<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::Rank };
      enum { M = DimMaxRank - N };
      m_rank = N ;
      DimHelper<N,N>::copy( m_stride , rhs.m_stride );
      DimHelper<M,M>::zero( m_stride + N );
    }

  //--------------------------------
  // 'DimFortran' and 'DimNatural' assignment operators:

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  DimBase & operator =
    ( const DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & rhs )
    {
      enum { N = DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::Rank };
      enum { M = DimMaxRank - N };
      m_rank = N ;
      DimHelper<N,N>::copy( m_stride , rhs.m_stride );
      DimHelper<M,M>::zero( m_stride + N );
      return *this ;
    }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  DimBase & operator =
    ( const DimNatural<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & rhs )
    {
      enum { N = DimNatural<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::Rank };
      enum { M = DimMaxRank - N };
      m_rank = N ;
      DimHelper<N,N>::copy( m_stride , rhs.m_stride );
      DimHelper<M,M>::zero( m_stride + N );
      return *this ;
    }

  //--------------------------------

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
    friend class DimFortran ;

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
    friend class DimNatural ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** \class DimNatural
 *  \brief Multi-index map with natural ordering.
 *
 *  Strongly typed multi-index map using Lexicographical semantics:
 *  right-most index has a stride of one and
 *  left-most index has the largest stride.
 */
template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class DimNatural {
private:
  typedef DimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> TagList ;

  typedef DimHelper< TagList::Rank , TagList::Rank > Helper ;

  unsigned m_stride[ TagList::Rank ? TagList::Rank : 1 ];
public:

  /** \brief Number of multi-index ordinates */
  enum { Rank = TagList::Rank };

  /** \brief Total size = upper bound of the multi-index range space. */ 
  unsigned size() const { return m_stride[ Rank - 1 ]; }

  //--------------------------------
  /** \brief Map for a Rank 8 natural multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 , unsigned i7 , unsigned i8 ) const
    {
      DimHelper<Rank,8>::good();
      return i8 + i7 * m_stride[0] + i6 * m_stride[1] + i5 * m_stride[2] +
                  i4 * m_stride[3] + i3 * m_stride[4] + i2 * m_stride[5] +
                  i1 * m_stride[6] ;
    }

  /** \brief Map for a Rank 7 natural multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 , unsigned i7 ) const
    {
      DimHelper<Rank,7>::good();
      return i7 + i6 * m_stride[0] + i5 * m_stride[1] + i4 * m_stride[2] +
                  i3 * m_stride[3] + i2 * m_stride[4] + i1 * m_stride[5] ;
    }

  /** \brief Map for a Rank 6 natural multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 ) const
    {
      DimHelper<Rank,6>::good();
      return i6 + i5 * m_stride[0] + i4 * m_stride[1] + i3 * m_stride[2] +
                  i2 * m_stride[3] + i1 * m_stride[4] ;
    }

  /** \brief Map for a Rank 5 natural multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 ) const
    {
      DimHelper<Rank,5>::good();
      return i5 + i4 * m_stride[0] + i3 * m_stride[1] + i2 * m_stride[2] +
                  i1 * m_stride[3] ;
    }

  /** \brief Map for a Rank 4 natural multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ) const
    {
      DimHelper<Rank,4>::good();
      return i4 + i3 * m_stride[0] + i2 * m_stride[1] + i1 * m_stride[2] ;
    }

  /** \brief Map for a Rank 3 natural multi-index */
  unsigned operator()( unsigned i1 , unsigned i2 , unsigned i3 ) const
    {
      DimHelper<Rank,3>::good();
      return i3 + i2 * m_stride[0] + i1 * m_stride[1] ;
    }

  /** \brief Map for a Rank 2 natural multi-index */
  unsigned operator()( unsigned i1 , unsigned i2 ) const
    {
      DimHelper<Rank,2>::good();
      return i2 + i1 * m_stride[0] ;
    }

  /** \brief Map for a Rank 1 natural multi-index */
  unsigned operator()( unsigned i1 ) const
    {
      DimHelper<Rank,1>::good();
      return i1 ;
    }

  //--------------------------------
  /** \brief Inverse map for a Rank 8 natural multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5, unsigned & i6, unsigned & i7, unsigned & i8) const
    {
      DimHelper<Rank,8>::good();
      i1 = offset / m_stride[6] ; offset %= m_stride[6] ;
      i2 = offset / m_stride[5] ; offset %= m_stride[5] ;
      i3 = offset / m_stride[4] ; offset %= m_stride[4] ;
      i4 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i5 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i6 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i7 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i8 = offset ;
    }

  /** \brief Inverse map for a Rank 7 natural multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5, unsigned & i6, unsigned & i7) const
    {
      DimHelper<Rank,7>::good();
      i1 = offset / m_stride[5] ; offset %= m_stride[5] ;
      i2 = offset / m_stride[4] ; offset %= m_stride[4] ;
      i3 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i4 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i5 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i6 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i7 = offset ;
    }

  /** \brief Inverse map for a Rank 6 natural multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5, unsigned & i6) const
    {
      DimHelper<Rank,6>::good();
      i1 = offset / m_stride[4] ; offset %= m_stride[4] ;
      i2 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i3 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i4 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i5 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i6 = offset ;
    }

  /** \brief Inverse map for a Rank 5 natural multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5) const
    {
      DimHelper<Rank,5>::good();
      i1 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i2 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i3 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i4 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i5 = offset ;
    }

  /** \brief Inverse map for a Rank 4 natural multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4) const
    {
      DimHelper<Rank,4>::good();
      i1 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i2 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i3 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i4 = offset ;
    }

  /** \brief Inverse map for a Rank 3 natural multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3) const
    {
      DimHelper<Rank,3>::good();
      i1 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i3 = offset ;
    }

  /** \brief Inverse map for a Rank 2 natural multi-index */
  void inv( unsigned offset , unsigned & i1, unsigned & i2) const
    {
      DimHelper<Rank,2>::good();
      i1 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i2 = offset ;
    }

  /** \brief Inverse map for a Rank 1 natural multi-index */
  void inv( unsigned offset , unsigned & i1) const
    {
      DimHelper<Rank,1>::good();
      i1 = offset ;
    }

  //--------------------------------
  /** \brief Map input validity check for Rank 8 natural multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 , unsigned i7 , unsigned i8 ) const
    {
      DimHelper<Rank,8>::good();
      return i8               < m_stride[0] &&
             i7 * m_stride[0] < m_stride[1] &&
             i6 * m_stride[1] < m_stride[2] &&
             i5 * m_stride[2] < m_stride[3] &&
             i4 * m_stride[3] < m_stride[4] &&
             i3 * m_stride[4] < m_stride[5] &&
             i2 * m_stride[5] < m_stride[6] &&
             i1 * m_stride[6] < m_stride[7] ;
    }

  /** \brief Map input validity check for Rank 7 natural multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 , unsigned i7 ) const
    {
      DimHelper<Rank,7>::good();
      return i7               < m_stride[0] &&
             i6 * m_stride[0] < m_stride[1] &&
             i5 * m_stride[1] < m_stride[2] &&
             i4 * m_stride[2] < m_stride[3] &&
             i3 * m_stride[3] < m_stride[4] &&
             i2 * m_stride[4] < m_stride[5] &&
             i1 * m_stride[5] < m_stride[6] ;
    }

  /** \brief Map input validity check for Rank 6 natural multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 ) const
    {
      DimHelper<Rank,6>::good();
      return i6               < m_stride[0] &&
             i5 * m_stride[0] < m_stride[1] &&
             i4 * m_stride[1] < m_stride[2] &&
             i3 * m_stride[2] < m_stride[3] &&
             i2 * m_stride[3] < m_stride[4] &&
             i1 * m_stride[4] < m_stride[5] ;
    }

  /** \brief Map input validity check for Rank 5 natural multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 ) const
    {
      DimHelper<Rank,5>::good();
      return i5               < m_stride[0] &&
             i4 * m_stride[0] < m_stride[1] &&
             i3 * m_stride[1] < m_stride[2] &&
             i2 * m_stride[2] < m_stride[3] &&
             i1 * m_stride[3] < m_stride[4] ;
    }

  /** \brief Map input validity check for Rank 4 natural multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ) const
    {
      DimHelper<Rank,4>::good();
      return i4               < m_stride[0] &&
             i3 * m_stride[0] < m_stride[1] &&
             i2 * m_stride[1] < m_stride[2] &&
             i1 * m_stride[2] < m_stride[3] ;
    }

  /** \brief Map input validity check for Rank 3 natural multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 ) const
    {
      DimHelper<Rank,3>::good();
      return i3               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] &&
             i1 * m_stride[1] < m_stride[2] ;
    }

  /** \brief Map input validity check for Rank 2 natural multi-index */
  bool valid( unsigned i1 , unsigned i2 ) const
    {
      DimHelper<Rank,2>::good();
      return i2               < m_stride[0] &&
             i1 * m_stride[0] < m_stride[1] ;
    }

  /** \brief Map input validity check for Rank 1 natural multi-index */
  bool valid( unsigned i1 ) const
    {
      DimHelper<Rank,1>::good();
      return i1 < m_stride[0] ;
    }

  //--------------------------------
  /** \brief Size of each ordinate for Rank 8 natural multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5, unsigned & n6, unsigned & n7, unsigned & n8) const
    {
      DimHelper<Rank,8>::good();
      n8 = m_stride[0] ;
      n7 = m_stride[1] / m_stride[0] ;
      n6 = m_stride[2] / m_stride[1] ;
      n5 = m_stride[3] / m_stride[2] ;
      n4 = m_stride[4] / m_stride[3] ;
      n3 = m_stride[5] / m_stride[4] ;
      n2 = m_stride[6] / m_stride[5] ;
      n1 = m_stride[7] / m_stride[6] ;
    }

  /** \brief Size of each ordinate for Rank 7 natural multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5, unsigned & n6, unsigned & n7) const
    {
      DimHelper<Rank,7>::good();
      n7 = m_stride[0] ;
      n6 = m_stride[1] / m_stride[0] ;
      n5 = m_stride[2] / m_stride[1] ;
      n4 = m_stride[3] / m_stride[2] ;
      n3 = m_stride[4] / m_stride[3] ;
      n2 = m_stride[5] / m_stride[4] ;
      n1 = m_stride[6] / m_stride[5] ;
    }

  /** \brief Size of each ordinate for Rank 6 natural multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5, unsigned & n6) const
    {
      DimHelper<Rank,6>::good();
      n6 = m_stride[0] ;
      n5 = m_stride[1] / m_stride[0] ;
      n4 = m_stride[2] / m_stride[1] ;
      n3 = m_stride[3] / m_stride[2] ;
      n2 = m_stride[4] / m_stride[3] ;
      n1 = m_stride[5] / m_stride[4] ;
    }

  /** \brief Size of each ordinate for Rank 5 natural multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5) const
    {
      DimHelper<Rank,5>::good();
      n5 = m_stride[0] ;
      n4 = m_stride[1] / m_stride[0] ;
      n3 = m_stride[2] / m_stride[1] ;
      n2 = m_stride[3] / m_stride[2] ;
      n1 = m_stride[4] / m_stride[3] ;
    }

  /** \brief Size of each ordinate for Rank 4 natural multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4) const
    {
      DimHelper<Rank,4>::good();
      n4 = m_stride[0] ;
      n3 = m_stride[1] / m_stride[0] ;
      n2 = m_stride[2] / m_stride[1] ;
      n1 = m_stride[3] / m_stride[2] ;
    }

  /** \brief Size of each ordinate for Rank 3 natural multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3) const
    {
      DimHelper<Rank,3>::good();
      n3 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
      n1 = m_stride[2] / m_stride[1] ;
    }

  /** \brief Size of each ordinate for Rank 2 natural multi-index */
  void size( unsigned & n1, unsigned & n2) const
    {
      DimHelper<Rank,2>::good();
      n2 = m_stride[0] ;
      n1 = m_stride[1] / m_stride[0] ;
    }

  /** \brief Size of each ordinate for Rank 1 natural multi-index */
  void size( unsigned & n1) const
    {
      DimHelper<Rank,1>::good();
      n1 = m_stride[0] ;
    }

  //--------------------------------
  /** \brief Constructor for Rank 8 natural multi-index */
  DimNatural( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 , unsigned n6 , unsigned n7 , unsigned n8 )
    {
      DimHelper<Rank,8>::good();
      m_stride[7] = n1 * (
      m_stride[6] = n2 * (
      m_stride[5] = n3 * (
      m_stride[4] = n4 * (
      m_stride[3] = n5 * (
      m_stride[2] = n6 * (
      m_stride[1] = n7 * (
      m_stride[0] = n8 )))))));
    }

  /** \brief Constructor for Rank 7 natural multi-index */
  DimNatural( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 , unsigned n6 , unsigned n7 )
    {
      DimHelper<Rank,7>::good();
      m_stride[6] = n1 * (
      m_stride[5] = n2 * (
      m_stride[4] = n3 * (
      m_stride[3] = n4 * (
      m_stride[2] = n5 * (
      m_stride[1] = n6 * (
      m_stride[0] = n7 ))))));
    }

  /** \brief Constructor for Rank 6 natural multi-index */
  DimNatural( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 , unsigned n6 )
    {
      DimHelper<Rank,6>::good();
      m_stride[5] = n1 * (
      m_stride[4] = n2 * (
      m_stride[3] = n3 * (
      m_stride[2] = n4 * (
      m_stride[1] = n5 * (
      m_stride[0] = n6 )))));
    }

  /** \brief Constructor for Rank 5 natural multi-index */
  DimNatural( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 )
    {
      DimHelper<Rank,5>::good();
      m_stride[4] = n1 * (
      m_stride[3] = n2 * (
      m_stride[2] = n3 * (
      m_stride[1] = n4 * (
      m_stride[0] = n5 ))));
    }

  /** \brief Constructor for Rank 4 natural multi-index */
  DimNatural( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 )
    {
      DimHelper<Rank,4>::good();
      m_stride[3] = n1 * (
      m_stride[2] = n2 * (
      m_stride[1] = n3 * (
      m_stride[0] = n4 )));
    }

  /** \brief Constructor for Rank 3 natural multi-index */
  DimNatural( unsigned n1 , unsigned n2 , unsigned n3 )
    {
      DimHelper<Rank,3>::good();
      m_stride[2] = n1 * (
      m_stride[1] = n2 * (
      m_stride[0] = n3 ));
    }

  /** \brief Constructor for Rank 2 natural multi-index */
  DimNatural( unsigned n1 , unsigned n2 )
    {
      DimHelper<Rank,2>::good();
      m_stride[1] = n1 * (
      m_stride[0] = n2 );
    }

  /** \brief Constructor for Rank 1 natural multi-index */
  explicit DimNatural( unsigned n1 )
    {
      DimHelper<Rank,1>::good();
      m_stride[0] = n1 ;
    }

  //--------------------------------

  DimNatural() { Helper::zero( m_stride ); }

  DimNatural( const DimNatural & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); }

  DimNatural & operator = ( const DimNatural & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); return *this ; }

  bool operator == ( const DimNatural & rhs ) const
    { return Helper::equal( m_stride , rhs.m_stride ); }

  bool operator != ( const DimNatural & rhs ) const
    { return ! Helper::equal( m_stride , rhs.m_stride ); }

  //--------------------------------
  /** \brief Compatible DimFortran type */
  typedef typename TagList::ReverseFortran DimFortranType ;

  /** \brief Implicit copy constructor for compatible DimFortran type */
  DimNatural( const DimFortranType & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); }

  /** \brief Assignment operator for compatible DimFortran type */
  DimNatural & operator = ( const DimFortranType & rhs )
      { Helper::copy( m_stride , rhs.m_stride ); return *this ; }

  bool operator == ( const DimFortranType & rhs ) const
    { return Helper::equal( m_stride , rhs.m_stride ); }

  bool operator != ( const DimFortranType & rhs ) const
    { return ! Helper::equal( m_stride , rhs.m_stride ); }

  //--------------------------------
  /** \brief Truncated dimension type */
  typedef typename TagList::TruncateNatural TruncatedType ;

  operator const TruncatedType & () const
    { return * reinterpret_cast<const TruncatedType*>(this); }

  /** \brief Offset for the greatest stride to support truncation. */
  unsigned operator[]( unsigned i ) const
    {
      DimHelper<Rank-1,Rank-1>::good();
      return i * m_stride[ Rank - 2 ];
    }
  //--------------------------------

  explicit
  DimNatural( const DimBase & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); }

  DimNatural( unsigned n , const DimBase & rhs )
    {
      DimHelper<Rank-1,Rank-1>::copy( m_stride , rhs.m_stride );
      m_stride[ Rank - 1 ] = n * m_stride[ Rank - 2 ] ;
    }

  //--------------------------------

  friend class DimBase ;

  template< class T1 , class T2 , class T3 , class T4 ,
            class T5 , class T6 , class T7 , class T8 >
  friend class DimNatural ;

  template< class T1 , class T2 , class T3 , class T4 ,
            class T5 , class T6 , class T7 , class T8 >
  friend class DimFortran ;
};

//----------------------------------------------------------------------
/** \class DimFortran
 *  \brief Multi-index map with Fortran ordering.
 *
 *  Strongly typed multi-index map using Fortran semantics:
 *  left-most index has a stride of one and
 *  right-most index has the largest stride.
 */
template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class DimFortran {
private:
  typedef DimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> TagList ;

  typedef DimHelper< TagList::Rank , TagList::Rank > Helper ;

  unsigned m_stride[ TagList::Rank ? TagList::Rank : 1 ];
public:

  /** \brief Number of multi-index ordinates */
  enum { Rank = TagList::Rank };

  /** \brief Total size = upper bound of the multi-index range space. */ 
  unsigned size() const { return m_stride[ Rank - 1 ]; }

  //--------------------------------
  /** \brief Map for a Rank 8 Fortran multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 , unsigned i7 , unsigned i8 ) const
    {
      DimHelper<Rank,8>::good();
      return i1 + i2 * m_stride[0] + i3 * m_stride[1] + i4 * m_stride[2] +
                  i5 * m_stride[3] + i6 * m_stride[4] + i7 * m_stride[5] +
                  i8 * m_stride[6] ;
    }

  /** \brief Map for a Rank 7 Fortran multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 , unsigned i7 ) const
    {
      DimHelper<Rank,7>::good();
      return i1 + i2 * m_stride[0] + i3 * m_stride[1] + i4 * m_stride[2] +
                  i5 * m_stride[3] + i6 * m_stride[4] + i7 * m_stride[5] ;
    }

  /** \brief Map for a Rank 6 Fortran multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 ) const
    {
      DimHelper<Rank,6>::good();
      return i1 + i2 * m_stride[0] + i3 * m_stride[1] + i4 * m_stride[2] +
                  i5 * m_stride[3] + i6 * m_stride[4] ;
    }

  /** \brief Map for a Rank 5 Fortran multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 ) const
    {
      DimHelper<Rank,5>::good();
      return i1 + i2 * m_stride[0] + i3 * m_stride[1] + i4 * m_stride[2] +
                  i5 * m_stride[3] ;
    }

  /** \brief Map for a Rank 4 Fortran multi-index */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ) const
    {
      DimHelper<Rank,4>::good();
      return i1 + i2 * m_stride[0] + i3 * m_stride[1] + i4 * m_stride[2] ;
    }

  /** \brief Map for a Rank 3 Fortran multi-index */
  unsigned operator()( unsigned i1 , unsigned i2 , unsigned i3 ) const
    {
      DimHelper<Rank,3>::good();
      return i1 + i2 * m_stride[0] + i3 * m_stride[1] ;
    }

  /** \brief Map for a Rank 2 Fortran multi-index */
  unsigned operator()( unsigned i1 , unsigned i2 ) const
    {
      DimHelper<Rank,2>::good();
      return i1 + i2 * m_stride[0] ;
    }

  /** \brief Map for a Rank 1 Fortran multi-index */
  unsigned operator()( unsigned i1 ) const
    {
      DimHelper<Rank,1>::good();
      return i1 ;
    }

  //--------------------------------
  /** \brief Inverse map for a Rank 8 Fortran multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5, unsigned & i6, unsigned & i7, unsigned & i8) const
    {
      DimHelper<Rank,8>::good();
      i8 = offset / m_stride[6] ; offset %= m_stride[6] ;
      i7 = offset / m_stride[5] ; offset %= m_stride[5] ;
      i6 = offset / m_stride[4] ; offset %= m_stride[4] ;
      i5 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i4 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i3 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i1 = offset ;
    }

  /** \brief Inverse map for a Rank 7 Fortran multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5, unsigned & i6, unsigned & i7) const
    {
      DimHelper<Rank,7>::good();
      i7 = offset / m_stride[5] ; offset %= m_stride[5] ;
      i6 = offset / m_stride[4] ; offset %= m_stride[4] ;
      i5 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i4 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i3 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i1 = offset ;
    }

  /** \brief Inverse map for a Rank 6 Fortran multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5, unsigned & i6) const
    {
      DimHelper<Rank,6>::good();
      i6 = offset / m_stride[4] ; offset %= m_stride[4] ;
      i5 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i4 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i3 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i1 = offset ;
    }

  /** \brief Inverse map for a Rank 5 Fortran multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4,
            unsigned & i5) const
    {
      DimHelper<Rank,5>::good();
      i5 = offset / m_stride[3] ; offset %= m_stride[3] ;
      i4 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i3 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i1 = offset ;
    }

  /** \brief Inverse map for a Rank 4 Fortran multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3, unsigned & i4) const
    {
      DimHelper<Rank,4>::good();
      i4 = offset / m_stride[2] ; offset %= m_stride[2] ;
      i3 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i1 = offset ;
    }

  /** \brief Inverse map for a Rank 3 Fortran multi-index */
  void inv( unsigned offset ,
            unsigned & i1, unsigned & i2, unsigned & i3) const
    {
      DimHelper<Rank,3>::good();
      i3 = offset / m_stride[1] ; offset %= m_stride[1] ;
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i1 = offset ;
    }

  /** \brief Inverse map for a Rank 2 Fortran multi-index */
  void inv( unsigned offset , unsigned & i1, unsigned & i2) const
    {
      DimHelper<Rank,2>::good();
      i2 = offset / m_stride[0] ; offset %= m_stride[0] ;
      i1 = offset ;
    }

  /** \brief Inverse map for a Rank 1 Fortran multi-index */
  void inv( unsigned offset , unsigned & i1) const
    {
      DimHelper<Rank,1>::good();
      i1 = offset ;
    }

  //--------------------------------
  /** \brief Map input validity check for Rank 8 Fortran multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 , unsigned i7 , unsigned i8 ) const
    {
      DimHelper<Rank,8>::good();
      return i1               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] &&
             i3 * m_stride[1] < m_stride[2] &&
             i4 * m_stride[2] < m_stride[3] &&
             i5 * m_stride[3] < m_stride[4] &&
             i6 * m_stride[4] < m_stride[5] &&
             i7 * m_stride[5] < m_stride[6] &&
             i8 * m_stride[6] < m_stride[7] ;
    }

  /** \brief Map input validity check for Rank 7 Fortran multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 , unsigned i7 ) const
    {
      DimHelper<Rank,7>::good();
      return i1               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] &&
             i3 * m_stride[1] < m_stride[2] &&
             i4 * m_stride[2] < m_stride[3] &&
             i5 * m_stride[3] < m_stride[4] &&
             i6 * m_stride[4] < m_stride[5] &&
             i7 * m_stride[5] < m_stride[6] ;
    }

  /** \brief Map input validity check for Rank 6 Fortran multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 ) const
    {
      DimHelper<Rank,6>::good();
      return i1               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] &&
             i3 * m_stride[1] < m_stride[2] &&
             i4 * m_stride[2] < m_stride[3] &&
             i5 * m_stride[3] < m_stride[4] &&
             i6 * m_stride[4] < m_stride[5] ;
    }

  /** \brief Map input validity check for Rank 5 Fortran multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 ) const
    {
      DimHelper<Rank,5>::good();
      return i1               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] &&
             i3 * m_stride[1] < m_stride[2] &&
             i4 * m_stride[2] < m_stride[3] &&
             i5 * m_stride[3] < m_stride[4] ;
    }

  /** \brief Map input validity check for Rank 4 Fortran multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ) const
    {
      DimHelper<Rank,4>::good();
      return i1               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] &&
             i3 * m_stride[1] < m_stride[2] &&
             i4 * m_stride[2] < m_stride[3] ;
    }

  /** \brief Map input validity check for Rank 3 Fortran multi-index */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 ) const
    {
      DimHelper<Rank,3>::good();
      return i1               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] &&
             i3 * m_stride[1] < m_stride[2] ;
    }

  /** \brief Map input validity check for Rank 2 Fortran multi-index */
  bool valid( unsigned i1 , unsigned i2 ) const
    {
      DimHelper<Rank,2>::good();
      return i1               < m_stride[0] &&
             i2 * m_stride[0] < m_stride[1] ;
    }

  /** \brief Map input validity check for Rank 1 Fortran multi-index */
  bool valid( unsigned i1 ) const
    {
      DimHelper<Rank,1>::good();
      return i1 < m_stride[0] ;
    }

  //--------------------------------
  /** \brief Size of each ordinate for Rank 8 Fortran multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5, unsigned & n6, unsigned & n7, unsigned & n8) const
    {
      DimHelper<Rank,8>::good();
      n1 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
      n3 = m_stride[2] / m_stride[1] ;
      n4 = m_stride[3] / m_stride[2] ;
      n5 = m_stride[4] / m_stride[3] ;
      n6 = m_stride[5] / m_stride[4] ;
      n7 = m_stride[6] / m_stride[5] ;
      n8 = m_stride[7] / m_stride[6] ;
    }

  /** \brief Size of each ordinate for Rank 7 Fortran multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5, unsigned & n6, unsigned & n7) const
    {
      DimHelper<Rank,7>::good();
      n1 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
      n3 = m_stride[2] / m_stride[1] ;
      n4 = m_stride[3] / m_stride[2] ;
      n5 = m_stride[4] / m_stride[3] ;
      n6 = m_stride[5] / m_stride[4] ;
      n7 = m_stride[6] / m_stride[5] ;
    }

  /** \brief Size of each ordinate for Rank 6 Fortran multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5, unsigned & n6) const
    {
      DimHelper<Rank,6>::good();
      n1 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
      n3 = m_stride[2] / m_stride[1] ;
      n4 = m_stride[3] / m_stride[2] ;
      n5 = m_stride[4] / m_stride[3] ;
      n6 = m_stride[5] / m_stride[4] ;
    }

  /** \brief Size of each ordinate for Rank 5 Fortran multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4,
             unsigned & n5) const
    {
      DimHelper<Rank,5>::good();
      n1 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
      n3 = m_stride[2] / m_stride[1] ;
      n4 = m_stride[3] / m_stride[2] ;
      n5 = m_stride[4] / m_stride[3] ;
    }

  /** \brief Size of each ordinate for Rank 4 Fortran multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3, unsigned & n4) const
    {
      DimHelper<Rank,4>::good();
      n1 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
      n3 = m_stride[2] / m_stride[1] ;
      n4 = m_stride[3] / m_stride[2] ;
    }

  /** \brief Size of each ordinate for Rank 3 Fortran multi-index */
  void size( unsigned & n1, unsigned & n2, unsigned & n3) const
    {
      DimHelper<Rank,3>::good();
      n1 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
      n3 = m_stride[2] / m_stride[1] ;
    }

  /** \brief Size of each ordinate for Rank 2 Fortran multi-index */
  void size( unsigned & n1, unsigned & n2) const
    {
      DimHelper<Rank,2>::good();
      n1 = m_stride[0] ;
      n2 = m_stride[1] / m_stride[0] ;
    }

  /** \brief Size of each ordinate for Rank 1 Fortran multi-index */
  void size( unsigned & n1) const
    {
      DimHelper<Rank,1>::good();
      n1 = m_stride[0] ;
    }

  //--------------------------------
  /** \brief Constructor for Rank 8 Fortran multi-index */
  DimFortran( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 , unsigned n6 , unsigned n7 , unsigned n8 )
    {
      DimHelper<Rank,8>::good();
      m_stride[7] = n8 * (
      m_stride[6] = n7 * (
      m_stride[5] = n6 * (
      m_stride[4] = n5 * (
      m_stride[3] = n4 * (
      m_stride[2] = n3 * (
      m_stride[1] = n2 * (
      m_stride[0] = n1 )))))));
    }

  /** \brief Constructor for Rank 7 Fortran multi-index */
  DimFortran( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 , unsigned n6 , unsigned n7 )
    {
      DimHelper<Rank,7>::good();
      m_stride[6] = n7 * (
      m_stride[5] = n6 * (
      m_stride[4] = n5 * (
      m_stride[3] = n4 * (
      m_stride[2] = n3 * (
      m_stride[1] = n2 * (
      m_stride[0] = n1 ))))));
    }

  /** \brief Constructor for Rank 6 Fortran multi-index */
  DimFortran( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 , unsigned n6 )
    {
      DimHelper<Rank,6>::good();
      m_stride[5] = n6 * (
      m_stride[4] = n5 * (
      m_stride[3] = n4 * (
      m_stride[2] = n3 * (
      m_stride[1] = n2 * (
      m_stride[0] = n1 )))));
    }

  /** \brief Constructor for Rank 5 Fortran multi-index */
  DimFortran( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
              unsigned n5 )
    {
      DimHelper<Rank,5>::good();
      m_stride[4] = n5 * (
      m_stride[3] = n4 * (
      m_stride[2] = n3 * (
      m_stride[1] = n2 * (
      m_stride[0] = n1 ))));
    }

  /** \brief Constructor for Rank 4 Fortran multi-index */
  DimFortran( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 )
    {
      DimHelper<Rank,4>::good();
      m_stride[3] = n4 * (
      m_stride[2] = n3 * (
      m_stride[1] = n2 * (
      m_stride[0] = n1 )));
    }

  /** \brief Constructor for Rank 3 Fortran multi-index */
  DimFortran( unsigned n1 , unsigned n2 , unsigned n3 )
    {
      DimHelper<Rank,3>::good();
      m_stride[2] = n3 * (
      m_stride[1] = n2 * (
      m_stride[0] = n1 ));
    }

  /** \brief Constructor for Rank 2 Fortran multi-index */
  DimFortran( unsigned n1 , unsigned n2 )
    {
      DimHelper<Rank,2>::good();
      m_stride[1] = n2 * (
      m_stride[0] = n1 );
    }

  /** \brief Constructor for Rank 1 Fortran multi-index */
  explicit DimFortran( unsigned n1 )
    {
      DimHelper<Rank,1>::good();
      m_stride[0] = n1 ;
    }

  //--------------------------------

  DimFortran() { Helper::zero( m_stride ); }

  DimFortran( const DimFortran & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); }

  DimFortran & operator = ( const DimFortran & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); return *this ; }

  bool operator == ( const DimFortran & rhs ) const
    { return Helper::equal( m_stride , rhs.m_stride ); }

  bool operator != ( const DimFortran & rhs ) const
    { return ! Helper::equal( m_stride , rhs.m_stride ); }

  //--------------------------------
  /** \brief Compatible DimNatural type */
  typedef typename TagList::ReverseNatural DimNaturalType ;

  DimFortran( const DimNaturalType & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); }

  DimFortran & operator = ( const DimNaturalType & rhs )
      { Helper::copy( m_stride , rhs.m_stride ); }

  bool operator == ( const DimNaturalType & rhs ) const
    { return Helper::equal( m_stride , rhs.m_stride ); }

  bool operator != ( const DimNaturalType & rhs ) const
    { return ! Helper::equal( m_stride , rhs.m_stride ); }

  //--------------------------------
  /** \brief Truncated dimension type */
  typedef typename TagList::TruncateFortran TruncatedType ;

  operator const TruncatedType & () const
    { return * reinterpret_cast<const TruncatedType*>(this); }

  /** \brief Offset for the greatest stride to support truncation. */
  unsigned operator[]( unsigned i ) const
    {
      DimHelper<Rank-1,Rank-1>::good();
      return i * m_stride[ Rank - 2 ];
    }
  //--------------------------------

  explicit
  DimFortran( const DimBase & rhs )
    { Helper::copy( m_stride , rhs.m_stride ); }

  DimFortran( const DimBase & rhs , unsigned n )
    {
      DimHelper<Rank-1,Rank-1>::copy( m_stride , rhs.m_stride );
      m_stride[Rank-1] = n * m_stride[Rank-2] ;
    }

  //--------------------------------

  friend class DimBase ;

  template< class T1 , class T2 , class T3 , class T4 ,
            class T5 , class T6 , class T7 , class T8 >
  friend class DimNatural ;

  template< class T1 , class T2 , class T3 , class T4 ,
            class T5 , class T6 , class T7 , class T8 >
  friend class DimFortran ;
};


//----------------------------------------------------------------------
//----------------------------------------------------------------------

void print( std::ostream & ,
            const DimBase & ,
            const DimTag ** , 
            const bool is_dim_natural );

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
inline
std::ostream & operator <<
  ( std::ostream & s ,
    const DimFortran<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & dim )
{
  DimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> tmp ;

  print( s, DimBase( dim ), tmp.tags , false );

  return s ;
}

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
inline
std::ostream & operator <<
  ( std::ostream & s ,
    const DimNatural<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & dim )
{
  DimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> tmp ;

  print( s, DimBase( dim ), tmp.tags , true );

  return s ;
}

//----------------------------------------------------------------------

}

#endif


