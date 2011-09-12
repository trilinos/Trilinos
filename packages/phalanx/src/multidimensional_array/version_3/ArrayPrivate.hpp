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

#ifndef util_ArrayPrivate_hpp
#define util_ArrayPrivate_hpp

namespace phdmesh {

//----------------------------------------------------------------------

unsigned array_stride_size(
  const unsigned  rank ,
  const unsigned * const stride );

void array_stride_to_natural_dimensions(
  const unsigned   rank ,
  const unsigned * const stride ,
        unsigned * const dim );

void array_stride_to_natural_indices(
  const unsigned   rank ,
  const unsigned * const stride ,
  const unsigned   offset ,
        unsigned * const indices );

//----------------------------------------------------------------------

void array_check_rank( const unsigned rank ,
                       const unsigned test_rank );

void array_check_ordinal( const unsigned rank ,
                          const unsigned test_ordinal );

void array_check_index( const unsigned dim ,
                        const unsigned index );

void array_check_offset(  const unsigned size ,
                          const unsigned test_offset );

void array_check_indices( const bool ,
                          const unsigned rank ,
                          const unsigned * const stride ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 );

//----------------------------------------------------------------------

template< unsigned , unsigned > struct array_check_ordinal_is_less ;

template<> struct array_check_ordinal_is_less<0,8> {};
template<> struct array_check_ordinal_is_less<1,8> {};
template<> struct array_check_ordinal_is_less<2,8> {};
template<> struct array_check_ordinal_is_less<3,8> {};
template<> struct array_check_ordinal_is_less<4,8> {};
template<> struct array_check_ordinal_is_less<5,8> {};
template<> struct array_check_ordinal_is_less<6,8> {};
template<> struct array_check_ordinal_is_less<7,8> {};

template<> struct array_check_ordinal_is_less<0,7> {};
template<> struct array_check_ordinal_is_less<1,7> {};
template<> struct array_check_ordinal_is_less<2,7> {};
template<> struct array_check_ordinal_is_less<3,7> {};
template<> struct array_check_ordinal_is_less<4,7> {};
template<> struct array_check_ordinal_is_less<5,7> {};
template<> struct array_check_ordinal_is_less<6,7> {};

template<> struct array_check_ordinal_is_less<0,6> {};
template<> struct array_check_ordinal_is_less<1,6> {};
template<> struct array_check_ordinal_is_less<2,6> {};
template<> struct array_check_ordinal_is_less<3,6> {};
template<> struct array_check_ordinal_is_less<4,6> {};
template<> struct array_check_ordinal_is_less<5,6> {};

template<> struct array_check_ordinal_is_less<0,5> {};
template<> struct array_check_ordinal_is_less<1,5> {};
template<> struct array_check_ordinal_is_less<2,5> {};
template<> struct array_check_ordinal_is_less<3,5> {};
template<> struct array_check_ordinal_is_less<4,5> {};

template<> struct array_check_ordinal_is_less<0,4> {};
template<> struct array_check_ordinal_is_less<1,4> {};
template<> struct array_check_ordinal_is_less<2,4> {};
template<> struct array_check_ordinal_is_less<3,4> {};

template<> struct array_check_ordinal_is_less<0,3> {};
template<> struct array_check_ordinal_is_less<1,3> {};
template<> struct array_check_ordinal_is_less<2,3> {};

template<> struct array_check_ordinal_is_less<0,2> {};
template<> struct array_check_ordinal_is_less<1,2> {};

template<> struct array_check_ordinal_is_less<0,1> {};

//----------------------------------------------------------------------

template< class , unsigned > struct ArrayTagAt ;

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,0>
{ typedef Tag1 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,1>
{ typedef Tag2 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,2>
{ typedef Tag3 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,3>
{ typedef Tag4 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,4>
{ typedef Tag5 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,5>
{ typedef Tag6 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,6>
{ typedef Tag7 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,7>
{ typedef Tag8 type ; };

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,TApp> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,TApp,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,Tag5,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,TApp,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,TApp,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,TApp,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,TApp,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,TApp,void,void,void,void,void,void> type ;
};

template< typename Scalar , class TApp >
struct ArrayAppend<
  Array<Scalar,RankZero,void,void,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,void,void,void,void,void,void,void,void> , TApp >
{ /* Cannot append to runtime-ranked array */ };

template< typename Scalar , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,void,void,void,void,void,void,void,void> , TApp >
{ /* Cannot append to runtime-ranked array */ };

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< class A > struct ArrayReverse ;

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> >
{
  typedef
    Array<Scalar,FortranOrder,Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1> type ;
};


template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void> type ;
};


template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void,void> type ;
};


template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag5,Tag4,Tag3,Tag2,Tag1,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag5,Tag4,Tag3,Tag2,Tag1,void,void,void> type ;
};


template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag4,Tag3,Tag2,Tag1,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag4,Tag3,Tag2,Tag1,void,void,void,void> type ;
};


template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag3,Tag2,Tag1,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag3,Tag2,Tag1,void,void,void,void,void> type ;
};


template< typename Scalar , class Tag1 , class Tag2 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,Tag2,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag2,Tag1,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag1,void,void,void,void,void,void> type ;
};


template< typename Scalar , class Tag1 >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 >
struct ArrayReverse<
  Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void> type ;
};


template< typename Scalar >
struct ArrayReverse<
  Array<Scalar,NaturalOrder,void,void,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,void,void,void,void,void,void,void,void> type ;
};

template< typename Scalar >
struct ArrayReverse<
  Array<Scalar,FortranOrder,void,void,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,NaturalOrder,void,void,void,void,void,void,void,void> type ;
};


//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< class A > struct ArrayTruncate ;

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTruncate<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8,void> type ;
};

template< typename Scalar , class Tag1 >
struct ArrayTruncate<
  Array<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void> >
{
  typedef Array<Scalar,RankZero,void,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> >
{
  typedef Array<Scalar,RankZero,void,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> type ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< ArrayOrder , unsigned Rank , unsigned Ordinal = 0 >
struct ArrayStrideDim ;

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<RankZero,Rank,Ordinal> {

  static unsigned dimension( const unsigned * )
    { return 0 ; }

  static unsigned dimension( const unsigned * , const unsigned & )
    { return 0 ; }
};

template< unsigned Rank >
struct ArrayStrideDim<FortranOrder,Rank,0> {
  static unsigned dimension( const unsigned * stride )
    { return stride[0]; }

  static unsigned dimension( const unsigned * stride ,
                             const unsigned & ordinal )
    { return ordinal ? stride[ordinal] / stride[ordinal-1] : stride[0] ; }
};

template< unsigned Rank >
struct ArrayStrideDim<NaturalOrder,Rank,0> {
  static unsigned dimension( const unsigned * stride ) { return stride[0]; }

  static unsigned dimension( const unsigned * stride ,
                             const unsigned & ordinal )
    {
      const unsigned i = ( Rank - 1 ) - ordinal ;
      return i ? stride[i] / stride[i-1] : stride[0] ;
    }
};

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<FortranOrder,Rank,Ordinal> {
  static unsigned dimension( const unsigned * stride )
    { return stride[Ordinal] / stride[Ordinal-1]; }
};

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<NaturalOrder,Rank,Ordinal> {
  static unsigned dimension( const unsigned * stride )
    {
      enum { I = ( Rank - 1 ) - Ordinal };
      return stride[I] / stride[I-1];
    }
};

//----------------------------------------------------------------------

namespace {

template< class Tag > const ArrayDimTag * array_dim_tag();

template<>
inline
const ArrayDimTag * array_dim_tag<void>() { return NULL ; }

template< class Tag >
inline
const ArrayDimTag * array_dim_tag() { return & Tag::tag(); }

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
const ArrayDimTag * const * array_dim_tags()
{
  static const ArrayDimTag * t[8] =
    {
      array_dim_tag< Tag1 >() ,
      array_dim_tag< Tag2 >() ,
      array_dim_tag< Tag3 >() ,
      array_dim_tag< Tag4 >() ,
      array_dim_tag< Tag5 >() ,
      array_dim_tag< Tag6 >() ,
      array_dim_tag< Tag7 >() ,
      array_dim_tag< Tag8 >()
    };

  return t ;
}

}

//----------------------------------------------------------------------
/** \cond */
//----------------------------------------------------------------------

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 , const unsigned & i8 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 );

//----------------------------------------------------------------------

template<> inline
unsigned array_offset<FortranOrder,8>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 , const unsigned & i8 )
{
  ARRAY_CHECK(array_check_indices(false,8,stride,i1,i2,i3,i4,i5,i6,i7,i8));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] + i6 * stride[4] +
         i7 * stride[5] + i8 * stride[6] ;
}

template<> inline
unsigned array_offset<FortranOrder,7>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 )
{
  ARRAY_CHECK(array_check_indices(false,7,stride,i1,i2,i3,i4,i5,i6,i7));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] + i6 * stride[4] +
         i7 * stride[5] ;
}

template<> inline
unsigned array_offset<FortranOrder,6>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 )
{
  ARRAY_CHECK(array_check_indices(false,6,stride,i1,i2,i3,i4,i5,i6));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] + i6 * stride[4] ;
}

template<> inline
unsigned array_offset<FortranOrder,5>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 )
{
  ARRAY_CHECK(array_check_indices(false,5,stride,i1,i2,i3,i4,i5));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] ;
}

template<> inline
unsigned array_offset<FortranOrder,4>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 )
{
  ARRAY_CHECK(array_check_indices(false,4,stride,i1,i2,i3,i4));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] ;
}

template<> inline
unsigned array_offset<FortranOrder,3>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 )
{
  ARRAY_CHECK(array_check_indices(false,3,stride,i1,i2,i3));
  return i1 + i2 * stride[0] + i3 * stride[1] ;
}

template<> inline
unsigned array_offset<FortranOrder,2>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 )
{
  ARRAY_CHECK(array_check_indices(false,2,stride,i1,i2));
  return i1 + i2 * stride[0] ;
}

template<> inline
unsigned array_offset<FortranOrder,1>(
  const unsigned * const ARRAY_CHECK( stride ) ,
  const unsigned & i1 )
{
  ARRAY_CHECK(array_check_indices(false,1,stride,i1));
  return i1 ;
}

//----------------------------------------------------------------------

template<> inline
unsigned array_offset<NaturalOrder,8>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 , const unsigned & i8 )
{
  ARRAY_CHECK(array_check_indices(true,8,stride,i1,i2,i3,i4,i5,i6,i7,i8));
  return i8             + i7 * stride[0] +
         i6 * stride[1] + i5 * stride[2] +
         i4 * stride[3] + i3 * stride[4] +
         i2 * stride[5] + i1 * stride[6] ;
}

template<> inline
unsigned array_offset<NaturalOrder,7>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 )
{
  ARRAY_CHECK(array_check_indices(true,7,stride,i1,i2,i3,i4,i5,i6,i7));
  return i7             + i6 * stride[0] +
         i5 * stride[1] + i4 * stride[2] +
         i3 * stride[3] + i2 * stride[4] +
         i1 * stride[5] ;
}

template<> inline
unsigned array_offset<NaturalOrder,6>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 )
{
  ARRAY_CHECK(array_check_indices(true,6,stride,i1,i2,i3,i4,i5,i6));
  return i6             + i5 * stride[0] +
         i4 * stride[1] + i3 * stride[2] +
         i2 * stride[3] + i1 * stride[4] ;
}

template<> inline
unsigned array_offset<NaturalOrder,5>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 )
{
  ARRAY_CHECK(array_check_indices(true,5,stride,i1,i2,i3,i4,i5));
  return i5             + i4 * stride[0] +
         i3 * stride[1] + i2 * stride[2] +
         i1 * stride[3] ;
}

template<> inline
unsigned array_offset<NaturalOrder,4>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 )
{
  ARRAY_CHECK(array_check_indices(true,4,stride,i1,i2,i3,i4));
  return i4             + i3 * stride[0] +
         i2 * stride[1] + i1 * stride[2] ;
}

template<> inline
unsigned array_offset<NaturalOrder,3>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 )
{
  ARRAY_CHECK(array_check_indices(true,3,stride,i1,i2,i3));
  return i3 + i2 * stride[0] + i1 * stride[1] ;
}

template<> inline
unsigned array_offset<NaturalOrder,2>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 )
{
  ARRAY_CHECK(array_check_indices(true,2,stride,i1,i2));
  return i2 + i1 * stride[0] ;
}

template<> inline
unsigned array_offset<NaturalOrder,1>(
  const unsigned * const ARRAY_CHECK( stride ) ,
  const unsigned & i1 )
{
  ARRAY_CHECK(array_check_indices(true,1,stride,i1));
  return i1 ;
}

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct Array<void,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
{
  enum { Rank = 8 };

  static void assign( unsigned * stride )
    {
        stride[7] = Tag8::Size * (
        stride[6] = Tag7::Size * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = Tag7::Size * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size  * (
        stride[1] = Tag2::Size  * (
        stride[0] = Tag1::Size  )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[7] = dims[7] * (
        stride[6] = dims[6] * (
        stride[5] = dims[5] * (
        stride[4] = dims[4] * (
        stride[3] = dims[3] * (
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] )))))));
    }
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct Array<void,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
{
  enum { Rank = 7 };

  static void assign( unsigned * stride )
    {
        stride[6] = Tag7::Size * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[6] = dims[6] * (
        stride[5] = dims[5] * (
        stride[4] = dims[4] * (
        stride[3] = dims[3] * (
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] ))))));
   }
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct Array<void,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
{
  enum { Rank = 6 };

  static void assign( unsigned * stride )
    {
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[5] = dims[5] * (
        stride[4] = dims[4] * (
        stride[3] = dims[3] * (
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] )))));
   }
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
struct Array<void,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
{
  enum { Rank = 5 };

  static void assign( unsigned * stride )
    {
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[4] = dims[4] * (
        stride[3] = dims[3] * (
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] ))));
   }
};


template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct Array<void,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
{
  enum { Rank = 4 };

  static void assign( unsigned * stride )
    {
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 )));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[3] = dims[3] * (
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] )));
   }
};

template< class Tag1 , class Tag2 , class Tag3 >
struct Array<void,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void>
{
  enum { Rank = 3 };

  static void assign( unsigned * stride )
    {
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 )
    {
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 ));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] ));
   }
};

template< class Tag1 , class Tag2 >
struct Array<void,FortranOrder,Tag1,Tag2,void,void,void,void,void,void>
{
  enum { Rank = 2 };

  static void assign( unsigned * stride )
    { stride[1] = Tag2::Size * ( stride[0] = Tag1::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n2 )
    { stride[1] = n2 * ( stride[0] = Tag1::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    { stride[1] = n2 * ( stride[0] = n1 ); }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    { stride[1] = dims[1] * ( stride[0] = dims[0] ); }
};

template< class Tag1 >
struct Array<void,FortranOrder,Tag1,void,void,void,void,void,void,void>
{
  enum { Rank = 1 };

  static void assign( unsigned * stride ) { stride[0] = Tag1::Size ; }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ) { stride[0] = n1 ; }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    { stride[0] = dims[0] ; }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct Array<void,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
{
  enum { Rank = 8 };

  static void assign( unsigned * stride )
    {
        stride[7] = Tag1::Size * (
        stride[6] = Tag2::Size * (
        stride[5] = Tag3::Size * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[7] = n1 * (
        stride[6] = Tag2::Size * (
        stride[5] = Tag3::Size * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = Tag3::Size * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = Tag6::Size  * (
        stride[1] = Tag7::Size  * (
        stride[0] = Tag8::Size  )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = n6 * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = n6 * (
        stride[1] = n7 * (
        stride[0] = Tag8::Size )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = n6 * (
        stride[1] = n7 * (
        stride[0] = n8 )))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[7] = dims[0] * (
        stride[6] = dims[1] * (
        stride[5] = dims[2] * (
        stride[4] = dims[3] * (
        stride[3] = dims[4] * (
        stride[2] = dims[5] * (
        stride[1] = dims[6] * (
        stride[0] = dims[7] )))))));
    }
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct Array<void,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
{
  enum { Rank = 7 };

  static void assign( unsigned * stride )
    {
        stride[6] = Tag1::Size * (
        stride[5] = Tag2::Size * (
        stride[4] = Tag3::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[6] = n1 * (
        stride[5] = Tag2::Size * (
        stride[4] = Tag3::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = Tag3::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = n5 * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = n5 * (
        stride[1] = n6 * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = n5 * (
        stride[1] = n6 * (
        stride[0] = n7 ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[6] = dims[0] * (
        stride[5] = dims[1] * (
        stride[4] = dims[2] * (
        stride[3] = dims[3] * (
        stride[2] = dims[4] * (
        stride[1] = dims[5] * (
        stride[0] = dims[6] ))))));
    }
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct Array<void,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
{
  enum { Rank = 6 };

  static void assign( unsigned * stride )
    {
        stride[5] = Tag1::Size * (
        stride[4] = Tag2::Size * (
        stride[3] = Tag3::Size * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[5] = n1 * (
        stride[4] = Tag2::Size * (
        stride[3] = Tag3::Size * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = Tag3::Size * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = n4 * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = n4 * (
        stride[1] = n5 * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = n4 * (
        stride[1] = n5 * (
        stride[0] = n6 )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[5] = dims[0] * (
        stride[4] = dims[1] * (
        stride[3] = dims[2] * (
        stride[2] = dims[3] * (
        stride[1] = dims[4] * (
        stride[0] = dims[5] )))));
    }
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
struct Array<void,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
{
  enum { Rank = 5 };

  static void assign( unsigned * stride )
    {
        stride[4] = Tag1::Size * (
        stride[3] = Tag2::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[4] = n1 * (
        stride[3] = Tag2::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = n3 * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = n3 * (
        stride[1] = n4 * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = n3 * (
        stride[1] = n4 * (
        stride[0] = n5 ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[4] = dims[0] * (
        stride[3] = dims[1] * (
        stride[2] = dims[2] * (
        stride[1] = dims[3] * (
        stride[0] = dims[4] ))));
    }
};


template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct Array<void,NaturalOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
{
  enum { Rank = 4 };

  static void assign( unsigned * stride )
    {
        stride[3] = Tag1::Size * (
        stride[2] = Tag2::Size * (
        stride[1] = Tag3::Size * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[3] = n1 * (
        stride[2] = Tag2::Size * (
        stride[1] = Tag3::Size * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[3] = n1 * (
        stride[2] = n2 * (
        stride[1] = Tag3::Size * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[3] = n1 * (
        stride[2] = n2 * (
        stride[1] = n3 * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n1 * (
        stride[2] = n2 * (
        stride[1] = n3 * (
        stride[0] = n4 )));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[3] = dims[0] * (
        stride[2] = dims[1] * (
        stride[1] = dims[2] * (
        stride[0] = dims[3] )));
    }
};

template< class Tag1 , class Tag2 , class Tag3 >
struct Array<void,NaturalOrder,Tag1,Tag2,Tag3,void,void,void,void,void>
{
  enum { Rank = 3 };

  static void assign( unsigned * stride )
    {
        stride[2] = Tag1::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag3::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[2] = n1 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag3::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[2] = n1 * (
        stride[1] = n2 * (
        stride[0] = Tag3::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[2] = n1 * (
        stride[1] = n2 * (
        stride[0] = n3 ));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[2] = dims[0] * (
        stride[1] = dims[1] * (
        stride[0] = dims[2] ));
    }
};

template< class Tag1 , class Tag2 >
struct Array<void,NaturalOrder,Tag1,Tag2,void,void,void,void,void,void>
{
  enum { Rank = 2 };

  static void assign( unsigned * stride )
    { stride[1] = Tag1::Size * ( stride[0] = Tag2::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    { stride[1] = n1 * ( stride[0] = Tag2::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    { stride[1] = n1 * ( stride[0] = n2 ); }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[1] = dims[0] * (
        stride[0] = dims[1] );
    }
};

template< class Tag1 >
struct Array<void,NaturalOrder,Tag1,void,void,void,void,void,void,void>
{
  enum { Rank = 1 };

  static void assign( unsigned * stride ) { stride[0] = Tag1::Size ; }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ) { stride[0] = n1 ; }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    { stride[0] = dims[0] ; }
};

//----------------------------------------------------------------------

template<>
struct Array<void,RankZero,void,void,void,void,void,void,void,void>
{
  enum { Rank = 0 };

  static void assign( unsigned * ) {}
};

//----------------------------------------------------------------------

}

#endif

