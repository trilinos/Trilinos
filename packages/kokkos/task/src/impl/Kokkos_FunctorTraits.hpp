/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_FUNCTORTRAITS_HPP
#define KOKKOS_IMPL_FUNCTORTRAITS_HPP

#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Functor must have a pattern_type and must not be the pattern_type

template< class FunctorType >
class FunctorTraits< FunctorType , void , void , void , void > {
public:
  typedef FunctorType                         functor_type ;
  typedef typename FunctorType::pattern_type  pattern_type ;
  typedef typename pattern_type::device_type  device_type ;
  typedef typename pattern_type::result_type  result_type ;
  typedef typename pattern_type::future_type  future_type ;
  typedef typename pattern_type::work_type    work_type ;
};

//----------------------------------------------------------------------------

template< class DeviceType , class Op1Type >
class FunctorTraits< task_serial< DeviceType , void > , Op1Type , void , void , void > {
public:
  typedef task_serial< DeviceType , void >    pattern_type ;
  typedef typename pattern_type::device_type  device_type ;
  typedef typename pattern_type::result_type  result_type ;
  typedef typename pattern_type::future_type  future_type ;
  typedef typename pattern_type::work_type    work_type ;

  class functor_type : public pattern_type {
  private:
    Op1Type m_op1 ;
  public:

    KOKKOS_INLINE_FUNCTION
    void apply() { m_op1(); }

    KOKKOS_INLINE_FUNCTION
    functor_type( const Op1Type & arg_op1 )
      : task_serial< DeviceType , void >()
      , m_op1( arg_op1 )
      {}
  };
};

template< class DeviceType , class ResultType , class Op1Type >
class FunctorTraits< task_serial< DeviceType , ResultType > , Op1Type , void , void , void > {
public:
  typedef task_serial< DeviceType , ResultType >  pattern_type ;
  typedef typename pattern_type::device_type      device_type ;
  typedef typename pattern_type::result_type      result_type ;
  typedef typename pattern_type::future_type      future_type ;
  typedef typename pattern_type::work_type        work_type ;

  class functor_type : public pattern_type {
  private:
    Op1Type m_op1 ;
  public:

    KOKKOS_INLINE_FUNCTION
    void apply( ResultType & result ) { m_op1(result); }

    KOKKOS_INLINE_FUNCTION
    functor_type( const Op1Type & arg_op1 )
      : task_serial< DeviceType , ResultType >()
      , m_op1( arg_op1 )
      {}
  };
};

//----------------------------------------------------------------------------

template< class DeviceType , class ResultType , class WorkType , class Op1Type >
class FunctorTraits< task_reduce< DeviceType , ResultType , WorkType >
                   , Op1Type , void , void , void >
{
  typedef task_reduce< DeviceType , ResultType , WorkType >  pattern_type ;
  typedef typename pattern_type::device_type  device_type ;
  typedef typename pattern_type::result_type  result_type ;
  typedef typename pattern_type::future_type  future_type ;
  typedef typename pattern_type::work_type    work_type ;

  class functor_type : public pattern_type {
  private:
    Op1Type  m_op1 ;
  public:

    KOKKOS_INLINE_FUNCTION
    void operator()( const work_type & work , result_type & update ) const
      { m_op1( work , update ); }

    KOKKOS_INLINE_FUNCTION
    functor_type( const WorkType & arg_work , const Op1Type & arg_op1 )
      : pattern_type( arg_work )
      , m_op1( arg_op1 )
      {}
  };
};


template< class DeviceType , class ResultType , class WorkType , class Op1Type , class Op2Type >
class FunctorTraits< task_reduce< DeviceType , ResultType , WorkType >
                   , Op1Type , Op2Type , void , void >
{
  typedef task_reduce< DeviceType , ResultType , WorkType >  pattern_type ;
  typedef typename pattern_type::device_type  device_type ;
  typedef typename pattern_type::result_type  result_type ;
  typedef typename pattern_type::future_type  future_type ;
  typedef typename pattern_type::work_type    work_type ;

  class functor_type : public pattern_type {
  private:
    Op1Type  m_op1 ;
    Op2Type  m_op2 ;
  public:

    KOKKOS_INLINE_FUNCTION
    void join( result_type volatile & update , result_type const volatile & input ) const
      { m_op2( update , input ); }

    KOKKOS_INLINE_FUNCTION
    void operator()( const work_type & work , result_type & update ) const
      { m_op1( work , update ); }

    KOKKOS_INLINE_FUNCTION
    functor_type( const WorkType & arg_work , const Op1Type & arg_op1 , const Op2Type & arg_op2 )
      : pattern_type( arg_work )
      , m_op1( arg_op1 )
      , m_op2( arg_op2 )
      {}
  };
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_FUNCTORTRAITS_HPP */

