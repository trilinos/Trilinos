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

#ifndef KOKKOS_SERIAL_TASK_HPP
#define KOKKOS_SERIAL_TASK_HPP

#include <string>
#include <typeinfo>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_Task.hpp>
#include <impl/Kokkos_IntPool.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Concrete base class for Kokkos::Serial tasks */
template<>
class Task< Kokkos::Serial , void > {
private:

  /**\brief  States of a task */
  enum { STATE_CONSTRUCTING = 0 , STATE_WAITING = 1 , STATE_EXECUTING = 2 , STATE_COMPLETE = 4 };

  /**\brief  Base dependence count when a task is allocated.
   *         A separate dependence array is allocated when the number
   *         of dependences exceeds this count.
   */
  enum { DEP_COUNT_BASE = 8 };

  typedef void (* function_type)( Task * );

public:

  std::type_info const & m_typeid ;

private:

  function_type    m_dealloc ;
  function_type    m_apply ;

  Task   * m_wait ;      ///< Linked list of tasks waiting on this task.
  Task   * m_next ;      ///< This task is a member of a linked list of
                         ///< tasks waiting on another task.
  Task  ** m_dep ;       ///< Dependences of this task
  int      m_dep_alloc ; ///< Allocation length of 'm_dep'
  int      m_dep_count ; ///< Length of 'm_dep'
  int      m_ref_count ; ///< Reference count
  int      m_state ;
  Task   * m_dep_base[ DEP_COUNT_BASE ];

  Task( const Task & );
  Task & operator = ( const Task & );

  //--------------------------------------------------------------------------

  void schedule();
  void schedule_dependence_reserve(int);

  //--------------------------------------------------------------------------

  template< class FunctorType >
  static
  void task_dealloc( Task * t )
    {
      FunctorType * f = static_cast<FunctorType*>(t);
      delete f ;
    }

  template< class Traits >
  KOKKOS_INLINE_FUNCTION
  void task_init()
    {
      typedef typename Traits::functor_type functor_type ;
      typedef typename Traits::pattern_type pattern_type ;

      m_dealloc = & task_dealloc<functor_type> ;
      m_apply   = & pattern_type::template task_apply< functor_type > ;
    }

  //--------------------------------------------------------------------------

protected:

  KOKKOS_INLINE_FUNCTION
  Task()
    : m_typeid( typeid(void) )
    , m_dealloc(0)
    , m_apply(0)
    , m_wait(0)
    , m_next(0)
    , m_dep(0)
    , m_dep_alloc(0)
    , m_dep_count(0)
    , m_ref_count(0)
    , m_state( STATE_CONSTRUCTING )
    { m_dep = m_dep_base ; }

  KOKKOS_INLINE_FUNCTION
  Task( std::type_info const & tid )
    : m_typeid(tid)
    , m_dealloc(0)
    , m_apply(0)
    , m_wait(0)
    , m_next(0)
    , m_dep(0)
    , m_dep_alloc(0)
    , m_dep_count(0)
    , m_ref_count(0)
    , m_state( STATE_CONSTRUCTING )
    { m_dep = m_dep_base ; }

public:

  ~Task();

  static void task_wait( Task * );
  static void increment( Task * );
  static void decrement( Task * );

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  Task * task_dependence( int i ) const
    { return ( STATE_EXECUTING == m_state && 0 <= i && i < m_dep_count ) ? m_dep[i] : (Task*) 0 ; }

  KOKKOS_INLINE_FUNCTION
  int task_dependence() const
    { return ( STATE_EXECUTING == m_state ) ? m_dep_count : 0 ; }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION static
  Task * verify_type( Task * t ) { return t ; }

  typedef TaskResultTypeIsVoid get_result_type ;

  get_result_type get() const { return get_result_type() ; }

  void set( get_result_type ) {}

  //--------------------------------------------------------------------------

  template< class FutureType >
  KOKKOS_INLINE_FUNCTION
  void task_respawn( FutureType const * const dep , const int count )
    {
      // Will throw if not STATE_CONSTRUCTING OR STATE_EXECUTING
      schedule_dependence_reserve( count );

      for ( int i = 0 ; i < count ; ++i ) {
        Task::increment( m_dep[i] = dep[i].m_task );
      }

      schedule();
    }

  //--------------------------------------------------------------------------

  template< class Traits , class Arg1Type >
  inline static
  Task * task_spawn( Arg1Type const & arg1 )
    {
      Task * const t = new typename Traits::functor_type(arg1);
      t->template task_init<Traits>();
      t->schedule();
      return t ;
    }

  template< class Traits
          , class Arg1Type
          , class Arg2Type
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2);
      t->template task_init<Traits>();
      t->schedule();
      return t ;
    }

  template< class Traits
          , class Arg1Type
          , class Arg2Type
          , class Arg3Type
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   , Arg3Type const & arg3
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2,arg3);
      t->template task_init<Traits>();
      t->schedule();
      return t ;
    }

  template< class Traits
          , class Arg1Type
          , class Arg2Type
          , class Arg3Type
          , class Arg4Type
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   , Arg3Type const & arg3
                   , Arg4Type const & arg4
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2,arg3,arg4);
      t->template task_init<Traits>();
      t->schedule();
      return t ;
    }

  template< class Traits
          , class Arg1Type
          , class Arg2Type
          , class Arg3Type
          , class Arg4Type
          , class Arg5Type
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   , Arg3Type const & arg3
                   , Arg4Type const & arg4
                   , Arg5Type const & arg5
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2,arg3,arg4,arg5);
      t->template task_init<Traits>();
      t->schedule();
      return t ;
    }

  template< class Traits
          , class Arg1Type
          , class FutureType
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , FutureType const * const ibegin
                   , FutureType const * const iend
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1);
      t->template init<Traits>();
      t->task_respawn( ibegin , iend - ibegin );
      return t ;
    }
 
  template< class Traits
          , class Arg1Type
          , class Arg2Type
          , class FutureType
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   , FutureType const * const ibegin
                   , FutureType const * const iend
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2);
      t->template init<Traits>();
      t->task_respawn( ibegin , iend - ibegin );
      return t ;
    }
 
  template< class Traits
          , class Arg1Type
          , class Arg2Type
          , class Arg3Type
          , class FutureType
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   , Arg3Type const & arg3
                   , FutureType const * const ibegin
                   , FutureType const * const iend
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2,arg3);
      t->template init<Traits>();
      t->task_respawn( ibegin , iend - ibegin );
      return t ;
    }
 
  template< class Traits
          , class Arg1Type
          , class Arg2Type
          , class Arg3Type
          , class Arg4Type
          , class FutureType
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   , Arg3Type const & arg3
                   , Arg4Type const & arg4
                   , FutureType const * const ibegin
                   , FutureType const * const iend
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2,arg3,arg4);
      t->template init<Traits>();
      t->task_respawn( ibegin , iend - ibegin );
      return t ;
    }
 
  template< class Traits
          , class Arg1Type
          , class Arg2Type
          , class Arg3Type
          , class Arg4Type
          , class Arg5Type
          , class FutureType
          >
  inline static
  Task * task_spawn( Arg1Type const & arg1
                   , Arg2Type const & arg2
                   , Arg3Type const & arg3
                   , Arg4Type const & arg4
                   , Arg5Type const & arg5
                   , FutureType const * const ibegin
                   , FutureType const * const iend
                   )
    {
      Task * const t = new typename Traits::functor_type(arg1,arg2,arg3,arg4,arg5);
      t->template init<Traits>();
      t->task_respawn( ibegin , iend - ibegin );
      return t ;
    }
};


template< class ResultType >
class Task< Kokkos::Serial , ResultType > : public Task< Kokkos::Serial > {
private:

  Task( const Task & );
  Task & operator = ( const Task & );

protected:

  ResultType m_result ;

  KOKKOS_INLINE_FUNCTION
  Task() : Task< Kokkos::Serial >( typeid(ResultType) ) , m_result() {}

public:

  KOKKOS_INLINE_FUNCTION static
  Task< Kokkos::Serial > * verify_type( Task< Kokkos::Serial > * t )
    {
      if ( t != 0 && t->m_typeid != typeid(ResultType) ) {
        throw std::runtime_error( std::string("Kokkos::Future bad cast for result type"));
      }
      return t ;
    }

  typedef const ResultType & get_result_type ;

  get_result_type get() const { return m_result ; }

  void set( get_result_type value ) { m_result = value ; }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template<>
class task_serial< Kokkos::Serial , void > : public Impl::Task< Kokkos::Serial > {
private:

  template< class , class > friend class Impl::Task ;

  template< class FunctorType >
  static
  void task_apply( Impl::Task< Kokkos::Serial > * t )
    {
      // Insure proper derivation chain:
      FunctorType & f = static_cast< FunctorType & >( * static_cast< task_serial * >( t ) );
      f.apply();
    }

protected:

  task_serial() : Impl::Task< Kokkos::Serial >() {}
  task_serial( const task_serial & ) : Impl::Task< Kokkos::Serial >() {}

public:

  typedef Kokkos::Serial  device_type ;
  typedef task_serial     pattern_type ;
  typedef void            result_type ;
  typedef void            work_type ;
  typedef Future< Kokkos::Serial >  future_type ;
};

template< class ResultType >
class task_serial< Kokkos::Serial , ResultType > : public Impl::Task< Kokkos::Serial , ResultType > {
private:

  template< class , class > friend class Impl::Task ;

  template< class FunctorType >
  static
  void task_apply( Impl::Task< Kokkos::Serial > * t )
    {
      // Insure proper derivation chain:
      FunctorType & f = static_cast< FunctorType & >( * static_cast< task_serial * >( t ) );
      f.apply( f.Impl::Task< Kokkos::Serial , ResultType >::m_result );
    }

protected:

  task_serial() : Impl::Task< Kokkos::Serial , ResultType >() {}
  task_serial( const task_serial & ) : Impl::Task< Kokkos::Serial , ResultType >() {}

public:

  typedef Kokkos::Serial  device_type ;
  typedef task_serial     pattern_type ;
  typedef void            result_type ;
  typedef void            work_type ;
  typedef Future< Kokkos::Serial , ResultType >  future_type ;
};

//----------------------------------------------------------------------------

template<>
class task_for< Kokkos::Serial , size_t > : public Impl::Task< Kokkos::Serial > {
private:

  template< class , class > friend class Impl::Task ;

  size_t m_count ;

  template< class FunctorType >
  static
  void task_apply( Impl::Task< Kokkos::Serial > * t )
    {
      // Insure proper derivation chain:
      FunctorType & f = static_cast< FunctorType & >( * static_cast< task_for * >( t ) );
      FunctorType const & cf = f ;

      for ( size_t i = 0 ; i < cf.task_for::m_count ; ++i ) { cf(i); }
      f.apply();
    }

protected:

  task_for( const size_t c ) : Impl::Task< Kokkos::Serial >() , m_count(c) {}
  task_for( const task_for & rhs ) : Impl::Task< Kokkos::Serial >() , m_count( rhs.m_count ) {}

public:

  KOKKOS_INLINE_FUNCTION
  void apply() {}

  typedef Kokkos::Serial  device_type ;
  typedef task_for        pattern_type ;
  typedef void            result_type ;
  typedef size_t          work_type ;
  typedef Future< Kokkos::Serial >  future_type ;
};

//----------------------------------------------------------------------------

template< class ResultType >
class task_reduce< Kokkos::Serial , ResultType , size_t >
  : public Impl::Task< Kokkos::Serial , ResultType > {
private:

  template< class , class > friend class Impl::Task ;

  size_t m_count ;

  template< class FunctorType >
  static
  void task_apply( Impl::Task< Kokkos::Serial > * t )
    {
      // Insure proper derivation chain:
      FunctorType & f = static_cast< FunctorType & >( * static_cast< task_reduce * >( t ) );
      FunctorType const & cf = f ;
      ResultType  & r = f.Impl::Task< Kokkos::Serial , ResultType >::m_result ;

      cf.init( r );
      for ( size_t i = 0 ; i < cf.task_reduce::m_count ; ++i ) { cf(i,r); }
      f.apply( r );
    }

protected:

  task_reduce( const size_t c ) : Impl::Task< Kokkos::Serial , ResultType >() , m_count(c) {}
  task_reduce( const task_reduce & rhs ) : Impl::Task< Kokkos::Serial , ResultType >() , m_count( rhs.m_count ) {}

public:

  void init( ResultType & update ) const { update = 0 ; }
  void join( ResultType & update , ResultType const & input ) const { update += input ; }
  void apply( ResultType & update ) {}

  typedef Kokkos::Serial  device_type ;
  typedef task_reduce     pattern_type ;
  typedef ResultType      result_type ;
  typedef size_t          work_type ;
  typedef Future< Kokkos::Serial , ResultType >  future_type ;
};

//----------------------------------------------------------------------------

template< class ResultType >
class task_scan< Kokkos::Serial , ResultType , size_t >
  : public Impl::Task< Kokkos::Serial , ResultType > {
private:

  template< class , class > friend class Impl::Task ;

  size_t m_count ;

  template< class FunctorType >
  static
  void task_apply( Impl::Task< Kokkos::Serial > * t )
    {
      // Insure proper derivation chain:
      FunctorType & f = static_cast< FunctorType & >( * static_cast< task_scan * >( t ) );
      FunctorType const & cf = f ;
      ResultType  & r = f.Impl::Task< Kokkos::Serial , ResultType >::m_result ;

      cf.init( r );
      for ( size_t i = 0 ; i < cf.task_scan::m_count ; ++i ) { cf(i,r,true); }
      f.apply( r );
    }

protected:

  task_scan( const size_t c ) : Impl::Task< Kokkos::Serial , ResultType >() , m_count(c) {}
  task_scan( const task_scan & rhs ) : Impl::Task< Kokkos::Serial , ResultType >() , m_count( rhs.m_count ) {}

public:

  void init( ResultType & update ) const { update = 0 ; }
  void join( ResultType & update , ResultType const & input ) const { update += input ; }
  void apply( ResultType & update ) {}

  typedef Kokkos::Serial  device_type ;
  typedef task_scan       pattern_type ;
  typedef ResultType      result_type ;
  typedef size_t          work_type ;
  typedef Future< Kokkos::Serial , ResultType >  future_type ;
};

//----------------------------------------------------------------------------

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class FunctorType , class PatternType >
class TaskPool< FunctorType , PatternType , Kokkos::Serial > {
public:

  typedef Kokkos::Serial device_type ;
  typedef Kokkos::Impl::IntPool<device_type>  int_pool_type ;

  class TaskType {
    FunctorType m_functor ;
    PatternType m_pattern ;
    int         m_ref_count ;
    int         m_state ;
    int         m_index ;

  };

  void apply( int itask ) const
    {
      m_task_pool[itask].apply();
    }

  template< class Arg1Type
          , class Arg2Type
          , class Arg3Type
          , class Arg4Type
          , class Arg5Type
          , class FutureType
          >
  KOKKOS_INLINE_FUNCTION
  typename PatternType::future_type
  spawn( Arg1Type const & arg1
       , Arg2Type const & arg2
       , Arg3Type const & arg3
       , Arg4Type const & arg4
       , Arg5Type const & arg5
       , FutureType const * const ibegin
       , FutureType const * const iend
       )
    {
      int value = 0 ;

      while ( int_pool_type::FAIL == m_task_wait.claim( value ) );

      TaskType * const t = & m_task_pool[ value ];

      new( & t->m_functor ) FunctorType(arg1,arg2,arg3,arg4,arg5);

      return typename PatternType::future_type(t);
    }


private:

  Kokkos::View<TaskType*,device_type> m_task_pool ; // Allocate without initializing
  Kokkos::Impl::IntPool<device_type>  m_task_wait ;
  Kokkos::Impl::IntPool<device_type>  m_task_depend ;

};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_SERIAL_TASK_HPP */

