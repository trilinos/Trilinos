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

#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

#include <Serial/Kokkos_Serial_Task.hpp>
#include <impl/Kokkos_IntPool.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

Task<Kokkos::Serial> * s_ready = 0 ;
Task<Kokkos::Serial> * const s_denied =
  reinterpret_cast< Task<Kokkos::Serial> * >( ~((unsigned long)(0)) );

}

/*
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
*/

Task< Kokkos::Serial , void >::~Task()
{
  if ( m_dep != 0 && m_dep != m_dep_base ) { free( m_dep ); }
  if ( m_ref_count ) {
    std::cerr << "ERROR Kokkos::Impl::Task destroyed with ref_count = " << m_ref_count << std::endl ;
  }
}

void Task< Kokkos::Serial , void >::schedule_dependence_reserve(int n)
{
  if ( STATE_CONSTRUCTING != m_state && STATE_EXECUTING != m_state ) {
    throw std::runtime_error(std::string("Kokkos::Impl::Task spawn or respawn state error"));
  }

  m_state = STATE_WAITING ;

  for ( int i = 0 ; i < m_dep_count ; ++i ) {
    decrement( m_dep[i] );
    m_dep[i] = 0 ;
  }

  if ( DEP_COUNT_BASE < n && m_dep_alloc < n ) {
    if ( m_dep_alloc ) {
      // Thread safety and scalability of 'free' ?
      free( m_dep );
    }
    m_dep_alloc = n ;
    // Thread safety and scalability of 'malloc' ?
    m_dep = reinterpret_cast<Task**>( malloc( n * sizeof(Task*) ) );

    for ( int i = 0 ; i < m_dep_alloc ; ++i ) { m_dep[i] = 0 ; }
  }

  m_dep_count = n ;
}

void Task< Kokkos::Serial , void >::increment( Task< Kokkos::Serial , void > * t )
{ if ( t ) { ++(t->m_ref_count); } }

void Task< Kokkos::Serial , void >::decrement( Task< Kokkos::Serial , void > * t )
{
  if ( t && 0 == --(t->m_ref_count) ) {

    while ( t->m_dep_count ) {
      const int i = --( t->m_dep_count );
      Task * const td = t->m_dep[i];
      t->m_dep[i] = 0 ;
      decrement( td );
    }

    t->m_dep = 0 ;
    t->m_dep_alloc = 0 ;

    function_type d = t->m_dealloc ;

    (*d)( t );
  }
}

void Task< Kokkos::Serial , void >::schedule()
{
  if ( 0 != m_next ) {
    throw std::runtime_error( std::string("ERROR in Kokkos::Impl::Task::schedule") );
  }

  // Insert this task into another dependence that is not complete

  int i = 0 ;
  for ( ; i < m_dep_count ; ++i ) {
    Task * const y = m_dep[i] ;
    if ( s_denied != ( m_next = y->m_wait ) ) {
      y->m_wait = this ; // CAS( & y->m_wait , m_next , this );
      break ;
    }
  }
  if ( i == m_dep_count ) {
    // All dependences are complete, insert into the ready list
    m_next = s_ready ;
    s_ready = this ; // CAS( & s_ready , m_next = s_ready , this );
  }
}

void Task< Kokkos::Serial , void >::task_wait( Task< Kokkos::Serial , void > * )
{
  while ( s_ready ) {

    // Remove this task from the ready list

    // Task * task ;
    // while ( ! CAS( & s_ready , task = s_ready , s_ready->m_next ) );

    Task * const task = s_ready ;
    s_ready = task->m_next ;

    task->m_next = 0 ;

    // precondition: task->m_state = STATE_WAITING
    // precondition: task->m_dep[i]->m_state == STATE_COMPLETE  for all i
    // precondition: does not exist T such that T->m_wait = task
    // precondition: does not exist T such that T->m_next = task

    task->m_state = STATE_EXECUTING ;

    (*task->m_apply)( task );

    if ( task->m_state == STATE_EXECUTING ) {
      // task did not respawn itself
      task->m_state = STATE_COMPLETE ;

      // release dependences:
      while ( task->m_dep_count ) {
        const int i = --( task->m_dep_count );
        decrement( task->m_dep[i] );
        task->m_dep[i] = 0 ;
      }

      // Stop other tasks from adding themselves to 'task->m_wait' ;

      Task * x ;
      // CAS( & task->m_wait , x = task->m_wait , s_denied );
      x = task->m_wait ; task->m_wait = (Task*) s_denied ;

      // update tasks waiting on this task
      while ( x ) {
        Task * const next = x->m_next ;

        x->m_next = 0 ;

        x->schedule();

        x = next ;
      }
    }
  }
}

} // namespace Impl
} // namespace Kokkos

