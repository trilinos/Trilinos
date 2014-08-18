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

#include <impl/Kokkos_Serial_TaskPolicy.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

typedef TaskMember<  Kokkos::Serial > Task ;
typedef TaskManager< Kokkos::Serial > Mgr ;

Mgr s_task_manager ;

Mgr::TaskManager()
  : m_ready(0)
  , m_denied( reinterpret_cast<Task*>( ~((unsigned long)0) ) )
{}

void Mgr::decrement( Task * t )
{
  if ( t && 0 == --(t->m_ref_count) ) {
    // Reference count at zero, delete it

    // Remove dependences:
    for ( int i = 0 ; i < MAX_DEPENDENCE ; ++i ) {
      Task * const td = t->m_dep[i];
      if ( td ) {
        t->m_dep[i] = 0 ;
        decrement( td );
      }
    }

    // Get deletion function and apply it
    Task::function_type d = t->m_dealloc ;

    (*d)( t );
  }
}

void Mgr::increment( Task * t )
{ if ( t ) { ++(t->m_ref_count); } }


void Mgr::set_dependence( Task * t , Task ** d )
{
  // Must be either constructing for original spawn or executing for a respawn.

  if ( Task::STATE_CONSTRUCTING != t->m_state &&
       Task::STATE_EXECUTING    != t->m_state ) {
    throw std::runtime_error(std::string("Kokkos::Impl::Task spawn or respawn state error"));
  }

  // Remove old dependences, for a respawn:
  for ( int i = 0 ; i < MAX_DEPENDENCE ; ++i ) {
    decrement( t->m_dep[i] );
    t->m_dep[i] = 0 ;
  }

  if ( d ) { // Assign new dependences:
    int i = 0 ;
    for ( ; i < MAX_DEPENDENCE && d[i] ; ++i ) { increment( t->m_dep[i] = d[i] ); }
    for ( ; i < MAX_DEPENDENCE ; ++i ) { t->m_dep[i] = 0 ; }
  }
}

void Mgr::schedule( Task * t )
{
  // Must not be in a dependence linked list:  0 == t->m_next

  if ( 0 != t->m_next ) {
    throw std::runtime_error(std::string("Kokkos::Impl::Task spawn or respawn state error"));
  }

  // Is waiting for execution

  t->m_state = Task::STATE_WAITING ;

  // Insert this task into another dependence that is not complete

  int i = 0 ;
  for ( ; i < MAX_DEPENDENCE ; ++i ) {
    Task * const y = t->m_dep[i] ;
    if ( y && m_denied != ( t->m_next = y->m_wait ) ) {
      y->m_wait = t ; // CAS( & y->m_wait , m_next , this );
      break ;
    }
  }
  if ( i == MAX_DEPENDENCE ) {
    // All dependences are complete, insert into the ready list
    t->m_next = m_ready ;
    m_ready = t ; // CAS( & s_ready , m_next = s_ready , this );
  }
}

void Mgr::wait()
{
  while ( m_ready ) {

    // Remove this task from the ready list

    // Task * task ;
    // while ( ! CAS( & s_ready , task = s_ready , s_ready->m_next ) );

    Task * const task = m_ready ;
    m_ready = task->m_next ;

    task->m_next = 0 ;

    // precondition: task->m_state = STATE_WAITING
    // precondition: task->m_dep[i]->m_state == STATE_COMPLETE  for all i
    // precondition: does not exist T such that T->m_wait = task
    // precondition: does not exist T such that T->m_next = task

    task->m_state = Task::STATE_EXECUTING ;

    (*task->m_apply)( task );

    if ( task->m_state == Task::STATE_EXECUTING ) {
      // task did not respawn itself
      task->m_state = Task::STATE_COMPLETE ;

      // release dependences:
      for ( int i = 0 ; i < MAX_DEPENDENCE ; ++i ) {
        decrement( task->m_dep[i] );
        task->m_dep[i] = 0 ;
      }

      // Stop other tasks from adding themselves to 'task->m_wait' ;

      Task * x ;
      // CAS( & task->m_wait , x = task->m_wait , s_denied );
      x = task->m_wait ; task->m_wait = (Task*) m_denied ;

      // update tasks waiting on this task
      while ( x ) {
        Task * const next = x->m_next ;

        x->m_next = 0 ;

        schedule( x );

        x = next ;
      }
    }
  }
}

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

TaskPolicy< Kokkos::Serial >::TaskPolicy()
  : m_task_manager( Impl::s_task_manager )
{}

} // namespace Kokkos

