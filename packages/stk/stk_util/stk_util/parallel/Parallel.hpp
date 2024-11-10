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

#ifndef stk_util_parallel_Parallel_hpp
#define stk_util_parallel_Parallel_hpp

// stk_config.h contains the #define macros that define
// which build-dependent features are enabled, such as 'STK_HAS_MPI'.

#include "stk_util/stk_config.h"

//----------------------------------------------------------------------
// Parallel machine

#if defined( STK_HAS_MPI )

#include "mpi.h"

namespace stk {

/** \addtogroup parallel_module
 * @{
 */

/// \todo REFACTOR: Figure out a better way to typedef for non-MPI builds

typedef MPI_Comm     ParallelMachine ;

/// \todo REFACTOR: Figure out a better way to typedef for non-MPI builds

typedef MPI_Datatype ParallelDatatype ;

/**
 * @brief <b>parallel_machine_null</b> returns MPI_COMM_NULL if MPI is enabled.
 *
 * @return			a <b>ParallelMachine</b> ...
 */
inline ParallelMachine parallel_machine_null() { return MPI_COMM_NULL ; }

/**
 * @brief <b>parallel_machine_init</b> calls MPI_Init.
 *
 * @return <b>ParallelMachine</b> (MPI_COMM_WORLD)
 */
inline ParallelMachine parallel_machine_init( int * argc , char *** argv )
{
  MPI_Init( argc , argv );
  return MPI_COMM_WORLD ; // CHECK: ALLOW MPI_COMM_WORLD
}

/**
 * @brief <b>parallel_machine_finalize</b> calls MPI_Finalize.
 *
 */
inline void parallel_machine_finalize()
{
  MPI_Finalize();
}

/** \} */

}

//----------------------------------------
// Other parallel communication machines go here
// as '#elif defined( STK_HAS_<name> )'

#else

//----------------------------------------
// Stub for non-parallel, no MPI

// Some needed stubs
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#define MPI_COMM_SELF 0
#define MPI_Barrier( a ) (static_cast<void>(a))

namespace stk {

typedef int ParallelMachine ;
typedef int ParallelDatatype ;

inline ParallelMachine parallel_machine_null() { return 0 ; }

inline ParallelMachine parallel_machine_init( int * , char *** )
{ return 0 ; }

inline void parallel_machine_finalize()
{}

}

#endif

//----------------------------------------------------------------------
// Common parallel machine needs.

namespace stk {

inline ParallelMachine parallel_machine_world()
{
  return MPI_COMM_WORLD; // CHECK: ALLOW MPI_COMM_WORLD
}

/**
 * @brief Member function <b>parallel_machine_size</b> ...
 *
 * @param parallel_machine	a <b>ParallelMachine</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
int parallel_machine_size( ParallelMachine parallel_machine );

/**
 * @brief Member function <b>parallel_machine_rank</b> ...
 *
 * @param parallel_machine	a <b>ParallelMachine</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
int parallel_machine_rank( ParallelMachine parallel_machine );

/**
 * @brief Member function <b>parallel_machine_barrier</b> ...
 */
void parallel_machine_barrier( ParallelMachine parallel_machine);

//----------------------------------------------------------------------

class Parallel {
//keeps a copy of the machine (really MPI_Comm), and the rank & size for fast access
//right now used in BulkData but could be used elsewhere
public: //methods
  Parallel(ParallelMachine parallel_in) :
    m_parallel_machine( parallel_in ),
    m_parallel_size( parallel_machine_size( parallel_in ) ),
    m_parallel_rank( parallel_machine_rank( parallel_in ) )
  {
    //no-op constructor right now
  }

  /** \brief  The parallel machine */
  inline ParallelMachine parallel() const { return m_parallel_machine ; }

  /** \brief  Size of the parallel machine */
  inline int parallel_size()   const { return m_parallel_size ; }

  /** \brief  Rank of the parallel machine's local processor */
  inline int parallel_rank()   const { return m_parallel_rank ; }

  ~Parallel() {}  //do nothing destructor

private: //data
  ParallelMachine m_parallel_machine;
  int m_parallel_size;
  int m_parallel_rank;

}; //class Parallel

} //namespace stk


#endif

