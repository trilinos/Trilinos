#ifndef stk_percept_ParallelUtil_hpp
#define stk_percept_ParallelUtil_hpp

#include <utility>

#include <stk_percept/Util.hpp>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif


namespace stk { 
  namespace percept { 


    /** Copied/Modeled after MPIX.[hc]pp in ~/code/stk/src/rsmesh and is an extension to the Sierra framework's
     *    Fmwk_global_min capability.
     */

    //========================================================================================================================
    // external interface
    namespace {

      template<class T>
      inline
      void stk_percept_global_lex_min(stk::ParallelMachine comm,  int n , T local_min[] , T global_min[] );
    }

  } // percept
}// stk

//========================================================================================================================
// implementation

#include <stk_percept/ParallelUtilDef.hpp>

#endif
