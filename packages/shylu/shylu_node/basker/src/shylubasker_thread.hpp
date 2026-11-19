// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_THREAD_HPP
#define SHYLUBASKER_THREAD_HPP

#include "shylubasker_types.hpp"
//#include "shylubasker_decl.hpp"

namespace BaskerNS
{
  
  template <class Int, class Entry, class Exe_Space>
  class BaskerPointBarrier
  {
  public:
    typedef Kokkos::TeamPolicy<Exe_Space>   TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;

    //Outer idx threads, Inner idx task
    Int ** volatile token;
    Int nthreads; 
    Int ntasks;
    Int nlvls;

    //An Exit
    BASKER_BOOL * volatile exit_token;

    BASKER_BOOL init_flg;
    
    inline
    BaskerPointBarrier()
    {
      init_flg    = BASKER_FALSE;
    }

    inline
    void init(Int _nthreads, Int _tasks, Int _lvls)
    {
      //printf("THREADS BARRIER INIT, %d  \n", _nthreads);

      Finalize();

      nthreads = _nthreads;
      ntasks   = _tasks;
      nlvls    = _lvls;
      Int msize = ntasks*nlvls;

      token = new Int*[_nthreads];
      for(Int i = 0; i < _nthreads;i++)
      {
        token[i] = new Int[msize];
        for(Int j = 0; j < msize; j++)
        {
          token[i][j] = BASKER_MAX_IDX;
        }
      }

      //init exit
      exit_token = new BASKER_BOOL[_nthreads];
      for(Int i = 0; i < _nthreads; i++)
      {
        exit_token[i] = BASKER_FALSE;
      }
      init_flg = BASKER_TRUE;
    }//end BaskerPointBarrier

    ~BaskerPointBarrier()
    {
      //JDB: Need to figure out why this does not work
      Finalize();
    }//end ~BaskerPointBarrier()
    
    inline
    void Finalize()
    {
      if(init_flg == BASKER_TRUE)
      {
        for(Int i = 0; i < nthreads; ++i)
        {
          delete [] token[i];
        }
        delete [] token;
        delete [] exit_token;
      }
      init_flg = BASKER_FALSE;
    }//end Finalize()
    
    inline
    void BarrierLeader(Int my_leader, Int my_id, Int task, Int k)
    {
      //printf("Entering barrier. leader: %d id: %d task: %d k: %d token: %d \n", 
      //    my_leader, my_id, task, token[my_leader][task]);
      //jdb: change from my_lead == my_id
      if(my_leader == my_id)
      {
        token[my_leader][task] = k;
      }
      else
      {
        while(token[my_leader][task] != k);
      }
      //printf("Leaving barrier. leader: %d id: %d task: %d k: %d token: %d \n",
      //    my_leader, my_id, task, token);

    }//end Barrier_Leader()

    inline
    void BarrierDomain
    (
     Int my_leader,
     Int my_id,
     Int barrier_id,
     Int num_threads_in_team,
     Int tug,
     Int lvl,
     Int max_thread,
     Int num_threads
    )
    {
      Int task_id = (lvl*ntasks) + barrier_id;
      Int last_thread = (my_leader+num_threads_in_team)-1;
      //printf("Enter Domain Barrier. leader=%d, num_threads_in_team=%d (%d:%d), my_id=%d, barrier_id=%d, tug=%d, lvl=%d -> task_id=%d, max_thread=%d\n",
      //        my_leader, num_threads_in_team, my_leader,last_thread, my_id, barrier_id, tug, lvl, task_id, max_thread); fflush(stdout);
      #if 0 // debug
      if (token[my_id][task_id] == tug) {
        printf( "\n BarrierDomain already tug ??\n " );
        exit(0);
      }
      #endif
      if (max_thread > 0) {
        // Partial Factorization
        if (my_leader >= max_thread) {
          //right-half of tree, skip
          //printf("> Enter Domain Barrier. leader=%d, num_threads_in_team=%d (%d:%d), my_id=%d, barrier_id=%d, tug=%d, lvl=%d -> task_id=%d, max_thread=%d, skip\n",
          //        my_leader, num_threads_in_team, my_leader,last_thread, my_id, barrier_id, tug, lvl, task_id, max_thread); fflush(stdout);
          return;
        }
        if (last_thread == num_threads-1) {
          // root of tree, wait for half of threads
          last_thread = (num_threads/2)-1;
        }
      }
      //printf("my_id=%d, my_leader=%d: set token[%d][%d] = %d, start checking %d .. %d ?\n", my_id,my_leader, my_id, task_id, tug, last_thread,my_leader); fflush(stdout);
      token[my_id][task_id] = tug;
      for(Int dp = last_thread; dp >= my_leader; dp--)
      {
        //printf(" > my_id=%d: checking location token[%d][%d] = %d ?\n", my_id, dp, task_id, tug); fflush(stdout);
        ///volatile Int ldp = token[dp][task];
        while(token[dp][task_id] != tug);
        //printf("kid: %d location: %d \n", my_id, token[dp][task_id]);

      }
      //printf("Leave Domain Barrier. leader=%d, num_threads_in_team=%d (%d:%d), my_id=%d, task=%d, tug=%d, lvl=%d -> task_id=%d\n",
      //        my_leader, num_threads_in_team, my_leader,last_thread, my_id, task, tug, lvl, task_id); fflush(stdout);
    }//end Barrier_Domain

    inline
    void ExitSet(Int my_leader, BASKER_BOOL flg)
    {
      exit_token[my_leader] = flg;
    }
    
    inline
    void ExitGet(Int my_leader, BASKER_BOOL &flg)
    {
      flg = exit_token[my_leader];
    }
   
  };//end BaskerPointBarrier


  //Old Atomic Count Barrier
  template <class Int,class Entry,class Exe_Space>
  class BaskerBarrier
  {
  public:

    typedef  Kokkos::TeamPolicy<Exe_Space>   TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;

    BaskerBarrier()
    {}

    //Kokkos thread
    BaskerBarrier(TeamMember &thread)
    {
      kokkos_barrier(thread);
    }

    //Atomic
    BaskerBarrier(volatile Int &value_in, volatile Int &value_out,
                  const Int l_size )
    {
      //jdb value ->value_in
      atomic_barrier(value_in,l_size);
    }

    BASKER_INLINE
    void Barrier(TeamMember &thread)
    {
      kokkos_barrier(thread);
    }

    BASKER_INLINE
    void Barrier
    (
     volatile Int &value_in, 
     volatile Int &value_out,
     const Int l_size
    )
    {
      atomic_barrier(value_in, value_out, l_size);
    }

    BASKER_INLINE
    void Barrier
    (
     TeamMember &thread,
     volatile Int &value_in, 
     volatile Int &value_out, 
     const Int l_size
    )
    {
      if(l_size <= thread.team_size())
      {
        kokkos_barrier(thread);
      }
      else
      {
        atomic_barrier(value_in, value_out,  l_size);
      }
    }//end Barrier()
  

  private:
    BASKER_INLINE
    void kokkos_barrier(TeamMember &thread)
    {
      thread.team_barrier();
    }

    BASKER_INLINE
    void atomic_barrier
    (
     volatile Int &value_in, 
     volatile Int &value_out,
     const Int l_size
    )
    {    
      atomic_barrier_fanin(value_in, l_size);
      //atomic_barrier_fanout(value_out, l_size);
    }

    BASKER_INLINE
    void atomic_barrier_fanin(volatile Int &value, const Int l_size)
    {
      //Note: need to comeback and makesure this is the best way
      /*
      Kokkos::atomic_fetch_add(&(value), 1);
      printf("fanin value: %d \n", value);
      while(value < l_size-1)
      {
        BASKER_NO_OP;
      }
      printf("done with fanin\n");
      */

      if(Kokkos::atomic_fetch_add(&(value), Int(1)) == (l_size-1))
      {
        value = 0;
      }
      else
      {
        while(value != 0);
      }
    }

    BASKER_INLINE
    void atomic_barrier_fanout(volatile Int &value, const Int l_size)
    {
      Kokkos::atomic_inc(&(value));
      while(value < l_size)
      {
        BASKER_NO_OP;
      }
    }
  }; //end BaskerBarrier
}//end namespace Basker

#endif //end ifndef BASKER_THREADS
