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
    #ifdef BASKER_KOKKOS
    typedef Kokkos::TeamPolicy<Exe_Space>   TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif
    
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
     Int task, 
     Int lsize, 
     Int k, 
     Int l
    )
    {
      Int ltask = (l*ntasks) + task;
      //printf("Enter Domain Barrier. leader=%d, lsize=%d (%d:%d), my_id=%d, task=%d, k=%d, l=%d -> ltask=%d\n",
      //        my_leader, lsize, my_leader,my_leader+lsize-1, my_id, task, k, l, ltask); fflush(stdout);
      token[my_id][ltask] = k;
      for(Int dp = (my_leader+lsize)-1; dp >= my_leader; dp--)
      {
        //printf("kid=%d: checking location token[%d][%d] = %d ?\n", my_id, dp, ltask, k);
        ///volatile Int ldp = token[dp][task];
        while(token[dp][ltask] != k);
        //{
        // printf("kid: %d location: %d \n", my_id, token[dp][ltask]);
        //}

      }
      //printf("Leave Domain Barrier. leader=%d, lsize=%d (%d:%d), my_id=%d, task=%d, k=%d, l=%d -> ltask=%d\n",
      //        my_leader, lsize, my_leader,my_leader+lsize-1, my_id, task, k, l, ltask); fflush(stdout);
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

    #ifdef BASKER_KOKKOS
    typedef  Kokkos::TeamPolicy<Exe_Space>   TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

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

      volatile BASKER_BOOL spin = BASKER_TRUE;

      if(Kokkos::atomic_fetch_add(&(value), Int(1)) == (l_size-1))
      {
        value = 0;
        //Kokkos::atomic_fetch_add(&(value), -1*l_size);
        //spin = BASKER_FALSE;
      }
      else
      {
        while(value != 0);
      }
    }

    BASKER_INLINE
    void atomic_barrier_fanout(volatile Int &value, const Int l_size)
    {
      Kokkos::atomic_increment(&(value));
      while(value < l_size)
      {
        BASKER_NO_OP;
      }
    }
			       
  }; //end BaskerBarrier


  /*   First attempt
  template <Int, Entry, Exe_Space>
  class BaskerBarrier
  {
  public:

#ifdef BASKER_KOKKOS
    typedef  Kokkos::TeampPolicy<Exe_Space>   TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
#endif

    BaskerBarrier()
    {}
    //Kokkos thread
    BaskerBarrier(TeamMember &thread)
    {
      kokkos_barrier(thread);
    }
    //Atomic
    BaskerBarrier(volatile Int &value, const Int l_size )
    {
      atomic_barrier(value,l_size);
    }
    BaskerBarrier(TeamMember &thread, 
        volatile Int &value, const Int l_size)
    {


    }
    BASKER_INLINE
      void Barrier(TeamMember &thread)
      {
        kokkos_barrier(thread);
      }
    BASKER_INLINE
      void Barrier(volatile Int &value, const Int l_size)
      {
        atomic_barrier(value, l_size);
      }
    BASKER_INLINE
      void Barrier(TeamMember &thread,
          volatile Int &value, const Int l_size)
      {
        if(l_size <= thread.team_size())
        {
          kokkos_barrier(thread);
        }
        else
        {
          atomic_barrier(value, l_size);
        }
      }//end Barrier()


  private:
    BASKER_INLINE
      void kokkos_barrier(TeamMember &thread)
      {
        thread.team_barrier();
      }
    BASKER_INLINE
      void 

      BASKER_INLINE
      void atomic_barrier(volatile Int &value, const Int l_size)
      {
        //Note: need to comeback and makesure this is the best way
        Kokkos::atomic_fetch_add(&(value), 1);
        while(value < l_size)
        {
          BASKER_NO_OP;
        }
      }
  }; //end BaskerBarrier
  */


}//end namespace Basker

#endif //end ifndef BASKER_THREADS
