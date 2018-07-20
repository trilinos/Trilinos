#ifndef SHYLUBASKER_UTIL_HPP
#define SHYLUBASKER_UTIL_HPP

/*Basker Includes*/
#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_types.hpp"

/*Kokkos Includes*/
#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif

/*System Includes*/
#include <iostream>
#include <string>
#include <fstream>

#include <assert.h>

using namespace std;

namespace BaskerNS
{

  //Kokkos struct for init 2D Structure of A
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_order_init_2D
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                    execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    #endif

    Basker<Int,Entry,Exe_Space>  *basker;
    BASKER_BOOL                  alloc;


    kokkos_order_init_2D()
    {}
    
    kokkos_order_init_2D(Basker<Int,Entry,Exe_Space> *_b)
    {
      basker = _b;
      alloc  = BASKER_TRUE;
    }//end kokkos_order_init_2D()

    kokkos_order_init_2D(Basker<Int,Entry,Exe_Space> *_b, BASKER_BOOL _alloc)
    {
      basker = _b;
      alloc  = _alloc;
    }
      

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    void operator()(const TeamMember &thread) const
    #else
    void operator()(Int kid) const
    #endif
    {
      #ifdef BASKER_KOKKOS
      Int kid = (Int)(thread.league_rank()*thread.team_size() + thread.team_rank());
      #endif
      {
        basker->t_init_2DA(kid, alloc);
      }
    }//end operator()

  };//end kokkos_order_init_2D
  
  //Kokkos struct for reinit for refactor
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_reset_factor
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                    execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    #endif

    Basker<Int,Entry,Exe_Space>  *basker;

    kokkos_reset_factor()
    {}
    
    kokkos_reset_factor(Basker<Int,Entry,Exe_Space> *_b)
    {
      basker = _b;
    }

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    void operator()(const TeamMember &thread) const
    #else
    void operator()(Int kid) const
    #endif
    {
      #ifdef BASKER_KOKKOS
      Int kid = (Int)(thread.league_rank()*thread.team_size() + thread.team_rank());
      #endif
      {
        if(basker->btf_tabs_offset!=0)
        {
          basker->t_reset_ND_factor(kid);
        }
        if(basker->btf_nblks > 1)
        {
          basker->t_reset_BTF_factor(kid);
        }
      }
    }//end operator()
  };//end kokkos_reset_factor


  //--------------------MEMORY RELATED UTIL ------------------//
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::init_value
  (
   INT_1DARRAY a, 
   Int size, 
   Int c
  )
  {
    for(Int i=0; i < size; i++)
    {
      a(i) = c;
    }
  }//end init_value 1d array host


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::init_value
  (
   INT_1DARRAY a, 
   Int size, 
   Int* c
  )
  {
    //Note, we will want to change this to a memcpy for the type
    for(Int i=0; i < size; i++)
    {
      a(i) = c[i];
    }
  }//end init_value 1d array host 


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::init_value
  (
   ENTRY_1DARRAY a, 
   Int size, 
   Entry c
  )
  {
    for(Int i=0; i < size; i++)
    {
      a(i) = c;
    }
  }//end init_value 1d array host


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::init_value
  (
   ENTRY_1DARRAY a, 
   Int size, Entry *c
  )
  {
    //Note: come back and make memcpy for this type
    for(Int i=0; i <size; i++)
    {
      a(i) = c[i];
    }
  }//end init_value 1d array host


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::init_value
  (
   BOOL_1DARRAY a, 
   Int size, BASKER_BOOL c
  )
  {
    for(Int i=0; i<size; i++)
    {
      a[i] = c;
    }
  }//end init_value 1d array host


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::init_value
  (
   BOOL_1DARRAY a,
   Int size, BASKER_BOOL* c
  )
  {
    for(Int i =0; i < size; i++)
    {
      a(i) = c[i];
    }
  }//end init_value 1d array host


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::init_value
  (
   INT_1DARRAY a, 
   Int size, Int c, Int kid
  )
  {
    printf("\n===SHOULD NOT BE CALLED\n");
    BASKER_ASSERT(0==1, "init_int_thread");

    #ifdef BASKER_KOKKOS
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    Kokkos::parallel_for(
			 TeamPolicy(Exe_Space::thread_pool_size(),1),
			 KOKKOS_LAMBDA(const TeamMember& thread)
    #else
    #pragma omp parallel
    #endif			 
    {
      #ifdef BASKER_KOKKOS
      if(kid == thread.league_rank())
      #else
      if(kid == omp_get_thread_num())
      #endif
      {
        for(Int i=0; i < size; i++)
        {
          a(i) = c;
        }
      }
    }
    #ifdef BASKER_KOKKOS
    );
    #endif
  }//end init_value int 1d 

  
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::init_value
  (
   ENTRY_1DARRAY a, 
   Int size, 
   Entry c, 
   Int kid
  )
  {
    printf("\n===SHOULD NOT BE CALLED===\n");
    BASKER_ASSERT(0==1, "INIT_VALUE_ENTRY_THREADS");

    #ifdef BASKER_KOKKOS
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    Kokkos::parallel_for(
			 TeamPolicy(Exe_Space::thread_pool_size(),1),
			 KOKKOS_LAMBDA(const TeamMember& thread)
    #else
    #pragma omp parallel
    #endif			 
    {
      #ifdef BASKER_KOKKOS
      if(kid == thread.league_rank())
      #else
      if(kid == omp_get_thread_num())
      #endif
      {
        for(Int i=0; i < size; i++)
        {
          a(i) = c;
        }
      }
    }
    #ifdef BASKER_KOKKOS
    );
    #endif
  }//end init_value entry 1d 


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_reset_BTF_factor(Int kid)
  {
    Int chunk_start = btf_schedule(kid);
    Int chunk_end   = btf_schedule(kid+1);
    Int chunk_size = chunk_end-chunk_start;
    
    if(chunk_size > 0)
    {
      for(Int b=chunk_start; b < chunk_end; b++)
      {
        BASKER_MATRIX &L = LBTF(b-btf_tabs_offset);
        L.clear_pend();
        L.nnz = L.mnnz;
      }//end-for over chunck
    }
  }//end t_reset_BTF_factor


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_reset_ND_factor(Int kid)
  {
    //L
    for(Int lvl = 0; lvl < tree.nlvls+1; lvl++)
    {
      if(kid%((Int)pow(2,lvl)) == 0)
      {
        Int b = S(lvl)(kid);

        for(Int row = 0; row < LL_size(b); row++)
        {
          #ifdef BASKER_DEBUG_INIT
          printf("L Factor Init: %d %d , kid: %d, nnz: %ld \n",
              b, row, kid, LL[b][row].nnz);
          #endif

          LL(b)(row).clear_pend();
          LL(b)(row).nnz = LL(b)(row).mnnz;

        }//end over all row
      }//end select which thread
    }//end for over all lvl

    for(Int lvl = 0; lvl < tree.nlvls+1; lvl++)
    {
      if(kid%((Int)pow(2,lvl)) == 0)
      {
        Int b = S(lvl)(kid);

        #ifdef BASKER_DEBUG_INIT
        printf("U Factor init: %d %d, nnz: %ld \n",
            b, LU_size[b]-1, 
            LU[b][LU_size[b]-1].nnz);
        #endif

        //LU(b)(LU_size(b)-1).nnz = 0;
        for(Int kk = 0; kk < LU(b)(LU_size(b)-1).ncol+1; kk++)
        {
          LU(b)(LU_size(b)-1).col_ptr(kk) = 0;
        }

        /*
           printf("flipU1 (%d,%d) %d %d \n",
           b, LU_size(b)-1,
           LU(b)(LU_size(b)-1).nnz,
           LU(b)(LU_size(b)-1).mnnz);
        */

        LU(b)(LU_size(b)-1).nnz = LU(b)(LU_size(b)-1).mnnz;

        for(Int l = lvl+1; l < tree.nlvls+1; l++)
        {
          Int U_col = S(l)(kid);

          Int my_row_leader = find_leader(kid, l-1);
          Int my_new_row = 
            b - S(0)(my_row_leader);

          Int U_row = (l==1)?(kid%2):S(lvl)(kid)%LU_size(U_col);

          //JDB TEST PASS
          U_row = my_new_row;

          #ifdef BASKER_DEBUG_INIT
          printf("Init U: %d %d lvl: %d l: %d kid: %d nnz: %ld \n",
              U_col, U_row, lvl, l, kid, 
              LU[U_col][U_row].nnz);
          #endif

          for(Int kk = 0; kk < LU(U_col)(U_row).ncol+1; kk++)
          {
            LU(U_col)(U_row).col_ptr(kk) = 0;
          }
          /*
             printf("flipU (%d,%d) %d %d \n",
             U_col, U_col,
             LU(U_col)(U_row).nnz,
             LU(U_col)(U_row).mnnz);
          */

          LU(U_col)(U_row).nnz = LU(U_col)(U_row).mnnz;

          //LU(U_col)(U_row).nnz = 0;
        }//over inner lvls
      }//if KID

    }//end over all lvls

  }//end t_reset_factor


  //--------------------  Workspace -------------------------//
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_init_factor(Int kid)
  {
    //L
    for(Int lvl = 0; lvl < tree.nlvls+1; lvl++)
    {
      if(kid%((Int)pow(2,lvl)) == 0)
      {
        Int b = S(lvl)(kid);

        for(Int row = 0; row < LL_size(b); row++)
        {
          #ifdef BASKER_DEBUG_INIT
          printf("L Factor Init: %d %d , kid: %d, nnz: %ld \n",
              b, row, kid, LL[b][row].nnz);
          #endif

          LL(b)(row).init_matrix("Loffdig",
              LL(b)(row).srow,
              LL(b)(row).nrow,
              LL(b)(row).scol,
              LL(b)(row).ncol,
              LL(b)(row).nnz);

          //Fix when this all happens in the future
          if(Options.incomplete == BASKER_TRUE)
          {
            LL(b)(row).init_inc_lvl();
          }
          LL(b)(row).fill();
          LL(b)(row).init_pend();

        }//end over all row
      }//end select which thread
    }//end for over all lvl
   
    //U
    for(Int lvl = 0; lvl < tree.nlvls+1; lvl++)
    {
      if(kid%((Int)pow(2,lvl)) == 0)
      {
        Int b = S(lvl)(kid);

        #ifdef BASKER_DEBUG_INIT
        printf("U Factor init: %d %d, nnz: %ld \n",
            b, LU_size[b]-1, 
            LU[b][LU_size[b]-1].nnz);
        #endif

        LU(b)(LU_size(b)-1).init_matrix("Udiag",
            LU(b)(LU_size(b)-1).srow,
            LU(b)(LU_size(b)-1).nrow,
            LU(b)(LU_size(b)-1).scol,
            LU(b)(LU_size(b)-1).ncol,
            LU(b)(LU_size(b)-1).nnz);

        LU(b)(LU_size(b)-1).fill();

        for(Int l = lvl+1; l < tree.nlvls+1; l++)
        {
          Int U_col = S(l)(kid);

          Int my_row_leader = find_leader(kid, l-1);
          Int my_new_row = 
            b - S(0)(my_row_leader);

          Int U_row = (l==1)?(kid%2):S(lvl)(kid)%LU_size(U_col);

          if( (b > 14) &&  // NDE: Why is 14 specifically used here?
              (b  > LU_size(U_col)) &&
              (l !=1)
            )
          {
            /*
               printf("ttp2 %d %d %d \n",
               kid,
               S(l)(kid), 
               LU_size(U_col));
            */

            Int tm = (b+1)/16;
            U_row = ((b+1)-(tm*16))%LU_size(U_col);
          }

          //printf("Init U kid: %d U: %d %d new: %d leader: %d %d \n",
          //  kid, U_col, U_row, my_new_row, S(0)(my_row_leader), b);

          //JDB TEST PASS
          U_row = my_new_row;

          #ifdef BASKER_DEBUG_INIT
          printf("Init U: %d %d lvl: %d l: %d kid: %d nnz: %ld \n",
              U_col, U_row, lvl, l, kid, 
              LU[U_col][U_row].nnz);
          #endif

          LU(U_col)(U_row).init_matrix("Uoffdiag",
              LU(U_col)(U_row).srow,
              LU(U_col)(U_row).nrow,
              LU(U_col)(U_row).scol, 
              LU(U_col)(U_row).ncol,
              LU(U_col)(U_row).nnz);

          LU(U_col)(U_row).fill();

          if(Options.incomplete == BASKER_TRUE)
          {
            LU(U_col)(U_row).init_inc_lvl();
          }

        }//over inner lvls
      }//if KID

    }//end over all lvls

  }//end t_init_factor()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::t_init_2DA
  (
   Int kid, BASKER_BOOL alloc
  )
  {
    //L
    for(Int lvl = 0; lvl < tree.nlvls+1; lvl++)
    {
      if(kid%((Int)pow(2,lvl)) == 0)
      {
        Int b = S(lvl)(kid);

        for(Int row = 0; row < LL_size(b); row++)
        {
          #ifdef BASKER_DEBUG_INIT
          printf("ALM Factor Init: %d %d , kid: %d, nnz: %d nrow: %d ncol: %d \n",
              b, row, kid, ALM(b)(row).nnz, 
              ALM(b)(row).nrow, 
              ALM(b)(row).ncol);
          #endif

          if(Options.btf == BASKER_FALSE)
          {
            ALM(b)(row).convert2D(A, alloc, kid);
          }
          else
          {
            //printf("Using BTF AL \n");
            #ifdef BASKER_DEBUG_INIT
            printf("ALM alloc: %d %d kid: %d \n",
                b, row, kid);
            #endif

            ALM(b)(row).convert2D(BTF_A, alloc, kid);
          }

        }//end over all row
      }//end select which thread
    }//end for over all lvl
   
    //U
    for(Int lvl = 0; lvl < tree.nlvls+1; lvl++)
    {
      if(kid%((Int)pow(2,lvl)) == 0)
      {
        Int b = S(lvl)(kid);

        #ifdef BASKER_DEBUG_INTI
        printf("AUM Factor init: %d %d, kid: %d nnz: %d nrow: %d ncol: %d \n",
            b, LU_size(b)-1, kid, 
            AVM(b)(LU_size(b)-1).nnz, 
            AVM(b)(LU_size(b)-1).nrow,
            AVM(b)(LU_size(b)-1).ncol);
        #endif

        if(Options.btf == BASKER_FALSE)
        {
          AVM(b)(LU_size(b)-1).convert2D(A, alloc, kid);
        }
        else
        {
          //printf("Using BTF AU\n");
          //printf("convert AVM: %d %d kid: %d  \n", 
          // b, LU_size(b)-1, kid);
          AVM(b)(LU_size(b)-1).convert2D(BTF_A, alloc, kid);
        }

        for(Int l = lvl+1; l < tree.nlvls+1; l++)
        {
          //MOVE LEFT TO RIGHT, FIX G-ROW

          //TEST
          Int my_leader = find_leader(kid,l-1);
          Int my_leader_row = S(0)(my_leader);
          //Int my_col_size  = pow(2,l); Not used
          Int my_new_row  = 
            (S(lvl)(kid) - my_leader_row);
          //my_new_row = my_new_row%my_col_size;

          /*
             printf("TEST lvl: %d l: %d leader: %d leader_r: %d my: %d col_size: %d new_row: %d \n",
             lvl, l,
             my_leader, my_leader_row, 
             S(lvl)(kid),
             my_col_size, my_new_row);
          */

          Int U_col = S(l)(kid);
          Int U_row = my_new_row;

          //Int U_row = (l==1)?(kid%2):S(lvl)(kid)%LU_size(U_col);
          //printf("U_col: %d U_row: %d lvl: %d l: %d \n",
          //   U_col, U_row, lvl, l);
          /*
             if( (S(lvl)(kid) > 14)&&
                 (S(lvl)(kid) > LU_size(U_col)) &&
                 (l!=1)
               )
             {
             printf("test point: %d %d %d  \n",
             S(lvl)(kid), LU_size(U_col),

             S(lvl)(kid)/16);

             Int tm = (S(lvl)(kid)+1)/16;
             U_row = ((S(lvl)(kid)+1)-(tm*16))%LU_size(U_col);

             }
          */

          #ifdef BASKER_DEBUG_INIT
          printf("Init AUM: %d %d lvl: %d l: %d kid: %d nnz: %d nrow: %d ncol: %d \n",
              U_col, U_row, lvl, l, kid, 
              AVM(U_col)(U_row).nnz, 
              AVM(U_col)(U_row).nrow,
              AVM(U_col)(U_row).ncol);
          #endif

          if(Options.btf == BASKER_FALSE)
          {
            BASKER_ASSERT(0==1, "SHOULD NOTH BE CALL\n");
            //AVM(U_col)(U_row).convert2D(A);
          }
          else
          {
            //printf("Using BTF AU\n");
            //printf("2nd convert AVM: %d %d size:%d kid: %d\n",
            //	   U_col, U_row, AVM(U_col)(U_row).nnz, 
            //	   kid);

            AVM(U_col)(U_row).convert2D(BTF_A, alloc, kid);
          }

        }//over inner lvls
      }//if KID

    }//end over all lvls

  }//end t_init_2DA

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_init_workspace(Int kid)
  {
    Int max_sep_size = 0;

    if(btf_tabs_offset != 0)
    {
      #ifdef BASKER_2DL
      Int            b  = S(0)(kid);
      //INT_1DARRAY    ws = LL[b][0].iws;
      //ENTRY_1DARRAY  X  = LL[b][0].ews;
      //Int      iws_size = LL[b][0].iws_size;
      //Int      ews_size = LL[b][0].ews_size;
      //Int      iws_mult = LL[b][0].iws_mult;
      //Int      ews_mult = LL[b][0].ews_mult;

      #else
      INT_1DARRAY  &ws = thread_array(kid).iws;
      ENTRY_1DARRAY &X = thread_array(kid).ews;
      Int iws_size     = thread_array(kid).iws_size;
      Int iws_mult     = thread_array(kid).iws_mult;
      Int ews_size     = thread_array(kid).ews_size;
      Int ews_mult     = thread_array(kid).ews_mult;
      #endif
      //Note: need to add a size array for all these

      #ifdef BASKER_2DL
      for(Int l = 0; l < LL_size(b); l++)
      {
        //defining here
        LL(b)(l).iws_size = LL(b)(l).nrow;
        //This can be made smaller, see notes in Sfactor_old
        LL(b)(l).iws_mult = 5;
        LL(b)(l).ews_size = LL(b)(l).nrow;
        //This can be made smaller, see notes in sfactor_old
        LL(b)(l).ews_mult = 2;

        Int iws_size = LL(b)(l).iws_size;
        Int iws_mult = LL(b)(l).iws_mult;
        Int ews_size = LL(b)(l).ews_size;
        Int ews_mult = LL(b)(l).ews_mult;

        if(iws_size > max_sep_size)
        {
          max_sep_size = iws_size;
        }

        if(iws_size == 0)
        {
          iws_size  = 1;
        }

        BASKER_ASSERT((iws_size*iws_mult)>0, "util iws");
        MALLOC_INT_1DARRAY(LL(b)(l).iws, iws_size*iws_mult);

        //TEST
        INT_1DARRAY att = LL(b)(l).iws; 
        if(ews_size == 0)
        {
          ews_size = 1;
        }

        BASKER_ASSERT((ews_size*ews_mult)>0, "util ews");
        MALLOC_ENTRY_1DARRAY(LL(b)(l).ews, ews_size*ews_mult);

        for(Int i=0; i<iws_mult*iws_size; i++)
        {
          LL(b)(l).iws(i) = 0;
        }

        for(Int i=0; i<ews_mult*ews_size; i++)
        {
          LL(b)(l).ews(i) = 0;
        }

        LL(b)(l).fill();

        if(l==0)
        {
          //Also workspace matrix 
          //This could be made smaller
          //printf("C: size: %d kid: %d \n",
          //	   iws_size, kid);

          //thread_array[kid].C.init_matrix("cwork", 
          //			     0, iws_size,
          //			     0, 2, 
          //			     iws_size*2);
        }
      } //end for l

      //Also workspace matrix 
      //This could be made smaller
      thread_array(kid).C.init_matrix("cwork", 0, max_sep_size,
          0, 2, max_sep_size*2);

    } //end if btf_tabs_offset != 0
    else
    {
      if(btf_nblks > 1)
      {
        if(Options.btf == BASKER_TRUE)
        {
          Int iws_mult = thread_array[kid].iws_mult;
          Int iws_size = thread_array[kid].iws_size;
          Int ews_mult = thread_array[kid].ews_mult;
          Int ews_size = thread_array[kid].ews_size;

          for(Int i=0; i < iws_mult*iws_size; i++)
          {
            thread_array[kid].iws[i] = 0;
          }

          for(Int i = 0; i < ews_mult*ews_size; i++)
          {
            thread_array[kid].ews[i] = 0.0;
          }
        }
      }

    }//else
   
    #else //ifdef basker_2dl
    printf("init_workspace 1d, kid: %d size: %d %d %d %d \n",
	    kid, iws_mult, iws_size, ews_mult, ews_size);
    for(Int i=0; i< iws_mult*iws_size; i++)
    {
      thread_array[kid].iws[i] = 0;
    }
    for(Int i = 0; i < ews_mult*ews_size; i++)
    {
      thread_array[kid].ews[i] = 0;
    }
    #endif  //endif def basker_2dl
    //return 0;
  }//end init_workspace
  
  //--------------------------PRINT RELATED UTIL---------------------------//

  //print a given submatrix
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::print_factor
  (
   BASKER_MATRIX &L, 
   BASKER_MATRIX &U
  )
  {
    printf("L.nnz: %d  L.ncol: %d \n", L.nnz, L.ncol);
    Int k;
    Int Lnnz = L.col_ptr[L.ncol];
    for(k = 0; k < Lnnz; k++)
    {
      printf("L[%d]=%f " ,k , L.val[k]); 
    }
    printf("\n");
    for(k = 0; k < Lnnz; k++)
    {
      printf("Li[%d]=%d ", k , L.row_idx[k]);
    }
    printf("\n");
    for(k = 0; k < L.ncol; k++)
    {
      printf("Lp[%d]=%d ", k, L.col_ptr[k]);
    }
    printf("\n");


    printf("U.nnz: %d  U.ncol: %d \n", U.nnz, U.ncol);
    Int Unnz = U.col_ptr[U.ncol];
    for(k = 0; k < Unnz; k++)
    {
      printf("U[%d]=%f ", k, U.val[k]);
    }
    printf("\n");
    for(k = 0; k < Unnz; k++)
    {
      printf("Ui[%d]=%d ", k, U.row_idx[k]);
    }
    cout << endl << endl;
    for(k =0; k < U.ncol; k++)
    {
      printf("Up[%d] = %d ", k, U.col_ptr[k]);
    }
    printf("\n\n");
  }//end print_factor()


  //Print L out to file L.txt
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::printL()
  {
    Int total_Lnnz = 0;

    FILE *fp;
    fp = fopen("L.txt", "w");

    for(Int l = 0; l < tree.nblks; l++)
    {
      BASKER_MATRIX &myL = LL[l][0];

      for(Int k = 0; k < myL.ncol; k++)
      {
        fprintf(fp, "k=%d \n", k+myL.scol);

        for(Int j = myL.col_ptr[k]; j < myL.col_ptr[k+1]; j++)
        {
          fprintf(fp, "(%li, %li, %li) %li %e , ",
              k+myL.scol, myL.row_idx[j],
              gperm[myL.row_idx[j]],
              j,
              myL.val[j]);

          total_Lnnz++;
        }//nnz in col

        fprintf(fp, " \n \n ");
      }// over all col
    }//over all Ls

    fclose(fp);
    printf("Lnnz: %d \n", total_Lnnz);
    return 0;
  }//end printL()


  //2D Print L
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::printL2D()
  {
    Int total_Lnnz = 0;
    FILE *fp;
    fp = fopen("L.txt", "w");
    //over each blk
    for(Int l = 0; l < tree.nblks; l++)
    {
      //over each column
      for(Int k=0; k < LL[l][0].ncol; k++)
      {
        fprintf(fp, "k=%ld \n", (long)k+LL[l][0].scol);

        for(Int r = 0; r < LL_size[l]; r++)
        {
          BASKER_MATRIX &myL = LL[l][r];
          for(Int j = myL.col_ptr[k]; j < myL.col_ptr[k+1]; j++)
          {
            fprintf(fp, "(%ld , %ld , %ld, %ld, %ld) %g , ",
                (long)k+myL.scol, (long)myL.row_idx[j], 
                (long)myL.row_idx[j]+myL.srow,
                (long)myL.srow,
                (long)gperm(myL.row_idx[j]+myL.srow),
                myL.val[j]);

            total_Lnnz++;
          }//end over each nnz in column (k) of local U              
        }//end over each matrix row

        fprintf(fp, " \n \n ");
      }//end over each column
    }//end over all nblks

    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myL = LBTF(i);

        for(Int k = 0; k < myL.ncol; k++)
        {
          fprintf(fp, "k=%ld \n", (long)myL.scol+k);

          for(Int j = myL.col_ptr(k); j< myL.col_ptr(k+1); j++)
          {
            fprintf(fp, "(%ld , %ld , %ld, %ld) %g , ", 
                (long)k+myL.scol, myL.row_idx[j], 
                (long)myL.row_idx[j]+myL.srow,
                (long)gperm(myL.row_idx[j]+myL.srow),
                myL.val[j]);

            total_Lnnz++;
          }//over all nnz

          fprintf(fp, " \n \n ");
        }//over all columns
      }//over all blks
    }//end option btf

    fclose(fp);

    printf("Total L nnz: %ld \n", (long)total_Lnnz);
    printf("-----Done Print L2d ---------\n");

    return 0;
  }//print L2d()


  template <class Int, class Entry, class Exe_Space>
  int Basker<Int, Entry, Exe_Space>::printLMTX()
  {
    Int nnz = get_Lnnz();
    FILE *fp;
    fp = fopen("L.mtx", "w");

    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%%Generated by **Basker** \n");
    fprintf(fp, "%ld %ld %ld \n", (long)gn, (long)gn, (long)nnz);

    //over each blk
    //for(Int l = 0; l < tree.nblks; l++)
    for(Int l = 0; l < tree.nblks; l++)
    {
      //over each column
      for(Int k=0; k < LL[l][0].ncol; k++)
      {
        //fprintf(fp, "k=%d \n", k+LL[l][0].scol);
        for(Int r = 0; r < LL_size[l]; r++)
        {
          BASKER_MATRIX &myL = LL[l][r];

          for(Int j = myL.col_ptr[k]; j < myL.col_ptr[k+1]; j++)
          {
            fprintf(fp, "%ld %ld %g \n",
                (long)gperm(myL.row_idx(j)+myL.srow)+1,
                (long)k+myL.scol+1,
                myL.val(j));
          }//end over each nnz in column (k) of local U              
        }//end over each matrix row
      }//end over each column
    }//end over all nblks

    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myL = LBTF(i);

        for(Int k = 0; k < myL.ncol; k++)
        {
          for(Int j = myL.col_ptr(k); j< myL.col_ptr(k+1); j++)
          {
            fprintf(fp, "%ld %ld %g \n",
                (long)gperm(myL.row_idx(j)+myL.srow)+1,
                (long)k+myL.scol+1,
                myL.val(j));
          }//over all nnz
        }//over all columns
      }//over all blks
    }//end option btf

    fclose(fp);
    
    return 0;
  }//end printLMTX
  
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry, Exe_Space>::printUMTX()
  {
    Int nnz_u = get_Unnz();

    FILE *fp;
    fp = fopen("U.mtx", "w");
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%%Generated by **Basker** \n");
    fprintf(fp, "%ld %ld %ld \n", (long)gn, (long)gn, (long)nnz_u);

    //over each blks
    for(Int l = 0; l < tree.nblks; l++)
    {
      //over each column
      for(Int k = 0; k < LU[l][0].ncol; k++)
      {
        //over each row of U
        for(Int r = 0; r < LU_size[l]; r++)
        {
          BASKER_MATRIX &myU = LU[l][r];

          //over each nnz in column (k) of local U
          for(Int j = myU.col_ptr[k]; j < myU.col_ptr[k+1]; j++)
          {

            BASKER_ASSERT((myU.row_idx(j)+myU.srow+1)>0, "location 1-1");
            BASKER_ASSERT((k+myU.scol+1)>0, "location 1-2");

            fprintf(fp, "%ld %ld %g \n", 
                (long)myU.row_idx(j)+myU.srow+1,
                (long)k+myU.scol+1,
                myU.val(j));

          }//end over each nnz in column (k) of local U              
        }//end over each matrix row
      }//end over each column
    }//end over nblks

    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myU = UBTF[i];
        for(Int k = 0; k < myU.ncol; k++)
        {
          for(Int j = myU.col_ptr[k]; j< myU.col_ptr[k+1]; j++)
          {
            BASKER_ASSERT((myU.row_idx(j)+myU.srow+1)>0, "location 2-1");
            BASKER_ASSERT((k+myU.scol+1)>0, "location 2-2");
            fprintf(fp, "%ld %ld %g \n",
                (long)myU.row_idx(j)+myU.srow+1,
                (long)k+myU.scol+1,
                myU.val(j));
          }//over all nnz
        }//over all columns
      }//over all blks
    }//end option btf

    fclose(fp);
    printf("---------Done Print U ------\n");

    return 0;
  }//end printUMTX


  //Print U out to file U.txt
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int, Entry, Exe_Space>::printU()
  {
    FILE *fp;
    fp = fopen("U.txt", "w");
    
    //over each blks
    for(Int l = 0; l < tree.nblks; l++)
    {
      //over each column
      for(Int k = 0; k < LU[l][0].ncol; k++)
      {
        fprintf(fp, "k=%ld \n", (long)k+LU[l][0].scol);

        //over each row of U
        for(Int r = 0; r < LU_size[l]; r++)
        {
          BASKER_MATRIX &myU = LU[l][r];

          //over each nnz in column (k) of local U
          for(Int j = myU.col_ptr[k]; j < myU.col_ptr[k+1]; j++)
          {

            fprintf(fp, "(%ld , %ld , %ld) %f , ",
                (long)k+myU.scol, 
                (long)myU.row_idx[j], 
                (long)myU.row_idx[j]+myU.srow,
                myU.val[j]);

          }//end over each nnz in column (k) of local U              
        }//end over each matrix row

        fprintf(fp, " \n \n ");
      }//end over each column
    }//end over nblks

    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myU = UBTF[i];
        for(Int k = 0; k < myU.ncol; k++)
        {
          fprintf(fp, "k=%ld \n", (long)k+myU.scol);

          for(Int j = myU.col_ptr[k]; j< myU.col_ptr[k+1]; j++)
          {
            fprintf(fp, "(%ld , %ld , %ld) %f , ", 
                (long)k+myU.scol, 
                (long)myU.row_idx[j], 
                (long)myU.row_idx[j]+myU.srow, 
                myU.val[j]);
          }//over all nnz

          fprintf(fp, " \n \n ");
        }//over all columns
      }//over all blks
    }//end option btf

    fclose(fp);

    printf("---------Done Print U ------\n");

    return 0;
  }//end printU()

  //Print MTX
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::printMTX
  (
   std::string fname,
   BASKER_MATRIX &M
  )
  {
    //Matrix has been initalized
    if(M.ncol == 0)
      return;

    FILE *fp;
    fp = fopen(fname.c_str(), "w");
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%%Generated by **Basker** \n");
    fprintf(fp, "%%Starting Row %ld  Starting Col %ld \n",
	    (long)M.srow, (long)M.scol);
    fprintf(fp, "%ld %ld %ld \n", (long)M.nrow, (long)M.ncol, (long)M.nnz);
   
    Int bcol=M.scol;
    for(Int k=M.scol; k < M.scol+M.ncol; k++)
    {
      for(Int j=M.col_ptr[k-bcol]; j<M.col_ptr[k-bcol+1]; j++)
      {
        fprintf(fp, "%ld %ld %e \n", 
            (long)M.row_idx[j]+1, (long)k-bcol+1, M.val[j]); 
      }//over nnz in each column
    }//over each column

    fclose(fp);
  }//end printMTX() 

   //Print MTX
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::printMTX
  (
   std::string fname,
   BASKER_MATRIX &M, 
   BASKER_BOOL off
  )
  {
    FILE *fp;
    fp = fopen(fname.c_str(), "w");
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%%Generated by **Basker** \n");
    fprintf(fp, "%%Starting Row %d  Starting Col %d \n",
	    M.srow, M.scol);
    fprintf(fp, "%ld %ld %ld \n", (long)M.nrow, (long)M.ncol, (long)M.nnz);
   
    Int bcol=M.scol;
    Int brow=M.srow;
    for(Int k=M.scol; k < M.scol+M.ncol; k++)
    {
      for(Int j=M.col_ptr[k-bcol]; j<M.col_ptr[k-bcol+1]; j++)
      {
        if(off == BASKER_FALSE)
        {
          fprintf(fp, "%ld %ld %e \n", 
              (long)M.row_idx[j]+1, (long)k-bcol+1, M.val[j]);
        }
        else
        {
          fprintf(fp, "%ld %ld %e \n", 
              (long)M.row_idx[j]+1-brow, (long)k-bcol+1,M.val[j]);
        }
      }//over nnz in each column
    }//over each column

    fclose(fp);

    printf("Done Writing Matrix \n");
  }//end printMTX() 

  template<class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::readMTX(std::string fname, BASKER_MATRIX &M)
  {
    //Note: Adapted from Siva's original bsk_util

    cout << "ReadingMTX " << fname << endl;
    ifstream inp_str;
    inp_str.open(fname, ios::in);
    Int i, j;
    Int nrows, ncols, nnz, true_nnz;
    Entry val;
    std::string s;
    size_t p1, p2, p3;
    Int ptype, sym_type;

    if (inp_str.is_open())
    {
      getline(inp_str, s);
      //cout << s << endl;

      // Check if matrix is pattern-only or symmetric
      p1 = s.find("pattern");
      if (p1 != string::npos)
      { ptype = 2; }
      else
      { ptype = 3; }

      p1 = s.find("symmetric");
      p2 = s.find("hermitian");
      p3 = s.find("skew-symmetric");

      if ((p1 != string::npos) || (p2 != string::npos) || (p3 != string::npos))
      { sym_type = 1; }
      else
      { sym_type = 0; }

      (void)ptype; //NDE silence warning
      (void)sym_type; //NDE silence warning

      while (inp_str.peek() == '%') // Skip the comments.
        getline(inp_str, s);

      // Find the # of rows, cols and nnzs.
      inp_str >> nrows;
      inp_str >> ncols;
      inp_str >> nnz;

      cout << nrows << " " << ncols  << " " << nnz << endl;
      M.ncol = ncols;
      M.nrow = nrows;
      M.nnz = nnz;

      BASKER_ASSERT((ncols+1)>0, "util ncols");
      MALLOC_INT_1DARRAY(M.col_ptr, ncols+1);
      init_value(M.col_ptr, ncols+1,(Int) 0);
      BASKER_ASSERT(nnz>0, "util nnz");
      MALLOC_INT_1DARRAY(M.row_idx, nnz);
      init_value(M.row_idx, nnz, (Int) 0);
      MALLOC_ENTRY_1DARRAY(M.val, nnz);
      init_value(M.val, nnz, (Entry) 0.0);
      Int innz = 0;

      while(nnz > 0)
      {
        inp_str >> i;
        //cout << "row: " << i-1 ;
        M.row_idx[innz] = i-1;
        inp_str >> j;
        //cout << " col: " << j-1;
        M.col_ptr[j] = M.col_ptr[j]+1;
        inp_str >> val;
        //cout << " val: " << val << endl;
        M.val[innz] = val;
        //Other type options..
        innz++;
        nnz--;
      }
      inp_str.close();

      //count col_sums
      for(Int k =1 ; k<(ncols+1); k++)
      {
        M.col_ptr[k] = M.col_ptr[k] +M.col_ptr[k-1];
      }

    }//end if open
  
    //M.print();

  }//end readMTX()


  //Print out RHS  RHS.txt
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::printRHS()
  {
    if(solve_flag == false)
    {return -1;}

    FILE *fp;
    fp = fopen("RHS.txt", "w");

    //over each row
    for(Int r = 0; r < A.nrow; r++)
    {
      //over each column NOTE: come back to
      //for(Int k = 0; k < rhs.size(); k++)
      for(Int k = 0; k < 1; k++)
      {
        fprintf(fp, "%ld %ld %f, ", (long)r, (long)gperm[r], rhs[k][r]);
      }//end over each column
      fprintf(fp, "\n");
    }//end over each row

    fclose(fp);

    return 0;
  }//end printRHS()

  //Print solution SOL.txt
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::printSOL()
  {
    if(solve_flag == false)
    {return -1;}

    FILE *fp;
    fp = fopen("SOL.txt", "w");
    
    //over each row
    for(Int r = 0; r < A.nrow; r++)
    {
      //over each column Note: come back to
      //for(Int k = 0; k < rhs.size(); k++)
      for(Int k = 0 ; k < 1; k++)
      {
        fprintf(fp, "%ld %ld %f, ", (long)r, (long)gperm[r], sol[k][r]);
      }//end over each column
      fprintf(fp, "\n");
    }//end over each row
    
    fclose(fp);

    return 0;
  }//end printSOL()

  //Prints the given tree into a file to analyze
  template<class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::printTree()
  {
    FILE *fp;
    fp = fopen("Tree_Info.txt", "w");

    fprintf(fp,"nblks: %ld \n", (long)tree.nblks);

    fprintf(fp,"row_tabs: ");

    for(Int i=0; i < tree.nblks+1; i++)
    {
      fprintf(fp, "%ld ", (long)tree.row_tabs[i]);
    }

    fprintf(fp, "\n");
    fprintf(fp, "col_tabs: ");

    for(Int i=0; i < tree.nblks+1; i++)
    {
      fprintf(fp, "%ld ", (long)tree.col_tabs[i]);
    }

    fprintf(fp, "\n");
    fprintf(fp, "tree_tabs: ");

    for(Int i=0; i< tree.nblks+1; i++)
    {
      fprintf(fp, "%ld ", (long)tree.treetab[i]);
    }

    fprintf(fp, "\n");
    fprintf(fp, "Sep Sizes: ");

    for(Int i=0; i < tree.nblks; i++)
    {
      fprintf(fp, "%ld %ld \n", 
          (long)i, (long)tree.col_tabs[i+1]-(long)tree.col_tabs[i]);
    }

    fprintf(fp, "\n");

    fclose(fp);
  }//end printTree()


  //just prints out the data for imbalance in last ll sep
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::print_sep_bal()
  {

    FILE *fp = fopen("ASEP.csv", "w");

    for(Int b = 0; b < LU_size[tree.nblks-1]; ++b)
    {
      BASKER_MATRIX &M = AVM(tree.nblks-1)(b);

      for(Int k = 0; k < M.ncol; ++k)
      {
        fprintf(fp, "%ld, ",(long)(M.col_ptr(k+1)-M.col_ptr(k)) );
      }//over all columns

      fprintf(fp, "\n");
    }//over all blks

    fclose(fp);

    FILE *fp2 = fopen("LSEP.csv", "w");

    for(Int b = 0; b < LU_size[tree.nblks-1]; ++b)
    {
      BASKER_MATRIX &M = LU(tree.nblks-1)(b);

      for(Int k = 0; k < M.ncol; ++k)
      {
        fprintf(fp2, "%ld, ", (long)(M.col_ptr(k+1)-M.col_ptr(k)) );
      }//over all columns

      fprintf(fp2, "\n");
    }//over all blks

    fclose(fp2);    
  }//end print_sep_bal()
  

  //This will be stored in the stats manager in future
  template <class Int, class Entry, class Exe_Space>
  Int Basker<Int,Entry,Exe_Space>::get_Lnnz()
  {
    //Add Check
    //A side note totoal LNNZ might be bigger than INT size
    //May want to use a LO and GO type in the future
    Int total_nnz = 0;

    //Get total from ND TREE
    for(Int l = 0; l < tree.nblks; ++l)
    {
      //over each Lower half  
      for(Int r = 0; r < LL_size(l); r++)
      {
        BASKER_MATRIX &myL = LL(l)(r);
        total_nnz += myL.col_ptr(myL.ncol);  
      }//end over each matrix row
    }//end over all nblks

    #ifdef BASKER_DEBUG_UTIL
    printf("nnz in ND: %d \n", total_nnz);
    #endif

    //Get nnz for BTF L
    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myL = LBTF(i);
        total_nnz += myL.col_ptr(myL.ncol);
      }//over all blks
    }//end option btf

    #ifdef BASKER_DEBUG_UTIL
    printf("Total L nnz: %d \n", total_nnz);
    #endif
    //printf("-----Done Print L2d ---------\n");

    return total_nnz;
  }//end get_Lnnz()


  template<class Int, class Entry, class Exe_Space>
  Int Basker<Int,Entry,Exe_Space>::get_Unnz()
  {
    Int total_nnz =0;

    //Get total from ND TREE
    for(Int l = 0; l < tree.nblks; ++l)
    {
      //over each Lower half  
      for(Int r = 0; r < LU_size(l); r++)
      {
        BASKER_MATRIX &myU = LU(l)(r);
        total_nnz += myU.col_ptr(myU.ncol);  
      }//end over each matrix row
    }//end over all nblks

    #ifdef BASKER_DEBUG_UTIL
    printf("nnz in ND: %d \n", total_nnz);
    #endif

    //Get nnz for BTF L
    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myU = UBTF(i);
        total_nnz += myU.col_ptr(myU.ncol);
      }//over all blks
    }//end option btf

    #ifdef BASKER_DEBUG_UTIL
    printf("Total U nnz: %d \n", total_nnz);
    #endif
    //printf("-----Done Print L2d ---------\n");

    return total_nnz;
  }//end get_Unnz()


  //Provides L to a set of arrays
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::get_L
  (
   Int      &n,
   Int      &nnz,
   Int**    col_ptr,
   Int**    row_idx,
   Entry**  val
  )
  {
    //Add Check
    //Kokkos::Impl::Timer timer;
    n   = gn;
    nnz = get_Lnnz();
    (*col_ptr) = new Int[gn+1]();
    (*row_idx) = new Int[nnz]();
    (*val)     = new Entry[nnz]();

    Int ptr = 0;
    Int total_Lnnz = 0;
    //over each blk
    for(Int l = 0; l < tree.nblks; l++)
    {
      //over each column
      for(Int k=0; k < LL(l)(0).ncol; ++k)
      {
        for(Int r = 0; r < LL_size(l); r++)
        {
          BASKER_MATRIX &myL = LL[l][r];
          Int brow = myL.srow;
          Int bcol = myL.scol;

          for(Int j = myL.col_ptr(k); j < myL.col_ptr(k+1); ++j)
          {
            (*row_idx)[ptr]=gperm(myL.row_idx(j)+myL.srow);
            (*val)[ptr]    = myL.val(j);
            ptr++;
            total_Lnnz++;
          }//end over each nnz in column (k) of local U              
        }//end over each matrix row
        (*col_ptr)[k+LL(l)(0).scol+1] = ptr; 
      }//end over each column
    }//end over all nblks
    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myL = LBTF(i);
        Int brow = myL.srow;
        Int bcol = myL.scol;

        for(Int k = 0; k < myL.ncol; k++)
        {
          for(Int j = myL.col_ptr(k); j< myL.col_ptr(k+1); j++)
          {
            (*row_idx)[ptr]= gperm(myL.row_idx(j)+myL.srow);
            (*val)[ptr]    = myL.val(j);
            ptr++;
            total_Lnnz++;
          }//over all nnz
          //printf( "\n \n");
          (*col_ptr)[k+myL.scol+1] = ptr;
        }//over all columns
      }//over all blks
    }//end option btf

    #ifdef BASKER_DEBUG_UTIL
    printf("Total L nnz: %d ptr: %d  \n", total_Lnnz, ptr);
    printf("-----Done Print GetL ---------\n");
    #endif

    return 0;
  }//end get_L();


  //Provides U to a set of arrays
  template<class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::get_U
  (
   Int     &n,
   Int     &nnz,
   Int**    col_ptr,
   Int**    row_idx,
   Entry**  val
  )
  {
    //Add Check
    n       = gn;
    nnz = get_Unnz();
    (*col_ptr) = new Int[gn+1]();
    (*row_idx) = new Int[nnz]();
    (*val)     = new Entry[nnz]();

    Int total_Unnz = 0;
    Int ptr = 0;
    //over each blks
    for(Int l = 0; l < tree.nblks; l++)
    {
      //over each column
      for(Int k = 0; k < LU(l)(0).ncol; k++)
      {
        for(Int r = 0; r < LU_size(l); r++)
        {
          BASKER_MATRIX &myU = LU(l)(r);

          for(Int j = myU.col_ptr(k); j < myU.col_ptr(k+1); j++)
          {
            (*row_idx)[ptr] = myU.row_idx(j)+myU.srow;
            (*val)[ptr]     = myU.val(j);
            ptr++;
            total_Unnz++;
          }//end over each nnz in column (k) of local U              
        }//end over each matrix row

        (*col_ptr)[k+LU(l)(0).scol+1] = ptr;
      }//end over each column
    }//end over nblks

    if(Options.btf == BASKER_TRUE)
    {
      Int nblks = btf_nblks-btf_tabs_offset;
      for(Int i =0; i < nblks; i++)
      {
        BASKER_MATRIX &myU = UBTF(i);
        Int brow = myU.srow;
        Int bcol = myU.scol;
        for(Int k = 0; k < myU.ncol; k++)
        {
          for(Int j = myU.col_ptr[k]; j< myU.col_ptr[k+1]; j++)
          {
            (*row_idx)[ptr] = myU.row_idx(j)+myU.srow;
            (*val)[ptr]     = myU.val(j);

            ptr++;
            total_Unnz++;
          }//over all nnz
          //printf("\n \n");
          (*col_ptr)[k+myU.scol+1] = ptr;
        }//over all columns
      }//over all blks
    }//end option btf

    #ifdef BASKER_DEBUG_UTIL
    printf("---------Done Get U : %d ------\n",  total_Unnz);
    #endif

    return 0;
  }//end get_U()


  //returns global p
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::get_p(Int** p)
  {
    //(*p) = new Int[gperm[0].size()];
    (*p) = new Int[A.nrow]; // NDE when is cleanup of this?

    for(Int i=0; i < A.nrow; i++)
    {
      (*p)[i] = gperm[i];
    }

    return 0;
  }//end get_p()
  
  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::matrix_transpose
  (
   BASKER_MATRIX &M, //BTF_A
   BASKER_MATRIX &AT //T, blank
  )
  {
    //const Int brow = M.srow; Not used
    //Setup what we do know
    AT.srow = M.srow;
    AT.nrow = M.ncol;
    AT.scol = M.scol;
    AT.ncol = M.nrow;
    AT.nnz  = M.nnz;

    BASKER_ASSERT((AT.ncol+1)>0, "util trans ncol");
    MALLOC_INT_1DARRAY(AT.col_ptr, AT.ncol+1);
    init_value(AT.col_ptr, AT.ncol+1, (Int)0);
    MALLOC_INT_1DARRAY(AT.row_idx, AT.nnz);
    BASKER_ASSERT((AT.nnz)>0, "util trans nnz");
    init_value(AT.row_idx, AT.nnz, (Int)0);
    MALLOC_ENTRY_1DARRAY(AT.val    , AT.nnz);
    init_value(AT.val,     AT.nnz, (Entry)1.0);

    //Setup a litte workspace
    const Int ws_size = M.nrow;
    INT_1DARRAY ws;
    BASKER_ASSERT(ws_size > 0, "util trans ws");
    MALLOC_INT_1DARRAY(ws, ws_size);

    for(Int j = 0; j < ws_size; ++j)
    {
      ws(j) = (Int) 0;
    }
    
    //Note could get number of nonzeros here inplace of nnz() for faster

    //get row counts
    Int total =0 ;
    Int maxv  =0;
    for(Int j = M.col_ptr(0); j < M.col_ptr(M.ncol); ++j)
    {
      if(M.row_idx(j) > maxv)
      {
        maxv = M.row_idx(j);
      }
      if(M.row_idx(j) > (ws_size-1))
      {
        printf("error type 1\n");
      }
      ws(M.row_idx(j)) = ws(M.row_idx(j)) + 1;

      total++;
    }
    //printf("got row counts, total: %d  %d \n", total, M.nnz);
    //printf("debug 0: %d 1: %d 2: %d 3: %d \n",
    //	   ws(0), ws(1), ws(2), ws(3));
    //printf("max idx: %d nrow: %d wssize: %d \n", 
    //	   maxv, M.nrow, ws_size); 

    //write stupid code!
    //add them all up
    total = 0;
    //Int total2 = 0; //Not used
    for(Int j = 0; j < ws_size; ++j)
    {
      total = total + ws(j);
      //total2 = total2 + ws_test[j];
    }
   
    for(Int j = 1; j < M.nrow; ++j)
    {
      ws(j)  = ws(j) + ws(j-1);
    }//for-j
  
    //copy over to AT
    AT.col_ptr(0) = (Int) 0;
    for(Int j = 1; j <= M.nrow; ++j)
    {
      AT.col_ptr(j) = ws(j-1);
    }

    //set ws
    for(Int j = 0; j < M.nrow; ++j)
    {
      ws(j) = AT.col_ptr(j);
    }

    for(Int k = 0; k < M.ncol; ++k)
    {
      for(Int j = M.col_ptr(k); j < M.col_ptr(k+1); ++j)
      {
        if(ws(M.row_idx(j)) >= AT.nnz)
        { 
          printf("error \n");
        }
        AT.row_idx(ws(M.row_idx(j))++) = k; 
        //starting at zero
      }
    }
    
    FREE_INT_1DARRAY(ws);
  }//end matrix_transpose

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::matrix_transpose
  (
   BASKER_MATRIX_VIEW &MV,
   BASKER_MATRIX &AT
  )
  {
    Int brow = MV.srow;
    //Setup what we do know
    AT.srow = MV.srow;
    AT.nrow = MV.nrow;
    AT.scol = MV.scol;
    AT.ncol = MV.ncol;
    AT.nnz  = MV.nnz();

    BASKER_ASSERT((AT.ncol+1)>0, "util trans ncol2");
    MALLOC_INT_1DARRAY(AT.col_ptr, AT.ncol+1);
    init_value(AT.col_ptr, AT.ncol+1, (Int)0);
    BASKER_ASSERT(AT.nnz > 0, "util trans nnz2");
    MALLOC_INT_1DARRAY(AT.row_idx, AT.nnz);
    init_value(AT.row_idx, AT.nnz, (Int)0);
    MALLOC_ENTRY_1DARRAY(AT.val    , AT.nnz);
    init_value(AT.val,     AT.nnz, (Entry)1.0);

    //Setup a litte workspace
    Int ws_size = MV.nrow;
    INT_1DARRAY ws;
    BASKER_ASSERT(ws_size > 0, "util trans ws2");
    MALLOC_INT_1DARRAY(ws, ws_size);
    init_value(ws, ws_size, (Int)0);
    
    //Note could get number of nonzeros here inplace of nnz() for faster

    //get row counts
    for(Int j = MV.col_ptr(MV.scol); j < MV.col_ptr(MV.scol+MV.ncol); ++j)
    {
      if(MV.good(j) != 0)
      {
        continue;
      }
      ws[MV.row_idx(j)-brow]++;
    }

    AT.col_ptr[1] = ws[0];
    for(Int j = 1; j < AT.nrow; j++)
    {
      ws[j] = ws[j]+ws[j-1];
      AT.col_ptr[j+1] = ws[j];
      ws[j-1] = AT.col_ptr[j-1];
    }
    ws[AT.nrow-1] = AT.col_ptr[AT.nrow-1];
    
    for(Int k = 0; k < AT.ncol; k++)
    {
      for(Int j = MV.col_ptr(MV.scol+k); j < MV.col_ptr(MV.scol+k+1); ++j)
      {
        if(MV.good(j) != 0)
        {
          continue;
        }

        AT.row_idx[ws[MV.row_idx(j)-brow]++] = k; //starting at zero
      }
    }

    FREE_INT_1DARRAY(ws);
  }//end matrix_transpose


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::matrix_transpose
  (
   const Int sm_ ,
   const Int m_,
   const Int sn_ ,
   const Int n_,
   const Int nnz_, 
   Int *col_ptr,
   Int *row_idx,
   Entry *val,
   BASKER_MATRIX &AT
  )
  {
    AT.srow  = sn_;
    AT.nrow  = n_;
    AT.scol  = sm_;
    AT.ncol  = m_;
    AT.nnz   = nnz_;

    if(matrix_flag == BASKER_FALSE)
    {
      BASKER_ASSERT((AT.ncol+1)>0, "util trans ncol");
      MALLOC_INT_1DARRAY(AT.col_ptr, AT.ncol+1);
      MALLOC_INT_1DARRAY(AT.row_idx, AT.nnz);
      BASKER_ASSERT((AT.nnz)>0, "util trans nnz");
      MALLOC_ENTRY_1DARRAY(AT.val    , AT.nnz);
    }    

    init_value(AT.col_ptr, AT.ncol+1, (Int)0);
    init_value(AT.row_idx, AT.nnz, (Int)0);
    init_value(AT.val,     AT.nnz, (Entry)1.0);

    //Setup a litte workspace
    const Int ws_size = m_;
    INT_1DARRAY ws;
    BASKER_ASSERT(ws_size > 0, "util trans ws");
    MALLOC_INT_1DARRAY(ws, ws_size);

    for(Int j = 0; j < ws_size; ++j)
    {
      ws(j) = (Int) 0;
    }

    Int total =0 ;
    Int maxv  =0;
    for(Int j = col_ptr[0]; j < col_ptr[n_]; ++j)
    {
      if(row_idx[j] > maxv)
      {
        maxv = row_idx[j];
      }
      if(row_idx[j] > (ws_size-1))
      {
        printf("error type 1\n");
      }
      ws(row_idx[j]) = ws(row_idx[j]) + 1;

      total++;
    }

    //write stupid code!
    //add them all up
    total = 0;
    //Int total2 = 0; //Not used
    for(Int j = 0; j < ws_size; ++j)
    {
      total = total + ws(j);
      //total2 = total2 + ws_test[j];
    }

    for(Int j = 1; j < n_; ++j)
    {
      ws(j)  = ws(j) + ws(j-1);

    }//for-j

    //copy over to AT
    AT.col_ptr(0) = (Int) 0;
    for(Int j = 1; j <= n_; ++j)
    {
      AT.col_ptr(j) = ws(j-1);
    }

    //set ws
    for(Int j = 0; j < n_; ++j)
    {
      ws(j) = AT.col_ptr(j);
    }

    for(Int k = 0; k < n_; ++k)
    {
      for(Int j = col_ptr[k]; j < col_ptr[k+1]; ++j)
      {
        if(ws(row_idx[j]) >= AT.nnz)
        { 
          printf("error \n");
        }
        //AT.row_idx(ws(row_idx[j])++) = k;
        AT.row_idx(ws(row_idx[j])) = k;
        AT.val(ws(row_idx[j])) = val[j]; // NDE: This will need to be tracked...
        ws(row_idx[j])++;
        //starting at zero
      }
    }

    sort_matrix(AT);
    //printMTX("A_TRANS.mtx", AT);

    FREE_INT_1DARRAY(ws);
  }//end matrix_transpos


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::matrix_transpose
  (
   const Int sm_ ,
   const Int m_,
   const Int sn_ ,
   const Int n_,
   const Int nnz_, 
   Int *col_ptr,
   Int *row_idx,
   Entry *val,
   BASKER_MATRIX &AT,
   INT_1DARRAY &vals_transpose_local
  )
  {
    AT.srow  = sn_;
    AT.nrow  = n_;
    AT.scol  = sm_;
    AT.ncol  = m_;
    AT.nnz   = nnz_;

    if(matrix_flag == BASKER_FALSE)
    {
      BASKER_ASSERT((AT.ncol+1)>0, "Basker matrix_transpose assert: ncol > 0 failed");
      MALLOC_INT_1DARRAY(AT.col_ptr, AT.ncol+1);
      BASKER_ASSERT((AT.nnz)>0, "Basker matrix_transpose assert: nnz > 0 failed");
      MALLOC_INT_1DARRAY(AT.row_idx, AT.nnz);
      MALLOC_ENTRY_1DARRAY(AT.val    , AT.nnz);
    }    

    init_value(AT.col_ptr, AT.ncol+1, (Int)0);
    init_value(AT.row_idx, AT.nnz, (Int)0);
    init_value(AT.val,     AT.nnz, (Entry)1.0);

    //Setup a litte workspace
    const Int ws_size = m_;
    INT_1DARRAY ws;
    BASKER_ASSERT(ws_size > 0, "util trans ws");
    MALLOC_INT_1DARRAY(ws, ws_size);

    for(Int j = 0; j < ws_size; ++j)
    {
      ws(j) = (Int) 0;
    }

    Int total =0 ;
    Int maxv  =0;
    for(Int j = col_ptr[0]; j < col_ptr[n_]; ++j)
    {
      if(row_idx[j] > maxv)
      {
        maxv = row_idx[j];
      }
      if(row_idx[j] > (ws_size-1))
      {
        printf("error type 1\n");
      }
      ws(row_idx[j]) = ws(row_idx[j]) + 1;

      total++;
    }

    //write stupid code!
    //add them all up
    total = 0;
    //Int total2 = 0; //Not used
    for(Int j = 0; j < ws_size; ++j)
    {
      total = total + ws(j);
      //total2 = total2 + ws_test[j];
    }

    for(Int j = 1; j < n_; ++j)
    {
      ws(j)  = ws(j) + ws(j-1);

    }//for-j

    //copy over to AT
    AT.col_ptr(0) = (Int) 0;
    for(Int j = 1; j <= n_; ++j)
    {
      AT.col_ptr(j) = ws(j-1);
    }

    //set ws
    for(Int j = 0; j < n_; ++j)
    {
      ws(j) = AT.col_ptr(j);
    }

    for(Int k = 0; k < n_; ++k)
    {
      for(Int j = col_ptr[k]; j < col_ptr[k+1]; ++j)
      {
        if(ws(row_idx[j]) >= AT.nnz)
        { 
          printf("error \n");
        }
        //AT.row_idx(ws(row_idx[j])++) = k;
        AT.row_idx(ws(row_idx[j])) = k;
        AT.val(ws(row_idx[j])) = val[j]; // NDE: This will need to be tracked...
          vals_transpose_local(j) = ws(row_idx[j]); // NDE: Added to track the transpose operation
          //vals_transpose_local(ws(row_idx[j])) = j; // NDE: Added to track the transpose operation
        ws(row_idx[j])++;
        //starting at zero
      }
    }

    sort_matrix_store_valperms(AT,vals_transpose_local);
   // sort_matrix(AT);
    //printMTX("A_TRANS.mtx", AT);

    FREE_INT_1DARRAY(ws);
  }//end matrix_transpos


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry,Exe_Space>::printVec
  (
   std::string fname, 
   INT_1DARRAY x, 
   Int n
  )
  {
    FILE *fp;
    fp = fopen(fname.c_str(), "w");

    for(Int i = 0; i < n; i++)
    {
      fprintf(fp, "%ld \n", (long)x(i));
    }

    fclose(fp);
  }//end printVec(file,Int);


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry,Exe_Space>::printVec
  (
   std::string fname, 
   ENTRY_1DARRAY x, 
   Int n
  )
  {
    FILE *fp;
    fp = fopen(fname.c_str(), "w");

    for(Int i = 0; i < n; i++)
    {
      fprintf(fp, "%f \n", x(i));
    }

    fclose(fp);
  }//end printVec(file,Int);


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::printVec
  (
   INT_1DARRAY x,
   Int n
  )
  {
    printf("---VECTOR: %d ----\n", n);

    for(Int i = 0 ; i < n;  i++)
    {
      printf("%ld %ld, \n", (long)i, (long)x(i));
    }

    printf("---END VECTOR %d --\n", n);
  }//end printVec(Int)


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::printVec
  (
   ENTRY_1DARRAY x,
   Int n
  )
  {
    printf("---VECTOR: %d ----\n", n);

    for(Int i = 0; i< n; i++)
    {
      printf("%ld %g, \n", (long)i, x[i]);
    }

    printf("---END VECTOR: %d ---\n", n);
  }//end printVec Entry


  template <class Int, class Entry, class Exe_Space>
  inline
  Int Basker<Int, Entry, Exe_Space>::t_get_kid
  (
   const TeamMember &thread
  )
  {
    return (Int)(thread.league_rank()*thread.team_size()+
		 thread.team_rank());
  }//end t_get_kid
 

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::get_total_perm
  (
   INT_1DARRAY outp_l,
   INT_1DARRAY outp_r
  )
  {
    INT_1DARRAY temp;
    MALLOC_INT_1DARRAY(temp, gn);

    //Fill with normal ordering
    for(Int i = 0; i < gm; ++i)
    {
      outp_l(i) = (Int) i;
    }

    for(Int i = 0; i < gn; ++i)
    {
      outp_r(i) = (Int) i;
    }

    //Add Match ordering
    if(match_flag == BASKER_TRUE)
    {
      permute_inv(outp_l, order_match_array,gn);
    }

    if(btf_flag   == BASKER_TRUE)
    {
      //Strong comp btf
      permute_inv(outp_l, order_btf_array,gn);
      permute_inv(outp_r, order_btf_array,gn);

      //AMD on blks
      permute_inv(outp_l, order_blk_amd_array,gn);
      permute_inv(outp_r, order_blk_amd_array,gn);
    }

    if(nd_flag == BASKER_TRUE)
    {
      //ND ordering
      for(Int i=0; i < BTF_A.ncol; i++)
      {
        temp(i) = part_tree.permtab(i);
      }
      for(Int i = BTF_A.ncol; i < gn; i++)
      {
        temp(i) = i;
      }

      permute_inv(outp_l,temp, gn);
      permute_inv(outp_r,temp, gn);
    }

    if(amd_flag == BASKER_TRUE)
    {
      for(Int i=0; i < BTF_A.ncol; i++)
      {
        temp(i) = order_csym_array(i);
      }

      for(Int i = BTF_A.ncol; i < gn; i++)
      {
        temp(i) = i;
      }

      //CAMDSYM
      permute_inv(outp_l, temp, gn);
      permute_inv(outp_r, temp, gn);
    }

    permute_inv(outp_l, gperm, gn);

    //Debug
    if(Options.verbose_matrix_out == BASKER_TRUE)
    {
      printVec("left_perm.csc", outp_l, gn);
      printVec("right_perm.csc", outp_r, gn);
    }

    FREE_INT_1DARRAY(temp);
  }//end get_total_perm


  //We need an easier and faster way to do this.  
  //Could get very big
  //We should use a dynamic build up
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int Basker<Int,Entry,Exe_Space>::find_leader(Int kid, Int l)
  {
    l = l+1;
    Int my_token = S(l)(kid);
    Int my_loc = kid;
    while((my_loc > 0))
    {
      my_loc--;
      if(S(l)(my_loc) != my_token)
      {
        my_loc++;
        break;
      }
    }

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("find_leader, kid: %d l: %d leader: %d \n", kid, l, my_loc);
    #endif

    return my_loc;
  }//end find_leader()


  //Added print function
  //I like printf because it is not a thread race dependend like 
  //c++ streams, however be may get compiler warnings
  //Come back and add own printf style calls
 
}//end namespace basker

#endif //end basker_util_hpp
