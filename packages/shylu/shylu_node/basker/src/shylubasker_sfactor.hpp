// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_SFACTOR_HPP
#define SHYLUBASKER_SFACTOR_HPP

//-------------NEED TO CONVERT to no view

//Note: need to finish the basker_max_idx trans form

/*This is the new outline for basker symbolic factoriziation*/
/*This framework assumes that A has been ordering/partitioned already*/
//Note: in future will want to use preallocated workspace

//Note: Options Used Here
//verbose
//symmetric
//AtA

//Include Files
#include "shylubasker_types.hpp"
#include "shylubasker_structs.hpp"
#include "shylubasker_util.hpp"

#include <iostream>
using namespace std;

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif
#include "Teuchos_OrdinalTraits.hpp"

//#define BASKER_TIMER 
//#define BASKER_DEBUG_SFACTOR

//Functor for Kokkos
namespace BaskerNS
{
  //Kokkos functor to init space for L and U
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_sfactor_init_workspace
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

    bool flag;
    Basker<Int,Entry,Exe_Space> *basker;
    
    kokkos_sfactor_init_workspace()
    {}

    kokkos_sfactor_init_workspace(bool _flag, Basker<Int,Entry,Exe_Space> *_b)
    {
      flag = _flag;
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
      //Int kid = (Int)(thread.league_rank()*thread.team_size()+
      //	      thread.team_rank());
      Int kid = basker->t_get_kid(thread);
      #endif

      basker->t_init_workspace(flag, kid);

    }//end operator()
  }; //end struct kokkos_sfactor_init_workspace

    //Kokkos functor to init space for L and U
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_sfactor_init_factor
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

    Basker<Int,Entry,Exe_Space> *basker;
    
    kokkos_sfactor_init_factor()
    {}

    kokkos_sfactor_init_factor(Basker<Int,Entry,Exe_Space> *_b)
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
      //Int kid = (Int)(thread.league_rank()*thread.team_size()+
      //	      thread.team_rank());
      Int kid = basker->t_get_kid(thread);
      #endif

      basker->t_init_factor(kid);

      //This needs to be done earlier in ordering now
      //basker->t_init_2DA(kid);

    }//end operator()
  }; //end struct kokkos_sfactor_init_factor

}//end namespace BaskerNS

//SFactor functions fo basker
namespace BaskerNS
{
 
template <class Int, class Entry, class Exe_Space>
BASKER_INLINE
int Basker<Int, Entry, Exe_Space>::sfactor()
{
  //check if already called
  if(symb_flag == true) 
  {return -1;}

  #ifdef BASKER_DEBUG_SFACTOR
  printf("Default Symbolic Called (btf_tabs_offset = %d)\n",btf_tabs_offset);
  #endif     

  // flag indicate if some symbolic & allocation are delayed
  //  because MWM & ND are applied to the big A block (in numeric Factor)
  bool allocate_nd_workspace = (Options.blk_matching == 0 && Options.static_delayed_pivot == 0);
  bool setup_flag = (allocate_nd_workspace || btf_tabs_offset == 0);
  if(btf_tabs_offset != 0 && allocate_nd_workspace)
  {
    symmetric_sfactor();

    if(Options.verbose == BASKER_TRUE)
    {
      printf("\n");
      printf("----------------------------------\n");
      printf("Total NNZ: %ld \n", (long)global_nnz);
      printf(" > blk_matching = %d\n", (int)Options.blk_matching );
      printf("----------------------------------\n");
      printf("\n");
    }
  }

  if(Options.btf == BASKER_TRUE)
  {
    btf_last_dense(setup_flag);
  }

  //Allocate Factorspace
  //printf(" >> kokkos_sfactor_init_factor( btf_tabs_offset = %d, allocate_nd_workspace = %d ) <<\n",
  //       btf_tabs_offset,allocate_nd_workspace);
  if(btf_tabs_offset != 0 && allocate_nd_workspace)
  {
  #ifdef BASKER_KOKKOS
    kokkos_sfactor_init_factor<Int,Entry,Exe_Space>
      iF(this);
    Kokkos::parallel_for(TeamPolicy(num_threads,1), iF);
    Kokkos::fence();
  #else
  #endif
  }

  //if(btf_tabs_offset != 0)
  {
    //Allocate workspace
  #ifdef BASKER_KOKKOS
    typedef Kokkos::TeamPolicy<Exe_Space>      TeamPolicy;
    kokkos_sfactor_init_workspace<Int,Entry,Exe_Space>
      iWS(setup_flag, this);
    Kokkos::parallel_for(TeamPolicy(num_threads,1), iWS);
    Kokkos::fence();
  #endif
  }

  BASKER_ASSERT(A.nrow > 0, "Basker sfactor assert: A.nrow > 0 failed");
  MALLOC_INT_1DARRAY(gperm, A.nrow);              // pivot
  init_value(gperm,A.nrow, BASKER_MAX_IDX);
  MALLOC_INT_1DARRAY(gpermi,A.nrow);              // ipiv
  init_value(gpermi, A.nrow, BASKER_MAX_IDX);
  MALLOC_INT_1DARRAY(gperm_array, A.nrow);        // perm
  MALLOC_INT_1DARRAY(gpermi_array, A.nrow);       // permi

  //Incomplete Factor Setup
  if(Options.incomplete == BASKER_TRUE)
  {
    Int lvl_nnz = 1.2*global_nnz;
    MALLOC_INT_1DARRAY(INC_LVL_ARRAY, lvl_nnz);
    init_value(INC_LVL_ARRAY, lvl_nnz, BASKER_MAX_IDX);
    MALLOC_INT_1DARRAY(INC_LVL_ARRAY_CNT, A.nrow);
    init_value(INC_LVL_ARRAY_CNT, A.nrow,(Int)0);
    MALLOC_INT_1DARRAY(INC_LVL_TEMP, A.nrow);
    init_value(INC_LVL_TEMP, A.nrow, BASKER_MAX_IDX);
  }

  return 0;
}//end default_symb()

  //Used with Matrix in ND-ORDERING-Assume Static-pivot
  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::symmetric_sfactor()
  {
    //Algorithm outline
    //Perform Normal col_etree on Domain
    //Perform simulated factorization on seps
    
    //Loop over all domain blocks
    //Note, that we could do this in parallel
    //(and should in future releases)
    //However parallel will require more memory
    //Need array of Strees

    //Need to make this sparse 
    //Without extra memory
    //-----OLD
    //INT_1DARRAY gScol;
    //MALLOC_INT_1DARRAY(gScol, A.ncol);
    //init_value(gScol, A.ncol, 0);
    //INT_1DARRAY gSrow;
    //MALLOC_INT_1DARRAY(gSrow, A.nrow);
    //init_value(gSrow, A.nrow, 0);

    #ifdef BASKER_TIMER 
    Kokkos::Timer timer;
    timer.reset();
    #endif

    //printf("N-lvls: %d \n", tree.nlvls);
    //printf("p/2: %d \n", num_threads/2);
    Int split_num = num_threads/2;
    if(split_num == 0)
    { 
      split_num = 1;
    }
    INT_2DARRAY gScol;
    INT_2DARRAY gSrow;
    //MALLOC_INT_2DARRAY(gScol, tree.nlvls);
    MALLOC_INT_2DARRAY(gScol, split_num);
    BASKER_ASSERT(split_num > 0, "Basker symmetric_sfactor assert: split_num > 0 failed");
    //MALLOC_INT_2DARRAY(gSrow, tree.nlvls);
    MALLOC_INT_2DARRAY(gSrow, split_num);
    split_num = num_threads/2;
    //for(Int ii=0; ii < tree.nlvls; ii++)
    for(Int ii=0; ii < split_num; ii++)
    {
      BASKER_ASSERT(A.ncol > 0, "Basker symmetric_sfactor assert: A.ncol malloc > 0 failed");
      MALLOC_INT_1DARRAY(gScol[ii], A.ncol);
      init_value(gScol[ii], A.ncol, (Int)0);
    }

 
    //INT_2DARRAY gSrow;
    //BASKER_ASSERT(split_num > 0, "sfactor split_num");
    //MALLOC_INT_2DARRAY(gSrow, tree.nlvls);
    //MALLOC_INT_2DARRAY(gSrow, split_num);
    //for(Int ii=0; ii< tree.nlvls; ii++)
    for(Int ii=0; ii < split_num; ii++)
    {
      BASKER_ASSERT(A.nrow > 0, "sfactor A.nrow malloc");
      MALLOC_INT_1DARRAY(gSrow[ii], A.nrow);
      init_value(gSrow[ii], A.nrow, (Int)0);
    }

    #ifdef BASKER_TIMER 
    std::cout << " >> symmetric_sfactor::init : " << timer.seconds() << " seconds" << std::endl;

    double time1 = 0.0;
    double time1_2 = 0.0;
    double time1_3 = 0.0;
    double time2 = 0.0;
    double time3 = 0.0;
    Kokkos::Timer timer1;
    timer.reset();
    #endif

    //split_num = num_threads/2;
    //for(Int p =0; p < 1; ++p)
    if(Options.verbose == BASKER_TRUE)
    {
      printf("\n");
      printf("\n --------------- OVER DOMS ---------------\n");
      printf("\n");
    }
    #define SHYLU_BASKER_STREE_LIST
    std::vector <BASKER_SYMBOLIC_TREE> stree_list (num_threads);
    #ifdef SHYLU_BASKER_STREE_LIST
    Kokkos::parallel_for(
      "permute_col", num_threads,
      KOKKOS_LAMBDA(const int p)
    #else
    for(Int p = 0; p < num_threads; ++p)
    #endif
    {
      Int blk = S(0)(p);
      if(Options.verbose == BASKER_TRUE)
      {
        printf(" ============= DOMAIN BLK (p=%d) ============\n",(int)p);
      }
      //printf("Sfactor blk: %d S[%d][0] \n", blk, p);

      //printf("\n\n STREE SIZE: %d \n", AL[blk][0].ncol);
      //printf("Here 0\n");
      //Find nnz_counts for leafs
      #ifdef BASKER_TIMER 
      timer1.reset();
      #endif
      #ifdef SHYLU_BASKER_STREE_LIST
      auto stree_p = stree_list[p];
      e_tree    (ALM(blk)(0), stree_p, 1);
      #else
      e_tree    (ALM(blk)(0), stree, 1);
      #endif
      #ifdef BASKER_TIMER 
      time1_2 += timer1.seconds();
      timer1.reset();
      #endif
      #ifdef SHYLU_BASKER_STREE_LIST
      post_order(ALM(blk)(0), stree_p);
      #else
      post_order(ALM(blk)(0), stree);
      #endif
      #ifdef BASKER_TIMER 
      time1_3 += timer1.seconds();
      timer1.reset();
      #endif
      #ifdef SHYLU_BASKER_STREE_LIST
      col_count (ALM(blk)(0), stree_p);
      #else
      col_count (ALM(blk)(0), stree);
      #endif
      #ifdef BASKER_TIMER 
      time1 += timer1.seconds();
      #endif

      //Assign nnz here
      //leaf_assign_nnz(LL[blk][0], stree, 0);
      //leaf_assign_nnz(LU[blk][LU_size[blk]-1], stree, 0);
      if(Options.verbose == BASKER_TRUE)
      {
        printf( " >> leaf_assign_nnz(LL(%d)(%d))\n",(int)blk,0);
        printf( " >> leaf_assign_nnz(LL(%d)(%d))\n",(int)blk,(int)LU_size(blk)-1);
      }
      #ifdef BASKER_TIMER 
      timer1.reset();
      #endif
      #ifdef SHYLU_BASKER_STREE_LIST
      leaf_assign_nnz(LL(blk)(0),              stree_p, 0);
      leaf_assign_nnz(LU(blk)(LU_size(blk)-1), stree_p, 0);
      #else
      leaf_assign_nnz(LL(blk)(0),              stree, 0);
      leaf_assign_nnz(LU(blk)(LU_size(blk)-1), stree, 0);
      #endif
      #ifdef BASKER_TIMER 
      time2 += timer1.seconds();
      #endif
    }
    #ifdef SHYLU_BASKER_STREE_LIST
    );
    Kokkos::fence();
    #endif

    for(Int p = 0; p < num_threads; ++p)
    {
      //Do off diag
      Int blk = S(0)(p);
      #ifdef SHYLU_BASKER_STREE_LIST
      auto stree_p = stree_list[p];
      #endif
      for(Int l =0; l < tree.nlvls; l++)
      {
        Int U_col = S(l+1)(p);
        //Note: Need to think more about this flow
        //Should be subtracted by how many times in the 
        //future

        Int my_row_leader = S(0)(find_leader(p,l));
        //Int my_new_row = 
        // blk - my_row_leader;
        Int U_row = blk-my_row_leader;
        Int glvl = p/2;

        #ifdef BASKER_DEBUG_SFACTOR
        printf("Proc off-diag block: %ld p: %ld loc: %ld %ld \n", 
            (long)l, (long)p, (long)U_col, (long)U_row);
        printf("BLK: %ld %ld col: %ld row: %ld \n",
            (long)U_col, (long)U_row, (long)l, (long)glvl);
        #endif

        Int off_diag = 1;
        //printf( " U_blk_sfactor(AVM(%d,%d))\n",U_col,U_row );
        //U_blk_sfactor(AV[U_col][U_row], stree,
        //		  gScol[l], gSrow[glvl],0);
        #ifdef BASKER_TIMER 
        timer1.reset();
        #endif
        #ifdef SHYLU_BASKER_STREE_LIST
        U_blk_sfactor(AVM(U_col)(U_row), stree_p,
                      gScol[l], gSrow[glvl], off_diag);
        #else
        U_blk_sfactor(AVM(U_col)(U_row), stree,
                      gScol[l], gSrow[glvl], off_diag);
        #endif
        #ifdef BASKER_TIMER 
        time3 += timer1.seconds();
        #endif

        //Determine lower blk nnz
        //Not need to be run in symmetric case
        //L_blk_factor(AL[?][?], stree);

        //Reduce all into global (don't need in serial)
        //S_sfactor_reduce(AV[U_col][U_row],
        //		     stree, gScol, gSrow);

        //Assign nnz counts for leaf off-diag
        //U_assign_nnz(LU[U_col][U_row], stree, 0);
        //L_assign_nnz(LL[blk][l+1], stree, 0);
        if(Options.verbose == BASKER_TRUE)
        {
          printf( "   ++ leaf_assign_nnz(LU(%d, %d))\n",(int)U_col,(int)U_row);
          printf( "   ++ leaf_assign_nnz(LL(%d, %d))\n",(int)blk,(int)l+1);
        }
        #ifdef BASKER_TIMER 
        timer1.reset();
        #endif
        //printf( " U_assign_nnz(LU(%d,%d))\n",U_col,U_row );
        double fill_factor = BASKER_DOM_NNZ_OVER+Options.user_fill; 
        #ifdef SHYLU_BASKER_STREE_LIST
        U_assign_nnz(LU(U_col)(U_row), stree_p, fill_factor, 0);
        L_assign_nnz(LL(blk)(l+1),     stree_p, fill_factor, 0);
        #else
        U_assign_nnz(LU(U_col)(U_row), stree, fill_factor, 0);
        L_assign_nnz(LL(blk)(l+1),     stree, fill_factor, 0);
        #endif
        #ifdef BASKER_TIMER 
        time2 += timer1.seconds();
        #endif

      }//end off diag
    }
    //over all domains
    #ifdef BASKER_TIMER 
    std::cout << " >> symmetric_sfactor::domain : " << timer.seconds() << " seconds" << std::endl;
    std::cout << "  ++ symmetric_sfactor::domain::postorder : " << time1_2 << " + " << time1_3 << " + " << time1 << " seconds" << std::endl;
    std::cout << "  ++ symmetric_sfactor::domain::init      : " << time2 << " seconds" << std::endl;
    std::cout << "  ++ symmetric_sfactor::domain::sfactor   : " << time3 << " seconds" << std::endl;
    timer.reset();
    #endif


    //do all the sep
    if(Options.verbose == BASKER_TRUE)
    {
      printf("\n");
      printf("\n --------------- OVER SEPS (nlvls=%d) ---------------\n", (int)tree.nlvls);
      printf("\n");
    }
    for(Int lvl=0; lvl < tree.nlvls; lvl++)
    {
      //Number of seps in the level
      Int p = pow(tree.nparts, tree.nlvls-lvl-1);

      //over all the seps in a lvle
      #ifdef SHYLU_BASKER_STREE_LIST
      Kokkos::parallel_for(
        "permute_col", p, KOKKOS_LAMBDA(const int pp)
      #else
      for(Int pp = 0; pp < p; pp++)
      #endif
      {
        //printf( " -- level = %d separator = %d --\n",lvl,pp );
        //S blks
        Int ppp;
        ppp =  pp*pow(tree.nparts, lvl+1);
        #ifdef BASKER_DEBUG_SFACTOR
        printf("p: %d pp: %d ppp: %d \n",
            p, pp, ppp);

        printf("\n\n");
        printf("Sep sfactor, lvl: %ld p: %ld U_col: %ld U_row: %ld \n",
            (long)lvl, (long)p, (long)U_col, (long)U_row);
        printf("BLK: %ld %ld Col: %ld Row: %ld \n",
            (long)U_col, (long)U_row, (long)lvl, (long)pp);
        #endif

        Int U_col = S(lvl+1)(ppp);
        Int U_row = 0;

        //S_blk_sfactor(AL[U_col][U_row], stree,
        //gScol[lvl], gSrow[pp]);

        //printf( " >>> S_blk_sfactor( ALM(%d)(%d) with %dx%d and nnz=%d) <<<\n",U_col,U_row, ALM(U_col)(U_row).nrow,ALM(U_col)(U_row).ncol,ALM(U_col)(U_row).nnz );
        #ifdef SHYLU_BASKER_STREE_LIST
        auto stree_p = stree_list[pp];
        S_blk_sfactor(ALM(U_col)(U_row), stree_p,
            gScol(lvl), gSrow(pp));
        #else
        S_blk_sfactor(ALM(U_col)(U_row), stree,
            gScol(lvl), gSrow(pp));
        #endif
        //printf( " >>> -> nnz = %d\n",ALM(U_col)(U_row).nnz );

        //S_assign_nnz(LL[U_col][U_row], stree, 0);
        if(Options.verbose == BASKER_TRUE)
        {
          printf( "   >>  S_assign_nnz( LL(%d,%d) )\n",(int)U_col,(int)U_row );
        }
        #ifdef SHYLU_BASKER_STREE_LIST
        S_assign_nnz(LL(U_col)(U_row), stree_p, 0);
        #else
        S_assign_nnz(LL(U_col)(U_row), stree, 0);
        #endif
        //S_assign_nnz(LU[U_col][LU_size[U_col]-1], stree,0);
        //printf( " >>>  S_assign_nnz( LU(%d,%d) )\n",U_col,LU_size(U_col)-1 );
        if(Options.verbose == BASKER_TRUE)
        {
          printf( "   ++ S_assign_nnz(LU(%d, %d))\n",(int)U_col,(int)LU_size(U_col)-1);
        }
        #ifdef SHYLU_BASKER_STREE_LIST
        S_assign_nnz(LU(U_col)(LU_size(U_col)-1), stree_p, 0);
        #else
        S_assign_nnz(LU(U_col)(LU_size(U_col)-1), stree, 0);
        #endif
      }
      #ifdef SHYLU_BASKER_STREE_LIST
      );
      Kokkos::fence();
      #endif

      for(Int pp = 0; pp < p; pp++)
      {
        #ifdef SHYLU_BASKER_STREE_LIST
        auto stree_p = stree_list[pp];
        #endif
        Int ppp;
        ppp =  pp*pow(tree.nparts, lvl+1);

        Int U_col = S(lvl+1)(ppp);
        Int U_row = 0;

        Int inner_blk = U_col;
        for(Int l = lvl+1; l < tree.nlvls; l++)
        {
          U_col = S(l+1)(ppp);
          U_row = S(lvl+1)(ppp)%LU_size(U_col);

          Int my_row_leader = S(0)(find_leader(ppp,l));
          //Int my_new_row = 
          // S(lvl+1)(ppp) - my_row_leader;
          U_row =  S(lvl+1)(ppp) - my_row_leader;

          #ifdef BASKER_DEBUG_SFACTOR
          printf("offida sep, lvl: %d l: %d U_col: %d U_row: %d \n", lvl, l, U_col, U_row);
          printf("BLK: %d %d Col: %d Row: %d \n", U_col, U_row, l, pp);
          #endif

          Int off_diag = 1;
          #ifdef SHYLU_BASKER_STREE_LIST
          U_blk_sfactor(AVM(U_col)(U_row), stree_p,
              gScol(l), gSrow(pp), off_diag);
          #else
          U_blk_sfactor(AVM(U_col)(U_row), stree,
              gScol(l), gSrow(pp), off_diag);
          #endif

          //In symmetric will not need
          //L_blk_factor(...)

          //Don't need in serial
          //S_sfactor_reduce(AV[U_col][U_row],
          //		     stree, gScol, gSrow);


          //Assign nnz
          if(Options.verbose == BASKER_TRUE)
          {
            printf( "   ++ leaf_assign_nnz(LU(%d, %d))\n",(int)U_col,(int)U_row);
            printf( "   ++ leaf_assign_nnz(LL(%d, %d))\n",(int)inner_blk,(int)(l-lvl));
          }
          double fill_factor = BASKER_SEP_NNZ_OVER+Options.user_fill; 
          #ifdef SHYLU_BASKER_STREE_LIST
          U_assign_nnz(LU(U_col)(U_row), stree_p, fill_factor, 0);
          L_assign_nnz(LL(inner_blk)(l-lvl), stree_p, fill_factor, 0);
          #else
          U_assign_nnz(LU(U_col)(U_row), stree, fill_factor, 0);
          L_assign_nnz(LL(inner_blk)(l-lvl), stree, fill_factor, 0);
          #endif
          //printf("Here 1 \n");
        }
      }
    } //end for lvl
    #ifdef BASKER_TIMER 
    std::cout << " >> symmetric_sfactor::separators : " << timer.seconds() << " seconds" << std::endl;
    #endif

    //free my memory
    for(Int ii = 0 ; ii < split_num; ++ii)
    {
      //printf("split\n");
      FREE(gScol[ii]);
      FREE(gSrow[ii]);
    }
    FREE(gScol);
    FREE(gSrow);

    return 0;
  }//end symmetric_symbolic()


  //Usd with matrix in HUND ordering
  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::unsymmetric_sfactor()
  {
    //Algorithm outline
    //Perform QR on Domains
    //Update sep 
    return -1;
  }//ed unsymmetric_symbolic
  

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::e_tree
  (
   BASKER_MATRIX &M,
   BASKER_SYMBOLIC_TREE &ST,
   Int ata_option
  )
  {
  #ifdef BASKER_DEBUG_SFACTOR
    printf("In E_tree\n");
  #endif
    BASKER_MATRIX *MV = &M;
  
    if((Options.symmetric == BASKER_TRUE))
    {
    #ifdef BASKER_DEBUG_SFACTOR
      printf("symmetric\n");
    #endif
      ata_option = 0;
    }
    else
    {
      //if(Options.AtA == BASKER_TRUE)
      {
    #ifdef BASKER_DEBUG_SFACTOR
        printf("AtransA option\n");
    #endif
        ata_option = 1;
      }
      /*
         else
         {
         ADD BACK IN ONCE YOU FIX A+A'
         printf("Aplus option \n");
         ata_option = 0;
         BASKER_MATRIX MVT;
         AplusAT(M,MVT);
         MV = &MVT; //need to make this a pointer
         }
       */
    }
  
    ST.init_parent(MV->ncol+1);

    // quick return
    if (MV->ncol <= 0) {
      return;
    }

    //Int brow = MV->srow; //Not used
    //Int bcol = MV->scol; //Not used
    INT_1DARRAY ws;
    Int ws_size = 2*MV->nrow;
  
    Int have_past = BASKER_MAX_IDX;
  
    BASKER_ASSERT(ws_size > 0 , "Basker unsymmetric_sfactor assert: sfactor ws_size > 0 failed");
    MALLOC_INT_1DARRAY(ws, ws_size);
  
    Int *past   = &(ws(0));
    Int *clique = &(ws(MV->ncol));
  
    //Zero out the cliques from ATA
    if(ata_option == 1)
    {
      for(Int ii = 0 ; ii < MV->nrow; ii++)
      {clique[ii] = BASKER_MAX_IDX;}
      //{clique[ii] = A.max_idx;}
    }
  
    // for each column
    for(Int k = 0; k < MV->ncol; ++k)
    {
      ST.parent[k] = BASKER_MAX_IDX;
      past[k]      = BASKER_MAX_IDX;
      for(Int j = MV->col_ptr(k); j < MV->col_ptr(k+1); ++j)
      {
        Int t = MV->row_idx(j);
  
        if(ata_option==1)
        {
          t = clique[MV->row_idx(j)];
        } //use clique instead
  
        for(; (t!=BASKER_MAX_IDX) && (t < k); t = have_past)
        {
          have_past  = past[t];
          past[t]    = k;

          if(have_past == BASKER_MAX_IDX)
          {
            ST.parent(t) = k;
          }
        }//climb up the tee until root and connect

        //connect others in clique
        if(ata_option == 1)
        {
          clique[MV->row_idx(j)] = k;
        }
      } 
    }//end over all columns
  
    FREE(ws);
  }//end e_tree()


  //Finds Elimination tree (Converted to use matrix_view so can be done ||)
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::e_tree
  (
   BASKER_MATRIX_VIEW &MV,
   BASKER_SYMBOLIC_TREE &ST,
   Int ata_option
  )
  {
    if((Options.symmetric == BASKER_TRUE))
    {
      ata_option = 0;
    }
    else
    {
      if(Options.AtA == BASKER_TRUE)
      {
        ata_option = 1;
      }
      else
      {
        ata_option = 0;
        //AplusAT(MV,MVT);
        //MV = MVT; //need to make this a pointer
      }
    }

    ST.init_parent(MV.ncol);

    Int brow = MV.srow;
    Int bcol = MV.scol;
    INT_1DARRAY ws;
    Int ws_size = MV.nrow;
    if(ata_option == 1)
    {ws_size += MV.nrow;}

    BASKER_ASSERT(ws_size > 0 , "Basker e_tree assert: sfactor ws_size > 0 failed");
    MALLOC_INT_1DARRAY(ws, ws_size);

    Int *past   = &(ws[0]);
    Int *clique = &(ws[MV.ncol]);

    //Zero out the cliques from ATA
    if(ata_option == 1)
    {
      for(Int ii = 0 ; ii < MV.nrow; ii++)
      {clique[ii] = BASKER_MAX_IDX;}
    }
    // for each column
    for(Int k = 0; k < MV.ncol; k++)
    {
      ST.parent[k] = BASKER_MAX_IDX;
      past[k]  = BASKER_MAX_IDX;
      for(Int j = MV.col_ptr(k+bcol); j < MV.col_ptr(k+1+bcol); j++)
      {
        if(MV.good(j) != 0)
        {continue;}

        Int t = MV.row_idx(j)-brow;
        if(ata_option==1)
        {t = clique[MV.row_idx(j)-brow];} //use clique instead
        Int have_past;
        for(; (t!=BASKER_MAX_IDX) && (t < k); t = have_past)
        {
          have_past = past[t];
          past[t] = k;
          if(have_past == BASKER_MAX_IDX)
          {
            ST.parent[t] = k;
          }
        }//climb up the tee until root and connect

        //connect others in clique
        if(ata_option == 1)
        { clique[MV.row_idx(j)-brow] = k; }
      } 
    }//end over all columns

    FREE(ws);
  }//end e_tree()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::post_order
  (
   BASKER_MATRIX &MV,
   BASKER_SYMBOLIC_TREE &ST
  )
  {
    const Int n = MV.ncol;
    if (n <= 0) {
      return;
    }

    INT_1DARRAY post;
    BASKER_ASSERT(n > 0, "Basker post_order assert: sfactor post_order n > 0 failed");
    MALLOC_INT_1DARRAY(post,n);
    init_value(post, n, (Int)0);
    Int *p = &(ST.parent(0)); //parent;
    const Int ws_size = 3*n;
    INT_1DARRAY ws;
    MALLOC_INT_1DARRAY(ws, ws_size);
    init_value(ws, ws_size, (Int)0);
    Int *head, *next, *stack;
    head  =  &(ws(0)); 
    next  =  &(ws(n)); 
    stack =  &(next[n]);
    for(Int j=0; j < n; j++)
      {head[j] = BASKER_MAX_IDX;}
  
    for(Int k = n; k > 0 ; k--)
    {
      Int j = k-1;
      if(p[j] == BASKER_MAX_IDX) 
      { continue; }

      next[j] = head[p[j]];
      head[p[j]] = j;
    }
    Int k = 0;
    for(Int j=0; j<MV.ncol; j++)
    {
      if(p[j] != BASKER_MAX_IDX)
      { continue; }
      k = post_dfs(j, k, head, next, post, stack);
    }
    FREE(ws);
    
    //Come back and make smaller
    //ST.init_post(n);
    ST.init_post(A.ncol);
    //Change quick fix
    for(Int i = 0; i < n; i++)
    {
      ST.post[i] = post[i];
    }

    FREE(post);
  }//end post_order()


  // Post orders a parent array representation of a tree
  //need to do something about workspace
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::post_order
  (
   BASKER_MATRIX_VIEW &MV,
   BASKER_SYMBOLIC_TREE &ST
  )
  {
    Int n = MV.ncol;
    INT_1DARRAY post;
    BASKER_ASSERT(n > 0, "Basker post_order assert: sfactor post_order n > 0 failed");
    MALLOC_INT_1DARRAY(post,n);
    init_value(post, n, (Int)0);
    Int *p = &(ST.parent[0]); //parent;
    Int ws_size = 3*n;
    INT_1DARRAY ws;
    MALLOC_INT_1DARRAY(ws, ws_size);
    init_value(ws, ws_size, (Int)0);
    Int *head, *next, *stack;
    head  =  &(ws[0]); 
    next  =  &(ws[n]); 
    stack =  &(next[n]);
    for(Int j=0; j < n; j++)
      {head[j] = BASKER_MAX_IDX;}
  
    for(Int k = n; k > 0 ; k--)
    {
      Int j = k-1;
      if(p[j] == BASKER_MAX_IDX)
      { continue; }
      next[j] = head[p[j]];
      head[p[j]] = j;
    }
    Int k = 0;
    for(Int j=0; j<MV.ncol; j++)
    {
      if(p[j] != BASKER_MAX_IDX)
      { continue; }

      k = post_dfs(j, k, head, next, post, stack);
    }
    FREE(ws);

    //Come back and make smaller
    //ST.init_post(n);
    ST.init_post(A.ncol);
    //Change quick fix
    for(Int i = 0; i < n; i++)
    {
      ST.post[i] = post[i];
    }
    FREE(post);
  }//end post_order()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int Basker<Int, Entry, Exe_Space>::post_dfs
  (
   Int j, Int k,
   Int *head, 
   Int *next,
   //Int *post,
   INT_1DARRAY post,
   Int *stack
  )
  {
    Int ii, jj, top;
    top = 0;
    //Start from j
    stack[0] = j;
    while(top >= 0)
    {
      jj = stack[top];
      ii = head[jj];
      if(ii == BASKER_MAX_IDX)
      {
        top--;
        post[k++] = jj;
      }
      else
      {
        head[jj] = next[ii];
        stack[++top] = ii;
      }
    }
    return k;
  }//end post_dfs()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::col_count
  (
   BASKER_MATRIX &MV,
   BASKER_SYMBOLIC_TREE &ST
   )
  {
    if (MV.ncol == 0) {
      return;
    }

    //Still like to find a way to do this without transpose
    BASKER_MATRIX  Mt;
    matrix_transpose(MV, Mt);
    Int *post   = &(ST.post(0));
    Int *parent = &(ST.parent(0));
   
    const Int ws_size = 4*MV.ncol+ (MV.ncol+MV.nrow+1);
    INT_1DARRAY ws;
    BASKER_ASSERT(ws_size > 0, "Basker col_count assert1: col_count ws_size > 0 failed");
    MALLOC_INT_1DARRAY(ws, ws_size);
    INT_1DARRAY delta;
    BASKER_ASSERT(MV.ncol > 0, "Basker col_count assert2: ncol > 0 failed");
    MALLOC_INT_1DARRAY(delta, MV.ncol);
    //init_value(delta, MV.ncol, (Int)0);
    //init_value(ws, ws_size, (Int)0);

    Int jleaf=0;
    Int *past, *mfirst, *pleaf, *first;
    past   = &(ws(0));
    mfirst = &(ws(MV.ncol));
    pleaf  = &(ws(MV.ncol+MV.ncol));
    first  = &(ws(MV.ncol+MV.ncol+MV.ncol));

    for(Int k = 0; k < ws_size; k++)
      {ws(k) = BASKER_MAX_IDX;}
   
    for(Int k = 0; k < MV.ncol; k++)
    {
      //Leaving post out
      Int j = post[k];

      if(first[j] == BASKER_MAX_IDX)
      {
        delta[j] = 1; // leaf node
      }
      else
      {
        delta[j] = 0;
      }
      for( ; 
          (j != BASKER_MAX_IDX) && 
          (first[j] == BASKER_MAX_IDX); 
          j = parent[j])
      {
        first[j] = k; // update with parent
      }
    }//initalize the delta counts for overlap

    // Create a linked list of the cliques
    //Cliques are need for nonsymmtrix A'A case
    //Int *head = ws+4*A.ncol;
    //Int *next = ws+5*A.ncol+1;
    Int *head, *next;
    head = &(ws(4*MV.ncol));
    next = &(ws(5*MV.ncol+1));

    if(Options.symmetric == BASKER_FALSE)
    {
    #ifdef BASKER_DEBUG_SFACTOR
      printf("\n\n\n\n");
      printf("NONSYM-SFACTOR");
      printf("\n\n\n\n");
    #endif

      for(Int k=0; k < MV.ncol; k++)
      {ws(post[k]) = k;} //overwriting past

      for(Int i = 0; i < MV.nrow; i++)
      {
        Int k=MV.ncol;
        for(Int p = Mt.col_ptr(i); p < Mt.col_ptr(i+1); ++p)
        {
          k = min(k,ws(Mt.row_idx(p)));
        }

        next[i] = head[k];
        head[k] = i;
      }
      // End create a linked list of the cliques

    }
    // reset past
    for(Int k = 0; k < MV.ncol; k++)
    {past[k] = k;}
    
    for(Int k = 0; k < MV.ncol; k++)
    {
      Int j = post[k];
      if(parent[j] != BASKER_MAX_IDX)
      {
        delta[parent[j]]--;
      }  // j is not a root node
      //loop over clique; In A this would only be done once

      /*
        for(Int J = ((Options.symmetric) ? j : head[k]);
      //J != A.max_idx;
      J != BASKER_MAX_IDX;
      //J = (Options.symmetric) ? A.max_idx : next[J])
      J = (Options.symmetric) ? BASKER_MAX_IDX: next[J])
      */
      for(Int J = ((Options.symmetric) ? j : head[k]);
              J != BASKER_MAX_IDX;
              J = ((Options.symmetric) ? BASKER_MAX_IDX : next[J]))
      {
        for(Int p = Mt.col_ptr(J); p < Mt.col_ptr(J+1); ++p)
        {
          Int i = Mt.row_idx(p);

          Int q = least_common(i, j, first, mfirst, pleaf, past, &jleaf);
          if(jleaf >=1) 
          {delta[j]++;}
          if(jleaf == 2)
          {delta[q]--;}
        }//for over row/col
      }//for over elements in clique
      //if(parent[j] != A.max_idx)

      if(parent[j] != BASKER_MAX_IDX)
      {past[j] = parent[j];}
    }//over all col/row

    for(Int k = 0; k < MV.ncol; k++)
    {
      if(parent[k] != BASKER_MAX_IDX)
      {
        delta[parent[k]] += delta[k];
      }
    }
    
    // Clean up AT
   
    //Comeback and make less later
    ST.init_col_counts(A.ncol);

    //copy column counts.
    for(Int i = 0; i < MV.ncol; i++)
    {
      ST.col_counts[i] = delta[i];
    }
   
    // Clean up workspace
    FREE(ws);
    FREE(delta);

  }//end col_count()
 

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::col_count
  ( 
   BASKER_MATRIX_VIEW &MV,
   BASKER_SYMBOLIC_TREE &ST
  )
  {
printf( " col_count:: view \n" );
    //Still like to find a way to do this without transpose
    BASKER_MATRIX  Mt;
    matrix_transpose(MV, Mt);
    Int *post = &(ST.post[0]);
    Int *parent = &(ST.parent[0]);
    Int ws_size = 4*A.ncol+ (A.ncol+A.nrow+1);
    INT_1DARRAY ws;
    BASKER_ASSERT(ws_size > 0, "Basker col_count assert1: ws_size > 0 failed");
    MALLOC_INT_1DARRAY(ws, ws_size);
    //init_value(ws, ws_size, (Int)0);
    //for(Int k = 0; k < ws_size; k++)
    //{ ws[k] = BASKER_MAX_IDX; }
    Kokkos::deep_copy(ws, BASKER_MAX_IDX);

    INT_1DARRAY delta;
    BASKER_ASSERT(MV.ncol > 0, "Basker col_count assert2: ncol > 0 failed");
    MALLOC_INT_1DARRAY(delta, MV.ncol);
    init_value(delta, MV.ncol, (Int)0);

    Int jleaf=0;
    Int *past, *mfirst, *pleaf, *first;
    past   = &(ws[0]);
    mfirst = &(ws[MV.ncol]);
    pleaf  = &(ws[MV.ncol+MV.ncol]);
    first  = &(ws[MV.ncol+MV.ncol+MV.ncol]);
    
    
    for(Int k = 0; k < MV.ncol; k++)
    {
      //Leaving post out
      Int j = post[k];

      if(first[j] == BASKER_MAX_IDX)
      {
        delta[j] = 1; // leaf node
      }
      else
      {
        delta[j] = 0;
      }

      for( ; 
          (j != BASKER_MAX_IDX) && 
          (first[j] == BASKER_MAX_IDX); 
          j = parent[j])
      {
        first[j] = k; // update with parent
      }
    } //initalize the delta counts for overlap
  
    // Create a linked list of the cliques
    // Cliques are need for nonsymmtrix A'A case
    Int *head = &(ws[4*MV.ncol]);
    Int *next = &(ws[5*MV.ncol+1]);
    
    if(Options.symmetric == BASKER_FALSE)
    {
      printf("\n\n\n\n");
      printf("NONSYM-SFACTOR");
      printf("\n\n\n\n");

      for(Int k=0; k < MV.ncol; k++)
      { ws[post[k]] = k; } //overwriting past
      for(Int i = 0; i < MV.nrow; i++)
      {
        Int k=MV.ncol;
        for(Int p = Mt.col_ptr[i]; p < Mt.col_ptr[i+1]; p++)
        {
          k = min(k,ws[Mt.row_idx[p]]);
        }

        next[i] = head[k];
        head[k] = i;
      }
      // End create a linked list of the cliques
    }
    // reset past
    for(Int k = 0; k < MV.ncol; k++)
    { past[k] = k; }
    
    for(Int k = 0; k < MV.ncol; k++)
    {
      Int j = post[k];
      if(parent[j] != BASKER_MAX_IDX)
      {
        delta[parent[j]]--;
      }  // j is not a root node
      // loop over clique; In A this would only be done once

      /*
       for(Int J = ((Options.symmetric) ? j : head[k]);
      //J != A.max_idx;
      J != BASKER_MAX_IDX;
      //J = (Options.symmetric) ? A.max_idx : next[J])
      J = (Options.symmetric) ? BASKER_MAX_IDX: next[J])
      */

      for(Int J = ((Options.symmetric) ? j : head[k]);
          J != BASKER_MAX_IDX;
          J = ((Options.symmetric) ? BASKER_MAX_IDX:next[J]))
      {
        for(Int p = Mt.col_ptr[J]; p < Mt.col_ptr[J+1]; p++)
        {
          Int i = Mt.row_idx[p];
          Int q = least_common(i, j, first, mfirst, pleaf, past, &jleaf);

          if(jleaf >= 1) 
          { delta[j]++; }
          if(jleaf == 2)
          { delta[q]--; }
        }//for over row/col
      }//for over elements in clique

      if(parent[j] != BASKER_MAX_IDX)
      { past[j] = parent[j]; }

    }//over all col/row

    for(Int k = 0; k < MV.ncol; k++)
    {
      if(parent[k] != BASKER_MAX_IDX)
      {
        delta[parent[k]] += delta[k];
      }
    }
    
    // Clean up AT
   
    //Comeback and make less later
    ST.init_col_counts(A.ncol);
    //copy column counts.
    for(Int i = 0; i < MV.ncol; i++)
    {
      ST.col_counts[i] = delta[i];
    }

    // Clean up workspace 
    FREE(ws);
    FREE(delta);
  }//end col_count()

  
  //Moved to basker_util
  //template <class Int, class Entry, class Exe_Space>
  //BASKER_INLINE
  //void Basker<Int,Entry,Exe_Space>::matrix_transpose(
  //                                  BASKER_MATRIX_VIEW &MV,
  //                                       BASKER_MATRIX &AT)
  //{
  //
  //Int brow = MV.srow;
  ////Setup what we do know
  //AT.srow = MV.srow;
  //AT.nrow = MV.nrow;
  //AT.scol = MV.scol;
  //AT.ncol = MV.ncol;
  //AT.nnz  = MV.nnz();
  //MALLOC_INT_1DARRAY(AT.col_ptr, AT.ncol+1);
  //init_value(AT.col_ptr, AT.ncol+1, (Int)0);
  //MALLOC_INT_1DARRAY(AT.row_idx, AT.nnz);
  //init_value(AT.row_idx, AT.nnz, (Int)0);
  //MALLOC_ENTRY_1DARRAY(AT.val    , AT.nnz);
  //init_value(AT.val,     AT.nnz, (Entry)0.0);
  //
  ////Setup a litte workspace
  //Int ws_size = MV.nrow;
  //INT_1DARRAY ws;
  //MALLOC_INT_1DARRAY(ws, ws_size);
  //init_value(ws, ws_size, (Int)0);
  //
  ////Note could get number of nonzeros here inplace of nnz() for fas//ter
  //
  ////get row counts
  //for(Int j = MV.col_ptr(MV.scol); j < MV.col_ptr(MV.scol+MV.ncol); j++)
  //  {
  //    if(MV.good(j) != 0)
  //      {
  //        continue;
  //      }
  //    ws[MV.row_idx(j)-brow]++;
  //  }
  //
  //AT.col_ptr[1] = ws[0];
  //for(Int j = 1; j < AT.nrow; j++)
  //  {
  //    ws[j] = ws[j]+ws[j-1];
  //    AT.col_ptr[j+1] = ws[j];
  //    ws[j-1] = AT.col_ptr[j-1];
  //  }
  //ws[AT.nrow-1] = AT.col_ptr[AT.nrow-1];
  //
  //for(Int k = 0; k < AT.ncol; k++)
  //  {
  //    for(Int j = MV.col_ptr(MV.scol+k); j < MV.col_ptr(MV.scol+k+1); j++)
  //      {
  //        if(MV.good(j) != 0)
  //          {
  //            continue;
  //          }
  //        
  //        AT.row_idx[ws[MV.row_idx(j)-brow]++] = k; //starting at zero
  //        //AT.row_idx[ws[MV.row_idx(j)-brow]++] = k+brow; //not o
  //      }
  //  }
  //FREE(ws);
  //}//end matrix_transpose

   template <class Int, class Entry, class Exe_Space>
  Int Basker<Int, Entry, Exe_Space>::least_common
  (
   Int i, Int j, Int *first, 
   Int *mfirst,Int *pleaf, 
   Int *past, Int *jleaf
  )
  {
    Int q, sparent;
    *jleaf = 0;
    //j not a leaf
    
    if((i <=j) || ((mfirst[i]!= BASKER_MAX_IDX)&& (first[j] <= mfirst[i])))
    {
      return BASKER_MAX_IDX;
    }

    mfirst[i] = first[j];
    Int jprev = pleaf[i];
    pleaf[i] = j;

    if(jprev == BASKER_MAX_IDX)
    { *jleaf = 1; }
    else
    { *jleaf = 2;}

    if(*jleaf == 1) 
    { return i; }

    for(q = jprev; q != past[q]; q = past[q]);

    for(Int k = jprev; k != q; k = sparent)
    {
      sparent = past[k];
      past[k] = q;
    }
    return q;      
  }//end least_common() 


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::U_blk_sfactor
  (
   BASKER_MATRIX &MV,
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol,
   INT_1DARRAY grow, 
   Int off_diag
  )
  {
    //printf("-----------UBLDSFACTO-------------\n");
    //Algorithm
    //Goal Reuse as much of LU etree as possible
    //Order A(:,k) in post order (Assume already)
    //Take smalled in post order, climb tree (Assume already)
    //Note----this post ordering can be done with MV.init_perm
    //Mark along tree.
    //repeat for all .... this should be slightly faster DFS

    //Looking forward.. we might want to save only the first
    //col pattern for Skyline model of S

    //Also could use some form of path compression

    //Set Column offset = MV.scol
    //Set Row-idx offset = 0;


    #ifdef BASKER_DEBUG_SFACTOR
    printf("\n\n");
    printf("BLK: %d %d \n", brow, bcol);
    printf("row: %d %d col: %d %d \n",
	   MV.srow, MV.nrow+MV.srow, 
	   MV.scol, MV.ncol+MV.scol);
    printf("\n\n");
    #endif

    //Start Algorithm
    #ifdef BASKER_TIMER 
    Kokkos::Timer timer;
    timer.reset();
    #endif
    //BASKER_ASSERT(MV.ncol > 0, "Basker U_blk_sfactor assert: ncol > 0 failed");
    INT_1DARRAY U_col_count;
    if(MV.ncol > 0)
    {
      MALLOC_INT_1DARRAY(U_col_count, MV.ncol);
      init_value(U_col_count, MV.ncol, (Int)0);
    }
    //May want to change color to bool
    //BASKER_ASSERT(MV.nrow > 0, "Basker U_blk_sfactor assert: nrow > 0 failed");
    INT_1DARRAY color;
    if(MV.nrow > 0)
    {
      MALLOC_INT_1DARRAY(color, MV.nrow);
      init_value(color, MV.nrow, (Int)0);
    }
    INT_1DARRAY pattern;
    if(MV.nrow > 0)
    {
      MALLOC_INT_1DARRAY(pattern, MV.nrow);
      init_value(pattern, MV.nrow, (Int)0);
    }

    //Waveform if symmetric
    INT_1DARRAY wave;
    if(MV.nrow > 0)
    {
      MALLOC_INT_1DARRAY(wave, MV.nrow);
      init_value(wave, MV.nrow, BASKER_MAX_IDX);
    }
    INT_1DARRAY wave_p;
    if(MV.nrow > 0)
    {
      MALLOC_INT_1DARRAY(wave_p, MV.nrow);
      init_value(wave_p, MV.nrow, BASKER_MAX_IDX);
    }

    //If symmetric (short cutting for S)
    //First row of S
    INT_1DARRAY first_color;
    if(MV.nrow > 0)
    {
      MALLOC_INT_1DARRAY(first_color, MV.nrow);
      init_value(first_color, MV.nrow, (Int)0);
    }
    INT_1DARRAY first_row;
    if(MV.ncol > 0)
    {
      MALLOC_INT_1DARRAY(first_row, MV.ncol);
      init_value(first_row, MV.ncol, (Int)0);
    }
    #ifdef BASKER_TIMER 
    std::cout << " >> U_blk_sfactor::init : " << timer.seconds() << " seconds" << std::endl;
    timer.reset();
    #endif
    //Could add first col if not symmetric

    // If symmetric  (short cutting for off-diag S)
#if 0
    INT_1DARRAY max_reach;
    //if(MV.ncol > 0)
    {
      MALLOC_INT_1DARRAY(max_reach, MV.ncol);
      init_value(max_reach, MV.ncol, BASKER_MAX_IDX);
    }
    //printf("-------------UBLK DONE ALLOC----");

    #ifdef BASKER_TIMER 
    Kokkos::Timer timer1;
    double time1 = 0.0;
    double time2 = 0.0;
    #endif

    //???Fix??
    // MV.init_offset(MV.scol,0);
    //Int brow = MV.srow; //Not used
    Int bcol = MV.scol;
 
    //Loop of all columns
    Int top = 0;
    for(Int k = 0; k < MV.ncol; ++k)
    {
      //Add any offdiag reach that might have already been
      if(off_diag == 1)
      {
        //printf("---OFF DIAG U-nnz called ----\n");
        Int t = grow(k+bcol);
        if((t != BASKER_MAX_IDX) && (t < MV.nrow))
        {
          while((t != BASKER_MAX_IDX) && (t <= MV.ncol) && (color(t) == 0))
          {
            U_col_count(k)++;
            color(t) = 1;
            pattern(top++) = t;
            t = ST.parent(t);
          }
        }
      }

      //Loop over rows
      #ifdef BASKER_TIMER 
      timer1.reset();
      #endif
      for(Int j = MV.col_ptr(k); j < MV.col_ptr(k+1); ++j)
      {
        //Climb tree
        //Will want to modify this to ST.post[row_idx(j)];
        Int t = MV.row_idx(j);

        // Processing element t
        while((t != BASKER_MAX_IDX) && (t <= MV.nrow) && (color(t) == 0))
        {
          U_col_count(k)++;
          color(t) = 1;
          pattern(top++) = t;
          t = ST.parent(t);
        }
      }// end for j
      #ifdef BASKER_TIMER 
      time1 += timer1.seconds();
      timer1.reset();
      #endif

      //clear color
      Int min_pos = max_reach(k);
      for(Int ii = 0; ii < top; ii++)
      {
        Int kk = pattern(ii);

        //if symmetric (shortcut for S)
        if(k == 0)
        {
          first_color(kk) = 1;
        }
        else
        {
          if(first_color(kk) == 1)
          {
            first_row(k) = 1;
          }
        }

        //Move into waveform
        if(wave_p(kk) == BASKER_MAX_IDX)
        {
          wave_p(kk) = k;
        }
        else
        {
          if(min_pos > wave_p(kk))
          {
            min_pos = wave_p(kk);
          }
        }

        //maybe own row
        if(kk < min_pos)
        {
          min_pos = k;
        }

        /*
         for(Int g=k-1; g>=0; g--)
         {
          //this is wrong, need to know pattern
          // cna use etree to fix cheaply?
          if(kk == max_reach[g])
          {
            min_pos = g;
          }
        }
        */

        pattern(ii) = 0;
        color(kk) = 0;
      }
      top = 0;
      max_reach(k) = min_pos;
      min_pos = MV.ncol;
      #ifdef BASKER_TIMER 
      time2 += timer1.seconds();
      #endif
    }//all columns

    //Copy into global 
    //col
    for(Int i = MV.scol; i < MV.scol+MV.ncol; i++)
    {
      gcol(i) += first_row(i-bcol);
    }

    //row
    //for(Int i = 0; i < MV.ncol; i++)
    //{
      //printf("%d %d \n", i, max_reach[i]);
    //}
    for(Int i = MV.scol; i < MV.scol+MV.ncol; i++)
    {
      Int l_min = max_reach(i-bcol);
      //printf("Test maxreach[%d]: %d grow: %d\n",
      //   i-bcol, l_min, grow[i]);
      //printf("Acessing: %d A.ncol: %d \n", i, A.ncol);
      if(l_min < grow(i))
      {
        grow(i) = l_min;
      }
    }
    
    //Copy column_counts
    ST.init_U_col_counts(MV.ncol);
    for(Int i = 0; i < MV.ncol; i++)
    {
      ST.U_col_counts(i) = U_col_count(i);
      //temp fix for elbow row
      //ST.U_col_counts[i] += (.4)*U_col_count[i];
    }

    //if symmetric
    ST.init_L_row_counts(MV.ncol);
    for(Int i = 0; i < MV.ncol; i++)
    {
      ST.L_row_counts(i) = U_col_count(i);
      //temp fix for elbow row
      //ST.L_row_counts[i] += (.4)*U_col_count[i];
    }
    #ifdef BASKER_TIMER 
    std::cout << " >> U_blk_sfactor::loop-columns               : " << timer.seconds() << " seconds" << std::endl;
    std::cout << " >> + U_blk_sfactor::loop-columns:loop-row    : " << time1           << " seconds" << std::endl;
    std::cout << " >> + U_blk_sfactor::loop-columns:clear-color : " << time2           << " seconds" << std::endl;
    timer.reset();
    #endif
#else
    // let me just count nnz in L and U
    // NOTE: a better way to avoid realloc
    ST.init_U_col_counts(MV.ncol, Options.verbose); // contain nnz / col in U
    ST.init_L_row_counts(MV.ncol, Options.verbose); // contain nnz / row in L
    if(off_diag != 1) {
      for(Int j = 0; j < MV.ncol; j++)
      {
        ST.L_row_counts(j) = 0;
        ST.U_col_counts(j) = 0;
        for(Int k = MV.col_ptr(j); k < MV.col_ptr(j+1); ++k) {
          Int i = MV.row_idx(k);
          if (i < j) {
            ST.U_col_counts(j) ++;
          } else if (i > j) {
            ST.L_row_counts(i) ++;
          } else {
            ST.U_col_counts(j) ++;
            ST.L_row_counts(i) ++;
          }
        }
      }
    }
    #ifdef BASKER_TIMER 
    std::cout << " >> U_blk_sfactor::loop-columns               : " << timer.seconds() << " seconds" << std::endl;
    #endif
#endif

    //Temp Patch fix
    //Note Comebaske
    if(off_diag == 1)
    {
      for(Int i = 0; i < MV.ncol; i++)
      {
        ST.U_col_counts(i) = MV.nrow;
        ST.L_row_counts(i) = MV.nrow;
      }
    }
    #ifdef BASKER_TIMER 
    std::cout << " >> U_blk_sfactor::copy : " << timer.seconds() << " seconds" << std::endl;
    timer.reset();
    #endif

    FREE(U_col_count);
    FREE(color);
    FREE(pattern);

    //if symmetric
    FREE(first_color);
    FREE(first_row);

  }//end U_blk_sfactor()


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::U_blk_sfactor
  (
   BASKER_MATRIX_VIEW &MV,
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol,
   INT_1DARRAY grow, 
   Int off_diag
  )
  {
    if(MV.ncol <= 0)
    {
      return;
    }
    
    //Algorithm
    //Goal Reuse as much of LU etree as possible
    //Order A(:,k) in post order (Assume already)
    //Take smalled in post order, climb tree (Assume already)
    //Note----this post ordering can be done with MV.init_perm
    //Mark along tree.
    //repeat for all .... this should be slightly faster DFS

    //Looking forward.. we might want to save only the first
    //col pattern for Skyline model of S

    //Also could use some form of path compression

    //Set Column offset = MV.scol
    //Set Row-idx offset = 0;
    MV.init_offset(MV.scol,0);
    Int brow = MV.srow;
    Int bcol = MV.scol;
 
    //Start Algorithm
    INT_1DARRAY U_col_count;
    BASKER_ASSERT(MV.ncol > 0, "Basker U_blk_sfactor assert: ncol > 0 failed");
    //if(MV.ncol > 0)
    {
      MALLOC_INT_1DARRAY(U_col_count, MV.ncol);
      init_value(U_col_count, MV.ncol, (Int)0);
    }
    //May want to change color to bool
    INT_1DARRAY color;
    BASKER_ASSERT(MV.nrow > 0, "Basker U_blk_sfactor assert: nrow > 0 failed");
    MALLOC_INT_1DARRAY(color, MV.nrow);
    init_value(color, MV.nrow, (Int)0);
    INT_1DARRAY pattern;
    MALLOC_INT_1DARRAY(pattern, MV.nrow);
    init_value(pattern, MV.nrow, (Int)0);
    Int top = 0;


    //Waveform if symmetric
    INT_1DARRAY wave;
    MALLOC_INT_1DARRAY(wave, MV.nrow);
    //init_value(wave, MV.nrow, (Int)A.max_idx);
    init_value(wave, MV.nrow, BASKER_MAX_IDX);
    INT_1DARRAY wave_p;
    MALLOC_INT_1DARRAY(wave_p, MV.nrow);
    //init_value(wave_p, MV.nrow, (Int)A.max_idx);
    init_value(wave_p, MV.nrow, BASKER_MAX_IDX);


    //If symmetric (short cutting for S)
    //First row of S
    INT_1DARRAY first_color;
    MALLOC_INT_1DARRAY(first_color, MV.nrow);
    init_value(first_color, MV.nrow, (Int)0);
    INT_1DARRAY first_row;
    MALLOC_INT_1DARRAY(first_row, MV.ncol);
    init_value(first_row, MV.ncol, (Int)0);
    //Could add first col if not symmetric


    // If symmetric  (short cutting for off-diag S)
    INT_1DARRAY max_reach;
    MALLOC_INT_1DARRAY(max_reach, MV.ncol);
    //init_value(max_reach, MV.ncol, A.max_idx);
    init_value(max_reach, MV.ncol, BASKER_MAX_IDX);


    //Loop of all columns
    for(Int k = 0; k < MV.ncol; k++)
    {
      //Add any offdiag reach that might have already been
      if(off_diag == 1)
      {
        //printf("---OFF DIAG U-nnz called ----\n");
        //Int t = gcol[k+bcol];
        Int t = grow[k+bcol];
        if(t < MV.nrow)
        {
          //t = t-brow;
          while((t <= MV.ncol)&&(color[t]==0))
          {
            U_col_count[k]++;
            color[t] = 1;
            pattern[top++] = t;
            t = ST.parent[t];
          }
        }
      }

      //Loop over rows
      for(Int j = MV.col_ptr(k+bcol); j < MV.col_ptr(k+1+bcol); j++)
      {
        if(MV.good(j) == BASKER_MAX_IDX)
        {
          continue;
        }

        //Climb tree
        //Will want to modify this to ST.post[row_idx(j)];
        Int t = MV.row_idx(j)-brow; 
        // Processing element t
        while((t<=MV.nrow)&&(color[t]==0))
        {
          U_col_count[k]++;
          color[t] = 1;
          pattern[top++] = t;
          t = ST.parent[t];
        }
      }// end for j

      //clear color
      Int min_pos = max_reach[k];
      for(Int ii = 0; ii < top; ii++)
      {
        Int kk = pattern[ii];

        //if symmetric (shortcut for S)
        if(k == 0)
        {
          first_color[kk] = 1;
        }
        else
        {
          if(first_color[kk] == 1)
          {
            first_row[k] = 1;
          }
        }

        //Move into waveform
        if(wave_p[kk] ==BASKER_MAX_IDX)
        {
          wave_p[kk] = k;
        }
        else
        {
          if(min_pos > wave_p[kk])
          {
            min_pos = wave_p[kk];
          }

        }

        //maybe own row
        if(kk < min_pos)
        {
          min_pos = k;
        }

        /*
         for(Int g=k-1; g>=0; g--)
         {
          //this is wrong, need to know pattern
          // cna use etree to fix cheaply?
            if(kk == max_reach[g])
            {
              min_pos = g;
            }
          }
        */

        pattern[ii] = 0;
        color[kk] = 0;
      }
      top = 0;
      max_reach[k] = min_pos;
      min_pos = MV.ncol;

    }//all columns

   
    //Copy into global 
    //col
    for(Int i = MV.scol; i < MV.scol+MV.ncol; i++)
    {
      gcol[i] += first_row[i-bcol];
    }

    //row
    //for(Int i = 0; i < MV.ncol; i++)
    //{
      //printf("%d %d \n", i, max_reach[i]);
    //}
    for(Int i = MV.scol; i < MV.scol+MV.ncol; i++)
    {
      Int l_min = max_reach[i-bcol];
      //printf("Test maxreach[%d]: %d grow: %d\n",
      //     i-bcol, l_min, grow[i]);
      //printf("Acessing: %d A.ncol: %d \n", i, A.ncol);
      if(l_min < grow[i])
      {
        grow[i] = l_min;
      }
    }
    
    //Copy column_counts
    //printf("MVncol: %d \n", MV.ncol);
    ST.init_U_col_counts(MV.ncol);
    for(Int i = 0; i < MV.ncol; i++)
    {
      ST.U_col_counts[i] = U_col_count[i];
      //temp fix for elbow row
      //ST.U_col_counts[i] += (.4)*U_col_count[i];
    }
    //if symmetric
    ST.init_L_row_counts(MV.ncol);
    for(Int i = 0; i < MV.ncol; i++)
    {
      ST.L_row_counts[i] = U_col_count[i];
      //temp fix for elbow row
      //ST.L_row_counts[i] += (.4)*U_col_count[i];
    }

    //Temp Patch fix
    //Note Comebaske
    if(off_diag ==1)
    {
      for(Int i = 0; i < MV.ncol; i++)
      {
        ST.U_col_counts[i] = MV.nrow;
        ST.L_row_counts[i] = MV.nrow;
      }
    }

    FREE(U_col_count);
    FREE(color);
    FREE(pattern);

    //if symmetric
    FREE(first_color);
    FREE(first_row);
  }//end U_blk_sfactor()


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::L_blk_sfactor
  (
   BASKER_MATRIX &MV,
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol,
   INT_1DARRAY grow
  )
  {
    printf("Basker: This L_blk_sfactor algorithm is not implemented\n");
    //Algorithm
    //You can either use the Row-count method or
    //Assume same as U_blk for symmtric case.
    //Note, Very unsymmtric and HUND will most likely not 
    //Need this called as we will use the QR on nxns 
  }//end L_blk_sfactor()


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::L_blk_sfactor
  (
   BASKER_MATRIX_VIEW &MV,
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol,
   INT_1DARRAY grow
  )
  {
    printf("Basker: This L_blk_sfactor algorithm is not implemented\n");
    //Algorithm
    //You can either use the Row-count method or
    //Assume same as U_blk for symmtric case.
    //Note, Very unsymmtric and HUND will most likely not 
    //Need this called as we will use the QR on nxns 
  }//end L_blk_sfactor()


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::S_sfactor_reduce
  (
   BASKER_MATRIX &MV, 
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol, 
   INT_1DARRAY grow
  )
  {
    printf("Basker: This S_sfactor_reduce algorithm is not implemented\n");
  }//end S_sfactor_reduce()



  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::S_sfactor_reduce
  (
   BASKER_MATRIX_VIEW &MV, 
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol, 
   INT_1DARRAY grow
  )
  {
    printf("Basker: This S_sfactor_reduce algorithm is not implemented\n");
  }//end S_sfactor_reduce()

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::S_blk_sfactor
  (
   BASKER_MATRIX &MV,
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol, 
   INT_1DARRAY grow
  )
  {

    //Algorithm
    //Most likely this will be handled in a dense way (future)
    //See comments in:
    //"Parallel Symbolic Factorization for Sparse LU with.."
    //By Grigori, Demmel, and Li
    //Section 4.1
    
    //Give a = nnz(L(:,1)) and b = nnz(U(1,:))
    //If a*b == (size-size)^2 .... adjust padding

    //Int brow = MV.srow; //Not used
    //Int bcol = MV.scol; //Not used
    
    //Find nnz L(:,1)
    Int nnz_c = 0;
    for(Int i = MV.srow; i < (MV.srow+MV.nrow); i++)
    {
      if(gcol(i) > 0)
      {
        nnz_c++;
      }
    }
    nnz_c += 1;

    #ifdef BASKER_DEBUG_SFACTOR
    printf("S - nnz(L(:,1)): %d \n", nnz_c);
    #endif

    //Int nnz_r = nnz_c; //Not used

    /*
    //Come back to, we know in symtric case same
    //Find nnz U(1,:)
    Int nnz_r = 0;
    for(Int i = MV.scol; i < (MV.scol+MV.ncol); i++)
    {
    Int l_min = grow[i];
    if(l_min == MV.scol)
    {
      nnz_r++;
    }
    }
    #ifdef BASKER_DEBUG_SFACTOR
    printf("S - nnz(U(1,:)): %d \n", nnz_r);
    #endif
    */
    
    //Int nnz_S = (nnz_c*nnz_r)/(MV.ncol);
    //Int nnz_S = (nnz_c*nnz_r);
    //Assumming dense
    Int nnz_S = ((MV.nrow*MV.ncol)+MV.nrow)/2;
    if (nnz_S < 0) {
      // overflow?
      nnz_S = Teuchos::OrdinalTraits<Int>::max()/2;
    }

    #ifdef BASKER_DEBUG_SFACTOR
    printf("Snnz: %d \n", nnz_S);
    #endif

    ST.init_S_col_counts(1);
    ST.S_col_counts[0] = nnz_S;
    
    //Build fake etree (dense) for next level
    ST.init_parent(MV.ncol); //make sure we have enough space
    
    for(Int i = 0; i < MV.ncol-1; i++)
    {
      ST.parent[i] = i+1;
    }

    if((MV.ncol-1)>=0)
    {
      ST.parent[MV.ncol-1] = BASKER_MAX_IDX;
    }
    else
    {
      ST.parent[0] = BASKER_MAX_IDX;
    }
    
  }//end S_blk_sfactor
  
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::S_blk_sfactor
  (
   BASKER_MATRIX_VIEW &MV,
   BASKER_SYMBOLIC_TREE &ST,
   INT_1DARRAY gcol, 
   INT_1DARRAY grow
  )
  {
    //Algorithm
    //Most likely this will be handled in a dense way (future)
    //See comments in:
    //"Parallel Symbolic Factorization for Sparse LU with.."
    //By Grigori, Demmel, and Li
    //Section 4.1
    
    //Give a = nnz(L(:,1)) and b = nnz(U(1,:))
    //If a*b == (size-size)^2 .... adjust padding

    Int brow = MV.srow;
    Int bcol = MV.scol;
    
    //Find nnz L(:,1)
    Int nnz_c = 0;
    for(Int i = MV.srow; i < (MV.srow+MV.nrow); i++)
    {
      if(gcol[i] > 0)
      {
        nnz_c++;
      }
    }
    nnz_c += 1;

    #ifdef BASKER_DEBUG_SFACTOR
    printf("S - nnz(L(:,1)): %d \n", nnz_c);
    #endif

    Int nnz_r = nnz_c;

    /*
    //Come back to, we know in symtric case samce
    //Find nnz U(1,:)
    Int nnz_r = 0;
    for(Int i = MV.scol; i < (MV.scol+MV.ncol); i++)
    {
      Int l_min = grow[i];
      if(l_min == MV.scol)
      {
        nnz_r++;
      }
    }
    #ifdef BASKER_DEBUG_SFACTOR
    printf("S - nnz(U(1,:)): %d \n", nnz_r);
    #endif
    */
    
    //Int nnz_S = (nnz_c*nnz_r)/(MV.ncol);
    //Int nnz_S = (nnz_c*nnz_r);
    //Assumming dense
    Int nnz_S = ((MV.nrow*MV.ncol)+MV.nrow)/2;
    if (nnz_S < 0) {
      // overflow?
      nnz_S = Teuchos::OrdinalTraits<Int>::max()/2;
    }

    #ifdef BASKER_DEBUG_SFACTOR
    printf("Snnz: %d \n", nnz_S);
    #endif

    ST.init_S_col_counts(1);
    ST.S_col_counts[0] = nnz_S;

    //Build fake etree (dense) for next level
    ST.init_parent(MV.ncol); //make sure we have enough space
    for(Int i = 0; i < MV.ncol-1; i++)
    {
      ST.parent[i] = i+1;
    }
    ST.parent[MV.ncol-1] = BASKER_MAX_IDX;

  }//end S_blk_sfactor


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::leaf_assign_nnz
  (
   BASKER_MATRIX &M,
   BASKER_SYMBOLIC_TREE &ST, 
   Int option
  )
  {
    if(option == 0)
    {
      const Int Int_MAX = std::numeric_limits<Int>::max();

      Int t_nnz = 0;
      for(Int i = 0; i < M.ncol; i++)
      {
        if (t_nnz <= Int_MAX - ST.col_counts[i]) {
          t_nnz += ST.col_counts[i];
        } else {
          // let's just hope it is enough, if overflow
          break;
        }
      }
      #ifdef BASKER_DEBUG_SFACTOR
      printf("leaf nnz: %ld \n", (long)t_nnz);
      #endif

      //double nnz_shoulder = 1.05;
      double fill_factor = BASKER_DOM_NNZ_OVER+Options.user_fill; // used to boost fill estimate
      Int temp = fill_factor*t_nnz;
      if (temp > t_nnz) {
        M.nnz = temp;
      } else {
        M.nnz = t_nnz;
      }
      #if 0
      printf( " debug: set leaf.nnz = 0 to force realloc\n" );
      M.nnz = 0;
      t_nnz = 0;
      #endif

      if (global_nnz <= Int_MAX-t_nnz) {
        // let's just hope it is enough, if overflow
        global_nnz += t_nnz;
      }
      if(Options.verbose == BASKER_TRUE)
      {
        printf("leaf with elbow-room global_nnz = %ld, t_nnz = %ld, M.nnz = %ld (%ld x %ld)\n",
               (long)global_nnz,(long)t_nnz,(long)M.nnz,(long)M.nrow,(long)M.ncol);
      }
    }
  }//end assign_leaf_nnz

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::U_assign_nnz
  (
   BASKER_MATRIX &M,
   BASKER_SYMBOLIC_TREE &ST,
   double fill_factor, // used to boost fill estimate
   Int option
  )
  {
    if(option == 0)
    {
      const Int Int_MAX = std::numeric_limits<Int>::max();

      Int t_nnz = 0; 
      for(Int i = 0; i < M.ncol; i++)
      {
        if (t_nnz <= Int_MAX-ST.U_col_counts[i]) {
          t_nnz += ST.U_col_counts[i];
        } else {
          // let's just hope it is enough, if overflow
          break;
        }
      }

      #ifdef BASKER_DEBUG_SFACTOR
      printf("U_assing_nnz: %ld \n", t_nnz);
      #endif

      //double fill_factor = 1.05;
      Int temp = fill_factor*t_nnz;
      if (temp >= t_nnz) {
        M.nnz = temp;
      } else {
        M.nnz = t_nnz;
      }
      if (global_nnz <= Int_MAX-t_nnz) {
        // let's just hope it is enough, if overflow
        global_nnz += t_nnz;
      }
      #if 0
      printf( " debug: set U.nnz = 0 to force realloc\n" );
      M.nnz = 0;
      t_nnz = 0;
      #endif
      if(Options.verbose == BASKER_TRUE)
      {
        printf("U_assing with elbow global_nnz = %ld, t_nnz = %ld (fill_factor = %f + %f = %f), M.nnz = %ld (%ld x %ld)\n",
               (long)global_nnz,(long)t_nnz, BASKER_DOM_NNZ_OVER,Options.user_fill,fill_factor, (long)M.nnz,(long)M.nrow,(long)M.ncol);
      }
    }
  }//end assign_upper_nnz

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::L_assign_nnz
  (
   BASKER_MATRIX &M,
   BASKER_SYMBOLIC_TREE &ST,
   double fill_factor, // used to boost fill estimate
   Int option
  )
  {
    if(option == 0)
    {
      const Int Int_MAX = std::numeric_limits<Int>::max();

      Int t_nnz = 0; 
      for(Int i = 0; i < M.nrow; i++)
      {
        if (t_nnz <= Int_MAX-ST.L_row_counts[i]) {
          t_nnz += ST.L_row_counts[i];
        } else {
          // let's just hope it is enough, if overflow
          break;
        }
      }

      #ifdef BASKER_DEBUG_SFACTOR
      printf("L_assign_nnz: %ld \n", t_nnz);
      #endif

      // double fill_factor = 2.05;
      double old_nnz = M.nnz;
      Int temp = fill_factor*t_nnz;
      if (temp >= t_nnz) {
        M.nnz = temp;
      } else {
        M.nnz = t_nnz;
      }
      #if 0
      printf( " debug: set L.nnz = 0 to force realloc\n" );
      M.nnz = 0;
      t_nnz = 0;
      #endif
      if (global_nnz <= Int_MAX-t_nnz) {
        // let's just hope it is enough, if overflow
        global_nnz += t_nnz;
      }
      if(Options.verbose == BASKER_TRUE)
      {
        printf("L_assign with elbow global_nnz = %ld, t_nnz = %ld (fill_factor = %e + %e = %e), M.nnz = %ld -> %ld (%ld x %ld)\n",
               (long)global_nnz,(long)t_nnz, BASKER_DOM_NNZ_OVER,Options.user_fill, fill_factor, (long)old_nnz,(long)M.nnz, (long)M.nrow,(long)M.ncol);
      }
    }
  }//end assign_lower_nnz


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::S_assign_nnz
  (
   BASKER_MATRIX &M,
   BASKER_SYMBOLIC_TREE &ST,
   Int option
  )
  {
    if(option == 0)
    {
      M.nnz = ST.S_col_counts(0);
      #ifdef BASKER_DEBUG_SFACTOR
      printf("S_assign_nnz: %ld  \n", M.nnz);
      #endif

      const Int Int_MAX = std::numeric_limits<Int>::max();
      if (global_nnz <= Int_MAX-M.nnz) {
        // let's just hope it is enough, if overflow
        global_nnz += M.nnz;
      }
      if(Options.verbose == BASKER_TRUE)
      {
        printf("S_assign elbow global_nnz = %ld, M.nnz = %ld + 2\n", (long)global_nnz, (long)M.nnz);
      }
      if (M.nnz <= Int_MAX - 2) {
        M.nnz += 2;
      }
    }
  }//end assign_sep_nnze


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::btf_last_dense(bool flag)
  {
    
    //Here we scan the lower part and malloc 
    #ifdef BASKER_DEBUG_SFACTOR
    //printf("Test btf_last_dense \n");
    //printf("btf_tabs_offset: %d col: %d \n",
    //	   btf_tabs_offset, btf_tabs[btf_tabs_offset]);
    //printf("number of blks: %d \n",
    //	   btf_nblks-btf_tabs_offset);
    #endif

    Int max_blk_size = 0;
#if defined(BASKER_SPLIT_A)
    if (btf_top_nblks > 0) {
      MALLOC_MATRIX_1DARRAY(L_D, btf_top_nblks);
      MALLOC_MATRIX_1DARRAY(U_D, btf_top_nblks);

      for(Int i = 0; i < btf_top_nblks; i++)
      {
        Int lblk_size = btf_tabs(i+1)-btf_tabs(i);

        if(lblk_size > max_blk_size)
        {
          max_blk_size = lblk_size;
        }

        Int nnz = (btf_blk_nnz(i)+lblk_size)*(BASKER_BTF_NNZ_OVER+Options.user_fill);
        if ((double)nnz > ((double)lblk_size)*((double)lblk_size)) {
            nnz = lblk_size*lblk_size;
        }
        //printf( " LBTF(%d, nnz = %d)\n",(int)(i-btf_tabs_offset), (int)nnz );
        L_D(i).init_matrix("LBFT",
          btf_tabs(i),
          lblk_size,
          btf_tabs(i),
          lblk_size,
          nnz);

        //For pruning
        L_D(i).init_pend();

        U_D(i).init_matrix("UBFT",
          btf_tabs(i),
          lblk_size,
          btf_tabs(i),
          lblk_size,
          nnz);
      }//over all blks
    }
#endif

    //Malloc L and U
    #ifdef BASKER_DEBUG_SFACTOR
    printf("btf_nblks %d btf_tabs_offset %d \n",
	   btf_nblks, btf_tabs_offset);
    #endif
    
    Int nblks_left = btf_nblks - btf_tabs_offset;
    if(nblks_left > 0) {

      BASKER_ASSERT(nblks_left > 0, "Basker btf_last_dense assert: nblks_left > 0 failed");
      MALLOC_MATRIX_1DARRAY(LBTF, nblks_left);
      MALLOC_MATRIX_1DARRAY(UBTF, nblks_left);

      for(Int i = btf_tabs_offset; i < btf_nblks; i++)
      {
        Int lblk_size = btf_tabs(i+1)-btf_tabs(i);

        if(lblk_size > max_blk_size)
        {
          max_blk_size = lblk_size;
        }

        Int nnz = (btf_blk_nnz(i)+lblk_size)*(BASKER_BTF_NNZ_OVER+Options.user_fill);
        if ((double)nnz > ((double)lblk_size)*((double)lblk_size)) {
          nnz = lblk_size*lblk_size;
        }
        //printf( " LBTF(%d, nnz = %d)\n",(int)(i-btf_tabs_offset), (int)nnz );
        LBTF(i-btf_tabs_offset).init_matrix("LBFT",
          btf_tabs(i),
          lblk_size,
          btf_tabs(i),
          lblk_size,
          nnz);

        //For pruning
        //printf( " LBTF(%d).init_pend()\n",(int)(i-btf_tabs_offset) );
        LBTF(i-btf_tabs_offset).init_pend();

        //printf( " UBTF(%d, nnz = %d)\n",(int)(i-btf_tabs_offset), (int)nnz );
        UBTF(i-btf_tabs_offset).init_matrix("UBFT",
          btf_tabs(i),
          lblk_size,
          btf_tabs(i),
          lblk_size,
          nnz);
        //(.5*lblk_size*lblk_size)+lblk_size);

        //will have to do the fill in nfactor
        //MALLOC workspace
      }//over all blks
    }

    //JDB: This needs to be fixed
    max_blk_size = BTF_D.nrow + BTF_C.nrow;

    if (flag) {
      // though offset=0, we may still have left-over BLK factorization
      //if(btf_tabs_offset == 0)
      {
        BASKER_ASSERT(num_threads > 0, "Basker btf_last_dense assert: num_threads > 0 failed");
        MALLOC_THREAD_1DARRAY(thread_array, num_threads);
      }

      for(Int i = 0 ; i < num_threads; i++)
      {
        thread_array(i).iws_size = max_blk_size;
        thread_array(i).ews_size = max_blk_size;

        //BASKER_ASSERT((thread_array(i).iws_size*thread_array(i).iws_mult) > 0, "Basker btf_last_dense assert: sfactor threads iws > 0 failed");
        //BASKER_ASSERT((thread_array(i).ews_size*thread_array(i).ews_mult) > 0, "Basker btf_last_dense assert: sfactor threads ews > 0 failed");
        if (max_blk_size > 0) {
          MALLOC_INT_1DARRAY(thread_array(i).iws, thread_array(i).iws_size*thread_array(i).iws_mult);
          MALLOC_ENTRY_1DARRAY(thread_array(i).ews, thread_array(i).ews_size*thread_array(i).ews_mult);
        }
        #ifdef BASKER_DEBUG_SFACTOR
        printf("Malloc Thread: %d iws: %d \n",
            i, (thread_array(i).iws_size*
              thread_array(i).iws_mult));
        printf("Malloc Thread: %d ews: %d \n", 
            i, (thread_array(i).ews_size*
              thread_array(i).ews_mult));
        #endif
      }
    }

  }//end btf_last_dense()
  
}//end namespace Bakser

#endif//endif BASKER_SFACTOR_NEWFRM_HPP
