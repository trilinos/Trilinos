// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_ORDER_SCOTCH_HPP
#define SHYLUBASKER_ORDER_SCOTCH_HPP

#include "shylubasker_types.hpp"
#include "shylubasker_sswrapper.hpp"

#ifdef HAVE_SHYLU_NODEBASKER_SCOTCH
  #include "scotch.h"
#endif
#ifdef HAVE_SHYLU_NODEBASKER_METIS
 #include "metis.h"
#endif

//#define BASKER_DEBUG_ORDER_SCOTCH
//#define BASKER_TIMER

//NOTE need to change all the max_idx here still

namespace BaskerNS
{
  template <typename iType>
  struct scotch_graph
  {
    static_assert( std::is_same<iType,int32_t>::value || std::is_same<iType,int64_t>::value
                 , "ShyLU Basker Error: scotch_graph members must be templated on type int32_t or int64_t only");
    scotch_graph()
    {};    
    iType m;
    iType nz;
    iType *Ap;
    iType *Ai;
    iType cblk;
    iType *permtab;
    iType *peritab;
    iType *rangtab;
    iType *treetab;
  };

  //A+A'
  //There is probably a much better and faster way to do this
  //Error here
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::AplusAT
  (
   BASKER_MATRIX &M,
   BASKER_MATRIX &C,
   BASKER_BOOL keep_zeros
  )
  {
    BASKER_MATRIX T;
    //get matrix transpose
    matrix_transpose(M, T, keep_zeros);
   
    C.set_shape(M.srow, M.nrow, M.scol, M.ncol);

    //BASKER_ASSERT((M.ncol+1) > 0, "scotch ncol"); 
    BASKER_ASSERT(M.nrow > 0, "scotch nrow");
    BASKER_ASSERT(M.nnz  > 0, "scotch nnz");

    MALLOC_INT_1DARRAY(C.col_ptr, M.ncol+1);
    MALLOC_INT_1DARRAY(C.row_idx, 2*M.nnz);

    //init_value(C.col_ptr, M.ncol+1, (Int) 0);
    //init_value(C.row_idx, 2*M.nnz, (Int) 0);
    
    INT_1DARRAY ws;
    MALLOC_INT_1DARRAY(ws, M.nrow);
    init_value(ws, M.nrow, BASKER_MAX_IDX);

    const Entry zero (0.0);
    Int c_nnz = 0;
    C.col_ptr(0) = 0;
    for(Int k = 0 ; k <  M.ncol; ++k)
    {
      //scatter M
      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        Int j = M.row_idx(i);
        if (ws(j) != k)
        {
          if (keep_zeros || M.val(i) != zero) {
            C.row_idx(c_nnz) = j;
            c_nnz++;
            ws(j) = k;
          }
        }
      }

      //scatter T
      for(Int i = T.col_ptr(k); i < T.col_ptr(k+1); ++i)
      {
        Int j = T.row_idx(i);
        if(ws(j) != k)
        {
          C.row_idx(c_nnz) = j;
          c_nnz++;
          ws(j) = k;
        }
      }

      C.col_ptr(k+1) = c_nnz;
    }
    
    //sort columns??

    #if 0
    //can be remove after debuggin
    BASKER_ASSERT(c_nnz > 0, "C.val scotch");
    MALLOC_ENTRY_1DARRAY(C.val, c_nnz);
    init_value(C.val,c_nnz, (Entry)1.0);
    #endif
    C.nnz = c_nnz;

    FREE(T);
    return 0;
  }//A+A'


  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::part_scotch
  (
   BASKER_MATRIX &M,
   BASKER_TREE &BT
  )
  {
    Int lvls = round(log(num_threads)/log(2));
    return part_scotch(M, BT, lvls);
  }//end part_scotch


  //Note that current part_scotch only works with symmetric matrix
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::part_scotch
  (
   BASKER_MATRIX &M,
   BASKER_TREE &BT,
   Int num_domains
  )
  {
    #if SHYLU_SCOTCH_64
    using scotch_integral_type = int64_t; //NDE: make this depend on the scotch type
    #else
    using scotch_integral_type = int32_t; //NDE: make this depend on the scotch type
    #endif
    Kokkos::Timer timer_scotch;

    Int num_levels = num_domains; 
    Int num_doms   = pow(2.0, (double)(num_levels+1)) - 1;
    scotch_graph<scotch_integral_type> sg;

    sg.m = M.nrow;
    sg.cblk = 0;
    sg.permtab = (scotch_integral_type *)malloc((sg.m)  *sizeof(scotch_integral_type));
    sg.peritab = (scotch_integral_type *)malloc((sg.m)  *sizeof(scotch_integral_type));
    sg.rangtab = (scotch_integral_type *)malloc((num_doms+1)*sizeof(scotch_integral_type));
    sg.treetab = (scotch_integral_type *)malloc((num_doms+1)*sizeof(scotch_integral_type));

    bool run_nd_on_leaves  = Options.run_nd_on_leaves;
    bool run_amd_on_leaves = Options.run_amd_on_leaves;
    if (num_levels == 0 && !run_nd_on_leaves && !run_amd_on_leaves) {
      if (Options.verbose == BASKER_TRUE) {
        std::cout << std::endl << " + No ND since one level +" << std::endl;
      }
      for(Int i = 0; i < sg.m; i++)
      {
        sg.permtab[i] = i;
        sg.peritab[i] = i;
      }
      // just one block
      sg.cblk = 1;
      sg.rangtab[0] = 0;
      sg.rangtab[1] = sg.m;
      // root
      sg.treetab[0] = -1;
    } else {
      #if defined (HAVE_SHYLU_NODEBASKER_METIS) || defined(HAVE_SHYLU_NODEBASKER_SCOTCH)
      if (Options.use_metis == BASKER_TRUE) {
        #if !defined(HAVE_SHYLU_NODEBASKER_METIS)
        BASKER_ASSERT(false, ">> BASKER ASSERT: METIS is not enabled <<");
        #else
        Kokkos::Timer timer_metis;
        double time_metis = 0.0;

        //idx_t  metis_offset = btf_tabs(btf_tabs_offset-1);
        idx_t  metis_offset = 0; // BTF_A now contains one big block
        idx_t  metis_size = M.nrow - metis_offset;
        idx_t  metis_nnz = M.col_ptr(M.nrow);

        using METIS_1DARRAY = Kokkos::View<idx_t*,  BASKER_EXE_SPACE>;
        METIS_1DARRAY metis_part   (BASKER_KOKKOS_NOINIT("metis_part"),   metis_size);
        METIS_1DARRAY metis_rowptr (BASKER_KOKKOS_NOINIT("metis_rowptr"), metis_size+1);
        METIS_1DARRAY metis_colidx (BASKER_KOKKOS_NOINIT("metis_colidx"), metis_nnz);

        METIS_1DARRAY metis_part_k  (BASKER_KOKKOS_NOINIT("metis_part_k"),  metis_size);
        METIS_1DARRAY metis_perm_k  (BASKER_KOKKOS_NOINIT("metis_perm_k"),  metis_size);
        METIS_1DARRAY metis_iperm_k (BASKER_KOKKOS_NOINIT("metis_iperm_k"), metis_size);

        // for calling AMD on leaves
        INT_1DARRAY amd_rowptr (BASKER_KOKKOS_NOINIT("amd_rowptr"), (run_amd_on_leaves ? metis_size+1 : 0));
        INT_1DARRAY amd_colidx (BASKER_KOKKOS_NOINIT("amd_colidx"), (run_amd_on_leaves ? metis_nnz : 0));
        INT_1DARRAY amd_perm_k (BASKER_KOKKOS_NOINIT("amd_perm_k"), (run_amd_on_leaves ? metis_size : 0));

        // to find vertex cover/separator
        using METIS_2DARRAY = Kokkos::View<idx_t**, BASKER_EXE_SPACE>;
        METIS_2DARRAY metis_vc (BASKER_KOKKOS_NOINIT("metis_vc"), metis_nnz, 3);
        METIS_1DARRAY metis_vc_score (BASKER_KOKKOS_NOINIT("metis_vc_score"),  metis_nnz);

        int info = 0;
        idx_t sepsize = 0;
        idx_t options[METIS_NOPTIONS];

        // --------------------------------------------------- 
        // compute post-order
        INT_1DARRAY metis_queue;
        INT_1DARRAY metis_check;
        INT_1DARRAY post_order;  // post(i) is the new ith domain id after ND
        INT_1DARRAY post_iorder; // the original ith domain becomes ipost(i) th domain after ND
        MALLOC_INT_1DARRAY(metis_queue, num_doms);
        MALLOC_INT_1DARRAY(metis_check, num_doms);
        MALLOC_INT_1DARRAY(post_order,  num_doms);
        MALLOC_INT_1DARRAY(post_iorder, num_doms);
        for (Int i = 0; i < num_doms; i++) {
          metis_check(i) = 0;
        }

        // id of the first leaf node (BF order, post_order maps from BF to ND)
        Int leaves_id = pow(2.0, (double)(num_levels)) - 1;
        //printf( " num_levels = %d, num_doms = %d, leves_id = %d\n",num_levels,num_doms,leaves_id );

        // > insert root
        Int num_queued = 0;
        metis_queue(num_queued) = 0;
        num_queued ++;

        num_doms = 0;
        while (num_queued > 0) {
          // pop a node from queue
          Int dom_id = metis_queue(num_queued-1);
          //printf( " > check (dom_id = %d) = %d\n",dom_id,metis_check(dom_id) );
          if (dom_id >= leaves_id ||     // leaf
              metis_check(dom_id) == 2)  // both childrend processed
          {
            post_order(num_doms) = dom_id;
            //printf( " pop queue(%d) = %d -> post(%d) = %d\n\n",num_queued-1,dom_id, num_doms,dom_id );
            num_doms ++;

            if (dom_id != 0) {
              // if not root, let the parent node know one of its children has been processed
              Int parent_id = (dom_id - 1)/2;
              metis_check(parent_id) ++;
            }

            num_queued --;
          } else {
            // push right child
            //printf( " push queue(%d) = %d\n",num_queued,2*dom_id+2 );
            metis_queue(num_queued) = 2*dom_id + 2;
            num_queued ++;
            // push left child
            //printf( " push queue(%d) = %d\n\n",num_queued,2*dom_id+1 );
            metis_queue(num_queued) = 2*dom_id + 1;
            num_queued ++;
          }
        }
        for (Int i = 0; i < num_doms; i++) {
          post_iorder(post_order(i)) = i;
        }
        /*for (Int i = 0; i < num_doms; i++) {
          printf( " post[%d] = %d, ipost[%d]=%d\n",i,post_order(i),i,post_iorder(i) );
        }*/
        //M.print_matrix("A.dat");

        // initial partition
        sg.cblk = 1;
        for (Int i = 0; i < M.nrow; i++) {
          metis_part[i] = post_iorder(0);
          sg.permtab[i] = i;
          sg.peritab[i] = i;
        }
        sg.treetab[post_iorder(0)] = -1;
        for (Int i = 0; i < num_doms; i++) {
          sg.rangtab[i] = 0;
        }
        sg.rangtab[num_doms] = M.nrow;
        Int last_level = num_levels - 1;
        if (run_nd_on_leaves || run_amd_on_leaves) {
          // level goes to num_leaves so that we can call ND on the final leaf nodes
          last_level = num_levels;
        }
        if (Options.verbose == BASKER_TRUE) {
          if (run_nd_on_leaves) {
            std::cout << std::endl << " + Using ND on leaves + " << std::endl;
          } else if (run_amd_on_leaves) {
            std::cout << std::endl << " + Using AMD on leaves + " << std::endl;
          }
        }
        // -------------------------------------------------- //
        Int level = 0;
        bool use_nodeNDP = Options.use_nodeNDP;
        if (use_nodeNDP && num_levels > 0) {
          if (METIS_OK != METIS_SetDefaultOptions(options)) {
            std::cout << std::endl << " > METIS_SetDefaultOptions failed < " << std::endl << std::endl;
            return BASKER_ERROR; // TODO: what to do here?
          }
          // remove diagonals
          idx_t nnz_k =0;
          metis_rowptr[0] = 0;
          for(Int j = 0; j < metis_size; j++) {
            for(Int k = M.col_ptr(j); k < M.col_ptr(j+1); k++)
            {
              Int i = M.row_idx(k);
              if(i != j) {
                metis_colidx[nnz_k] = i;
                nnz_k ++;
              }
            }
            metis_rowptr[j+1] = nnz_k;
          }
          idx_t *vwgt = nullptr;    // contraints (n * num_constraints)
          Int num_leaves = pow(2.0, (double)(num_levels));
          METIS_1DARRAY metis_sep_sizes (BASKER_KOKKOS_NOINIT("metis_isizes_k"), 2*num_doms-1);
          if (Options.verbose == BASKER_TRUE) {
            std::cout << std::endl << " > calling METIS_NodeNDP ( n = " << metis_size
                      << ", num_leaves = " << num_leaves << " ) << " << std::endl;
          }
          timer_metis.reset();
          METIS_NodeNDP(metis_size, &(metis_rowptr(0)), &(metis_colidx(0)), vwgt,
                        num_leaves, options, &(metis_perm_k(0)), &(metis_iperm_k(0)), &(metis_sep_sizes(0)));
          time_metis += timer_metis.seconds();
          #if 0
          // debug: merging all the subdomains into one domain
          //metis_sep_sizes(0) = metis_size;
          //for (int i=1; i < 2*num_leaves-1; i++) metis_sep_sizes(i) = 0;

          // debug: merging the first two subdomains
          //metis_sep_sizes(0) += metis_sep_sizes(1);
          //metis_sep_sizes(1) = 0;
          #endif
          //for (int i=0; i < 2*num_leaves-1; i++) printf( " > size[%d] = %d\n",i,metis_sep_sizes(i) );
          for(Int i = 0; i < metis_size; i++) {
            sg.peritab[i] = metis_perm_k[i];
            sg.permtab[i] = metis_iperm_k[i];
          }

          // construct ND tree
          INT_1DARRAY revert_tree; // invert ND node the breadth-first order from top-to-bottom to bottom-to-top
          MALLOC_INT_1DARRAY(revert_tree, num_doms);
          for (Int level_k = 0; level_k <= num_levels; level_k++) {
            Int first_sep1  = pow(2.0, (double)(  level_k)) - 1; // id of the first leaf at this level         (top-to-bottom)
            Int first_leaf1 = pow(2.0, (double)(1+level_k)) - 1; // id of the first new leaf at the next level (top-to-bottom)

            Int first_sep2  = num_doms - first_leaf1; // id of the first leaf at this level (bottom-to-top)
            for (Int j = first_sep1; j < first_leaf1; j++) {
              //printf( " %d: revert_trree(%d) = %d + %d\n",level_k, j,first_sep2,j-first_sep1 );
              revert_tree(j) = first_sep2+(j-first_sep1);
            }

            if (level_k < num_levels) {
              Int num_leaves_k = pow(2.0, (double)(level_k));       // number of leaves at this level
              //printf(" num_leaves(%d) = %d\n",level_k,num_leaves_k );
              for (Int leaf_id = 0; leaf_id < num_leaves_k; leaf_id++) {
                Int dom_id = 2 * leaf_id + first_leaf1;
                Int dom_id1 = dom_id;     // id of left-child after bisection
                Int dom_id2 = dom_id + 1; // id of right-child after bisection
                Int sep_id = (dom_id1)/2; // id of this domain before bisection (becomes separator after bisection)

                //printf( " + level=%d: post(%d)=%d,post(%d)=%d -> post(%d)=%d\n",level_k, dom_id1,post_iorder(dom_id1),dom_id2,post_iorder(dom_id2),sep_id,post_iorder(sep_id) );
                dom_id1 = post_iorder(dom_id1);  // id of left-child after bisection
                dom_id2 = post_iorder(dom_id2);  // id of right-child after bisection
                sep_id  = post_iorder(sep_id);   // id of this domain before bisection (becomes separator after bisection)
                sg.treetab[dom_id1] = sep_id;
                sg.treetab[dom_id2] = sep_id;
              }
            } else {
              sg.treetab[num_doms-1] = -1;
            }
          }
          sg.treetab[num_doms] = 0;
          sg.cblk = num_doms;

          // compute row offset to separator
          sg.rangtab[0] = 0;
          for(Int i = 0; i < sg.cblk; i++) {
            sg.rangtab[i+1] = sg.rangtab[i] + metis_sep_sizes[revert_tree[post_order(i)]];
            for (Int j = sg.rangtab[i]; j < sg.rangtab[i+1]; j++) {
              Int col = sg.peritab[j]; // j is after ND, col is original
              metis_part[col] = i;
            }
            if (Options.verbose == BASKER_TRUE) {
              std::cout << " + metis_sep_sizes[" << i << " -> " << revert_tree[post_order(i)] << "] = "
                        << metis_sep_sizes[revert_tree[post_order(i)]] << std::endl;
            }
          }
          //for(Int i = 0; i <= sg.cblk; i++) {
          //  printf( " + row_tabs[%d] = %d, col_tabs[%d] = %d, treetab[%d] = %d\n",i,sg.rangtab[i], i,sg.rangtab[i], i,sg.treetab[i] );
          //}
          // may run ND on leaves
          level = num_levels;
        }
        // calling bisection each domain, or NodeND/AMD on the leaf
        {
          // -------------------------------------------------- //
          for (; level <= last_level; level++) {
            //printf( "\n ================== (level = %d) =========================\n",level );
            Int num_leaves = pow(2.0, (double)(level));       // number of leaves at this level
            Int first_sep  = pow(2.0, (double)(level)) - 1;   // id of the first leaf at this level
            Int first_leaf = pow(2.0, (double)(1+level)) - 1; // id of the first new leaf at the next level
            //printf( "  num_leaves = %d, first_sep = %d, first_leaf = %d\n",num_leaves,first_sep,first_leaf );

            for (Int leaf_id = 0; leaf_id < num_leaves; leaf_id++) {
              // extract k-th interoior
              Int dom_id = 2 * leaf_id + first_leaf;
              Int dom_id1 = dom_id;     // id of left-child after bisection
              Int dom_id2 = dom_id + 1; // id of right-child after bisection
              Int sep_id = (dom_id1)/2; // id of this domain before bisection (becomes separator after bisection)
              //printf( "\n > level = %d, dom_id = 2*%d + %d = %d, dom_id1 = %d, dom_id2 = %d, sep_id = %d\n",level,leaf_id,first_leaf,dom_id,dom_id1,dom_id2,sep_id );

              // trying to figure out which block comes right before me
              Int sep_left_sibling = 0;
              if (sep_id == first_sep) {
                // this is the first leaf
                sep_left_sibling = 0;
              } else if (leaf_id%2 == 1) {
                // this is right-child, so assign to left-sibling
                //printf( " set_left_sibling = 1 + iorder(%d - 1) = 1 + %d\n",sep_id,post_iorder(sep_id - 1) );
                sep_left_sibling = 1+post_iorder(sep_id - 1);
              } else {
                // this is left-child, so assign to left-sibling of its separator
                //sep_left_sibling = 1+(post_iorder(sep_id/2 - 1));
                Int last_left_child = pow(2.0, (num_levels - level))*(1+sep_id) - 1;
                //printf( " set_left_sibling = 1 + iorder(2^(%d-%d) * (1+%d)) - 1)-1 = 1 + iorder(%d)-1 = 1+%d-1\n",num_levels,level,sep_id,last_left_child,post_iorder(last_left_child) );
                sep_left_sibling = 1+(post_iorder(last_left_child)-1);
              }
              //printf( "  -> level = %d, sep_left_sibling = %d (sep_id = %d, first_sep = %d)\n",level,sep_left_sibling,sep_id,first_sep );

              if (level < num_levels) 
              {
                // post order
                //printf( " * level = %d, dom_id1 = %d, dom_id2 = %d, sep_id = %d\n",level,dom_id1,dom_id2,sep_id );
                dom_id1 = post_iorder(dom_id1);
                dom_id2 = post_iorder(dom_id2);
                //printf( "  -> level = %d, dom_id1 = %d, dom_id2 = %d, sep_id = %d\n",level,dom_id1,dom_id2,sep_id );
                //printf( "  -> level = %d, rantab[%d] = %d, rangtab[%d] = %d, rangtab[%d] = %d\n",level,dom_id1,sg.rangtab[dom_id1],dom_id2,sg.rangtab[dom_id2],sep_id,sg.rangtab[sep_id] );
                //printf( "  -> level = %d, dom = %d -> %d (nrow = %d)\n",level,dom_id,sep_id,M.nrow );
              }
              sep_id = post_iorder(sep_id);
              Int frow = sg.rangtab[sep_left_sibling];
              Int metis_offset_k = frow + metis_offset;
              //printf( "frow = rangtab[%d] = %d\n",sep_left_sibling,frow );
              //printf( " -> metis_offset_k = %d + %d\n",frow,metis_offset );
              //for(Int i = 0; i < M.nrow; i++) {
              //  printf( " %d %d %d %d\n",i, metis_part(i), sg.permtab[i], sg.peritab[i] );
              //}
              idx_t metis_size_k = 0;
              idx_t nnz_k = 0;
              metis_rowptr[0] = 0;
              for(Int j = metis_offset; j < M.ncol; j++) {
                Int col = sg.peritab[j]; // j is after ND, col is original
                if (metis_part(col) == sep_id) {
                  for(Int k = M.col_ptr(col); k < M.col_ptr(col+1); k++)
                  {
                    Int row = M.row_idx(k);
                    Int i = sg.permtab[row]; // i is after ND, row is original
                    if(row != col && i >= metis_offset_k) {
                      if(metis_part(row) == sep_id)
                      {
                        metis_colidx[nnz_k] = i - frow;
                        nnz_k ++;
                      }
                    }
                  }
                  metis_rowptr[metis_size_k+1] = nnz_k;
                  metis_size_k ++;
                }
              }

              //printf( " metis_size = %d (sizeof = %d), n = %d (sizeof = %d)\n",metis_size_k,sizeof(idx_t),M.nrow,sizeof(Int) );
              /*printf( " Ke = [\n" );
              for(Int i = metis_offset; i < M.nrow; i++)
              {
                Int col = sg.peritab[i];
                for(Int k = M.col_ptr(col); k <M.col_ptr(col+1); k++)
                {
                  Int j = M.row_idx(k);
                  Int row = sg.permtab[j];
                  printf( "%d %d %d %d\n",j,i,row,col );
                }
              }
              printf( "];\n" );*/
              /*printf( " Me = [\n" );
              for(Int i = 0; i < metis_size_k; i++) {
                for(Int k = metis_rowptr(i); k < metis_rowptr(i+1); k++) printf( "%d %d %d\n",i,metis_colidx(k),k );
              }
              printf( "];\n" );*/

              sepsize = 0;
              info = 0;
              if (METIS_OK != METIS_SetDefaultOptions(options)) {
                std::cout << std::endl << " > METIS_SetDefaultOptions failed < " << std::endl << std::endl;
                return BASKER_ERROR; // TODO: what to do here?
              }
              if (level == num_levels) {
                // ===============================================
                // calling ND or AMD on final leaves (fill-reducing)
                if (Options.verbose == BASKER_TRUE) {
                  std::cout << std::endl << " > calling " << (run_amd_on_leaves ? "Basker_AMD" : "METIS_NodeND" ) << " on leaf " << leaf_id
                            << " at final level " << level << ", size = " << metis_size_k << ", sep-id = " << sep_id
                            << " < " << std::endl;
                }
                if (run_nd_on_leaves) {
                  // calling ND
                  // perm(i) of the original matrix is i-th row in the new matrix
                  timer_metis.reset();
                  info = METIS_NodeND(&metis_size_k,
                                      &(metis_rowptr(0)),
                                      &(metis_colidx(0)),
                                       nullptr,
                                       options,
                                      &(metis_perm_k(0)),
                                      &(metis_iperm_k(0)));
                  time_metis += timer_metis.seconds();
                  if (info != METIS_OK) {
                    std::cout << std::endl << " > METIS_NodeND failed < " << std::endl << std::endl;
                    return BASKER_ERROR; // TODO: what to do here?
                  }
                }
                else if (run_amd_on_leaves) {
                  // calling AMD
                  // just copying for now, in case of Int and idx_t type-mismatch
                  for (Int i = 0; i <= metis_size_k; i++) amd_rowptr(i) = metis_rowptr(i);
                  for (Int i = 0; i <  nnz_k; i++) amd_colidx(i) = metis_colidx(i);
                  double l_nnz, lu_work;
                  info = BaskerSSWrapper<Int>::amd_order(metis_size_k,
                                                         &(amd_rowptr(0)),
                                                         &(amd_colidx(0)),
                                                         &(amd_perm_k(0)),
                                                         l_nnz, lu_work, Options.verbose);
                  for (Int i = 0; i < metis_size_k; i++) metis_perm_k(i) = amd_perm_k(i);

                  if (info != TRILINOS_AMD_OUT_OF_MEMORY && info != TRILINOS_AMD_INVALID) {
                    for(Int i = 0; i < metis_size_k; i++) {
                      metis_iperm_k(metis_perm_k(i)) = i;
                    }
                    info = METIS_OK;
                  } else {
                    std::cout << std::endl << " > Basker AMD failed < " << std::endl << std::endl;
                    return BASKER_ERROR; // TODO: what to do here?
                  }
                }

                // update perm/
                for(Int i = 0; i < metis_size_k; i++) {
                  // metis_part_k[i] is the original global index before ND
                  metis_part_k[i] = sg.peritab[frow + metis_perm_k[i]];
                }
                for(Int i = 0; i < metis_size_k; i++) {
                  sg.peritab[frow + i] = metis_part_k[i];
                  sg.permtab[sg.peritab[frow + i]] = frow + i;
                }
              } else {
                // ===============================================
                // finding separator on the current leaves (ND)
                bool use_metis_kway = false;
                timer_metis.reset();
                if (use_metis_kway) {
                  // compute two-way partition
                  if (Options.verbose == BASKER_TRUE) {
                    std::cout << std::endl << " > calling METIS_PartGraphKway on leaf " << leaf_id
                              << " at level " << level << " < " << std::endl;
                  }
                  idx_t num_constraints = 1;
                  idx_t num_parts = 2;
                  idx_t objval = -1;

                  idx_t *vwgt = nullptr;    // contraints (n * num_constraints)
                  idx_t *vsize = nullptr;   // for total comm vol
                  idx_t *adjwgt = nullptr;  // for reducing cut

                  real_t *tpwgts = nullptr; // weights for eachh partition              (num_parts * num_constraints)
                  real_t *ubvec = nullptr;  // imbalance tolerance for each constraints (num_constraints)
                  info = METIS_PartGraphKway(&metis_size_k,
                                             &num_constraints,
                                             &(metis_rowptr(0)),
                                             &(metis_colidx(0)),
                                              vwgt,
                                              vsize,
                                              adjwgt,
                                             &num_parts,
                                              tpwgts,
                                              ubvec,
                                              options,
                                             &objval,
                                             &(metis_part_k(0)));
                  // look for edge separator
                  for(Int i = 0; i < metis_size_k; i++) {
                    Int con1 = 0;
                    Int con2 = 0;
                    for (idx_t k = metis_rowptr(i); k < metis_rowptr(i+1); k++) {
                      if (metis_part_k(metis_colidx(k)) == 0 || metis_part_k(metis_colidx(k)) == -2) {
                        con1 ++;
                      } else {
                        con2 ++;
                      }
                    }
                    if (con1 > 0 && con2) {
                      if (metis_part_k(i) == 0) {
                        metis_part_k(i) = -2;
                      } else {
                        metis_part_k(i) = 2;
                      }
                    }
                  }
                  #if 0
                  for(Int i = 0; i < metis_size_k; i++) {
                    if (metis_part_k(i) == -2) {
                      // put it to edge separator
                      metis_part_k(i) = 2;
                    }
                  }
                  #else
                  // look for vertex cover/separator
                  // > look for edges within edge separator
                  Int num_edges = 0;
                  for(Int i = 0; i < metis_size_k; i++) {
                    if (metis_part_k(i) == 2 || metis_part_k(i) == -2) {
                      metis_vc_score(i) = 0;
                    }
                  }
                  for(Int i = 0; i < metis_size_k; i++) {
                    if (metis_part_k(i) == 2 || metis_part_k(i) == -2) {
                      for (idx_t k = metis_rowptr(i); k < metis_rowptr(i+1); k++) {
                        if (metis_part_k(metis_colidx(k)) == 2 || metis_part_k(metis_colidx(k)) == -2) {
                          metis_vc(num_edges, 0) = i;
                          metis_vc(num_edges, 1) = metis_colidx(k);
                          metis_vc(num_edges, 2) = 0;
                          num_edges ++;

                          metis_vc_score(i) ++;
                          metis_vc_score(metis_colidx(k)) ++;
                        }
                      }
                    }
                  }
                  //for (Int i = 0; i < num_edges; i++) std::cout << " edge(" << metis_vc(i, 0) << ", " << metis_vc(i, 1) << ")" << std::endl;
                  // > look for vertex cover/separator
                  Int num_edges_left = num_edges;
                  while (num_edges_left > 0) {
                    // look for un-processed edge
                    Int next_edge = 0;
                    #if 0
                    while (metis_vc(next_edge, 2) != 0) {
                      next_edge ++;
                    }
                    #else
                    Int max_score = 0;
                    for (Int i = 0; i < num_edges; i++) {
                      Int v1 = metis_vc(i, 0);
                      Int v2 = metis_vc(i, 1);
                      Int score = metis_vc_score(v1) + metis_vc_score(v2);
                      if (metis_vc(i, 2) == 0 && score > max_score) {
                        next_edge = i;
                      }
                    }
                    #endif
                    num_edges_left --;
                    Int v1 = metis_vc(next_edge, 0);
                    Int v2 = metis_vc(next_edge, 1);
                    metis_part_k(v1) = 3;
                    metis_part_k(v2) = 3;
                    metis_vc_score(v1) --;
                    metis_vc_score(v2) --;
                    //std::cout << " >> next edge = " << next_edge << ": edges left = " << num_edges_left << std::endl;
                    // > remove all the incidental edges
                    metis_vc(next_edge, 2) = 2;
                    for (Int i = 0; i < num_edges; i++) {
                      bool edge_removed = false;
                      if (metis_vc(i, 0) == v1 || metis_vc(i, 0) == v2) {
                        if (metis_vc(i, 2) == 0) {
                          num_edges_left --;
                        }
                        //std::cout << "   ++ incident edge = " << i << std::endl;
                        metis_vc(i, 2) ++;
                        edge_removed = true;
                      }
                      if (metis_vc(i, 1) == v1 || metis_vc(i, 1) == v2) {
                        if (metis_vc(i, 2) == 0) {
                          num_edges_left --;
                        }
                        //std::cout << "   -- incident edge = " << i << std::endl;
                        metis_vc(i, 2) ++;
                        edge_removed = true;
                      }
                      if (edge_removed) {
                        metis_vc_score(metis_vc(i, 0)) --;
                        metis_vc_score(metis_vc(i, 1)) --;
                      }
                    }
                  }
                  for(Int i = 0; i < metis_size_k; i++) {
                    if (metis_part_k(i) == -2) {
                      // put it back to dom-0
                      metis_part_k(i) = 0;
                    }
                    if (metis_part_k(i) == 2) {
                      // put it back to dom-1
                      metis_part_k(i) = 1;
                    }
                    if (metis_part_k(i) == 3) {
                      // put it to vertex separator
                      metis_part_k(i) = 2;
                    }
                  }
                  #endif
                } else {
                  // find vertex separator
                  if (Options.verbose == BASKER_TRUE) {
                    std::cout << std::endl << " > calling METIS_ComputeVertexSeparator on leaf " << leaf_id
                              << " at level " << level << " < " << std::endl;
                  }
                  info = METIS_ComputeVertexSeparator(&metis_size_k,
                                                      &(metis_rowptr(0)),
                                                      &(metis_colidx(0)),
                                                       nullptr,
                                                       options,
                                                      &sepsize,
                                                      &(metis_part_k(0)));
                  if (info != METIS_OK) {
                    std::cout << std::endl << " > METIS_ComputeVertexSeparator failed < " << std::endl << std::endl;
                    return BASKER_ERROR; // TODO: what to do here?
                  }
                }
                time_metis += timer_metis.seconds();

                Int dom1 = 0;
                Int dom2 = 0;
                Int sep  = 0;
                for(Int i = 0; i < metis_size_k; i++)
                {
                  if (metis_part_k[i] == 0) {
                    dom1 ++;
                  } else if (metis_part_k[i] == 1) {
                    dom2 ++;
                  } else {
                    sep ++;
                  }
                }
                if(Options.verbose == BASKER_TRUE) {
                  std::cout << " METIS: info = " << info << "(okay = " << METIS_OK << ")"
                            << " size = " << metis_size_k << " sepsize = " << sepsize << std::endl;
                  std::cout << " dom1=(id=" << dom_id1 << ", size=" << dom1 << "),"
                            << " dom2=(id=" << dom_id2 << ", size=" << dom2 << "),"
                            << "  sep=(id=" << sep_id  << ", size=" << sep  << ")" << std::endl;
                  if (dom1 == 0 || dom2 == 0) {
                    std::cout << std::endl << " > METIS returned an empty domain  "
                                           << dom1 << " + " << dom2 << " + " << sep
                                           << " (n = " << M.nrow << ")"
                                           << std::endl << std::endl;
                  }
                }

                // update num doms
                if (dom1 > 0) {
                  sg.cblk ++;
                }
                if (dom2 > 0) {
                  sg.cblk ++;
                }
                if (sep == 0) {
                  sg.cblk --;
                }
                // update permtab, permtab maps original to ND order (permtab[i] is the new row index after nd)
                sep  = frow + dom1 + dom2;
                dom2 = frow + dom1;
                dom1 = frow + 0;
                for(Int i = 0; i < metis_size_k; i++)
                {
                  // i_g is the original global index before ND
                  Int i_g = sg.peritab[frow + i];
                  if (metis_part_k[i] == 0) {
                    sg.permtab[i_g] = dom1;
                    dom1 ++;
                  } else if (metis_part_k[i] == 1) {
                    sg.permtab[i_g] = dom2;
                    dom2 ++;
                  } else {
                    sg.permtab[i_g] = sep;
                    sep ++;
                  }
                }

                // update global part
                for(Int i = 0; i < metis_size_k; i++) {
                  Int i_g = sg.peritab[frow + i];
                  if (metis_part_k[i] == 0) {
                    metis_part[i_g] = dom_id1;
                  } else if (metis_part_k[i] == 1) {
                    metis_part[i_g] = dom_id2;
                  } else {
                    metis_part[i_g] = sep_id;
                  }
                }

                // update peritab
                for(Int i = 0; i < sg.m; i++)
                {
                  sg.peritab[sg.permtab[i]] = i;
                }

                // update rangtab
                sg.rangtab[dom_id1+1] = dom1;
                sg.rangtab[dom_id2+1] = dom2;
                sg.rangtab[sep_id+1] = sep;

                // update postorder tree
                sg.treetab[dom_id1] = sep_id;
                sg.treetab[dom_id2] = sep_id;
                //printf( " -> rangtab[%d] = %d, rangtab[%d] = %d, rangtab[%d]=%d\n",dom_id1+1,dom1,dom_id2+1,dom2,sep_id+1,sep);
              } // if level == num_levels
            } // for each leaves at this level
          } // for each level

          sg.cblk = num_doms;
        }
        // --------------------------------------------------- 
        // done with ND using METIS
        if(Options.verbose == BASKER_TRUE) {
          std::cout << std::endl << " > Time to call METIS : " << time_metis << std::endl;
        }
        #if 0
        printf( " sg.cblk = %d\n",sg.cblk );
        for (Int i = 0; i < num_doms; i++) {
          printf( " - post_order(%d) = %d, post_iorder(%d) = %d\n",i,post_order(i), i,post_iorder(i) );
        }
        for(Int i = 0; i <= sg.cblk; i++)
        {
          printf( " + row_tabs[%d] = %d, col_tabs[%d] = %d, treetab[%d] = %d",i,sg.rangtab[i], i,sg.rangtab[i], i,sg.treetab[i] );
          if (i > 0) {
            printf( " (%d)", sg.rangtab[i]-sg.rangtab[i-1] );
          }
          printf("\n");
        }
        /*for (Int i = 0; i < M.nrow; i++) {
          printf( " > permtab[%d] = %d, peritab[%d] = %d\n",i,sg.permtab[i],i,sg.peritab[i] );
        }*/
        #endif
        #endif
      } else
      { // using SCOTCH
        #if !defined(HAVE_SHYLU_NODEBASKER_SCOTCH)
        BASKER_ASSERT(false, ">> BASKER ASSERT: Scotch is not enabled <<");
        #else
        //----------------------INIT Scotch Graph------------//
        sg.Ap = (scotch_integral_type *)malloc((sg.m+1)     *sizeof(scotch_integral_type));
        sg.Ai = (scotch_integral_type *)malloc((M.nnz)      *sizeof(scotch_integral_type));

        sg.Ap[0] = 0;
        Int sj;
        Int sptr = 0;
        Int self_edge = 0; //If we do not have them, the matrix order will be bad
        for(Int i = 0; i < sg.m; i++)
        {
          sj=0;
          for(Int k = M.col_ptr(i); k <M.col_ptr(i+1); k++)
          {
            if(M.row_idx(k) != i)
            {
              BASKER_ASSERT(sptr < M.nnz);
              sg.Ai[sptr++] = M.row_idx(k);
              sj++;
            }
            else
            {
              self_edge++;
            }
          }
          sg.Ap[i+1] = sg.Ap[i]+sj;
        }
        sg.nz = sg.Ap[sg.m];

        //printf("num self_edge: %d sg.m: %d \n",
        //	   self_edge, sg.m);
        if(self_edge != (sg.m))
        {
          BASKER_ASSERT(self_edge == (sg.m-1), 
              "ZERO ON DIAGONAL, SCOTCH FAIL\n");
          //JDB: comeback need to have a better way to exit
          exit(0);
          //Need to clean up this 
        }

        for(Int i =0; i < sg.m; i++)
        {
          sg.permtab[i] = 0;
          sg.peritab[i] = 0;
        }
        for(Int i =0; i < num_doms; i++)
        {
          sg.rangtab[i] = 0;
          sg.treetab[i] = 0;
        }
        sg.rangtab[num_doms] = 0;

        SCOTCH_Strat strdat;
        SCOTCH_Graph cgrafptr;
        int err;
        scotch_integral_type *vwgts; vwgts = NULL;
    
        if(SCOTCH_graphInit(&cgrafptr) != 0)
        {
          printf("Scotch: error initalizing graph \n");
          return BASKER_ERROR; // TODO: what to do here?
        }

        if( SCOTCH_graphBuild(&cgrafptr, 0, sg.m, sg.Ap, sg.Ap+1, vwgts, NULL,
                              sg.nz, sg.Ai, NULL) !=0 )
        {
          printf("Scotch: failed to build scotch graph \n");
          return BASKER_ERROR; // TODO: what to do here?
        }

        //Need to come back to this so update based on nthreads
        Int flagval = SCOTCH_STRATLEVELMAX | SCOTCH_STRATLEVELMIN | SCOTCH_STRATLEAFSIMPLE | SCOTCH_STRATSEPASIMPLE;
        double balrat = 0.2;
        SCOTCH_stratInit(&strdat);
        err = SCOTCH_stratGraphOrderBuild(&strdat, flagval, num_levels, balrat);

        /*
        err = SCOTCH_stratGraphOrderBuild(&strdat, 
                                          SCOTCH_STRATLEVELMAX | SCOTCH_STRATLEVELMIN,
                                          num_levels, 0.2);
        */

        if(err != 0)
        {
          printf("Scotch: cannot build strategy \n");
          return BASKER_ERROR; // TODO: what to do here?
        }

        // permtab[i] = dom id
        if(SCOTCH_graphOrder(&cgrafptr, &strdat, sg.permtab, sg.peritab, 
                             &sg.cblk, sg.rangtab, sg.treetab) != 0)
        {
          printf("Scotch: cannot compute ordering \n");
          return BASKER_ERROR; // TODO: what to do here?
        }

        if(Options.verbose == BASKER_TRUE) {
          std::cout << " calling SCOTCH_graphOrder(" << M.nrow << " x " << M.ncol 
                    << ", num_levels = " << num_domains << ")" << std::endl;
        }
        free(sg.Ap);
        free(sg.Ai);

        SCOTCH_stratExit(&strdat);
        SCOTCH_graphFree(&cgrafptr);
        #endif
      } // end of SCOTCH
      #else
      BASKER_ASSERT(false, ">> BASKER ASSERT: needs Metis or Scotch to run with multiple threasds <<");
      #endif
    }
    if(Options.verbose == BASKER_TRUE) {
      double time_scotch = timer_scotch.seconds();
      std::cout << " > Time to compute ND : " << time_scotch << std::endl << std::endl;
      for(Int i = 0; i < sg.cblk; i++) {
        printf( " dom-%d : size = %d (%d:%d), tab = %d\n",(int)i,(int)(sg.rangtab[i+1]-sg.rangtab[i]),(int)sg.rangtab[i],(int)(sg.rangtab[i+1]-1),(int)(sg.treetab[i]) );
      }
    }

    //Scan see how many -1
    Int num_trees = 0;
    for(Int i = 0; i < sg.cblk; i++)
    {
      if(sg.treetab[i] == -1)
      {
        num_trees++;
      }
    }
    
    #ifdef BASKER_DEBUG_ORDER_SCOTCH
     printf("FIX SCOTCH PRINT OUT\n");
     printf("SCOTCH: NUM_LEVELS ASKED = %d,  NUM DOMS GOT = %d, NUM TREES = %d \n",
	    num_levels, sg.cblk, num_trees);
     printf("\n");
     printf("%d %d should blks: %f \n",
	    2, ((Int)num_levels+1),
	    pow(2.0,((double)num_levels+1))-1);
    #endif
     
     if(((sg.cblk) != pow(2.0,((double)num_levels+1))-1) || (num_trees != 1))
     {
       //printf("ERROR:  SCOTCH DID NOT PROVIDE A SET BASED ON BISECTION \n");

       Int iblks = pow(2, num_levels+1)-1;

       #ifdef BASKER_DEBUG_ORDER_SCOTCH
       printf("lvl: %d iblks: %d \n", num_levels, iblks);
       #endif

       INT_1DARRAY ttree;
       BASKER_ASSERT((iblks+1) > 0, "scotch iblks");
       MALLOC_INT_1DARRAY(ttree, iblks+1);
       init_value(ttree, iblks+1,(Int) -1);
       INT_1DARRAY ttabs;
       MALLOC_INT_1DARRAY(ttabs, iblks+1);
       init_value(ttabs, iblks+1, (Int) M.ncol);

       for(Int i = 0; i < sg.cblk; i++)
       {
         ttree(i) = sg.treetab[i];
       }

       for(Int i = 0; i < sg.cblk+1; i++)
       {
         ttabs(i) = sg.rangtab[i];
       }

       #ifdef BASKER_DEBUG_ORDER_SCOTCH
       printf("\n\n Starting DEBUG COMPLETE OUT \n\n");
       printf("Tree: ");
       `	for(Int i = 0; i < iblks+1; i++)
       {
         printf("%d, ", ttree(i));
       }
       printf("\n");
       printf("Tabs: ");
       for(Int i = 0; i < iblks+1; i++)
       {
         printf("%d, ", ttabs(i));
       }
       printf("\n");
       printf("\n");
       #endif

       if(Options.verbose == BASKER_TRUE) {
         printf(" > calling to_complete_tree (cblk = %d) <\n",sg.cblk );
       }
       to_complete_tree( num_levels, iblks, sg.cblk,
                         ttabs, ttree );


       #ifdef BASKER_DEBUG_ORDER_SCOTCH
       printf("\n\n DEBUG COMPLETE OUT \n\n");
       printf("Tree: ");
       for(Int i = 0; i < iblks+1; i++)
       {
         printf("%d, ", ttree(i));
       }
       printf("\n");
       printf("Tabs: ");
       for(Int i = 0; i < iblks+1; i++)
       {
         printf("%d, ", ttabs(i));
       }
       printf("\n");
       #endif

       //copy back into scotch
       sg.cblk = iblks;
       for(Int i =0; i < iblks; i++)
       {
         sg.treetab[i] = ttree(i);
         sg.rangtab[i] = ttabs(i);
       }
       sg.rangtab[iblks] = ttabs(iblks);
     }
     #ifdef BASKER_DEBUG_ORDER_SCOTCH
     printf("SCOTCH: NUM LEVELS ASKED = %d,  NUM DOMS GOT = %d \n",
             num_levels, sg.cblk);
     #endif

    //Find the leaf nad non-leaf nodes
    //Int is_nonleaf[sg.cblk];
    Int *is_nonleaf = new Int[sg.cblk];
    
    for(Int i = 0; i < sg.cblk; i++)
    { is_nonleaf[i] = 0; }
    
    for(Int i = 0; i < sg.cblk; i++)
    {
      if(sg.treetab[i] != -1)
      {
        is_nonleaf[sg.treetab[i]] = 1;
      }
    }

    #ifdef BASKER_DEBUG_SCOTCH
    printf("\n\n");
    printf("-----------------------SCOTCH--------------------\n");
    printf("\n\n");
    printf("Scotch Nodes: %d \n", sg.cblk);
    printf("Scotch tree: \n");
    for(Int i = 0; i < sg.cblk; i++)
    {
      printf("%d, ", sg.treetab[i]);
    }
    printf("\n");
    printf("Scotch rangtab: \n");
    for(Int i=0; i < sg.cblk+1; i++)
    {
      printf("%d, ", sg.rangtab[i]);
    }
    printf("\n");
    printf("Scotch Parts \n");
    Int p;
    Int part = 0;
    for(Int i = 0 ; i < sg.cblk; i++)
    {
      printf("Column in ");
      if(is_nonleaf[i])
        printf("interior node (ID: %d ): ", i);
      else
        printf("leaf node (ID: %d ): ", i);

      for(Int j = sg.rangtab[i]; j < sg.rangtab[i+1]; j++)
      {
        printf("%d, ", sg.peritab[j]);
      }
      for(Int j = sg.rangtab[i]; j < sg.rangtab[i+1]; j++)
      {
        printf(" (%d , %d) , ", sg.peritab[j]+1, p ); 
      }
      printf("\n");
    }
    #endif

    //Copy into a temp basker tree
    BT.nblks = sg.cblk;
    BASKER_ASSERT((BT.nblks+1) > 0, "scotch bt.nblks+1");
    MALLOC_INT_1DARRAY(BT.row_tabs, BT.nblks+1);
    init_value(BT.row_tabs, BT.nblks+1, M.max_idx);
    MALLOC_INT_1DARRAY(BT.col_tabs, BT.nblks+1);
    init_value(BT.col_tabs, BT.nblks+1, M.max_idx);
    BASKER_ASSERT(M.nrow > 0, "scotch M.nrow");
    MALLOC_INT_1DARRAY(BT.permtab, M.nrow);
    init_value(BT.permtab, M.nrow, M.max_idx);
    MALLOC_INT_1DARRAY(BT.ipermtab, M.nrow);
    init_value(BT.ipermtab, M.nrow, M.max_idx);
    MALLOC_INT_1DARRAY(BT.treetab, BT.nblks+1);
    init_value(BT.treetab, BT.nblks+1, M.max_idx);
    
    for(Int i = 0; i < BT.nblks+1; i++)
    {
      BT.row_tabs[i] = sg.rangtab[i];
      BT.col_tabs[i] = sg.rangtab[i];
      BT.treetab[i]  = sg.treetab[i];
    }

    //printf( " + permtab, peritab\n" );
    //FILE *fp = fopen("permtab_p.dat","w");
    for(Int i = 0; i < M.nrow; i++)
    {
      BT.permtab[i] = sg.permtab[i];
      BT.ipermtab[i] = sg.peritab[i];
      //fprintf(fp, " %d, %d\n",BT.permtab[i],BT.ipermtab[i] );
    }
    //fclose(fp);

    //Used for recursing easier
    BT.treetab[BT.nblks-1]   = BT.nblks;
    BT.treetab[BT.nblks]     = -1;

    //Right now defaulting parts to two
    BT.nparts = 2;

    free(sg.permtab);
    free(sg.peritab);
    free(sg.rangtab);
    free(sg.treetab);
    delete [] is_nonleaf;

    return BASKER_SUCCESS;
  }//end part_scotch()


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::to_complete_tree
  (
   Int lvl,
   Int iblks,
   Int nblks,
   INT_1DARRAY  tabs,
   INT_1DARRAY _tree
  )
  {
    //Goal is to turn the incomplete tree 
    //Scotch gives to a complete tree we need

    //This might not be the best way
    //This scans in a post-order and adds to the tree 
    //as need

    //1. Prealloc output
    Int onblks;
    INT_1DARRAY otabs;
    INT_1DARRAY otree;
    onblks = iblks;
    BASKER_ASSERT((onblks+1) > 0, "scotch onblks");
    MALLOC_INT_1DARRAY(otabs, onblks+1);
    init_value(otabs, onblks+1, (Int) -1);
    MALLOC_INT_1DARRAY(otree, onblks+1);
    init_value(otree, onblks+1, (Int) -1);

    
    Int lpos = 0;
    Int rpos = iblks-1;
    Int mynum = iblks-1;
    otree(iblks) = -1;
    rec_build_tree(lvl, 
		   lpos,rpos, 
		   mynum,
		   otree);

		 
    INT_1DARRAY ws;
    BASKER_ASSERT((iblks+1)>0, "scotch iblks 2");
    MALLOC_INT_1DARRAY(ws, iblks+1);
    init_value(ws, iblks+1, (Int) 0);


    //DEBUG
    if(Options.verbose == BASKER_TRUE)
    {
      printf("test - otree\n");
      for(Int t_blk = 1; t_blk < iblks+1; t_blk++)
      {
        printf("%ld,", (long)otree(t_blk));
      }
      printf("\n");
      printf("nblks %ld \n", (long)nblks);
    }


    //test if enough domain
    //this is a huge error that we need to take care of
    Int ndomains   = 0;
    Int indomains  = 0;
    //Int c_treenode = tree(0);
    //scan over all and count set of pairs
    
    for(Int t_blk = 1; t_blk < nblks; t_blk++)
    {
      if(_tree(t_blk-1) == _tree(t_blk))
      {
        ndomains++;
      } 
    }
    for(Int t_blk = 1; t_blk < iblks+1; t_blk++)
    {
      if(otree(t_blk) == -1)
      {
        break;
      }
      if(otree(t_blk-1) == otree(t_blk))
      {
        indomains++;
      }
    }

    Int bdomains = pow(2, lvl+1)-1; // I am hoping for binary tree..
    if(Options.verbose == BASKER_TRUE)
    {
      printf("Domains Found: %ld \n", (long)ndomains);
      printf("Domains Ideal: %ld \n", (long)indomains);
      printf("Domains Binary: %ld \n", (long)ndomains);
    } 

    if(ndomains != indomains || ndomains != bdomains)
    {
      printf(" ShyLU Basker Error: do_complete_tree routine \n");
      printf(" > ERROR: NOT ENOUGH DOMAINS FOR THREADS\n");
      printf(" > REDUCE THREAD COUNT AND TRY AGAIN\n\n");
      printf("  This error occurs when the matrix to be solved is too small for given number of threads\n \
          The number of threads must match the number of leaves in the tree\n \
          To resolve this, rerun with fewer threads\n");
      // Make this throw exception instead
      exit(EXIT_FAILURE);
    }

    //scan correct
    Int s_tree_p = 0;
    Int m_tree_p = 0;
    Int s_tab_p = 1;
    Int m_tab_p = 1;
    
    otabs(0) = 0;

    #ifdef BASKER_DEBUG_ORDER_SCOTCH
    printf("\n Start Debug Print, %d\n", m_tree_p);
    printf("WS: ");
    for(Int i=0; i< iblks; i++)
    {
      printf("%d, ", ws(i));
    }
    printf("\n");
    printf("IN Tree: ");
    for(Int i=0; i < iblks+1; i++)
    {
      printf("%d, ", _tree(i));
    }
    printf("\n");
    printf("Out Tree: ");
    for(Int i=0; i < iblks+1; i++)
    {
      printf("%d, ", otree(i));
    }
    printf("\n");
    printf("Tabs: ");
    for(Int i =0; i < iblks+1; i++)
    {
      printf("%d, ", otabs(i));
    }
    printf("\n");
    #endif

    for(m_tree_p = 0; m_tree_p < iblks; m_tree_p++)
    {

    #ifdef BASKER_DEBUG_ORDER_SCOTCH
      printf("top of loop: %d \n",
          _tree(s_tree_p));
    #endif

      if(ws(m_tree_p) == 0)
      {
        //Not assigned yet
        if((_tree(s_tree_p) == otree(m_tree_p)))
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("same case\n");
        #endif
          otabs(m_tab_p) = tabs(s_tab_p);
        }
        else if(_tree(s_tree_p) == -1)
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("-1 case \n");
        #endif
          _tree(s_tree_p) = otree(m_tree_p);
          otabs(m_tab_p) = tabs(s_tab_p);
        }
        else if(_tree(s_tree_p) == 0)
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("end case\n");
        #endif
          _tree(s_tree_p) = otree(m_tree_p);
          otabs(m_tab_p) = tabs(s_tab_p-1);
          s_tab_p--;
        }
        else
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("no sep \n");
          printf("tree: %d %d otree: %d %d \n",
              s_tree_p, _tree(s_tree_p),
              m_tree_p, otree(m_tree_p));
        #endif
          if(ws(otree(m_tree_p)) == 0)
          {
            //need to offset to make space
            for(Int jj = iblks; jj > otree(m_tree_p); jj--)
            {
              if((_tree(jj-1) == 0) || (_tree(jj-1) ==-1))
              {
                _tree(jj) = _tree(jj-1);
                #ifdef BASKER_DEBUG_ORDER_SCOTCH
                printf("sliding1: %d %d %d \n",
                    jj, jj-1, _tree(jj-1));
                #endif
              }
              else
              {
                _tree(jj) = _tree(jj-1)+1;
                #ifdef BASKER_DEBUG_ORDER_SCOTCH
                printf("sliding2: %d %d %d \n",
                    jj, jj-1, _tree(jj-1));
                #endif
              }
            }//end for-over upper

            _tree(otree(m_tree_p)) = otree(otree(m_tree_p));
            ws(otree(m_tree_p)) = 1;

          }//ws == 0;

          _tree(s_tree_p) = otree(m_tree_p);
          otabs(m_tab_p) = tabs(s_tab_p);
        }//eles

        ws(m_tree_p) = -1;
        s_tree_p++;
        m_tab_p++;
        s_tab_p++;

      }//ws == 0
      else if(ws(m_tree_p) == 1)
      {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
        printf("Fix sep\n");
        #endif
        //Note assigned sep
        otabs(m_tab_p) = tabs(s_tab_p-1);
        _tree(s_tree_p) = otree(m_tree_p);
        s_tree_p++;
        m_tab_p++;
      }//ws == 1
      else
      {
        //Should not go here
        BASKER_ASSERT(0==1, "ERROR");
      }

      #ifdef BASKER_DEBUG_ORDER_SCOTCH
      printf("\n Debug Print, %d\n", m_tree_p);
      printf("WS: ");
      for(Int i=0; i< iblks; i++)
      {
        printf("%d, ", ws(i));
      }
      printf("\n");
      printf("Tree: ");
      for(Int i=0; i < iblks+1; i++)
      {
        printf("%d, ", _tree(i));
      }
      printf("\n");
      printf("Tabs: ");
      for(Int i =0; i < iblks+1; i++)
      {
        printf("%d, ", otabs(i));
      }
      printf("\n");
      #endif

    }//for over all

    #ifdef BASKER_DEBUG_ORDER_SCOTCH
    printf("\n Debug Print\n");
    printf("WS: ");
    for(Int i=0; i< iblks; i++)
    {
      printf("%d, ", ws(i));
    }
    printf("\n");
    printf("Tree: ");
    for(Int i=0; i < iblks+1; i++)
    {
      printf("%d, ", _tree(i));
    }
    printf("\n");
    printf("Tabs: ");
    for(Int i =0; i < iblks+1; i++)
    {
      printf("%d, ", otabs(i));
    }
    printf("\n");
    #endif

    //copy back
    for(Int i = 0; i < iblks+1; i++)
    {
      tabs(i) = otabs(i);
    }
    
  }//end to_complete_tree


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::rec_build_tree
  (
   Int lvl,
   Int &lpos, Int &rpos, 
   Int &mynum,
   INT_1DARRAY _tree
  )
  {
    //printf("assign, lpos: %d rpos: %d  number: %d\n",
    //	   lpos, rpos, mynum);
    
    if(lvl > 0)
    {
      Int rightc = rpos -1;
      //Int leftc  = rightc - pow(2,lvl-1)-1;
      Int leftc = (rpos+lpos-1)/2;

      #ifdef BASKER_DEBUG_ORDER_SCOTCH
      printf("Left Child: %d Right Child: %d \n",
          leftc, rightc);
      printf("lpos: %d rpos: %d  lvl: %d \n",
          lpos, rpos, lvl);
      #endif

      _tree(leftc)  = mynum;
      _tree(rightc) = mynum;

      #ifdef BASKER_DEBUG_ORDER_SCOTCH
      printf("assign: %d %d \n", leftc, mynum);
      printf("assign: %d %d \n", rightc,mynum);
      #endif

      mynum = rightc;
      rec_build_tree(lvl-1, 
          leftc, rightc,
          mynum,
          _tree);

      mynum = leftc;
      rec_build_tree(lvl-1,
          lpos, leftc, 
          mynum,
          _tree);

    } // end if lvl > 0

  }//end rec_build_tree

}//end namespace Basker

#undef BASKER_TIMER
#endif //ifndef basker_order_scotch_hpp
