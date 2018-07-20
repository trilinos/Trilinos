#ifndef SHYLUBASKER_ORDER_SCOTCH_HPP
#define SHYLUBASKER_ORDER_SCOTCH_HPP

#include "shylubasker_types.hpp"
#include "scotch.h"

//#define BASKER_DEBUG_ORDER_SCOTCH

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
   BASKER_MATRIX &C
  )
  {
    BASKER_MATRIX T;
    //get matrix transpose
    matrix_transpose(M,T);
   
    C.set_shape(M.srow, M.nrow, M.scol, M.ncol);

    BASKER_ASSERT((M.ncol+1) > 0, "scotch ncol"); 
    MALLOC_INT_1DARRAY(C.col_ptr, M.ncol+1);
    init_value(C.col_ptr, M.ncol+1, (Int) 0);
    BASKER_ASSERT(M.nnz > 0, "scotch nnz");
    MALLOC_INT_1DARRAY(C.row_idx, 2*M.nnz);
    init_value(C.row_idx, 2*M.nnz, (Int) 0);

    
    INT_1DARRAY ws;
    BASKER_ASSERT(M.nrow > 0, "scotch nrow");
    MALLOC_INT_1DARRAY(ws, M.nrow);
    init_value(ws, M.nrow, BASKER_MAX_IDX);
    

    Int c_nnz = 0;
    for(Int k = 0 ; k <  M.ncol; ++k)
    {
      //scatter M
      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        Int j = M.row_idx(i);
        if(ws(j) != k)
        {
          C.row_idx(c_nnz) = j;
          c_nnz++;
          ws(j) = k;
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

    //can be remove after debuggin
    BASKER_ASSERT(c_nnz > 0, "C.val scotch");
    MALLOC_ENTRY_1DARRAY(C.val,c_nnz);
    init_value(C.val,c_nnz, (Entry)1.0);
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
    part_scotch(M, BT, lvls);
    return 0;
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
  //----------------------INIT Scotch Graph------------//
  #if SHYLU_SCOTCH_64
    using scotch_integral_type = int64_t; //NDE: make this depend on the scotch type
  #else
    using scotch_integral_type = int32_t; //NDE: make this depend on the scotch type
  #endif
    scotch_graph<scotch_integral_type> sg;
    sg.m = M.nrow;
    sg.Ap = (scotch_integral_type *)malloc((sg.m+1)     *sizeof(scotch_integral_type));
    sg.Ai = (scotch_integral_type *)malloc((M.nnz)      *sizeof(scotch_integral_type));
    sg.permtab = (scotch_integral_type *)malloc((sg.m)  *sizeof(scotch_integral_type));
    sg.peritab = (scotch_integral_type *)malloc((sg.m)  *sizeof(scotch_integral_type));
    sg.rangtab = (scotch_integral_type *)malloc((sg.m+1)*sizeof(scotch_integral_type));
    sg.treetab = (scotch_integral_type *)malloc((sg.m)  *sizeof(scotch_integral_type));
    
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
          ASSERT(sptr < M.nnz);
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
      sg.rangtab[i] = 0;
      sg.treetab[i] = 0;
    }
    sg.cblk = 0;

    SCOTCH_Strat strdat;
    SCOTCH_Graph cgrafptr;
    int err;
    scotch_integral_type *vwgts; vwgts = NULL;
    
    if(SCOTCH_graphInit(&cgrafptr) != 0)
    {
      printf("Scotch: error initalizing graph \n");
      return -1;
    }

    if( SCOTCH_graphBuild(&cgrafptr, 0, sg.m, sg.Ap, sg.Ap+1, vwgts, NULL,
          sg.nz, sg.Ai, NULL) !=0 )
    {
      printf("Scotch: failed to build scotch graph \n");
      return -1;
    }

    //Need to come back to this so update based on nthreads
    
    Int num_levels = num_domains; 
    SCOTCH_stratInit(&strdat);
    
    err = SCOTCH_stratGraphOrderBuild(&strdat, 
				      SCOTCH_STRATLEVELMAX | SCOTCH_STRATLEVELMIN | SCOTCH_STRATLEAFSIMPLE | SCOTCH_STRATSEPASIMPLE,
		  num_levels, 0.2);
    

    /*
    err = SCOTCH_stratGraphOrderBuild(&strdat, 
				      SCOTCH_STRATLEVELMAX | SCOTCH_STRATLEVELMIN,
		  num_levels, 0.2);
    */

    if(err != 0)
    {
      printf("Scotch: cannot build strategy \n");
      return -1;
    }

    if(SCOTCH_graphOrder(&cgrafptr, &strdat, sg.permtab, sg.peritab, 
          &sg.cblk, sg.rangtab, sg.treetab) != 0)
    {
      printf("Scotch: cannot compute ordering \n");
      return -1;
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
     printf("SCOTCH: ASKED: %d  GOT : %d TREES: %d \n",
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

       to_complete_tree( num_levels,iblks, sg.cblk,
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
    printf("SCOTCH: ASKED: %d  GOT : %d \n",
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

    for(Int i = 0; i < M.nrow; i++)
    {
      BT.permtab[i] = sg.permtab[i];
      BT.ipermtab[i] = sg.peritab[i];
    }

    //Used for recursing easier
    BT.treetab[BT.nblks-1]   = BT.nblks;
    BT.treetab[BT.nblks]     = -1;

    //Right now defaulting parts to two
    BT.nparts = 2;

    free(sg.Ap);
    free(sg.Ai);
    free(sg.permtab);
    free(sg.peritab);
    free(sg.rangtab);
    free(sg.treetab);
    free(is_nonleaf);

    SCOTCH_stratExit(&strdat);
    SCOTCH_graphFree(&cgrafptr);

    return 0;
  }//end part_scotch()


  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::to_complete_tree
  (
   Int lvl,
   Int iblks,
   Int nblks,
   INT_1DARRAY tabs,
   INT_1DARRAY tree
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
      if(tree(t_blk-1) == tree(t_blk))
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

    if(Options.verbose == BASKER_TRUE)
    {
      printf("Domains Found: %ld \n", (long)ndomains);
      printf("Domains Ideal: %ld \n", (long)indomains);
    }  

    if(ndomains != indomains)
    {
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Domains Found: %ld \n", (long)ndomains);
        printf("Domains Ideal: %ld \n", (long)indomains);
        printf("ERROR: NOT ENOUGH DOMAINS FOR THREADS\n");
        printf("REDUCE THREAD COUNT AND TRY AGAIN\n");
      }
      printf(" ShyLU Basker Error: do_complete_tree routine \n");
      printf("   num domains != ideal num domains\n");
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
      printf("%d, ", tree(i));
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
          tree(s_tree_p));
    #endif

      if(ws(m_tree_p) == 0)
      {
        //Not assigned yet
        if((tree(s_tree_p) == otree(m_tree_p)))
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("same case\n");
        #endif
          otabs(m_tab_p) = tabs(s_tab_p);
        }
        else if(tree(s_tree_p) == -1)
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("-1 case \n");
        #endif
          tree(s_tree_p) = otree(m_tree_p);
          otabs(m_tab_p) = tabs(s_tab_p);
        }
        else if(tree(s_tree_p) == 0)
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("end case\n");
        #endif
          tree(s_tree_p) = otree(m_tree_p);
          otabs(m_tab_p) = tabs(s_tab_p-1);
          s_tab_p--;
        }
        else
        {
        #ifdef BASKER_DEBUG_ORDER_SCOTCH
          printf("no sep \n");
          printf("tree: %d %d otree: %d %d \n",
              s_tree_p, tree(s_tree_p),
              m_tree_p, otree(m_tree_p));
        #endif
          if(ws(otree(m_tree_p)) == 0)
          {
            //need to offset to make space
            for(Int jj = iblks; jj > otree(m_tree_p); jj--)
            {
              if((tree(jj-1) == 0) || (tree(jj-1) ==-1))
              {
                tree(jj) = tree(jj-1);
                #ifdef BASKER_DEBUG_ORDER_SCOTCH
                printf("sliding1: %d %d %d \n",
                    jj, jj-1, tree(jj-1));
                #endif
              }
              else
              {
                tree(jj) = tree(jj-1)+1;
                #ifdef BASKER_DEBUG_ORDER_SCOTCH
                printf("sliding2: %d %d %d \n",
                    jj, jj-1, tree(jj-1));
                #endif
              }
            }//end for-over upper

            tree(otree(m_tree_p)) = otree(otree(m_tree_p));
            ws(otree(m_tree_p)) = 1;

          }//ws == 0;

          tree(s_tree_p) = otree(m_tree_p);
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
        tree(s_tree_p) = otree(m_tree_p);
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
        printf("%d, ", tree(i));
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
      printf("%d, ", tree(i));
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
   INT_1DARRAY tree
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

      tree(leftc)  = mynum;
      tree(rightc) = mynum;

      #ifdef BASKER_DEBUG_ORDER_SCOTCH
      printf("assign: %d %d \n", leftc, mynum);
      printf("assign: %d %d \n", rightc,mynum);
      #endif

      mynum = rightc;
      rec_build_tree(lvl-1, 
          leftc, rightc,
          mynum,
          tree);

      mynum = leftc;
      rec_build_tree(lvl-1,
          lpos, leftc, 
          mynum,
          tree);

    } // end if lvl > 0

  }//end rec_build_tree

}//end namespace Basker

#endif //ifndef basker_order_scotch_hpp
