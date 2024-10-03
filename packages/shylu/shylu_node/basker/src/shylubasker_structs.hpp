// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_STRUCTS_HPP
#define SHYLUBASKER_STRUCTS_HPP

#include "shylubasker_types.hpp"
#include "shylubasker_vector.hpp"
#include "shylubasker_scalar_traits.hpp"

namespace BaskerNS
{

  //-------------basker - thread------------//
  template <class Int, class Entry, class Exe_Space>
  struct basker_thread
  {
    BASKER_INLINE
    basker_thread()
    {
      iws_size = 0;
      ews_size = 0;
      iws_mult = 5; //Note:come back and reduce
      ews_mult = 2; //Note:come back and reduce


      error_type   = BASKER_ERROR_NOERROR;
      error_blk    = BASKER_MAX_IDX;
      error_subblk = BASKER_MAX_IDX;
      error_info   = BASKER_MAX_IDX;
      

      init_ops();
    }//end basker_thread
    
    BASKER_INLINE
    ~basker_thread()
    {
      Finalize();
    }//end ~basker_thread

    BASKER_INLINE
    void Finalize()
    {
      //Note,
      //basker thread is always used in a array so not need.
      //Added for completeness/OpenMP
      #ifndef BASKER_KOKKOS
      FREE_INT_1DARRAY(iws);
      FREE_ENTRY_1DARRAY(ews);
      C.Finalize();
      #endif
    }

     
    //arrays
    Int iws_size;
    Int ews_size;
    Int iws_mult;
    Int ews_mult;
    INT_1DARRAY iws;
    ENTRY_1DARRAY ews;

    //Each thread gets a column matrix
    BASKER_MATRIX C;
    
    //Manage realloc 
    Int error_type; //Remalloc, pivot, etc
    Int error_blk;
    Int error_subblk;
    Int error_info; //number 
    

    //----------------Depracted----------
    //volatile int **token;
    Int ** volatile token;

    BASKER_INLINE
    void init_token(const Int size)
    {
      printf("token size: %ld \n", (long)size);
      token = (Int **)malloc(size*sizeof(Int*));
      for(Int i = 0; i < size; i++)
      {
        token[i] = (Int *)malloc(8*sizeof(Int));
        token[i][0] = 0;
        token[i][1] = 0;
        token[i][2] = 0;
        token[i][3] = 0;
        token[i][4] = 0;
        token[i][5] = 0;
        token[i][6] = 0;
        token[i][7] = 0;
      }
    }


    Int ops_counts[10][10];

    void init_ops()
    {
      for(Int i=0; i < 10; i++) {
        for(Int j = 0; j < 10; j++) {
          ops_counts[i][j] = 0;
        }
      }
    }

  };

  //------------------------Basker-Tree-------------------------------//
  //Used to store information about the tree
  template <class Int, class Entry, class Exe_Space >
  struct  basker_tree
  {  
    BASKER_INLINE
    basker_tree()
    {
      nroots = 0;
      nparts = 2;
      nblks = 0;
    }
    BASKER_INLINE
    ~basker_tree()
    {
      //Finalize();
    }//end ~basker_tree

    BASKER_INLINE
    void Finalize()
    {
      //printf("basker_tree Finalize todo \n");
      if(nroots > 0)
      {
        FREE_INT_1DARRAY(roots);
      }
      if(nblks > 0)
      {
        FREE_INT_1DARRAY(permtab);
        FREE_INT_1DARRAY(ipermtab);
        FREE_INT_1DARRAY(row_tabs);
        FREE_INT_1DARRAY(col_tabs);
        FREE_INT_1DARRAY(lvlset);
        FREE_INT_1DARRAY(treetab);
        FREE_INT_1DARRAY(lvltreetabs);
        FREE_INT_1DARRAY(treeptr);
        FREE_INT_1DARRAY(rowptr);
        FREE_INT_1DARRAY(child);
        FREE_INT_1DARRAY(sibling);
      }
    }//end Finalize();

    BASKER_INLINE
    void basic_convert
    (
     Int _m, 
     Int *_perm, 
     Int _nblks,
     Int _parts, 
     Int *_row_tabs, 
     Int *_col_tabs,
     Int *_treetab
    )
    {
      nblks = _nblks;
      //For now defaulting parts to 2
      //nparts = _parts;
      nparts = 2;
      //For now do nothing with perm
      //Note: cannt use init_value because circle depend basker
      BASKER_ASSERT((nblks+1)>0, "struct basic convert");
      MALLOC_INT_1DARRAY(row_tabs, nblks+1);
      for(Int i =0; i < nblks+1; i++)
      {
        row_tabs[i] = _row_tabs[i];
      }
      MALLOC_INT_1DARRAY(col_tabs, nblks+1);
      for(Int i=0; i < nblks+1; i++)
      {
        col_tabs[i] = _col_tabs[i];
      }
      MALLOC_INT_1DARRAY(treetab, nblks+1);
      for(Int i=0; i < nblks+1; i++)
      {
        treetab[i] = _treetab[i];
      }
    }//end basic_convert

    BASKER_INLINE
    void print()
    {
      std::cout << std::endl;
      std::cout << "basker_tree info: " << std::endl;
      //std::cout << "nparts: " << nparts << " nblks: " << nblks 
      //        <<  " nlvls: " << nlvls <<std::endl;
      std::cout << "row_tabs: < ";
      for(Int i=0 ; i < nblks+1; i++)
        {std::cout << row_tabs[i] << " " ;}
      std::cout << " > " << std::endl;
      std::cout << "col_tabs: < " ;
      for(Int i=0; i < nblks+1; i++)
        {std::cout << col_tabs[i] << " " ;}
      std::cout << " > " << std::endl;
      std::cout << "treetab: < ";
      for(Int i=0; i < nblks+1; i++)
        {std::cout << treetab[i] << " " ;}
      std::cout << " > " << std::endl;
      std::cout << std::endl;
    }
    
    BASKER_INLINE
    void info()
    {
      print();
    }

    Int nroots;
    INT_1DARRAY roots;
    Int nparts;
    INT_1DARRAY permtab;
    INT_1DARRAY ipermtab;
    Int nblks;
    INT_1DARRAY row_tabs;
    INT_1DARRAY col_tabs;
    Int nlvls;
    INT_1DARRAY lvlset;
    INT_1DARRAY treetab;
    INT_1DARRAY lvltreetabs;
    INT_1DARRAY treeptr;
    BOOL_1DARRAY pivot;
    INT_1DARRAY  rowptr;
    INT_1DARRAY  child;
    INT_1DARRAY  sibling;
  };//end basker_tree


  //---------------------------BASKER-S-TREE-----------------------------------//
  template <class Int, class Entry, class Exe_Space>
  struct basker_symbolic_tree
  {
    typedef basker_vector<Int,Entry,Exe_Space> BV;

    basker_symbolic_tree()
    {
      parent_flg       = 0;
      post_flg         = 0;
      col_ptr_flg      = 0;
      left_most_flg    = 0;
      row_counts_flg   = 0;
      col_counts_flg   = 0;
      WS_flg           = 0;
      S_col_flg        = 0;
      S_row_flg        = 0;
      L_row_counts_flg = 0;
      L_col_counts_flg = 0;
      U_row_counts_flg = 0;
      U_col_counts_flg = 0;
      S_col_counts_flg = 0;
      S_row_counts_flg = 0;
    }

    ~basker_symbolic_tree()
    {
      //Finalize();
    }//end ~basker_symbolic_tree

    BASKER_INLINE
    void Finalize()
    {
      //This can be removed to use only in sfactor stage than delete
      //printf("basker_symbolic_tree finalize todo \n");
      if(parent_flg > 0)
      {
        FREE_INT_1DARRAY(parent);
      }
      if(post_flg > 0)
      {
        FREE_INT_1DARRAY(post);
      }
      if(col_ptr_flg > 0)
      {
        FREE_INT_1DARRAY(col_ptr);
      }
      if(left_most_flg > 0)
      {
        FREE_INT_1DARRAY(left_most);
      }
      if(row_counts_flg > 0)
      {
        FREE_INT_1DARRAY(row_counts);
      }
      if(col_counts_flg > 0)
      {
        FREE_INT_1DARRAY(col_counts);
      }
      if(WS_flg > 0)
      {
        FREE_INT_1DARRAY(WS);
      }
      if(S_col_flg > 0)
      {
        FREE_INT_1DARRAY(S_col);
      }
      if(S_row_flg > 0)
      {
        FREE_INT_1DARRAY(S_row);
      }
      if(L_row_counts_flg > 0)
      {
        FREE_INT_1DARRAY(L_row_counts);
      }
      if(L_col_counts_flg > 0)
      {
        FREE_INT_1DARRAY(L_col_counts);
      }
      if(U_row_counts_flg > 0)
      {
        FREE_INT_1DARRAY(U_row_counts);
      }
      if(U_col_counts_flg > 0)
      {
        FREE_INT_1DARRAY(U_col_counts);
      }
      if(S_col_counts_flg > 0)
      {
        FREE_INT_1DARRAY(S_col_counts);
      }
      if(S_row_counts_flg > 0)
      {
        FREE_INT_1DARRAY(S_row_counts);
      }

    }//end Finalize()
    
    BASKER_INLINE
    void init_parent(Int size)
    {
      if(size <=0)
        //{return;} JDB: Changed because of kokkos
      {size = 1;} //defaul to 1
      if(parent_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct parent size");
        MALLOC_INT_1DARRAY(parent, size);
        parent_flg =size;
      }
      else if(size > parent_flg)
      {
        BASKER_ASSERT(size > 0, "struct parent size2");
        REALLOC_1DARRAY(parent, parent_flg, size);
        parent_flg = size;
      }
      //BV::init_value(parent, parent_flg, 0);
      for(Int j = 0; j < parent_flg; j++)
      {
        parent(j) = (Int) 0;
      }

    }//end init_parent

    BASKER_INLINE
    void init_post(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_post: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(post_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct post size");
        MALLOC_INT_1DARRAY(post, size);
        post_flg = size;
      }
      else if(size > parent_flg)
      {
        BASKER_ASSERT(size > 0, "struct post size2");
        REALLOC_1DARRAY(post, post_flg, size);
        post_flg = size;
      }
      BV::init_value(post, post_flg, 0);
    }//end init_post


    BASKER_INLINE
    void init_col_ptr(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_col_ptr: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(col_ptr_flg == 0)
      {
        BASKER_ASSERT(size > 0, "structs col size");
        MALLOC_INT_1DARRAY(col_ptr, size);
        col_ptr_flg =size;
      }
      else if(size > parent_flg)
      {
        BASKER_ASSERT(size > 0, "structs col size2");
        REALLOC_1DARRAY(col_ptr, col_ptr_flg, size);
        col_ptr_flg = size;
      }
      BV::init_value(col_ptr, col_ptr_flg, 0);
    }//end init_col_ptr


    BASKER_INLINE
    void init_left_most(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_left_most: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(left_most_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct left size");
        MALLOC_INT_1DARRAY(left_most, size);
        left_most_flg =size;
      }
      else if(size > left_most_flg)
      {
        BASKER_ASSERT(size > 0, "struct left size2");
        REALLOC_1DARRAY(left_most, left_most_flg, size);
        left_most_flg = size;
      }
      BV::init_value(left_most, left_most_flg, 0);
    }//end init_left_most


    BASKER_INLINE
    void init_row_counts(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_row_counts: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(row_counts_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct row size");
        MALLOC_INT_1DARRAY(row_counts, size);
        row_counts_flg =size;
      }
      else if(size > row_counts_flg)
      {
        BASKER_ASSERT(size > 0, "struct row size2");
        REALLOC_1DARRAY(row_counts, row_counts_flg, size);
        row_counts_flg = size;
      }
      BV::init_value(row_counts, row_counts_flg, 0);

    }//end init_row_counts
    

    BASKER_INLINE
    void init_col_counts(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_col_counts: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(col_counts_flg == 0)
      {
        //printf("first call init_col_counts\n");
        BASKER_ASSERT(size > 0, "struct col count size");
        MALLOC_INT_1DARRAY(col_counts, size);
        col_counts_flg =size;
      }
      else if(size > col_counts_flg)
      {
        printf("realloc call init_col_counts\n");
        printf("oldsize: %ld newsize: %ld \n",
            (long)col_counts_flg, (long)size);
        BASKER_ASSERT(size > 0, "struct col count size2");
        REALLOC_1DARRAY(col_counts, col_counts_flg, size);
        col_counts_flg = size;
      }
      //printf("zero out col_counts\n");
      BV::init_value(col_counts, col_counts_flg, 0);
    }//end init_col_counts
	

    BASKER_INLINE
    void init_WS(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_WS: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(WS_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct WS size");
        MALLOC_INT_1DARRAY(WS, size);
        WS_flg =size;
      }
      else if(size > WS_flg)
      {
        BASKER_ASSERT(size > 0, "struct WS size2");
        REALLOC_1DARRAY(WS, WS_flg, size);
        WS_flg = size;
      }
      BV::init_value(WS, WS_flg, 0);
    }//end init_WS()


    BASKER_INLINE
    void init_S_col(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_S_col: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(S_col_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct scol size");
        MALLOC_INT_1DARRAY(S_col, size);
        S_col_flg =size;
      }
      else if(size > S_col_flg)
      {
        BASKER_ASSERT(size > 0, "struct scol size2");
        REALLOC_1DARRAY(WS, S_col_flg, size);
        S_col_flg = size;
      }
      BV::init_value(WS, S_col_flg, 0);
    }//end init_S_col()
    

    BASKER_INLINE
    void init_S_row(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_S_row: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(S_row_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct s row size");
        MALLOC_INT_1DARRAY(S_row, size);
        S_row_flg =size;
      }
      else if(size > S_row_flg)
      {
        BASKER_ASSERT(size > 0, "struct s row size2");
        REALLOC_1DARRAY(S_row, S_row_flg, size);
        S_row_flg = size;
      }
      BV::init_value(S_row, S_row_flg, 0);
    }//end init_S_row()


    BASKER_INLINE
    void init_L_row_counts(Int size, BASKER_BOOL verbose)
    {
      if(size <=0)
      {
        if (verbose) {
          std::cout << " init_L_row_counts: size <= 0 no alloc " << std::endl;
        }
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(L_row_counts_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct lrowcount size");
        MALLOC_INT_1DARRAY(L_row_counts, size);
        L_row_counts_flg =size;
      }
      else if(size > L_row_counts_flg)
      {
        BASKER_ASSERT(size > 0, "struct lrowcount size2");
        REALLOC_1DARRAY(L_row_counts, L_row_counts_flg, size);
        L_row_counts_flg = size;
      }
      BV::init_value(L_row_counts, L_row_counts_flg, 0);

    }//end init_L_row_counts

    BASKER_INLINE
    void init_L_col_counts(Int size, BASKER_BOOL verbose)
    {
      if(size <=0)
      {
        if (verbose) {
          std::cout << " init_L_col_counts: size <= 0 no alloc " << std::endl;
        }
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(L_col_counts_flg == 0)
      {
        BASKER_ASSERT(size > 0, "strcut lcolcountsize");
        MALLOC_INT_1DARRAY(L_col_counts, size);
        L_col_counts_flg =size;
      }
      else if(size > L_col_counts_flg)
      {
        BASKER_ASSERT(size > 0, "strcut lcolcountsize2");
        REALLOC_1DARRAY(L_col_counts, L_col_counts_flg, size);
        L_col_counts_flg = size;
      }
      BV::init_value(L_col_counts, L_col_counts_flg, 0);

    }//end init_L_col_counts

    BASKER_INLINE
    void init_U_row_counts(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_U_row_counts: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(U_row_counts_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct urowcount size");
        MALLOC_INT_1DARRAY(U_row_counts, size);
        U_row_counts_flg =size;
      }
      else if(size > U_row_counts_flg)
      {
        BASKER_ASSERT(size > 0, "strcut urowcount size2");
        REALLOC_1DARRAY(U_row_counts, U_row_counts_flg, size);
        U_row_counts_flg = size;
      }
      BV::init_value(U_row_counts, U_row_counts_flg, 0);

    }//end init_U_row_counts

    BASKER_INLINE
    void init_U_col_counts(Int size, BASKER_BOOL verbose)
    {
      if(size <=0)
      {
        if (verbose) {
          std::cout << " init_U_col_counts: size <= 0 no alloc " << std::endl;
        }
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(U_col_counts_flg == 0)
      {
        BASKER_ASSERT(size > 0, "strcut ucolcount size");
        MALLOC_INT_1DARRAY(U_col_counts, size);
        U_col_counts_flg =size;
      }
      else if(size > U_col_counts_flg)
      {
        BASKER_ASSERT(size > 0, "struct ucolcont size2");
        REALLOC_1DARRAY(U_col_counts, U_col_counts_flg, size);
        U_col_counts_flg = size;
      }
      BV::init_value(U_col_counts, U_col_counts_flg, 0);

    }//end init_U_col_counts

    BASKER_INLINE
    void init_S_row_counts(Int size, BASKER_BOOL verbose)
    {
      if(size <=0)
      {
        if (verbose) {
          std::cout << " init_S_row_counts: size <= 0 no alloc " << std::endl;
        }
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(S_row_counts_flg == 0)
      {
        BASKER_ASSERT(size > 0, "struct srowcount size");
        MALLOC_INT_1DARRAY(S_row_counts, size);
        S_row_counts_flg =size;
      }
      else if(size > S_row_counts_flg)
      {
        BASKER_ASSERT(size > 0, "struct srowcount size2");
        REALLOC_1DARRAY(S_row_counts, S_row_counts_flg, size);
        S_row_counts_flg = size;
      }
      BV::init_value(S_row_counts, S_row_counts_flg, 0);

    }//end init_S_row_count

    BASKER_INLINE
    void init_S_col_counts(Int size)
    {
      if(size <=0)
      {
        std::cout << " init_S_col_counts: size <= 0 no alloc " << std::endl;
        return;  //NDE - info of failure to allocate should be returned...
      }
      if(S_col_counts_flg == 0)
      {
        BASKER_ASSERT(size > 0, "strcut scolcount size");
        MALLOC_INT_1DARRAY(S_col_counts, size);
        S_col_counts_flg =size;
      }
      else if(size > S_col_counts_flg)
      {
        BASKER_ASSERT(size > 0, "strcut scolcount size2");
        REALLOC_1DARRAY(S_col_counts, S_col_counts_flg, size);
        S_col_counts_flg = size;
      }
      BV::init_value(S_col_counts, S_col_counts_flg, 0);

    }//end init_S_col_counts


    INT_1DARRAY parent;       BASKER_INT parent_flg;
    INT_1DARRAY post;         BASKER_INT post_flg;
    INT_1DARRAY col_ptr;      BASKER_INT col_ptr_flg;
    INT_1DARRAY left_most;    BASKER_INT left_most_flg;
    INT_1DARRAY row_counts;   BASKER_INT row_counts_flg;
    INT_1DARRAY col_counts;   BASKER_INT col_counts_flg;
    INT_1DARRAY WS;           BASKER_INT WS_flg;
    INT_1DARRAY S_col;        BASKER_INT S_col_flg;
    INT_1DARRAY S_row;        BASKER_INT S_row_flg;
    INT_1DARRAY L_row_counts; BASKER_INT L_row_counts_flg;
    INT_1DARRAY L_col_counts; BASKER_INT L_col_counts_flg;
    INT_1DARRAY U_row_counts; BASKER_INT U_row_counts_flg;
    INT_1DARRAY U_col_counts; BASKER_INT U_col_counts_flg;
    INT_1DARRAY S_col_counts; BASKER_INT S_col_counts_flg;
    INT_1DARRAY S_row_counts; BASKER_INT S_row_counts_flg;
    Int lnnz;
    Int unnz;
  };//end basker_symbolic_tree


  //---------------BASKER--Options------------------//
  template <class Int, class Entry, class Exe_Space>
  struct basker_options
  {
    using STS = Basker_ScalarTraits<Entry>;
    using BASKER_MAGNITUDE = typename STS::magnitudeType;

    basker_options()
    {
      //Operation Options
      same_pattern = BASKER_FALSE;

      //Note verbose will not give debug
      //Debug is a compile time option
      verbose    = BASKER_FALSE;
      verbose_matrix_out = BASKER_FALSE;
     
      //Memory Options
      realloc    = BASKER_FALSE;

      //Symmetric Options
      symmetric  = BASKER_TRUE;

      //Nonsymmetric Options
      AtA        = BASKER_TRUE;
      A_plus_At  = BASKER_FALSE;  //Experimental status
      
      //Transpose Option
      transpose  = BASKER_FALSE;

      //Matching Ordering Options
      //Default is on using bottle-neck
      matching      = BASKER_TRUE;
      matching_type = BASKER_MATCHING_BN;
      // matching before BTF in symbolic 
      // > 0 = none, 1 = ShyLUBasker::mwm, 2 = trilinos_btf_maxtrans, or 3 = MC64 if enabled
      btf_matching = 2;
      min_block_size = 0; // no merging blocks

      // Matrix scaling
      matrix_scaling = 0;

      //BTF Ordering Options
      btf             = BASKER_TRUE;
      btf_max_percent = BASKER_BTF_MAX_PERCENT;
      btf_large       = BASKER_BTF_LARGE;
      use_sequential_diag_facto = BASKER_FALSE;
      // MWM matching before AMD 
      //  0: no matching, 1: ShyLUBasker::mwm, 2: MC63 if enabled
      blk_matching = 1; // if 0, then ND & AMD are applied in symbolic

      // ND Ordering Options (Should METIS optional?)
      use_metis = true;
      run_nd_on_leaves = false;
      run_amd_on_leaves = false;
      use_nodeNDP = true;

      // AMD Option
      amd_dom = true;

      //Pivot
      no_pivot   = BASKER_FALSE;
      pivot_tol  = BASKER_PIVOT_TOL; 
      static_delayed_pivot = 0;

      //Prune (if not pruned, check for numerical cancelatin)
      prune = BASKER_FALSE;
      replace_zero_pivot = BASKER_TRUE;
      replace_tiny_pivot = BASKER_FALSE;

      //BTF Options
      btf_prune_size = (Int)BASKER_BTF_PRUNE_SIZE;

      //Incomplete Factorization Options
      //incomplete = (Int) BASKER_INCOMPLETE_LVL;
      incomplete = BASKER_FALSE;
      incomplete_type = BASKER_INCOMPLETE_LVL;
      //incomplete_type = BASKER_INCOMPLETE_RLVL_LIMITED;
      inc_lvl    = BASKER_INC_LVL_VALUE;
      inc_tol    = BASKER_INC_TOL_VALUE;
      user_fill  = BASKER_FILL_USER;
    }

    //Reuse Pattern (Save time if same pattern can be used)
    BASKER_BOOL same_pattern;

    //Operation Options
    BASKER_BOOL verbose; 
    BASKER_BOOL verbose_matrix_out;

    //Memory Options
    BASKER_BOOL  realloc;
    BASKER_BOOL  transpose;
    
    //Symmetric Options
    BASKER_BOOL  symmetric;
    
    //NonSymmetri Option
    BASKER_BOOL AtA;
    BASKER_BOOL A_plus_At;
    
    //Matching Ordering
    BASKER_BOOL matching;
    Int         matching_type;

    //BTF Ordering Options
    BASKER_BOOL  btf;
    BASKER_MAGNITUDE btf_max_percent;
    BASKER_MAGNITUDE btf_large;
    BASKER_BOOL  use_sequential_diag_facto;
    // TODO: remove this (matching during symbolic)
    int btf_matching; // carbinality matching before BTF (Symbolic):                0 = none, 1 = Basker, or 2 = Trilinos (default) (3 = MC64 if enable)
    int blk_matching; // max weight matching on each of diagonal blocks (Numeric):  0 = none (default), or 1 = Basker (2 = MC64 if enable)
    Int min_block_size; // min size of blocks
    // matrix scaling
    int matrix_scaling; // 0 = no scaling, 1 = symmetric diagonal scaling, 2 = non-symm row-max, and then col-max scaling

    // Replace tiny pivott
    BASKER_BOOL replace_zero_pivot;
    BASKER_BOOL replace_tiny_pivot;

    // Option to apply MWM at numeric through "delayed" pivot
    int static_delayed_pivot;

    //AMD Ordering Options
    BASKER_BOOL  amd_dom;

    // ND Ordering Options
    BASKER_BOOL use_metis;
    BASKER_BOOL run_nd_on_leaves;
    BASKER_BOOL run_amd_on_leaves;
    BASKER_BOOL use_nodeNDP;

    //Pivot Options
    BASKER_BOOL  no_pivot;
    BASKER_MAGNITUDE pivot_tol;

    //Prune (if not pruned, check for numerical cancelatin)
    BASKER_BOOL prune;

    //BTF Options
    BASKER_INT   btf_prune_size;
  
    //Incomplete Factorization Options
    BASKER_BOOL  incomplete;
    BASKER_INT   incomplete_type;
    BASKER_INT   inc_lvl;
    BASKER_MAGNITUDE inc_tol;    //Not Used
    BASKER_MAGNITUDE user_fill;
    
    /* ---- todo add more ----*/
  }; // end bask_options

}//end namespace Basker

#endif //end ifdef basker_structs.hpp
