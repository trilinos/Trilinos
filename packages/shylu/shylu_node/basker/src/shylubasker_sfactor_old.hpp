#ifndef SHYLUBASKER_SFACTOR_OLD_HPP
#define SHYLUBASKER_SFACTOR_OLD_HPP

#include "shylubasker_decl.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_def.hpp"
#include "shylubasker_structs.hpp"

#include <iostream>
using namespace std;

namespace BaskerNS
{

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
      Int kid = (Int)(thread.league_rank()*thread.team_size()+
		      thread.team_rank());
      #endif

      basker->t_init_factor(kid);

    }//end operator()
  }; //end struct kokkos_sfactor_init_factor


  //Allocate space for LU when sybolic factorization has not been called
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::default_symb()
   {
     //check if already called
       if(symb_flag == true) 
       {return -1;}

     #ifdef BASKER_DEBUG_SFACTOR
     printf("Default Symbolic Called \n");
     #endif     

     return 0;
   }//end default_symb()

  //Note: Removed allocate_workspace_old() //

  //allocate workspace for each thread
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::allocate_workspace()
  {
    for(Int t = 0; t < num_threads; t++)
    {
      #ifdef BASKER_2DL
      Int b = S[0][t];
      MATRIX_1DARRAY &Ltemp = LL[b];

      for(Int j = 0; j < LL_size[b]; j++)
      {
        Ltemp[j].ews_size = Ltemp[j].nrow;
        Ltemp[j].iws_size = Ltemp[j].nrow;
        Ltemp[j].ews_mult = 2; //Note: comeback and make smaller
        Ltemp[j].iws_mult = 5; //Note: comeback and make smaller
        MALLOC_ENTRY_1DARRAY(LL[b][j].ews, 
            Ltemp[j].ews_size*Ltemp[j].ews_mult);
        MALLOC_INT_1DARRAY(LL[b][j].iws,
            Ltemp[j].iws_size*Ltemp[j].iws_mult);
      }

      #else  //else if BASKER_2DL
      //Alloc workspace arrays
      #ifdef BASKER_ATOMIC 
      printf("allocate_workspace, size: %d %d %d %d \n",
          5, A.nrow, 2, A.nrow);
      MALLOC_INT_1DARRAY(thread_array[t].iws, 5*A.nrow);
      //note could be made maller
      MALLOC_ENTRY_1DARRAY(thread_array[t].ews, 2*A.nrow);
      //note could be made smaller
      thread_array[t].iws_size = A.nrow;
      thread_array[t].iws_mult = 5;
      thread_array[t].ews_size = A.nrow;
      thread_array[t].ews_mult = 2;

      #else
      MALLOC_INT_1DARRAY(thread_array[t].iws,
          thread_array[t].iws_size*thread_array[t].iws_mult);
      MALLOC_ENTRY_1DARRAY(thread_array[t].ews, 
          thread_array[t].ews_size*thread_array[t].ews_mult);
      #endif //Endif BASKER_ATOMIC

      #endif //Endif BASKER_2DL
      //Alloc C matrix
      //Note: could be made smaller
      thread_array[t].C.init_matrix("column matrix",
          thread_array[t].ews_size, 2,
          thread_array[t].ews_size);
      thread_array[t].C.fill();
      thread_array[t].C.malloc_union_bit(); //Note: may be removed
      thread_array[t].C.init_union_bit(); //Note: may be removed
    }//for every thread
    workspace_flag = true;
    return 0;
  }//end allocate_workspace()


  //Note: comeback and fix
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::init_stree()
  {
    MALLOC_INT_1DARRAY(stree.parent, A.ncol);
    init_value(stree.parent, A.ncol, 0);

    MALLOC_INT_1DARRAY(stree.post, A.ncol);
    init_value(stree.post, A.ncol, 0);

    MALLOC_INT_1DARRAY(stree.col_ptr, A.ncol);
    init_value(stree.col_ptr, A.ncol, 0);

    MALLOC_INT_1DARRAY(stree.left_most, A.ncol);
    init_value(stree.left_most, A.ncol, 0);

    MALLOC_INT_1DARRAY(stree.row_counts, A.nrow);
    init_value(stree.row_counts, A.nrow, 0);

    MALLOC_INT_1DARRAY(stree.col_counts, A.ncol);
    init_value(stree.col_counts, A.ncol, 0);

    MALLOC_INT_1DARRAY(stree.WS, 10*A.ncol);
    init_value(stree.WS, 10*A.ncol, 0);
    
    stree.lnnz = 0;
    stree.unnz = 0;

    return 0;
  }//end init_stree()


  //Note: lots of work to fix!
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::nontree_QR_symb()
  {
     init_stree();
  
     //Get nontree QR counts 
     //not this can be made much better using the tree structure
     //plus can help remove structure singluar blocks
 
     // e_tree of A'A , stored in etree
     e_tree(1);
     
     #ifdef BASKER_DEBUG_SFACTOR
     printf("etree\n");
     #endif

     //get post order of etree
     post_order(stree.parent, A.ncol);

     row_count(stree.parent, stree.post);
   
     #ifdef BASKER_DEBUG_SFACTOR
     cout << "Row Count: ";
     for(Int i = 0; i < 10; i++)
       {
         cout << stree.row_counts[i] << " , " ;
       }
     cout << endl;
     #endif

     
     //We can only assure this works with hund ordering
     //because of pivots 

     householder_count(stree.parent);

     for(Int i = 0; i < A.ncol; i++)
     {
       stree.unnz += stree.row_counts[i];
       stree.lnnz += stree.col_counts[i];
     }
     
     MALLOC_MATRIX_2DARRAY(LL, tree.nblks);
     MALLOC_MATRIX_2DARRAY(LU, tree.nblks);
     

     INT_1DARRAY LU_count;
     MALLOC_INT_1DARRAY(LU_count, tree.nblks);
     init_value(LU_count, tree.nblks, 0);

     
     #ifdef BASKER_2DL
     INT_1DARRAY LL_count;
     MALLOC_INT_1DARRAY(LL_count, tree.nblks);
     init_value(LL_count, tree.nblks, 0);
     #endif

     for(Int i = 0; i < tree.nblks; i++)
     {
       //Malloc Inner L
    	 #ifdef BASKER_2DL
       MALLOC_MATRIX_1DARRAY(LL[i], LL_size[i]);
       printf("L size %d \n", LL_size[i]);
    	 #else
       MALLOC_MATRIX_1DARRAY(LL[i], 1);
    	 #endif

       //Malloc Inner U
       MALLOC_MATRIX_1DARRAY(LU[i], LU_size[i]);

       //Debug
       #ifdef BASKER_DEBUG_SFACTOR
       printf("AV.size(%d) = %d \n", i, LU_size[i]);
       printf("AL.size(%d) = %d \n", i, LL_size[i]);
       #endif
     }
     

     //Init submatrix variables (come back and fix up)
     #ifdef BASKER_2DL
     for(Int i=0; i < tree.nblks; i++)
     {

       //----------------Get count for L------------------//
       MATRIX_1DARRAY &Ltemp = LL[i];
       MATRIX_VIEW_1DARRAY &ALtemp = AL[i];
       Int local_nnz = 0;

       //Get local count
       for(Int k = tree.col_tabs[i]; k < tree.col_tabs[i+1]; k++)
       {
         local_nnz += stree.col_counts[k];
       }

       //for each L in the column
       #ifdef BASKER_DEBUG
       printf("Access blk: %d LL_size: %d LU_size: %d \n",
           i, LL_size[i], LU_size[i]);
       #endif


       for(Int j = 0; j < LL_size[i]; j++)
       {

         printf("%d %d %d %d %d \n",
             ALtemp[j].srow,
             ALtemp[j].nrow,
             ALtemp[j].scol,
             ALtemp[j].ncol,
             local_nnz+1);

         Ltemp[j].init_matrix("L:", 
             ALtemp[j].srow,
             ALtemp[j].nrow,
             ALtemp[j].scol,
             ALtemp[j].ncol,
             local_nnz+1);	     

         #ifdef BASKER_DEBUG
         printf("Init L(%d,%d), sr: %d m: %d sc: %d n: %d nnz: %d \n",
             i, j,
             Ltemp[j].srow, Ltemp[j].nrow,
             Ltemp[j].scol, Ltemp[j].ncol,
             Ltemp[j].nnz);
         #endif
       }//end over all L submatrices in column


       //-----------Get Count for U--------------------//
       //Iterate over all columns in a row
       for(Int j=i; j!=A.ncol; j=tree.treetab[j])
       {
         MATRIX_1DARRAY &Utemp = LU[j];
         MATRIX_VIEW_1DARRAY &AVtemp = AV[j];
         local_nnz = 0;

         //go back and change this to the intersection of col and row counts
         //need a smarter way to do this
         for(Int k = tree.col_tabs[j]; k < tree.col_tabs[j+1]; k++)
         {
           local_nnz += stree.row_counts[k];
         }
         Utemp[LU_count[j]].init_matrix("U:", 
             AVtemp[LU_count[j]].srow,
             AVtemp[LU_count[j]].nrow,
             AVtemp[LU_count[j]].scol,
             AVtemp[LU_count[j]].ncol,
             local_nnz+1);

         #ifdef BASKER_DEBUG
         printf("Init U(%d,%d), sr: %d m: %d sc: %d n: %d nnz: %d \n",
             j, LU_count[j],
             Utemp[LU_count[j]].srow, Utemp[LU_count[j]].nrow,
             Utemp[LU_count[j]].scol, Utemp[LU_count[j]].ncol,
             Utemp[LU_count[j]].nnz);
         #endif

         LU_count[j] = LU_count[j] + 1;

       }//end over all row

     }//end over all nblks

     #else //Elseif BAKSER_2DL

     for(Int i=0; i < tree.nblks; i++)
     {
       MATRIX_1DARRAY temp = LL[i];
       Int local_nnz = 0;

       //count nnz
       for(Int k = tree.col_tabs[i]; k < tree.col_tabs[i+1]; k++)
       {
         local_nnz += stree.col_counts[k];
       }

       temp[0].init_matrix("L:" , tree.col_tabs[i], 
           (A.nrow - tree.col_tabs[i]),
           tree.col_tabs[i],
           (tree.col_tabs[i+1] - tree.col_tabs[i]),
           local_nnz+1);

       //Note, we can get rid of this
       temp[0].malloc_perm(A.nrow);

       #ifdef BASKER_DEBUG_SFACTOR
       cout << "Alloced L(" << i << ") with m: "<< temp[0].nrow
         << " n: " << temp[0].ncol 
         << " nnz: " << temp[0].nnz << endl;
       #endif

       for(Int j=i; j!=A.ncol; j=tree.treetab[j])
       {
         //Init LU
         temp = LU[j];
         local_nnz = 0;

         //go back and change this to the intersection of col and row counts
         //need a smarter way to do this
         for(Int k = tree.col_tabs[j]; k < tree.col_tabs[j+1]; k++)
         {
           local_nnz += stree.row_counts[k];
         }

         temp[LU_count[j]].init_matrix("U: " ,//tree.row_tabs[i],
             //(tree.row_tabs[i+1] - tree.row_tabs[i]),
             tree.col_tabs[i],
             (tree.col_tabs[i+1] - tree.col_tabs[i]),
             tree.col_tabs[j],
             (tree.col_tabs[j+1] - tree.col_tabs[j]), 
             local_nnz+1);

       #ifdef BASKER_DEBUG_SFACTOR              
         cout << "Allocating U " << j << " , " << i << " , " << LU_count[j] 
           << "m: " << (tree.col_tabs[i+1] - tree.col_tabs[i])
           << "n: " <<  (tree.col_tabs[j+1] - tree.col_tabs[j])
           << "nnz: " << local_nnz << endl;
       #endif

         LU_count[j] = LU_count[j] + 1;
       }//over all leaves

     }//over all blks
    #endif //Endif BASKER_2DL

     FREE(LU_count);
     FREE(LL_count);


     //test what is going on
     if(workspace_flag == false)
     {allocate_workspace();}

     //init factor
     #ifdef BASKER_KOKKOS
     typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
     #ifdef BASKER_NO_LAMBDA
     kokkos_sfactor_init_factor<Int,Entry,Exe_Space> iF(this);
     Kokkos::parallel_for(TeamPolicy(num_threads,1), iF);
     Kokkos::fence();
     #else //else if BASKER__NO_LAMBDA

     #endif //endif BASKER_NO_LAMBDA
     #else //else if BASKER_KOKKOS
     
     #endif //endif BASKER_KOKKOS
     
     MALLOC_INT_1DARRAY(gperm, A.nrow);
     init_value(gperm, A.nrow, A.max_idx);
     MALLOC_INT_1DARRAY(gpermi, A.nrow);
     init_value(gpermi, A.nrow, A.max_idx);


     symb_flag = true;
     
     //Note: tree STREEE ???

   
    return 0;
  }//end nontree_QR_symb()


  //Creates an e-tree of matrix
  //need to change workspace
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::e_tree(Int ata)
  {
    Int *ws;
    Int ws_size = A.ncol;

    if(ata == 1)
    {ws_size += A.nrow;}

    ws = new Int[ws_size]();
    Int *past = ws;
    Int *clique = ws+A.ncol;

    //Zero out the cliques from ATA
    if(ata == 1)
    {
      for(Int ii = 0 ; ii < A.nrow; ii++)
      {clique[ii] = A.ncol;}
    }
    // for each column
    for(Int k = 0; k < A.ncol; k++)
    {
      stree.parent[k] = A.ncol;
      past[k] = A.ncol;
      for(Int j = A.col_ptr[k]; j < A.col_ptr[k+1]; j++)
      {
        Int t = A.row_idx[j];
        if(ata==1)
        {t = clique[A.row_idx[j]];} //use clique instead

        Int have_past;
        for(; (t!=A.ncol) && (t<k); t = have_past)
        {
          have_past = past[t];
          past[t] = k;
          if(have_past == A.ncol)
          {stree.parent[t] = k;}
        }//climb up the tee until root and connect
        //connect others in clique
        if(ata == 1)
        {clique[A.row_idx[j]] = k;}
      } 
    }//end over all columns
    delete [] ws;
    return 0;
  }//end e_tree()


  //Post orders a parent array representation of a tree
  //need to do something about workspace
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::post_order(INT_1DARRAY &p, Int n)
  {
    //stree.post.malloc("stree:post", n); stree.post.init();
    MALLOC_INT_1DARRAY(stree.post, n);
    init_value(stree.post, n, 0);

    Int *post = new Int[n]();
    Int ws_size = 3*n;
    Int *ws = new Int[ws_size]();
    Int *head, *next, *stack;
    head = ws; next = ws+n; stack = ws+2*n;

    for(Int j=0; j < n; j++)
    {head[j] = A.ncol;}

    for(Int k = n; k > 0 ; k--)
    {
      Int j = k-1;
      if(p[j] == A.ncol)
      { continue; }
      next[j] = head[p[j]];
      head[p[j]] = j;
    }

    Int k = 0;
    for(Int j=0; j<A.ncol; j++)
    {
      if(p[j]!=A.ncol) 
      { continue; }
      k = post_dfs(j, k, head, next, post, stack);
    }
    delete [] ws;

    //Change quick fix
    for(Int i = 0; i < n; i++)
    {
      stree.post[i] = post[i];
    }
    delete [] post;

    //return post;
  }//end post_order()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  Int Basker<Int, Entry, Exe_Space>::post_dfs
  (
   Int j,
   Int k,
   Int *head,
   Int *next,
   Int *post,
   Int *stack
  )
  {
    Int ii, jj, top;
    top = 0;
    //Start from j
    stack[0] = j;
    while(top != A.ncol)
    {
      jj = stack[top];
      ii = head[jj];
      if(ii == A.ncol)
      {
        if(top == 0)
        {top = A.ncol;}
        else
        {top--;}
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


  //Like to find a way to do this without needed the transpose
  //don't need the transpose if the matrix is symmtric so use cheaper function call!
  template <class Int, class Entry, class Exe_Space>
  Int* Basker<Int, Entry, Exe_Space>::col_counts_ata(Int *parent, Int *post)
  {
    transpose(); //create a transpose so can find count of ATA \approx R \approx U
    Int ws_size = 4*A.ncol+ (A.ncol+A.nrow+1);
    Int *ws =   new Int[ws_size]();
    Int *delta = new Int[A.ncol]();
    
    Int jleaf;
    Int *past, *mfirst, *pleaf, *first;
    past = ws;  mfirst = ws+A.ncol; pleaf = ws+2*A.ncol; first = ws+ 3*A.ncol;
    
    for(Int k = 0; k < ws_size; k++) 
    {ws[k] = A.ncol;}
    
    for(Int k = 0; k < A.ncol; k++)
    {
      Int j = post[k];
      if(first[j] == A.ncol)
      {
        delta[j] = 1; // leaf node
      }
      else
      {
        delta[j] = 0;
      }

      for( ; (j != A.ncol) && (first[j] == A.ncol); j = parent[j])
      {
        first[j] = k; // update with parent
      }

    }//initalize the delta counts for overlap
  
    //Create a linked list of the cliques
    Int *head = ws+4*A.ncol;
    Int *next = ws+5*A.ncol+1;

    for(Int k=0; k < A.ncol; k++)
    {ws[post[k]] = k;} //overwriting past

    for(Int i = 0; i < A.nrow; i++)
    {
      Int k=A.ncol;
      for(Int p = At.col_ptr[i]; p < At.col_ptr[i+1]; p++)
      { 
        k = min(k, ws[At.row_idx[p]]);
      }
      next[i] = head[k];
      head[k] = i;
    }
    //End create a linked list of the cliques

    // reset past
    for(Int k = 0; k < A.ncol; k++)
    {past[k] = k;}
    
    for(Int k = 0; k < A.ncol; k++)
    {
      Int j = post[k];
      if(parent[j] != A.ncol)
      {
        delta[parent[j]]--;
      }  //j is not a root node

      //loop over clique; In A this would only be done once
      for(Int J = head[k]; J != A.ncol; J = next[J])
      {
        for(Int p = At.col_ptr[J]; p < At.col_ptr[J+1]; p++)
        {
          Int i = At.row_idx[p];
          Int q = least_common(i, j, first, mfirst, pleaf, past, &jleaf);

          if(jleaf >=1) 
          {delta[j]++;}

          if(jleaf == 2)
          {delta[q]--;}

        }//for over row/col
      }//for over elements in clique

      if(parent[j] != A.ncol)
      {past[j] = parent[j];}

    }//over all col/row

    for(Int k = 0; k < A.ncol; k++)
    {
      if(parent[k] != A.ncol)
      {
        delta[parent[k]] += delta[k];
      }
    }
    
    // Clean up AT 
       
    // Clean up workspace
    delete [] ws;

    return delta ;
  }//end col_count_ata()


  template <class Int, class Entry, class Exe_Space>
  Int Basker<Int, Entry, Exe_Space>::least_common
  (
   Int i, 
   Int j, 
   Int *first
   Int *mfirst,
   Int *pleaf, 
   Int *past,
   Int *jleafi
  )
  {
    Int q, sparent;
    *jleaf = 0;
    mfirst[i] = first[j];
    Int jprev = pleaf[i];
    pleaf[i] = j;

    if(jprev == A.ncol)
    {*jleaf = 1;}
    else
    {*jleaf = 2;}

    if(*jleaf == 1) 
    {return i;}

    for(q = jprev; q != past[q]; q = past[q]);

    for(Int k = jprev; k != q; k = sparent)
    {
      sparent = past[k];
      past[k] = q;
    }

    return q;      
  }//end least_common() 


  //malloc of matrix need to change
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::transpose()
  {
    //Malloc the space for the transpose
    At.init_matrix("tranpose" , A.ncol, A.nrow, A.nnz);
   
    //All done in serial so can fill here
    //init here
    for(int i = 0; i < A.ncol+1; i++)
    {
      At.col_ptr[i] = 0;
    }
    for(int i = 0; i < A.nnz; i++)
    {
      At.row_idx[i] = 0;
      At.val[i] = 0;
    }

    Int ws_size = A.nrow;
    Int *ws = new Int[ws_size]();
    //get row counts
    for(Int j = 0; j < A.col_ptr[A.ncol]; j++)
    {ws[A.row_idx[j]]++;}
    
    At.col_ptr[1] = ws[0];
    for(Int j = 1; j < A.nrow; j++)
    {
      ws[j] = ws[j]+ws[j-1];
      At.col_ptr[j+1] = ws[j];
      ws[j-1] = At.col_ptr[j-1];
    }

    ws[A.nrow-1] = At.col_ptr[A.nrow-1];

    for(Int k = 0; k < A.ncol; k++)
    {
      for(Int j = A.col_ptr[k]; j < A.col_ptr[k+1]; j++)
      {
        At.row_idx[ws[A.row_idx[j]]++] = k;
      }
    } 
    
    delete [] ws;

    return 0;
  }//end transpose()


  //need to fix work spacce
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int, Entry, Exe_Space>::householder_count(INT_1DARRAY &parent)
  {
    //we are not going to worry about structure singluar (yet!)
    Int *leftmost = new Int[A.nrow]();
    Int ws_size = A.nrow + 3 *A.ncol;
    Int *ws = new Int[ws_size]();
    Int *nq  = new Int[A.nrow]();
    Int *lpinv = new Int[A.ncol + A.nrow];
 
    //Divide the workspace
    Int *next,*head, *tail;
    next = ws; head = ws + A.nrow; tail = ws + A.nrow + A.ncol; 
    //nq = ws + A->nrow + 2*A->ncol; //This seems wrong
    for(Int k=0; k < A.ncol; k++) head[k] = A.ncol;

    for(Int k=0; k < A.ncol; k++) tail[k] = A.ncol;

    for(Int k=0; k < A.ncol; k++) nq[k] = 0;

    for(Int i=0; i < A.nrow; i++) leftmost[i] = A.ncol;

    //Find leftmost
    for(Int c = A.ncol; c > 0 ; c--)
    {
      Int k = c-1;
      for(Int p = A.col_ptr[k]; p < A.col_ptr[k+1]; p++)
      {
        leftmost[A.row_idx[p]] = k;
      }
    }

    //Build list of connections
    for(Int ii=A.nrow; ii > 0; ii--)
    {
      Int i = ii -1;
      Int k = leftmost[i];

      if(k == A.ncol)
      {continue;}

      if(nq[k]++ == 0)
      {
        tail[k] = i;
      }
      next[i] = head[k];
      head[k] = i;
    }

    //Count
    Int m2 = A.nrow;
    Int count = 0;
    for(Int k = 0; k < A.ncol; k++)
    {
      Int i = head[k];
      count++;

      if(i < 0)
      {
        i = m2++;
      }

      lpinv[i] = k;

      if(--nq[k] <=0)
      { continue; }

      count += nq[k];
      Int pa = parent[k];
      if(pa != A.ncol)
      {
        if(nq[pa] == 0)
        {
          tail[pa] = tail[k];
        }
        next[tail[k]] = head[pa];
        head[pa] = next[i];
        nq[pa] += nq[k];
      }
    }

    for(Int k = 0; k < A.ncol; k++)
    { nq[k]++;}
    
    delete [] leftmost;
    delete [] lpinv;
    delete [] ws;

    for(Int i = 0; i < A.ncol; i++)
    {
      stree.col_counts[i] = nq[i];
    }
    
    delete [] nq;
  }//end householder_count()


  //need to do something about ws malloc
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry, Exe_Space>::row_count(INT_1DARRAY &parent, INT_1DARRAY &post)
  {
    //matrix transpose
    transpose();
    Int *head = new Int [A.nrow+1]();
    Int *next = new Int [A.nrow+1]();
    Int *first = new Int[A.nrow+1](); // first[i] = k is the least postordering of a
                                      //  any descendent i
    Int *lvl = new Int[A.nrow](); //length of path from i to root
    
    Int *ws = new Int[2*A.ncol]();

    for(Int i = 0; i < A.nrow+1; i++)
    {
      head[i] = A.ncol;
      first[i] = A.ncol;
      next[i] = A.ncol;
    }

    //post order traveral to find first and lvl
    for(Int k = 0; k < A.nrow; k++)
    {
      Int j = post[k];
      Int r;
      Int len = 0; 

      for(r =j ; (r!=A.ncol) && (first[r] == A.ncol); r = parent[r])
      {
        first[r] = k;
        len++;
      }

      if(r == A.ncol)
      {
        len--;
      }
      else
      {
        len += lvl[r];
      }

      for(Int s = j; s != r; s = parent[s])
      {
        lvl[s] = len--;
      }
    }
    
    //Build link list representing the cliques of A'A
    for(Int k = 0 ; k < A.ncol; k++)
    { 
      ws[post[k]] = k;
    }

    for(Int i =0 ; i < A.nrow; i++)
    {
      Int k, p;
      k = A.ncol;

      for( p = At.col_ptr[i]; p < At.col_ptr[i+1]; p++)
      {
        k = min(k, ws[At.row_idx[p]]);
      }

      next[i] = head[k];
      head[k] = i;
    }

    Int *rowcount = new Int[A.nrow]();
    Int *pleaf = new Int[A.ncol]();
    Int *mfirst = new Int[A.ncol]();
    Int *past = new Int[A.ncol]();
    Int jleaf;

    for(Int i = 0 ; i < A.nrow; i++)
    {
      rowcount[i] = 1;
    }

    for(Int i = 0; i < A.ncol; i++)
    {
      pleaf[i] = A.ncol;
      mfirst[i] = A.ncol;
      past[i] = i;
    }

    for(Int l = 0; l < A.ncol; l++)
    {
      Int j = post[l];

      //loop over the linked list
      for(Int J = head[l]; J!= A.ncol; J = next[J])
      {
        for(Int p = At.col_ptr[J]; p < At.col_ptr[J+1]; p++)
        {
          Int i = At.row_idx[p];
          Int q = least_common(i, j, first, mfirst, pleaf, past, &jleaf);

          if(jleaf)
          {
            rowcount[i] += (lvl[j] - lvl[q]);
          }
        }
      }

      if(parent[j] != A.ncol)
      {
        past[j] = parent[j];
      }
    }

    delete [] head;
    delete [] next;
    delete [] first;
    delete [] lvl;  
    delete [] ws ;
    delete [] pleaf;
    delete [] mfirst;
    delete [] past;
    
    for(Int i = 0; i < A.nrow; i++)
    {
      stree.row_counts[i] = rowcount[i];
    }

    delete [] rowcount;
  }//end row_count()

}//end namespace baskerNS

#endif //end basker_sfactor_hpp
