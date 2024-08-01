// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_TREE_HPP
#define SHYLUBASKER_TREE_HPP

//Basker includes
//#include "shylubasker_decl.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_structs.hpp"

//C++ includes
#include <cstdlib>
#include <iostream>
#include <math.h>

using namespace std;

namespace BaskerNS
{
  ///wrapper init tree function for selection
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::init_tree(
  Int *perm, Int nblks, Int parts, 
  Int *row_tabs, Int *col_tabs, Int *treetab, Int options)       
  {
    if((tree_flag == true) || (matrix_flag == false))
    {return -1;}

    //Copy into our tree_struct
    part_tree.basic_convert(A.nrow, perm, nblks, parts,
                            row_tabs, col_tabs, treetab);

    int result = 0;
    if(options == 0) //default case (based on threads)
    {
      /*
         result = init_tree_thread(perm, nblks, parts, row_tabs, 
         col_tabs, treetab);
      */
      init_tree_thread();

    }
    else if(options > 0)
    {
      /*
         result = init_tree_lvl(perm, nblks, parts, row_tabs, col_tabs, 
         treetab, options);
      */
    }
    
    return result;
  }//init_tree


  //Used to provide memory
  //Do we still need this??
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::malloc_tree(Int parts, Int nblks, Int nlvls)
  {
    tree.nroots = 1;
    MALLOC_INT_1DARRAY(tree.roots, 1);
    tree.nparts = parts;  
    tree.nblks = nblks;
    BASKER_ASSERT((nblks+1)>0, "tree nblks malloc");
    BASKER_ASSERT(A.ncol > 0, "tree ncol malloc");
    BASKER_ASSERT(A.nrow > 0, "tree nrow malloc");
    MALLOC_INT_1DARRAY(tree.row_tabs, nblks+1);
    MALLOC_INT_1DARRAY(tree.col_tabs, nblks+1);
    MALLOC_INT_1DARRAY(tree.treetab, nblks+1);
    MALLOC_INT_1DARRAY(tree.treeptr, A.ncol);
    MALLOC_BOOL_1DARRAY(tree.pivot, A.nrow);
    MALLOC_INT_1DARRAY(tree.rowptr, nblks+1);
    tree.nlvls = nlvls;
    MALLOC_INT_1DARRAY(tree.lvlset, nblks+1);
    MALLOC_INT_1DARRAY(tree.child, nblks+1);
    MALLOC_INT_1DARRAY(tree.sibling, nblks+1);
    MALLOC_INT_2DARRAY(S, nblks+1);
    MALLOC_INT_1DARRAY(LL_size, nblks+1);
    MALLOC_INT_1DARRAY(LU_size, nblks+1);

    init_tree_struc();
  }//end malloc_tree


  //Used to init values
  //Can be done with functor call in future
  //Do we still need this
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::init_tree_struc()
  {   
    //Ok because will only be accessed by master thrd
    for(Int i=0; i < tree.nblks+1; i++)
    {
      tree.row_tabs[i] = 0;
      tree.col_tabs[i] = 0;
      tree.treetab[i]  = 0;
      tree.rowptr[i]   = 0;
      //tree.child[i]    = A.max_idx;
      tree.child[i]     = BASKER_MAX_IDX;
      //tree.sibling[i]  = A.max_idx;
      tree.sibling[i]   = BASKER_MAX_IDX;
    }

    //S will be accessed by all threads (but is first created using one thread
    //Should we use first touch?
    for(Int i =0; i < tree.nblks+1; i++)
    {
      BASKER_ASSERT(num_threads > 0, "tree num_threads");
      MALLOC_INT_1DARRAY(S[i], num_threads);
    }

    //this will want to be across all threads 
    //(but is first created using on thread)
    for(Int i=0; i < A.nrow; i++)
    {
      tree.pivot[i] = true;
    }

  }// end init_tree()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::init_tree_thread()
  {
    if (Options.verbose == BASKER_TRUE)
    {
      printf("----init_tree_thread()----\n");
      printf("---Starting Tree----\n");
      part_tree.info();
    }

    //Note: Convert this to use these values everywhere
    INT_1DARRAY ttemp;
    Int parts = part_tree.nparts;
    Int nblks = part_tree.nblks;
    Int *row_tabs = &(part_tree.row_tabs(0));
    Int *col_tabs = &(part_tree.col_tabs(0));
    Int *treetab  = &(part_tree.treetab(0));
    
    #ifdef BASKER_DEBUG_TREE
    printf("Start, nparts: %d nblks: %d \n",
           parts, nblks);
    #endif


    //Calculate maximum number of levels
    // parts^x = threads
    #ifdef BASKER_DEBUG_TREE
    printf("nblks: %d parts: %d \n", nblks, parts);
    printf("log: %f %f %f \n", 
           log(nblks+1), log(parts), log(num_threads));
    #endif

    Int maxlvls = round(
                        log((double)num_threads)/
                        log((double)parts));

    Int treelvls = round(
                         log((double)(nblks+1))/
                         log((double)parts)) - 1;


    //Int treelvls = ((((double)log(nblks+1))/((double)log(parts))) - 1.0);
    //Int treelvls = (((double)log((double)(nblks+1)))/
    //                    ((double)log((double)parts))) ;
    //Int treelvls = crazy_temp -1;

    //Number of blocks now is parts^lvl+1 - 1
    Int newnblks = round(pow(parts, maxlvls+1)) - 1;
    
    #ifdef BASKER_DEBUG_TREE
    printf("maxlvl: %d treelvl: %d newnblks: %d \n",
           maxlvls, treelvls, newnblks);
    #endif

    if(newnblks < nblks)
    {
      cout << "Reducing nblks to: " << newnblks << 
        " in order to fit threads" <<endl;
    }
    else if(newnblks > nblks)
    {
      cout << "Error: need more blocks to fit threads" << endl;
      return -1;
    }

    malloc_tree(parts, newnblks, maxlvls);
    //thread array
    //MALLOC_THREAD_1DARRAY(thread_array, num_threads);

    //Step 1
    //Build new row and col tabs (DFS)
    Int old_pos = nblks;
    Int new_pos = tree.nblks;
    
    INT_1DARRAY temp_tree;
    BASKER_ASSERT((tree.nblks+1)>0, "tree blks here");
    MALLOC_INT_1DARRAY(temp_tree, tree.nblks+1);
    init_value(temp_tree, tree.nblks+1, (Int)0);
    INT_1DARRAY temp_col_tabs;
    MALLOC_INT_1DARRAY(temp_col_tabs, tree.nblks+1);
    init_value(temp_col_tabs, tree.nblks+1, (Int)0);
    INT_1DARRAY temp_row_tabs;
    MALLOC_INT_1DARRAY(temp_row_tabs, tree.nblks+1);
    init_value(temp_row_tabs, tree.nblks+1, (Int)0);
    
    rec_tabs(0, tree.nlvls, treelvls, tree.nblks, old_pos,
             &new_pos, col_tabs, row_tabs, treetab,
       temp_col_tabs, temp_row_tabs, temp_tree);

    for(Int i=0; i < tree.nblks+1; i++)
    {
      tree.row_tabs[i] = temp_row_tabs[i];
      tree.col_tabs[i] = temp_col_tabs[i];
      tree.treetab[i] = temp_tree[i];
    }

    FREE(temp_tree);
    FREE(temp_col_tabs);
    FREE(temp_row_tabs);

    #ifdef BASKER_DEBUG_TREE
    cout << "row_tabs: " << endl;
    for(Int i = 0; i < tree.nblks+1; i++)
      {cout << tree.row_tabs[i] << " " ;}
    cout << endl;
    cout << "col_tabs: " << endl;
    for(Int i = 0; i < tree.nblks+1; i++)
      {cout << tree.col_tabs[i] << " ";}
    cout << endl;
    cout << "tree_tabs: " << endl;
    for(Int i = 0; i < tree.nblks+1; i++)
      {cout << tree.treetab[i] << " ";}
    cout << endl;
    cout << "SIZES: " << endl;
    for(Int i = 0; i < tree.nblks; i++)
      {cout << i << "  " << tree.col_tabs[i+1]-tree.col_tabs[i] << endl;}
    #endif

    BASKER_ASSERT(num_threads > 0, "tree num_threas");
    MALLOC_THREAD_1DARRAY(thread_array, num_threads);

    //Step 2
    //Build treeptr list
    for(Int i=0; i < tree.nblks; i++)
    {
      for(Int j = tree.col_tabs[i]; j < tree.col_tabs[i+1]; j++)
      {tree.treeptr[j] = i;}
    }
    
    //Step 3
    //Get level sets  (BFS) 
    //Note: needs refactored for nparts
    INT_1DARRAY temptree;
    BASKER_ASSERT(tree.nblks > 0, "tree treenblks");
    MALLOC_INT_1DARRAY(temptree, tree.nblks);
    init_value(temptree, tree.nblks, (Int)0);

    new_pos = tree.nblks-1;
    temptree[new_pos] = tree.treetab[new_pos];
    tree.lvlset[new_pos] = new_pos;
    
    new_pos--;
    for(Int l = 1; l < tree.nlvls+1; l++)
    {
      Int nparents = pow(tree.nparts, l-1);
      Int pp = new_pos;
      //ttemp = lvl_task[l];
      //Int task_pos = 0; //Not used
      //Int task_offset = pow(tree.nparts, tree.nlvls-l); //NU
      for(Int p = 0 ; p < nparents; p++)
      {
        Int parentptr = pp+nparents-p;
        Int leftc     = tree.lvlset[parentptr]-1;
        Int rightc    = tree.lvlset[parentptr] - 
          pow(tree.nparts, tree.nlvls-l+1);

        tree.lvlset[new_pos] = leftc;
        temptree[new_pos] = tree.treetab[leftc];
        new_pos--;

        tree.lvlset[new_pos] = rightc;
        temptree[new_pos] = tree.treetab[rightc];
        tree.child[tree.lvlset[parentptr]] = rightc;
        tree.sibling[rightc] = leftc;
        new_pos--; 
      }//over all parents
    }// over all levels


    #ifdef BASKER_DEBUG_TREE
    cout << "Treetabs: " << endl;
    for(Int b = 0 ; b < tree.nblks; b++)
      { cout << temptree[b] << " " ;}
    cout << endl;
    
    cout << "Level Set: " << endl;
    for(Int b = 0; b < tree.nblks; b++)
      { cout << tree.lvlset[b] << " " ;}
    cout << endl;

    cout << "Children: " << endl;
    for(Int b = 0; b < tree.nblks; b++)
      { cout << tree.child[b] << " " ;}
    cout << endl;
   
    cout << "Siblings: " << endl;
    for(Int b = 0; b < tree.nblks; b++)
      { cout << tree.sibling[b] << " " ;} 
    #endif

    //Step 4 
    //Setup schedule (S)
    Int lvl_idx = 0;
    for(Int l=0; l < tree.nlvls+1; l++)
    {
      Int lvl_counter = 0;
      for(Int t=0; t < num_threads; t++)
      {
        #ifdef BASKER_DEBUG_TREE
        printf("l %d t %d lvl_counter %d lvl_idx %d size: %d \n",
            l, t, lvl_counter ,lvl_idx, tree.nblks);
        #endif

        S[l][t] = tree.lvlset[lvl_idx];
        if(lvl_counter >= (pow(tree.nparts,l)-1))
        {
          lvl_idx++;
          lvl_counter = 0;
        }
        else
        {
          lvl_counter++;
        }
      }//end over threads
    }//end over nlvls

    //printf("===TEST==\n");

    #ifdef BASKER_DEBUG_TREE
    cout << "Schedule: " << endl; 
    for(Int l=0; l < tree.nlvls+1; l++)
    {
      for(Int t=0; t < num_threads; t++)
      {
        cout << S[l][t] << " , " ; 
      }//end over nhreads
      cout << endl;
    }//end over nlvls
    #endif
    
    //Step 5
    //Get workspace size
    for(Int l=0; l < tree.nlvls+1; l++)
    {
      for(Int t=0; t < num_threads; t++)
      {
        Int s_element = S[l][t];
        Int row_size = (tree.row_tabs[s_element+1] - 
            tree.row_tabs[s_element]);
        thread_array[t].iws_size += row_size;
        thread_array[t].ews_size += row_size;
      }//end over threads
    }//end over lvls
    
    //Step 6
    //Generate the do not pivot list
    //Don't need anymore

    /*
    Int leafStart = pow(tree.nparts, tree.nlvls);
    for(Int i = (tree.nblks-1); i > leafStart-1; i--)
      {
        for(Int j = tree.row_tabs[tree.lvlset[i]]; 
            j < tree.row_tabs[tree.lvlset[i]+1]; j++)
          {
            #ifdef BASKER_DEBUG_TREE
            //cout << "DENSE ROW AT: " << j << endl;
            #endif
            tree.pivot[j] = false;
          }
      }
    */

    //Comeback and clean this up!!!
    if(Options.btf == BASKER_FALSE)
    {
      //printf("Given A End\n");
      tree.treetab[tree.nblks-1] = A.ncol;
    }
    else
    {
      //printf("Given BTF_A End\n");
      tree.treetab[tree.nblks-1] = BTF_A.ncol;
    }
    FREE(temptree);

    tree_flag = true;
    return 0;

  }//end init_tree_thread()


  //function init tree using omp threads to determine tree size
  //OLD DEFUNCT
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::init_tree_thread
  (
   Int *perm, Int nblks, 
   Int parts,Int *row_tabs, Int *col_tabs, Int *treetab
  )
  {

    #ifdef BASKER_DEBUG_TREE
    printf("----init_tree_thread----\n");
    #endif

    //test
    //printf("---------------test -----------------\n");
    //basker_thread<Int,Entry,Exe_Space> tt();
    //printf("-----------after test--------------\n");


    INT_1DARRAY ttemp;

    //Calculate maximum number of levels
    // parts^x = threads
    Int maxlvls = log(num_threads)/log(parts);
    Int treelvls = log(nblks+1)/log(parts) - 1;
    //Number of blocks now is parts^lvl+1 - 1
    Int newnblks = pow(parts, maxlvls+1) - 1;

    if(newnblks < nblks)
    {
      cout << "Reducing nblks to: " << newnblks << 
        " in order to fit threads" <<endl;
    }
    else if(newnblks > nblks)
    {
      cout << "Error: need more blocks to fit threads" << endl;
      return -1;
    }

    malloc_tree(parts, newnblks, maxlvls);

    //thread array
    //MALLOC_THREAD_1DARRAY(thread_array, num_threads);
    //cout << "Done alloc thread_array" << endl;
    //Step 1
    //Build new row and col tabs (DFS)
    Int old_pos = nblks;
    Int new_pos = tree.nblks;

    BASKER_ASSERT((tree.nblks+1)>0, "tree nblks2");
    INT_1DARRAY temp_tree;
    MALLOC_INT_1DARRAY(temp_tree, tree.nblks+1);
    init_value(temp_tree, tree.nblks+1, (Int)0);
    INT_1DARRAY temp_col_tabs;
    MALLOC_INT_1DARRAY(temp_col_tabs, tree.nblks+1);
    init_value(temp_col_tabs, tree.nblks+1, (Int)0);
    INT_1DARRAY temp_row_tabs;
    MALLOC_INT_1DARRAY(temp_row_tabs, tree.nblks+1);
    init_value(temp_row_tabs, tree.nblks+1, (Int)0);
    
    rec_tabs(0, tree.nlvls, treelvls, tree.nblks, old_pos, &new_pos,
             col_tabs, row_tabs, treetab,
             temp_col_tabs, temp_row_tabs, temp_tree);

    for(Int i=0; i < tree.nblks+1; i++)
    {
      tree.row_tabs[i] = temp_row_tabs[i];
      tree.col_tabs[i] = temp_col_tabs[i];
      tree.treetab[i] = temp_tree[i];
    }

    FREE(temp_tree);
    FREE(temp_col_tabs);
    FREE(temp_row_tabs);

    #ifdef BASKER_DEBUG_TREE
    cout << "row_tabs: " << endl;
    for(Int i = 0; i < tree.nblks+1; i++)
      {cout << tree.row_tabs[i] << " " ;}
    cout << endl;
    cout << "col_tabs: " << endl;
    for(Int i = 0; i < tree.nblks+1; i++)
      {cout << tree.col_tabs[i] << " ";}
    cout << endl;
    cout << "tree_tabs: " << endl;
    for(Int i = 0; i < tree.nblks+1; i++)
      {cout << tree.treetab[i] << " ";}
    cout << endl;
    cout << "SIZES: " << endl;
    for(Int i = 0; i < tree.nblks; i++)
      {cout << i << "  " << tree.col_tabs[i+1]-tree.col_tabs[i] << endl;}
    #endif

    BASKER_ASSERT(num_threads > 0, "tree num_threads2");
    MALLOC_THREAD_1DARRAY(thread_array, num_threads);

    //Step 2
    //Build treeptr list
    for(Int i=0; i < tree.nblks; i++)
    {
      for(Int j = tree.col_tabs[i]; j < tree.col_tabs[i+1]; j++)
      {tree.treeptr[j] = i;}
    }
    
    //Step 3
    //Get level sets  (BFS) 
    //Note: needs refactored for nparts
    INT_1DARRAY temptree;
    BASKER_ASSERT(tree.nblks > 0, "tree nblks 3");
    MALLOC_INT_1DARRAY(temptree, tree.nblks);
    init_value(temptree, tree.nblks, (Int)0);

    new_pos = tree.nblks-1;
    temptree[new_pos] = tree.treetab[new_pos];
    tree.lvlset[new_pos] = new_pos;
    
    new_pos--;
    for(Int l = 1; l < tree.nlvls+1; l++)
    {
      Int nparents = pow(tree.nparts, l-1);
      Int pp = new_pos;
      //ttemp = lvl_task[l];
      Int task_pos = 0;
      Int task_offset = pow(tree.nparts, tree.nlvls-l);
      for(Int p = 0 ; p < nparents; p++)
      {
        Int parentptr = pp+nparents-p;
        Int leftc     = tree.lvlset[parentptr]-1;
        Int rightc    = tree.lvlset[parentptr] - 
          pow(tree.nparts, tree.nlvls-l+1);
        tree.lvlset[new_pos] = leftc;
        temptree[new_pos] = tree.treetab[leftc];
        new_pos--;
        tree.lvlset[new_pos] = rightc;
        temptree[new_pos] = tree.treetab[rightc];
        tree.child[tree.lvlset[parentptr]] = rightc;
        tree.sibling[rightc] = leftc;
        new_pos--; 
      }//over all parents
    }// over all levels

    #ifdef BASKER_DEBUG_TREE
    cout << "Treetabs: " << endl;
    for(Int b = 0 ; b < tree.nblks; b++)
      { cout << temptree[b] << " " ;}
    cout << endl;
    
    cout << "Level Set: " << endl;
    for(Int b = 0; b < tree.nblks; b++)
      { cout << tree.lvlset[b] << " " ;}
    cout << endl;

    cout << "Children: " << endl;
    for(Int b = 0; b < tree.nblks; b++)
      { cout << tree.child[b] << " " ;}
    cout << endl;
   
    cout << "Siblings: " << endl;
    for(Int b = 0; b < tree.nblks; b++)
      { cout << tree.sibling[b] << " " ;} 
    cout << endl;
    #endif

    //Step 4 
    //Setup schedule (S)
    Int lvl_idx = 0;
    for(Int l=0; l < tree.nlvls+1; l++)
    {
      Int lvl_counter = 0;
      for(Int t=0; t < num_threads; t++)
      {
        #ifdef BASKER_DEBUG_TREE
        printf("l %d t %d lvl_counter %d lvl_idx %d size: %d \n",
            l, t, lvl_counter ,lvl_idx, tree.nblks);
        #endif

        S[l][t] = tree.lvlset[lvl_idx];
        if(lvl_counter >= (pow(tree.nparts,l)-1))
        {
          lvl_idx++;
          lvl_counter = 0;
        }
        else
        {
          lvl_counter++;
        }
      }//end over threads
    }//end over nlvls

    #ifdef BASKER_DEBUG_TREE
    cout << "Schedule: " << endl; 
    for(Int l=0; l < tree.nlvls+1; l++)
    {
      for(Int t=0; t < num_threads; t++)
      {
        cout << S[l][t] << " , " ; 
      }//end over nhreads
      cout << endl;
    }//end over nlvls
    #endif
    
    
    //Step 5
    //Get workspace size
    for(Int l=0; l < tree.nlvls+1; l++)
    {
      for(Int t=0; t < num_threads; t++)
      {
        Int s_element = S[l][t];
        Int row_size = (tree.row_tabs[s_element+1] - tree.row_tabs[s_element]);
        thread_array[t].iws_size += row_size;
        thread_array[t].ews_size += row_size;
      }//end over threads
    }//end over lvls
    
    //Step 6
    //Generate the do not pivot list
    Int leafStart = pow(tree.nparts, tree.nlvls);
    for(Int i = (tree.nblks-1); i > leafStart-1; i--)
    {
      for(Int j = tree.row_tabs[tree.lvlset[i]]; 
          j < tree.row_tabs[tree.lvlset[i]+1]; j++)
      {
        #ifdef BASKER_DEBUG_TREE
        cout << "DENSE ROW AT: " << j << endl;
        #endif
        tree.pivot[j] = false;
      }
    }
    
    FREE(temptree);

    tree_flag = true;
    return 0;
  }//init_tree_omp()


  //recursive function to find col and row pos
  //come back and write as an iterative function
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::rec_tabs
  (
   Int lvl, Int mlvl, Int tlvl, Int mynum,
   Int old_pos, Int *new_pos,
   Int* old_col, Int* old_row, 
   Int* old_tree,
   INT_1DARRAY new_col, INT_1DARRAY new_row, 
   INT_1DARRAY new_tree
  )
  {
    //Set Values
    #ifdef BASKER_DEBUG_TREE
    printf("Rec Tabs, new_pos: %d old_pos %d\n", 
               *new_pos, old_pos);
    printf("Rec Tabs, oldc: %d oldr: %d mynum %d \n",
               old_col[old_pos], old_row[old_pos], mynum);
    #endif

    new_col[*new_pos]  = old_col[old_pos];
    new_row[*new_pos]  = old_row[old_pos];
    new_tree[(*new_pos)-1] = mynum;
    //Note: needs to be refactored for nparts
    if(lvl < mlvl)
    {
      Int rightc = old_pos - 1;
      Int leftc  = old_pos - pow(tree.nparts, tlvl-(lvl));
      #ifdef BASKER_DEBUG_TREE
      printf("Left Child, old_pos: %d tlvl: %d lvl: %d parts: %d \n", old_pos, tlvl, lvl, tree.nparts);
      printf("Left child: %d \n", leftc);
      #endif

      mynum = (*new_pos)-1;
      lvl++;
      (*new_pos)--;
      rec_tabs(lvl, mlvl, tlvl, mynum, rightc, new_pos, 
          old_col, old_row, old_tree,
          new_col, new_row, new_tree);

      (*new_pos)--;


      rec_tabs(lvl, mlvl, tlvl, mynum, leftc, new_pos,
          old_col, old_row, old_tree,
          new_col, new_row, new_tree);
    }

  }//rec_tabs()


  //Still need to add
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::init_tree_lvl
  (
   Int *perm, 
   Int nblks, 
   Int parts,
   Int *row_tabs, 
   Int *col_tabs, 
   Int *treetab, 
   Int lvls
  )
  {
    return -1;
  }//init_tree_lvl()


  //function to update the do not pivot list as we move up the factor tree
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::update_pivot_tree(Int start, Int end)
  {
    Int b;
    for(b = start; b < end; b++)
    {
      Int bb = tree.lvlset[b];
      for(Int i = tree.row_tabs[bb]; i < tree.row_tabs[bb+1]; i++)
      {
        tree.pivot[i] = true;
      }//Allow rows in this level to pivot
    }//end for each block
  }//end update_pivot_tree()


  //function to update the local pinv into paren -- reformat
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::update_lpinv_tree(Int start, Int end)
  {
    return -1;
  }//end update_lpinv_tree


  // ALM : strictly-lower triangular blocks (no diagonal blocks)
  //  LL : -> L-factor
  // AVM :          upper triangular blocks (with diagonal blocks)
  //  LU : -> U-factor
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry, Exe_Space>::matrix_to_views_2D
  (
   BASKER_MATRIX &flat
  )
  {
    #ifdef BASKER_DEBUG_TREE
    printf("Called matrix_to_views_2d \n");
    #endif
    
    //----------Count number of off-diagonal blks--------//
    Int nnzblks = 0 ;
    Int nblks = tree.nblks;
    INT_1DARRAY U_view_count;
    INT_1DARRAY L_view_count;
    BASKER_ASSERT(nblks> 0, "tree nblks 12");
    MALLOC_INT_1DARRAY(U_view_count, nblks);
    MALLOC_INT_1DARRAY(L_view_count, nblks);
    init_value(U_view_count, nblks, (Int)0);
    init_value(L_view_count, nblks, (Int)0);

    //Should also do LL and LU while it
    #ifdef BASKER_DEBUG_TREE
    flat.info();
    printf("nblks: %d \n", nblks);
    #endif

    // ----- kind of hack -----
    // In some rare-cases (nblk == ncol), 
    // which messes up the termination of the following j-loop
    // to avoid it, tree[nblks-1] is set to be..
    Int last_tree_tab = tree.treetab[nblks-1];
    tree.treetab[nblks-1] = -flat.ncol;
    //#define MY_DEBUG
    #ifdef MY_DEBUG
    printf( "\n ** nblks = %d, flat.ncol = %d **\n",nblks,flat.ncol );   
    for(Int i=0; i < nblks; i++) printf( " treetab[%d] = %d\n",i,tree.treetab[i] );
    #endif
    // ------------------------

    // i <= j: (i,j) (lowertri) and (j,i) (uppertri) are block index pairs of partitioned 2D blocks
    // Note block structure is symmetric
    for(Int i=0; i < nblks; i++)
    {
      for(Int j=i; j != -flat.ncol; nnzblks++, j=(tree.treetab[j]))
      {
        U_view_count[j] = U_view_count[j] + 1;
        L_view_count[i] = L_view_count[i] + 1;
      }
    }

    #ifdef BASKER_DEBUG_TREE
    printf("Make Hier View of size: %d %d \n", nblks, nnzblks);
    #endif    

    //Malloc  AU and AL views, and needed space
    BASKER_ASSERT(nblks > 0, "tree nblks 99");
    //Malloc AU and AL, matrix
    MALLOC_MATRIX_2DARRAY(AVM, nblks);
    MALLOC_MATRIX_2DARRAY(ALM, nblks);
    //Malloc LU and LL matrix ARRAYs
    MALLOC_MATRIX_2DARRAY(LU, nblks);
    MALLOC_MATRIX_2DARRAY(LL, nblks);

    for(Int i=0; i < nblks; i++)
    {

      // Malloc AU subarray
      // NOTE: size at least one to allow empty block
      Int U_view_size = (U_view_count(i) > 0 ? U_view_count(i) : 1);
      if (U_view_size > 0)
      {
        MALLOC_MATRIX_1DARRAY(AVM(i),     U_view_size);
        MALLOC_MATRIX_1DARRAY(LU(i),      U_view_size);
      }
      //Malloc AL subarray
      // NOTE: size at least one to allow empty block
      Int L_view_size = (L_view_count(i) > 0 ? L_view_count(i): 1);
      if (L_view_size > 0) 
      {
        MALLOC_MATRIX_1DARRAY(ALM(i),     L_view_size);
        MALLOC_MATRIX_1DARRAY(LL(i),      L_view_size);
      }

      LU_size(i) = U_view_count(i);
      LL_size(i) = L_view_count(i);
      U_view_count(i) = 0;
      L_view_count(i) = 0;
    }
    
    //Loop over again and fill
    nnzblks = 0;

    for(Int i=0; i < nblks; i++)
    {
      #ifdef MY_DEBUG
      printf( "\n >> i = %d <<\n",i );
      #endif
      for(Int j=i; j != -flat.ncol; j=tree.treetab[j])
      {
        MATRIX_1DARRAY &UMtemp = AVM[j];
        MATRIX_1DARRAY &LMtemp = ALM[i];

        MATRIX_1DARRAY &LUtemp = LU[j];
        MATRIX_1DARRAY &LLtemp = LL[i];

        #ifdef MY_DEBUG
        printf( " AVM(%d)(%d).set_shape(%dx%d)\n",j,U_view_count[j], tree.col_tabs[i+1]-tree.col_tabs[i],tree.col_tabs[j+1]-tree.col_tabs[j] );
        #endif
        UMtemp[U_view_count[j]].set_shape(
            tree.col_tabs[i], 
            (tree.col_tabs[i+1]-tree.col_tabs[i]),
            tree.col_tabs[j], 
            (tree.col_tabs[j+1]-tree.col_tabs[j])
            );
        if (tree.col_tabs[i+1] >= tree.col_tabs[i]) {
          UMtemp[U_view_count[j]].init_col();
        }

        //Int LU
        #ifdef MY_DEBUG
        printf( " LU (%d)(%d).set_shape(%dx%d)\n",j,U_view_count[j], tree.col_tabs[i+1]-tree.col_tabs[i],tree.col_tabs[j+1]-tree.col_tabs[j] );
        #endif
        LUtemp[U_view_count[j]].set_shape(
            tree.col_tabs[i], 
            (tree.col_tabs[i+1]-tree.col_tabs[i]),
            tree.col_tabs[j], 
            (tree.col_tabs[j+1]-tree.col_tabs[j])
            );        

        U_view_count[j] = U_view_count[j]+1;

        #ifdef MY_DEBUG
        printf( " ALM(%d)(%d).set_shape(%dx%d)\n",i,L_view_count[i], tree.col_tabs[j+1]-tree.col_tabs[j],tree.col_tabs[i+1]-tree.col_tabs[i] );
        #endif
        LMtemp[L_view_count[i]].set_shape(
            tree.col_tabs[j],
            (tree.col_tabs[j+1]-tree.col_tabs[j]),
            tree.col_tabs[i],
            (tree.col_tabs[i+1] - tree.col_tabs[i])
            );
        if (tree.col_tabs[i+1] >= tree.col_tabs[i]) {
          LMtemp[L_view_count[i]].init_col();
        }

        //Init LL
        #ifdef MY_DEBUG
        printf( " LL (%d)(%d).set_shape(%dx%d)\n",i,L_view_count[i], tree.col_tabs[j+1]-tree.col_tabs[j],tree.col_tabs[i+1]-tree.col_tabs[i] );
        #endif
        LLtemp[L_view_count[i]].set_shape(
            tree.col_tabs[j],
            (tree.col_tabs[j+1]-tree.col_tabs[j]),
            tree.col_tabs[i],
            (tree.col_tabs[i+1] - tree.col_tabs[i])
            );

        L_view_count[i] = L_view_count[i]+1;

        nnzblks++;
        #ifdef MY_DEBUG
        printf( " -> j = tree.treetab[%d] = %d\n",j,tree.treetab[j] );
        #endif
      }//end over all inner
    }//end over all blks

    // ----- kind of hack -----
    // revert the last tree tab
    tree.treetab[nblks-1] = last_tree_tab;
    // ------------------------

    FREE(U_view_count);
    FREE(L_view_count); 
  }//end matrix_to_view_2D()
  

  //Translate A -> B
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int, Entry, Exe_Space>::matrix_to_views
  (
   BASKER_MATRIX &flat, 
   MATRIX_VIEW_2DARRAY &view
  )
  {

    #ifdef BASKER_DEBUG_TREE
    printf("Called matrix_to_view \n");
    #endif

    Int nnzblks = 0 ;
    Int nblks = tree.nblks;
    INT_1DARRAY view_count;
    BASKER_ASSERT(nblks > 0, "tree mtv nblks");
    MALLOC_INT_1DARRAY(view_count, nblks);
    init_value(view_count, nblks, (Int)0);
    
    MALLOC_MATRIX_VIEW_2DARRAY(view, nblks);

    for(Int i=0; i < nblks; i++)
    {
      for(Int j=i; j != flat.ncol; nnzblks++, j=(tree.treetab[j]))
      {
        view_count[j] = view_count[j] + 1;
      }
    }
   
    #ifdef BASKER_DEBUG_TREE
    printf("Make Hier View of size: %d %d \n", nblks, nnzblks);
    #endif    
  
    for(Int i=0; i < nblks; i++)
    {
      BASKER_ASSERT(view_count[i] > 0 , "tree vc matrix2");
      MALLOC_MATRIX_VIEW_1DARRAY(view[i], view_count[i]);
      LU_size[i] = view_count[i];
      view_count[i] = 0;
    }     

    nnzblks = 0;
    for(Int i=0; i < nblks; i++)
    {
      for(Int j=i; j != flat.ncol; j=tree.treetab[j])
      {
        //M_V temp = view[j];
        MATRIX_VIEW_1DARRAY temp = view[j];
        temp[view_count[j]].init(&(flat), 
            tree.col_tabs[i], 
            (tree.col_tabs[i+1] - tree.col_tabs[i]),
            tree.col_tabs[j], 
            (tree.col_tabs[j+1]-tree.col_tabs[j])
            );
        view_count[j] = view_count[j]+1;
        nnzblks++;
      }
    }

    FREE(view_count);
  }//end matrix_to_views()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::find_2D_convert
  (
   BASKER_MATRIX &M
  )
  {
    //This will scan the whole matrix and find the nnz and cuts
    //In parallel the values will be copied over

    Int L_col = 0; //Lower Blk col
    Int L_row = 0; //Lower Blk row
    Int U_col = 0; //Upper Blk col
    Int U_row = 0; //Upper Blk row
    Int c_idx = 0; //Col tab offset; used to iterate through tree.col_tabs

    #ifdef MY_DEBUG
    printf( "\n >> find_2D_convert (%d x %d) <<\n\n",M.nrow,M.ncol );
    for (Int k = 0; k <= tree.nblks; k++) printf( "  row_tabs[%d] = %d\n",k,tree.row_tabs(k));
    for (Int k = 0; k <  tree.nblks; k++) printf( "   LU_size[%d] = %d\n",k,LU_size(k));
    M.print_matrix("M.dat");
    #endif
    for(Int k = 0; k < M.ncol; ++k)
    {
      // rest row indexes
      L_row = 0;
      U_row = 0;
      Int r_idx = 0; //used to iterate through tree.row_tabs

      //Fast-forward to first column tab in this column
      while((k >= tree.col_tabs(c_idx+1)) ||
            (c_idx < tree.nblks && tree.row_tabs(c_idx+1) == tree.row_tabs(c_idx))) // skip empty blocks
      {
        c_idx++;
        L_col++;
        U_col++;

        //r_idx++;
        //L_row++;
        //U_row++;
      }

      //Get the first blks
      BASKER_BOOL start_col = BASKER_TRUE;

      #ifdef MY_DEBUG
      printf( "\n >> k = %d, (c_idx=%d, L_col=%d, U_col=%d) <<\n",k,c_idx,L_col,U_col );
      #endif
      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) //offsets to row_idx - will yield range of row values to col k
      {
        Int j = M.row_idx(i); //first row id entry in col k
        #ifdef MY_DEBUG
        printf( "   + row_idx(%d) = %d vs. row_tabs(%d) = %d <<\n", i,j, r_idx+1,tree.row_tabs(r_idx+1) );
        #endif

        //Get right blk
        if(j > k) { //lower
          while((j >= tree.row_tabs(r_idx+1)) ||
                (r_idx < tree.nblks && tree.row_tabs(r_idx+1) == tree.row_tabs(r_idx))) // skip empty blocks
          {
            if((L_row+1 < LL_size(L_col)) &&
               (tree.row_tabs(r_idx+1) == ALM(L_col)(L_row+1).srow))
            {
              //printf( " > ALM(%d)(%d).srow = %d, row_tab(%d) = %d\n",L_col,L_row+1,ALM(L_col)(L_row+1).srow, r_idx+1,tree.row_tabs(r_idx+1) );
              L_row++;
              BASKER_ASSERT(L_row < LL_size(L_col), " Wrong L in A to 2d");
              start_col = BASKER_TRUE;
            }
            r_idx++;
          }
          //printf( "   -> lower(r_idx=%d)\n",r_idx );
        } else if(j <= k) { //upper
          while((j >= tree.row_tabs(r_idx+1)) ||
                (r_idx < tree.nblks && tree.row_tabs(r_idx+1) == tree.row_tabs(r_idx))) // skip empty blocks
          {
            if((U_row+1 < LU_size(U_col)) &&
               (tree.row_tabs(r_idx+1) == AVM(U_col)(U_row+1).srow))
            {
              //printf( " + AVM(%d)(%d).srow = %d, row_tab(%d) = %d\n",U_col,U_row+1,AVM(U_col)(U_row+1).srow, r_idx+1,tree.row_tabs(r_idx+1) );
              U_row++;
              BASKER_ASSERT(U_row < LU_size(U_col), " Wrong U in A to 2d");
              start_col = BASKER_TRUE;
            }
            r_idx++;
          }
          //printf( "   -> upper (r_idx=%d)\n",r_idx );
        }
        #ifdef MY_DEBUG
        std::cout << " > " << (j > k ? " lower" : " upper" ) << " --> "
                  << "  L(" << L_col << ", " << L_row << ") "
                  << "  U(" << U_col << ", " << U_row << ") "
                  << " with j = " << j << " and k = " << k
                  << ", LU_size( " << U_col << " ) = " << LU_size(U_col)
                  << " r_idx = " << r_idx
                  << " start_col = " << start_col
                  << std::endl;
        #endif


        //Get Matrix Ref
        BASKER_MATRIX &Ltemp = ALM(L_col)(L_row);
        BASKER_MATRIX &Utemp = AVM(U_col)(U_row);
        Int bcol  = Ltemp.scol;

        //diag blk
        if(L_row == 0 && U_row == LU_size(U_col)-1)
        {
          if(start_col == BASKER_TRUE)
          {
            Ltemp.col_ptr(k-bcol) = i;
            #ifdef MY_DEBUG
            std::cout << " > L ( " << L_col << ", " << L_row << " ).col_ptr( " << k - bcol << ") = " << i << " with k = " << k << " and j = " << j << std::endl;
            #endif
          }
          Ltemp.nnz = Ltemp.nnz+1;
        }
        else //offdig
        {
          if (j > k)
          {
            if(start_col == BASKER_TRUE)
            {
              Ltemp.col_ptr(k-bcol) = i;
              #ifdef MY_DEBUG
              std::cout << " + L ( " << L_col << ", " << L_row << " ).col_ptr( " << k - bcol << ") = " << i << std::endl;
              #endif
            }
            Ltemp.nnz = Ltemp.nnz+1;
          } else if (j < k)
          {
            if(start_col == BASKER_TRUE)
            {
              Utemp.col_ptr(k-bcol) = i;
              #ifdef MY_DEBUG
              std::cout << " + U ( " << U_col << ", " << U_row << " ).col_ptr( " << k - bcol << ") = " << i << std::endl;
              #endif
            }
            Utemp.nnz = Utemp.nnz+1;
          } else //if(j == k)
          {
            // printf("Error: L: %d %d U: %d %d Usize: %d \n",L_col, L_row, U_col, U_row, LU_size[U_col]);
            std::cout << std::endl
                      << "Error: L(" << L_col << ", " << L_row << ") "
                            << " U(" << U_col << ", " << U_row << ") "
                      << " with j = " << j << ", k = " << k << ", and " 
                      << " Usize = " << LU_size[U_col] 
                      << std::endl << std::endl;
            BASKER_ASSERT(0==1, "A2D offdiag blk with diag element");
          }
        }
        start_col = BASKER_FALSE;
      }//over each row        
    }//over each colunm

  }//end find_2d_convert()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::clean_2d()
  {

    for(Int b = 0 ; b < tree.nblks; ++b)
    {
      for(Int sb = 0; sb < LL_size(b); ++sb)
      {
        //printf( " ALM(%d)(%d).clean_col()\n",b,sb );
        ALM(b)(sb).clean_col();
      }
      for(Int sb = 0; sb < LU_size(b); ++sb)
      {
        AVM(b)(sb).clean_col();
      }
    }//for - over all blks

    return 0;
  }//end clean_2d


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::sfactor_copy()
  {
    //Reorder A;
    //Match order
    if(match_flag == BASKER_TRUE)
    {
      permute_row(A,order_match_array);
    }
    sort_matrix(A);

    //BTF order
    if(btf_flag == BASKER_TRUE)
    {
      permute_col(A,order_btf_array);
      permute_row(A,order_btf_array);
      break_into_parts(A, btf_nblks, btf_tabs);
    }

    //ND order
    if(nd_flag == BASKER_TRUE)
    {
      if(btf_tabs_offset != 0)
      {
        sort_matrix(BTF_A);

        //Permute the A
        permute_row(BTF_A, part_tree.permtab);
        permute_col(BTF_A, part_tree.permtab);

        //Permute the B
        if(btf_nblks > 1)
        {
          permute_row(BTF_B, part_tree.permtab);
        }

        sort_matrix(BTF_A);
        if(btf_nblks > 1)
        {
          sort_matrix(BTF_B);
          sort_matrix(BTF_C);
        }
      }
    }

    //AMD
    if(amd_flag == BASKER_TRUE)
    {
      if(btf_tabs_offset != 0)
      {
        //Permute A
        permute_col(BTF_A, order_csym_array);
        sort_matrix(BTF_A);
        permute_row(BTF_A, order_csym_array);
        sort_matrix(BTF_A);

        //Permute B
        if(btf_nblks > 1)
        {
          permute_row(BTF_B, order_csym_array);
          sort_matrix(BTF_B);
        }
      }
    }

    if(btf_tabs_offset != 0)
    {
      //=====Move into 2D ND-Structure====/
      //Find submatices view shapes
      clean_2d();

      //matrix_to_views_2D(BTF_A);
      //Find starting point
      find_2D_convert(BTF_A);

      //Fill 2D structure
      #ifdef BASKER_KOKKOS
      BASKER_BOOL keep_zeros = BASKER_FALSE;
      kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this, BASKER_FALSE, keep_zeros);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
      Kokkos::fence();
      #else
      //Comeback
      #endif
    }

    // Initialize C & B blocks
    {
      sort_matrix(BTF_C);
      permute_col(BTF_C, order_c_csym_array);
      sort_matrix(BTF_C);
      permute_row(BTF_C, order_c_csym_array);
      sort_matrix(BTF_C);

      if(btf_tabs_offset != 0)
      {
        permute_col(BTF_B, order_c_csym_array);        
        sort_matrix(BTF_B);
      }
    }

    //test
    //printMTX("A_BTF.mtx", BTF_A);

    return 0;
  }//sfactor_copy()


  // NDE: sfactor_copy2 is now only responsible for mapping blocks to 2D blocks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::sfactor_copy2(bool alloc_BTFA, bool copy_BTFA)
  {
    //Timers
    #ifdef BASKER_TIMER_FINE
    std::ios::fmtflags old_settings = cout.flags();
    int old_precision = std::cout.precision();
    std::cout.setf(ios::fixed, ios::floatfield);
    std::cout.precision(8);
    double sort_time = 0.0;
    double tmp_time = 0.0;
    Kokkos::Timer timer_permute;
    Kokkos::Timer timer_sort;
    #endif

    if(copy_BTFA && btf_tabs_offset != 0)
    {
      //=====Move into 2D ND-Structure====/
      //Find submatices view shapes
      if(Options.verbose == BASKER_TRUE)
      { printf("btf tabs reorder\n"); }

      #ifdef BASKER_TIMER_FINE 
      double twod_time = 0.0;
      Kokkos::Timer timer_twod;
      #endif

      clean_2d(); // clear vals from ALM, AVM - views of views that store the local 2D block CCS reordered matrix info

      //matrix_to_views_2D(BTF_A);
      //Find starting point
      find_2D_convert(BTF_A); //prepare CCS 'sizes' of each ALM(i)(j), AVM(i)(j) (nnz, col_idx, )

      //Fill 2D structure
      #ifdef BASKER_KOKKOS
      BASKER_BOOL keep_zeros = BASKER_FALSE;
      BASKER_BOOL alloc      = alloc_BTFA; //BASKER_FALSE;
      kokkos_order_init_2D<Int,Entry,Exe_Space> iO(this, alloc, keep_zeros); // t_init_2DA; fill row_idx, vals into ALM, AVM calling convert2D
      Kokkos::parallel_for(TeamPolicy(num_threads,1), iO);
      Kokkos::fence();
      #else
      //Comeback
      #endif

      #ifdef BASKER_TIMER_FINE
      tmp_time = timer_twod.seconds();
      twod_time += tmp_time;
      std::cout << "    Basker move into 2D ND reorder time: " << tmp_time << std::endl;
      #endif
    }

    if(Options.verbose_matrix_out == BASKER_TRUE)
    {
      printMTX("C_Factor.mtx", BTF_C);
    }

    //If same pattern, permute using pivot, and reset
    #ifdef BASKER_TIMER_FINE 
    double gperm_time = 0.0;
    Kokkos::Timer timer_gperm;
    #endif
    if((Options.same_pattern == BASKER_TRUE))
    {
      if(same_pattern_flag == BASKER_FALSE)
      {
        MALLOC_INT_1DARRAY(gperm_same, gn);
        for(Int i = 0; i < gn; i++)
        {
          gperm_same(i) = gperm(i);
          gperm(i) = BASKER_MAX_IDX;
        }
        same_pattern_flag = BASKER_TRUE;
      }
    }
    else
    {
      for(Int i = 0; i < gn; ++i)
      {
        gperm(i) = BASKER_MAX_IDX;
      }
    }
    #ifdef BASKER_TIMER_FINE
    tmp_time = timer_gperm.seconds();
    gperm_time += tmp_time;
    std::cout << "    Basker gperm (pivot) reset time: " << tmp_time << std::endl;
    timer_gperm.reset();
    #endif

    // reset all the time (factor may have failed, and some variables may not be cleared)
    //if(factor_flag == BASKER_TRUE)
    {
      typedef Kokkos::TeamPolicy<Exe_Space> TeamPolicy;
      kokkos_reset_factor<Int,Entry,Exe_Space> reset_factors(this); //t_reset_ND_factor, BTF; reset LL and LU for factorization factor_notoken step
      Kokkos::parallel_for(TeamPolicy(num_threads,1), reset_factors);
      Kokkos::fence();
    }
    #ifdef BASKER_TIMER_FINE
    tmp_time = timer_gperm.seconds();
    std::cout << "    Basker 2D reset_factors time: " << tmp_time << std::endl;
    std::cout << "    Basker sorts total time: " << sort_time << std::endl;
    std::cout.precision(old_precision);
    std::cout.flags(old_settings);
    #endif

    return 0;
  }//sfactor_copy2()

}//end namespace basker

#endif //end ifndefbasker_tree_hpp
