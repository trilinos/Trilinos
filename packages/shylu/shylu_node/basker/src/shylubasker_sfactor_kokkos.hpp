// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_SFACTOR_KOKKOS_HPP
#define SHYLUBASKER_SFACTOR_KOKKOS_HPP

#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_types.hpp"

#include <Kokkos_Core.hpp>
//#include "shylubasker_sfactor_kokkos.hpp" //is this correct??

// BUG: This file contains types that do not exist in current Basker impl
// Likely carryover from earlier experimental code before being ported to Trilinos

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  struct workspace_init_funct
  {
    //Kokkos Typedefs

    //Basker Typedefs
    //Matrix data structures
    typedef  BaskerMatrix<Int, Entry, Exe_Space>                  MATRIX;
    typedef  BaskerMatrixView<Int, Entry, Exe_Space>              MATRIX_VIEW;
    typedef  basker_2Darray<Int, MATRIX_VIEW, Exe_Space>          MATRIX_VIEW_ARRAY;
    typedef  basker_2Darray<Int, MATRIX, Exe_Space>               MATRIX_ARRAY;

    //Arrays
    typedef  basker_2Darray<Int, Int,   Exe_Space>                INT_2DARRAY;
    typedef  basker_2Darray<Int, Entry, Exe_Space>                ENTRY_2DARRAY;
    typedef  basker_1Darray<Int, Int,   Exe_Space>                INT_1DARRAY;
    typedef  basker_1Darray<Int, Entry, Exe_Space>                ENTRY_1DARRAY;
   
    //Tree
    typedef  basker_tree<Int, Entry, Exe_Space>                   TREE;

    
    INT_2DARRAY   ws;
    ENTRY_2DARRAY X;
    TREE          tree;
    INT_2DARRAY   S;
    Int           lvl;

    KOKKOS_INLINE_FUNCTION
    workspace_init_funct()
    {}

    KOKKOS_INLINE_FUNCTION
    workspace_init_funct
    (
     INT_2DARRAY   _ws,
     ENTRY_2DARRAY _X, 
		 TREE          _tree, 
     INT_2DARRAY   _S,
     Int           _lvl
    )
    {
      ws = _ws;
      X  = _X;
      tree = _tree;
      S = _S;
      lvl  = _lvl;
    }// workspace_init_funct()

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const
    {
      printf("Kokkos: workspace: %d lvl: %d \n", i, lvl);
    } // operator()
   
  };

}//end Basker namespace

#endif //end ifndef
