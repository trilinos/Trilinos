// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_TYPES_HPP
#define SHYLUBASKER_TYPES_HPP

#include <exception>
#include "Teuchos_TestForException.hpp"

//#define BASKER_DEBUG

//MACRO TURN ON FUCNTIONS
#define BASKER_ATOMIC         //Use Atomics (OLD)
#define BASKER_ATOMIC_2       //Use Atomics (OLD)
#define BASKER_NO_LAMBDA      //Do not use lambda
#define BASKER_2DL            //Use a 2D L
#define BASKER_MULTIPLE_UPPER //Use multiple threads for upper
#define BASKER_MULTIPLE_LOWER //Use multiple threads for lower
#define HAVE_AMESOS           //Use Amesos orderings
#define BASKER_SCOTCH         //Use Scotch

//MACRO TYPE
#define BASKER_INT            Int
#define BASKER_ENTRY          Entry
#define BASKER_BOOL           bool
#define BASKER_TRUE           true
#define BASKER_FALSE          false
#define BASKER_NO_OP          ((void)0)
#define BASKER_MAX_IDX        -1         //What we check against
#define BASKER_ERROR          -1         //Return error no recov
#define BASKER_ERROR_RETRY    -2
#define BASKER_SUCCESS        0

#define BASKER_RESTART       gn

#define BASKER_EPSILON       1e-6

#define BASKER_PIVOT_TOL     0.001

//Error Codes
enum BASKER_ERROR_CODE 
{
  BASKER_ERROR_NOERROR,     //No error
  BASKER_ERROR_SINGULAR,    //Singular during factorization
  BASKER_ERROR_NAN,         //NaN during factorization
  BASKER_ERROR_REMALLOC,    //Need to be realloc
  BASKER_ERROR_NOMALLOC,    //Failed with nomalloc option
  BASKER_ERROR_OTHER
};

//Matching Types
enum BASKER_MATCHING_CODE
{
  BASKER_MATCHING_COMMON,
  BASKER_MATCHING_BN,
  BASKER_MATCHING_SUM,
  BASKER_MATCHING_PRODUCT,
  BASKER_MATCHING_EXP
};

//MACRO BTF METHOD
#define BASKER_BTF_MAX_PERCENT   1.00
#define BASKER_BTF_LARGE         500  //Made smaller for unit test
#define BASKER_BTF_IMBALANCE     0.10
#define BASKER_BTF_SMALL         100
#define BASKER_BTF_NNZ_OVER      2.0 //Upped from 1.20
#define BASKER_BTF_PRUNE_SIZE    100

#define BASKER_DOM_NNZ_OVER      1.0 //Added to control estimate for DOM blocks
#define BASKER_SEP_NNZ_OVER      2.0 //Added to control estimate for SEP blocks

enum BASKER_INCOMPLETE_CODE
{
  BASKER_INCOMPLETE_LVL,           //ilu(k) 
  BASKER_INCOMPLETE_RLVL,          //milu(k)
  BASKER_INCOMPLETE_RLVL_LIMITED,  //milu(k) -- no float on offd
  BASKER_INCOMPLETE_TOL,
  BASKER_INCOMPLETE_LVL_TOL,
  BASKER_INCOMPLETE_EXP
};

//MACRO DEFAULT INC SETTINGS
#define BASKER_INC_LVL_VALUE        0
#define BASKER_INC_TOL_VALUE      0.0001

//MACRO INC FILL (this will become dynamic in the future)
#define BASKER_FILL_USER           0.00
#define BASKER_FILL_LESTIMATE      1.50
#define BASKER_FILL_UESTIMATE      1.50
#define BASKER_FILL_LLOWERESTIMATE 2.00
#define BASKER_FILL_UUPPERESTIMATE 2.00
#define BASKER_FILL_LSEPESTIMATE   2.00
#define BASKER_FILL_USEPESTIMATE   2.00

//MACRO SYSTEM FUNCTIONS
/*#ifdef BASKER_DEBUG
 #include <assert.h>
 #define ASSERT(a)             assert(a)
#else
 //#define ASSERT(a)             BASKER_NO_OP
 #include <assert.h>
 #define ASSERT(a)           assert(a)
#endif*/

#ifdef BASKER_DEBUG
#define BASKER_ASSERT(a,s)       \
  {                              \
    if(!(a))                     \
      {printf("\n\n%s \nLINE: %d \nFILE: %s\n\n", s, __LINE__, __FILE__);} \
    assert(a);	                 \
    if(!(a))                     \
      exit(0);                   \
  }
#else
#define BASKER_ASSERT(a,s)      \
  {                             \
    TEUCHOS_TEST_FOR_EXCEPTION(!(a), \
      std::runtime_error, " ShyLUBasker:: error "+std::string(s)+"."); \
  }
#endif


//Note:  Should see if Kokkos has a fast memory cpy in place of for-loop
//MACRO ARRAY FUNCTIONS
//Execution Space
#include <Kokkos_Core.hpp>
#define BASKER_EXE_SPACE     Kokkos::DefaultHostExecutionSpace
//ReMacro Basker Classes
#define BASKER_SOLVER        Basker<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
#define BASKER_MATRIX        BaskerMatrix<BASKER_INT, BASKER_ENTRY, BASKER_EXE_SPACE>
#define BASKER_MATRIX_VIEW   BaskerMatrixView<BASKER_INT, BASKER_ENTRY, BASKER_EXE_SPACE>
#define BASKER_STATS         BaskerStats<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
//ReMacro Basker Structs
#define BASKER_TREE          basker_tree<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
#define BASKER_SYMBOLIC_TREE basker_symbolic_tree<BASKER_INT, BASKER_ENTRY, BASKER_EXE_SPACE>
#define BASKER_THREAD        basker_thread<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
//Macro Arrays
#define BASKER_KOKKOS_NOINIT      Kokkos::ViewAllocateWithoutInitializing
#define INT_RANK2DARRAY           Kokkos::View<BASKER_INT**,         BASKER_EXE_SPACE>
#define INT_1DARRAY               Kokkos::View<BASKER_INT*,          BASKER_EXE_SPACE>
#define ENTRY_1DARRAY             Kokkos::View<BASKER_ENTRY*,        BASKER_EXE_SPACE>
#define BOOL_1DARRAY              Kokkos::View<BASKER_BOOL*,         BASKER_EXE_SPACE>
#define BOOL_2DARRAY              Kokkos::View<BOOL_1DARRAY*,        BASKER_EXE_SPACE>

#define INT_2DARRAY               Kokkos::View<INT_1DARRAY*,          BASKER_EXE_SPACE>
#define ENTRY_2DARRAY             Kokkos::View<ENTRY_1DARRAY*,        BASKER_EXE_SPACE>
#define MATRIX_1DARRAY            Kokkos::View<BASKER_MATRIX*,        BASKER_EXE_SPACE>
#define MATRIX_2DARRAY            Kokkos::View<MATRIX_1DARRAY*,       BASKER_EXE_SPACE>
#define MATRIX_VIEW_1DARRAY       Kokkos::View<BASKER_MATRIX_VIEW*,   BASKER_EXE_SPACE>
#define MATRIX_VIEW_2DARRAY       Kokkos::View<MATRIX_VIEW_1DARRAY*,  BASKER_EXE_SPACE>
#define THREAD_1DARRAY            Kokkos::View<BASKER_THREAD*,        BASKER_EXE_SPACE>

#define INT_1DARRAY_PAIRS        Kokkos::View<std::pair<Int,Int>*,  BASKER_EXE_SPACE>
//Macro Memory Calls
//MALLOC
#define MALLOC_INT_1DARRAY_PAIRS(a,s)   \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC malloc_pairs_1d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {                 \
      Kokkos::resize(a, s);                               \
      if(a.data() == NULL)                                \
        throw std::bad_alloc();                           \
    }                                                     \
  }
#define MALLOC_INT_1DARRAY(a,s)   \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC int_1d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {         \
      Kokkos::resize(a, s);                       \
      if(a.data() == NULL)                        \
        throw std::bad_alloc();                   \
    }                                             \
  }
#define MALLOC_INT_RANK2DARRAY(a,s0,s1)   \
  { \
    BASKER_ASSERT(s0>0, "BASKER ASSERT MALLOC int_rank2d: size to alloc > 0 fails"); \
    BASKER_ASSERT(s1>0, "BASKER ASSERT MALLOC int_rank2d: size to alloc > 0 fails"); \
    if (Int(a.extent(0)) != s0 || Int(a.extent(1)) != s1) { \
      Kokkos::resize(a, s0,s1);                             \
      if(a.data() == NULL)                                  \
        throw std::bad_alloc();                             \
    }                                                       \
  }
#define MALLOC_INT_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0,"BASKER ASSERT MALLOC int_2d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {                                          \
      a = INT_2DARRAY(Kokkos::view_alloc("int_2d", Kokkos::SequentialHostInit),s); \
      if(a.data() == NULL)                                                         \
        throw std::bad_alloc();                                                    \
    }                                                                              \
  }
#define MALLOC_ENTRY_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC entry_1d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {             \
      Kokkos::resize(a, s);                           \
      if(a.data() == NULL)                            \
        throw std::bad_alloc();                       \
    }                                                 \
  }
#define MALLOC_ENTRY_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC entry_2d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {   \
      a = ENTRY_2DARRAY(Kokkos::view_alloc("matrix_2d", Kokkos::SequentialHostInit),s); \
      if(a.data() == NULL)                  \
        throw std::bad_alloc();             \
    }                                       \
  }
#define MALLOC_BOOL_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC bool_1d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {            \
      Kokkos::resize(a, s);                          \
      if(a.data() == NULL)                           \
        throw std::bad_alloc();                      \
    }                                                \
  }
#define MALLOC_BOOL_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC bool_2d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {  \
      Kokkos::resize(a, s);                \
      if(a.data() == NULL)                 \
        throw std::bad_alloc();            \
    }                                      \
  }
#define MALLOC_MATRIX_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC matrix_1d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {  \
      a = MATRIX_1DARRAY(Kokkos::view_alloc("matrix_1d", Kokkos::SequentialHostInit),s); \
      if(a.data() == NULL)                 \
        throw std::bad_alloc();            \
    }                                      \
  }
#define MALLOC_MATRIX_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC matrix_2d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {  \
      a = MATRIX_2DARRAY(Kokkos::view_alloc("matrix_2d", Kokkos::SequentialHostInit),s); \
      if(a.data() == NULL)                 \
        throw std::bad_alloc();            \
    }                                      \
  }
#define MALLOC_MATRIX_VIEW_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC matrix_view_1d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {  \
      a = MATRIX_VIEW_1DARRAY(Kokkos::view_alloc("matrix_view_1d", Kokkos::SequentialHostInit),s); \
      if(a.data() == NULL)                 \
        throw std::bad_alloc();            \
    }                                      \
  }
#define MALLOC_MATRIX_VIEW_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC matrix_view_2d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {  \
      a = MATRIX_VIEW_2DARRAY(Kokkos::view_alloc("matrix_view_2d", Kokkos::SequentialHostInit),s); \
      if(a.data() == NULL)                 \
        throw std::bad_alloc();            \
    }                                      \
  }
#define MALLOC_THREAD_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT MALLOC thread_1d: size to alloc >= 0 fails"); \
    if (s > 0 && Int(a.extent(0)) != s) {  \
      a = THREAD_1DARRAY(Kokkos::view_alloc("thread_1d", Kokkos::SequentialHostInit),s); \
      if(a.data() == NULL)                 \
        throw std::bad_alloc();            \
    }                                      \
  }

//REALLOC (no copy)
#define REALLOC_1DARRAY(a,os,s)          \
  { \
    BASKER_ASSERT(s >= 0, "BASKER ASSERT REALLOC 1D ARRAY: size to alloc >= 0 fails"); \
    Kokkos::realloc(a,s);                 \
  }
#define REALLOC_2DRRAAY(a,os1,os2,s1,s2) \
  { \
    BASKER_ASSERT(s1 >= 0 && s2 >= 0, "BASKER ASSERT REALLOC 2D ARRAY: size to alloc >= 0 fails"); \
    Kokkos::realloc(a,s1,s2);            \
  }
#define REALLOC_INT_1DARRAY(a,os,s)      REALLOC_1DARRAY(a,os,s)
#define REALLOC_ENTRY_1DARRAY(a,os,s)    REALLOC_1DARRAY(a,os,s)

//Set values
#define SET_INT_1DARRAY(a, b, s)    \
  {                                 \
  MALLOC_INT_1DARRAY(a,s);          \
  for(BASKER_INT MACRO_I = 0; MACRO_I< s; MACRO_I++)\
    a[MACRO_I] = b[MACRO_I];        \
  }
       
#define SET_ENTRY_1DARRAY(a, b,s)   \
  {                                 \
  MALLOC_ENTRY_1DARRAY(a,s);        \
  for(BASKER_INT MACRO_I =0; MACRO_I <s; MACRO_I++)\
    a[MACRO_I] = b[MACRO_I];        \
  }
  
#define SET_BOOL_1DARRAY(a, b, s)   \
  {                                 \
  MALLOC_BOOL_1DARRY(a,s);          \
  for(BASKER_INT MACRO_I =0; MACRO_I <s; MACRO_I++)\
    a[MACRO_I] = b[MACRO_I];        \
  }

#define FREE(a)                        BASKER_NO_OP

#define FREE_INT_1DARRAY_PAIRS(a) \
  {                               \
    Kokkos::resize(a,0);          \
  }

#define FREE_INT_1DARRAY(a) \
  {                         \
    Kokkos::resize(a,0);    \
  }

#define FREE_INT_RANK2DARRAY(a) \
  {                             \
    Kokkos::resize(a,0);        \
  }

#define FREE_INT_2DARRAY(a,n) \
  {                           \
    Kokkos::resize(a,0);      \
  }

#define FREE_ENTRY_1DARRAY(a) \
  {                           \
    Kokkos::resize(a,0);      \
  }

#define FREE_ENTRY_2DARRAY(a,n) \
  {                             \
    Kokkos::resize(a,0);        \
  }

#define FREE_BOOL_1DARRAY(a) \
  {                          \
    Kokkos::resize(a,0);     \
  }

#define FREE_BOOL_2DARRAY(a,n) \
  {                            \
    Kokkos::resize(a,0);       \
  }

#define FREE_MATRIX_1DARRAY(a) \
  {                            \
    Kokkos::resize(a,0);       \
  }

#define FREE_MATRIX_2DARRAY(a,n) \
  {                              \
    Kokkos::resize(a,0);         \
  }

#define FREE_MATRIX_VIEW_1DARRAY(a) \
  {                                 \
    Kokkos::resize(a,0);            \
  }

#define FREE_MATRIX_VIEW_2DARRAY(a,n) \
  {                                   \
    Kokkos::resize(a,0);              \
  }

#define FREE_THREAD_1DARRAY(a) \
  {                            \
    Kokkos::resize(a,0);       \
  }

//Inline command
#define BASKER_INLINE   inline
#define BASKER_LAMBDA   [&]

//Time Macro
#ifdef BASKER_TIME
#define BASKER_TIMER
#define BASKER_TIMER_FINE
#endif

#endif //end basker_types_hpp
