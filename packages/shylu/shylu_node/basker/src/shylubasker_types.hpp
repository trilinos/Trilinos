#ifndef SHYLUBASKER_TYPES_HPP
#define SHYLUBASKER_TYPES_HPP

#include <exception>

#define BASKER_DEBUG

//MACRO TURN ON FUCNTIONS
#define BASKER_KOKKOS         //Use Kokkos
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

#define BASKER_PIVOT_TOL     .0001
#define BASKER_PIVOT_BIAS    1.001

//Error Codes
enum BASKER_ERROR_CODE 
{
  BASKER_ERROR_NOERROR,     //No error
  BASKER_ERROR_SINGULAR,    //Singular during factorization
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
#define BASKER_BTF_MAX_PERCENT  1.00
#define BASKER_BTF_LARGE        500  //Made smaller for unit test
#define BASKER_BTF_IMBALANCE     0.10
#define BASKER_BTF_SMALL         100
#define BASKER_BTF_NNZ_OVER      2.0 //Upped from 1.20
#define BASKER_BTF_PRUNE_SIZE    100

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
#define BASKER_FILL_USER           1.00
#define BASKER_FILL_LESTIMATE      1.50
#define BASKER_FILL_UESTIMATE      1.50
#define BASKER_FILL_LLOWERESTIMATE 2.00
#define BASKER_FILL_UUPPERESTIMATE 2.00
#define BASKER_FILL_LSEPESTIMATE   2.00
#define BASKER_FILL_USEPESTIMATE   2.00

//MACRO SYSTEM FUNCTIONS
#ifdef BASKER_DEBUG
#include <assert.h>
#define ASSERT(a)             assert(a)
#else
//#define ASSERT(a)             BASKER_NO_OP
#include <assert.h>
#define ASSERT(a)           assert(a)
#endif

#ifdef BASKER_DEBUG
#define BASKER_ASSERT(a,s)       \
  {                              \
    if(!(a))                     \
      {printf("\n\n%s \nLINE: %d \nFILE: %s\n\n", s, __LINE__, __FILE__);} \
    ASSERT(a);                   \
    assert(a);			 \
    if(!(a))                     \
      exit(0);                  \
  }
#else
#define BASKER_ASSERT(a,s)      \
  {                             \
    BASKER_NO_OP;               \
  }
#endif


//Note:  Should see if Kokkos has a fast memory cpy in place of for-loop
//MACRO ARRAY FUNCTIONS
#ifdef BASKER_KOKKOS
//Execution Space
#include <Kokkos_Core.hpp>
#define BASKER_EXE_SPACE     Kokkos::OpenMP
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
#define KOKKOS_NOINIT             Kokkos::ViewAllocateWithoutInitializing
#define INT_1DARRAY               Kokkos::View<BASKER_INT*,          BASKER_EXE_SPACE>
#define INT_2DARRAY               Kokkos::View<INT_1DARRAY*,         BASKER_EXE_SPACE> 
#define ENTRY_1DARRAY             Kokkos::View<BASKER_ENTRY*,        BASKER_EXE_SPACE>
#define ENTRY_2DARRAY             Kokkos::View<ENTRY_1DARRAY*,       BASKER_EXE_SPACE>  
#define BOOL_1DARRAY              Kokkos::View<BASKER_BOOL*,         BASKER_EXE_SPACE>
#define BOOL_2DARRAY              Kokkos::View<BOOL_1DARRAY*,        BASKER_EXE_SPACE>
#define MATRIX_1DARRAY            Kokkos::View<BASKER_MATRIX*,       BASKER_EXE_SPACE>
#define MATRIX_2DARRAY            Kokkos::View<MATRIX_1DARRAY*,      BASKER_EXE_SPACE>
#define MATRIX_VIEW_1DARRAY       Kokkos::View<BASKER_MATRIX_VIEW*,  BASKER_EXE_SPACE>
#define MATRIX_VIEW_2DARRAY       Kokkos::View<MATRIX_VIEW_1DARRAY*, BASKER_EXE_SPACE>
#define THREAD_1DARRAY            Kokkos::View<BASKER_THREAD*,       BASKER_EXE_SPACE>
#define THREAD_2DARRAY            Kokkos::View<THREAD_1DARRAY*,      BASKER_EXE_SPACE>

#define INT_1DARRAY_PAIRS        Kokkos::View<std::pair<Int,Int>*,  BASKER_EXE_SPACE>
//Macro Memory Calls
//MALLOC
#define MALLOC_INT_1DARRAY_PAIRS(a,s)   \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC malloc_pairs_1d: size to alloc > 0 fails"); \
    a = INT_1DARRAY_PAIRS(KOKKOS_NOINIT("pairs_1d"),s); \
    if(a.data() == NULL)           \
      throw std::bad_alloc();	   \
  }
#define MALLOC_INT_1DARRAY(a,s)   \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC int_1d: size to alloc > 0 fails"); \
    a = INT_1DARRAY(KOKKOS_NOINIT("int_1d"),s); \
    if(a.data() == NULL)           \
      throw std::bad_alloc();	   \
  }
#define MALLOC_INT_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0,"BASKER ASSERT MALLOC int_2d: size to alloc > 0 fails"); \
    a = INT_2DARRAY("int_2d",s); \
    if(a.data() == NULL)         \
      throw std::bad_alloc();    \
  }
#define MALLOC_ENTRY_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC entry_1d: size to alloc > 0 fails"); \
    a = ENTRY_1DARRAY(KOKKOS_NOINIT("entry_1d"),s); \
    if(a.data() == NULL)           \
      throw std::bad_alloc();      \
  }
#define MALLOC_ENTRY_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC entry_2d: size to alloc > 0 fails"); \
    a = ENTRY_2DARRAY("entry_2d",s); \
    if(a.data() == NULL)             \
      throw std::bad_alloc();        \
  }
#define MALLOC_BOOL_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC bool_1d: size to alloc > 0 fails"); \
    a = BOOL_1DARRAY(KOKKOS_NOINIT("bool_1d"), s); \
    if(a.data() == NULL)           \
      throw std::bad_alloc();      \
  }
#define MALLOC_BOOL_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC bool_2d: size to alloc > 0 fails"); \
    a = BOOL_2DARRAY("bool_2d", s); \
    if(a.data() == NULL)            \
      throw std::bad_alloc();       \
  }
#define MALLOC_MATRIX_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC matrix_1d: size to alloc > 0 fails"); \
    a = MATRIX_1DARRAY("matrix_1d",s); \
    if(a.data() == NULL)              \
      throw std::bad_alloc();         \
  }
#define MALLOC_MATRIX_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC matrix_2d: size to alloc > 0 fails"); \
    a = MATRIX_2DARRAY("matrix_2d",s); \
    if(a.data() == NULL)              \
      throw std::bad_alloc();         \
  }
#define MALLOC_MATRIX_VIEW_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC matrix_view_1d: size to alloc > 0 fails"); \
    a = MATRIX_VIEW_1DARRAY("matrix_view_1d",s); \
    if(a.data() == NULL)                 \
      throw std::bad_alloc();            \
  }
#define MALLOC_MATRIX_VIEW_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC matrix_view_2d: size to alloc > 0 fails"); \
    a = MATRIX_VIEW_2DARRAY("matrix_view_2d",s); \
    if(a.data() == NULL)                 \
      throw std::bad_alloc();            \
  }
#define MALLOC_THREAD_1DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC thread_1d: size to alloc > 0 fails"); \
    a = THREAD_1DARRAY("thread_1d",s); \
    if(a.data() == NULL)              \
      throw std::bad_alloc();         \
  }
#define MALLOC_THREAD_2DARRAY(a,s) \
  { \
    BASKER_ASSERT(s>0, "BASKER ASSERT MALLOC thread_2d: size to alloc > 0 fails"); \
    a = THREAD_2DARRAY("thread_2d",s); \
    if(a.data() == NULL)             \
      throw std::bakd_alloc();       \
  }
//RESIZE (with copy)
#define RESIZE_1DARRAY(a,os,s)           Kokkos::resize(a,s)
#define RESIZE_2DARRAY(a,os1,os2,s1,s2)   Kokkos::resize(a,s1,s2)
#define RESIZE_INT_1DARRAY(a,os,s)       RESIZE_1DARRAY(a,os,s)
#define RESIZE_ENTRY_1DARRAY(a,os,s)     RESIZE_1DARRAY(a,os,s)
//REALLOC (no copy)
#define REALLOC_1DARRAY(a,os,s)          Kokkos::realloc(a,s)
#define REALLOC_2DRRAAY(a,os1,os2,s1,s2) Kokkos::realloc(a,s1,s2)
#define REALLOC_INT_1DARRAY(a,os,s)      REALLOC_1DARRAY(a,os,s)
#define REALLOC_ENTRY_1DARRAY(a,os,s)    REALLOC_1DARRAY(a,os,s)
//Set values
#define SET_INT_1DARRAY(a, b, s)    	\
  {						\
  MALLOC_INT_1DARRAY(a,s);			     \
  for(BASKER_INT MACRO_I = 0; MACRO_I< s; MACRO_I++)	\
    a[MACRO_I] = b[MACRO_I];				\
  }
       
#define SET_ENTRY_1DARRAY(a, b,s)		\
  {						\
  MALLOC_ENTRY_1DARRAY(a,s);			    \
  for(BASKER_INT MACRO_I =0; MACRO_I <s; MACRO_I++)	\
    a[MACRO_I] = b[MACRO_I];				\
  }
  
#define SET_BOOL_1DARRAY(a, b, s)		\
  {						\
  MALLOC_BOOL_1DARRY(a,s);			    \
  for(BASKER_INT MACRO_I =0; MACRO_I <s; MACRO_I++)	\
    a[MACRO_I] = b[MACRO_I];				\
  }

#define FREE(a)                        BASKER_NO_OP

#define FREE_INT_1DARRAY_PAIRS(a)      \
  { \
    a = INT_1DARRAY_PAIRS(); \
  }

#define FREE_INT_1DARRAY(a)      \
  { \
    a = INT_1DARRAY(); \
  }

#define FREE_INT_2DARRAY(a,n)                    \
  { \
    a = INT_2DARRAY(); \
  }

#define FREE_ENTRY_1DARRAY(a)    \
  { \
    a = ENTRY_1DARRAY(); \
  }

#define FREE_ENTRY_2DARRAY(a,n)                  \
  { \
    a = ENTRY_2DARRAY(); \
  }

#define FREE_BOOL_1DARRAY(a)    \
  { \
    a = BOOL_1DARRAY(); \
  }

#define FREE_BOOL_2DARRAY(a,n)                   \
  { \
    a = BOOL_2DARRAY(); \
  }

#define FREE_MATRIX_1DARRAY(a)  \
  { \
    a = MATRIX_1DARRAY(); \
  }

#define FREE_MATRIX_2DARRAY(a,n)                 \
  { \
    a = MATRIX_2DARRAY(); \
  }

#define FREE_MATRIX_VIEW_1DARRAY(a) \
  { \
    a = MATRIX_VIEW_1DARRAY(); \
  }

#define FREE_MATRIX_VIEW_2DARRAY(a,n)            \
  { \
    a = MATRIX_VIEW_2DARRAY(); \
  }

#define FREE_THREAD_1DARRAY(a) \
  { \
    a = THREAD_1DARRAY(); \
  }

#define FREE_THREAD_2DARRAY(a,n)                 \
  { \
    a = TRHEAD_2DARRAY(); \
  }

#else
//Execution Space
#define BASKER_EXE_SPACE     void*
//ReMacro Basker Classes
#define BASKER_SOLVER        Basker<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
#define BASKER_MATRIX        BaskerMatrix<BASKER_INT, BASKER_ENTRY, BASKER_EXE_SPACE>
#define BASKER_MATRIX_VIEW   BaskerMatrixView<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
#define BASKER_STATS         BaskerStats<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
//ReMacor Basker Structs
#define BASKER_TREE          basker_tree<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
#define BASKER_SYMBOLIC_TREE basker_symbolic_tree<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
#define BASKER_THREAD        basker_thread<BASKER_INT,BASKER_ENTRY,BASKER_EXE_SPACE>
//Array Types
#define INT_1DARRAY          BASKER_INT*
#define INT_2DARRAY          BASKER_INT**
#define ENTRY_1DARRAY        BASKER_ENTRY*
#define ENTRY_2DARRAY        BASKER_ENTRY**
#define BOOL_1DARRAY         BASKER_BOOL*
#define BOOL_2DARRAY         BASKER_BOOL**
#define MATRIX_1DARRAY       BASKER_MATRIX*
#define MATRIX_2DARRAY       BASKER_MATRIX**
#define MATRIX_VIEW_1DARRAY  BASKER_MATRIX_VIEW*
#define MATRIX_VIEW_2DARRAY  BASKER_MATRIX_VIEW**
#define THREAD_1DARRAY       BASKER_THREAD*
#define THREAD_2DARRAY       BASKER_THREAD**

//Macro Memory Calls
//Malloc
#define MALLOC_INT_1DARRAY(a,s)          a = new BASKER_INT         [s]
#define MALLOC_INT_2DARRAY(a,s)          a = new INT_1DARRAY        [s]
#define MALLOC_ENTRY_1DARRAY(a,s)        a = new BASKER_ENTRY       [s]
#define MALLOC_ENTRY_2DARRAY(a,s)        a = new ENTRY_1DARRAY      [s]
#define MALLOC_BOOL_1DARRAY(a,s)         a = new BASKER_BOOL        [s]
#define MALLOC_BOOL_2DARRAY(a,s)         a = new BOOL_1DARRAY       [s]
#define MALLOC_MATRIX_1DARRAY(a,s)       a = new BASKER_MATRIX      [s]
#define MALLOC_MATRIX_2DARRAY(a,s)       a = new MATRIX_1DARRAY     [s]
#define MALLOC_MATRIX_VIEW_1DARRAY(a,s)  a = new BASKER_MATRIX_VIEW [s]
#define MALLOC_MATRIX_VIEW_2DARRAY(a,s)  a = new MATRIX_VIEW_1DARRAY[s]
#define MALLOC_THREAD_1DARRAY(a,s)       a = new BASKER_THREAD      [s]
#define MALLOC_THREAD_2DARRAY(a,s)       a = new THREAD_1DARRAY     [s]
//Resize (copy old data) (come back and add)
#define RESIZE_1DARRAY(a,os,s)               BASKER_NO_OP
#define RESIZE_2DARRAY(a,os1,os2,s1,s2)      BASKER_NO_OP
#define RESIZE_INT_1DARRAY(a,os,s)           BASKER_NO_OP
#define RESIZE_ENTRY_1DARRAY(a,os,s)         BASKER_NO_OP
//Realloc (dont copy old data)
#define REALLOC_1DARRAY(a,os,s)              BASKER_NO_OP
#define REALLOC_2DARRAY(a,os1,os2,s1,s2)     BASKER_NO_OP
#define REALLOC_INT_1DARRAY(a,os,s)          BASKER_NO_OP
#define REALLOC_ENTRY_1DARRAY(a,os,s)        BASKER_NO_OP
//Set functions
#define SET_INT_1DARRAY(a,b,s)           a = b
#define SET_ENTRY_1DARRAY(a,b,s)         a = b
#define SET_ENTRY_1DARRAY(a,b,s)         a = b  
#define FREE(a)                    delete [] a

#define FREE_INT_1DARRAY(a)      \
  { \
    FREE(a); \
  }

#define FREE_INT_2DARRAY(a,s)                    \
  { \
    for(BASKER_INT MACRO_I = 0; MACRO_I < s; MACRO_I++) \
      FREE(a[MACRO_I]); \
    FREE(a); \
  }

#define FREE_ENTRY_1DARRAY(a)    \
  { \
    FREE(a); \
  }

#define FREE_ENTRY_2DARRAY(a,s)                  \
  { \
    for(BASKER_INT MACRO_I = 0; MACRO_I < s; MACRO_I++) \
        FREE(a[MACRO_I]); \
    FREE(a); \
  }

#define FREE_BOOL_1DARRAY(a)    \
  { \
    FREE(a); \
  }

#define FREE_BOOL_2DARRAY(a,n)    \
  { \
    for(BASKER_INT MACRO_I = 0; MACRO_I < s; MACRO_I++) \
      FREE(a[MACRO_I]); \
    FREE(a); \
  }

#define FREE_MATRIX_1DARRAY(a)  \
  { \
    FREE(a); \
  }

#define FREE_MATRIX_2DARRAY(a,s)  \
  { \
    for(BASKER_INT MACRO_I = 0; MACRO_I < s; MACRO_I++) \
      FREE(a[MARCO_I]); \
    FREE(a); \
  }

#define FREE_MATRIX_VIEW_1DARRAY(a) \
  { \
    FREE(a); \
  }

#define FREE_MATRIX_VIEW_2DARRAY(a,s)            \
  { \
    for(BASKER_INT MACRO_I = 0; MACRO_I < s; MACRO_I++) \
      FREE(a[MACRO_I]); \
    FREE(a); \
  }

#define FREE_THREAD_1DARRAY(a) \
  { \
    FREE(a);  \
  }

#define FREE_THREAD_2DARRAY(a,n)                 \
  { \
    for(BASKER_INT MACRO_I = 0; MACRO_I < s; MACRO_I++) \
      FREE(a[MACRO_I]); \
    FREE(a); \
  }

#endif //end ifdef BASKER_KOKKOS

//Inline command
#ifdef BASKER_KOKKOS
#define BASKER_INLINE   KOKKOS_INLINE_FUNCTION
#else
#define BASKER_INLINE    inline
#endif

#define BASKER_FINLINE  inline

//Time Macro
#ifdef BASKER_TIME
#ifdef BASKER_KOKKOS
#define BASKER_TIMER
#define BASKER_TIMER_FINE
#else
#define BASKER_OMP_TIME
#endif
#endif

#endif //end basker_types_hpp
