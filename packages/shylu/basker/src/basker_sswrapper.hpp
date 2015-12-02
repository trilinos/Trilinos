#ifndef BASKER_SSWRAPER_HPP
#define BASKER_SSWRAPPER_HPP

#include "basker_types.hpp"

#ifdef HAVE_AMESOS
#include "amesos_btf_decl.h"
#endif


namespace BaskerNS
{

  template<class Int>
  class BaskerSSWrapper
  {
  public:
    static 
    inline
    int my_strong_component
    (
     Int           &n,
     Int           *col_ptr,
     Int           *row_idx,
     Int           &nblks,
     Int           *perm,
     Int           *perm_in,
     Int           *CC
     )
    {
      return -1;
    }//end strong_component




  }; //end BaskerSSWrapper template <Int>

  template <>
  class BaskerSSWrapper <int>
  {
  public:
    static 
    inline
    int my_strong_component 
    (
     int           &n,
     int           *col_ptr,
     int           *row_idx,
     int           &nblks,
     int           *perm,
     int           *perm_in,
     int           *CC
     )
    {
      //typedef long int  l_Int;
      typedef int l_Int;
      
      //l_Int p[M.nrow]; //output row_per
      l_Int p[n]; //output row_per
      //l_Int r[M.nrow+1]; //comp_tabs
      l_Int r[n+1]; //comp_tabs
      //We will want to add option to use q in the future
      
      l_Int work[n*4];
      //l_Int work[M.nrow*4];
      
      //printf("before amesos call \n");
      /*
        nblks = amesos_btf_strongcomp(M.ncol,&(M.col_ptr[0]),
        &(M.row_idx[0]), 
        &(perm_in[0]), p, r, work);
      */
      nblks = amesos_btf_strongcomp(n, col_ptr,
				    row_idx, 
				    perm_in, p, r, work);
      //printf("after amesos call \n");
      
      
      
#ifdef BASKER_DEBUG_ORDER_BTF
      
      printf("\nBTF perm: \n");
      //for(Int i=0; i <M.nrow; i++)
      for(Int i=0; i < n; i++)
        {
          printf("%d, ", p[i]);
        }
      
      printf("\n\nBTF tabs: <right> \n");
      for(Int i=0; i < nblks+1; i++)
        {
          printf("%d, ", r[i]);
        }
      printf("\n");
#endif
      
      //BASKER_ASSERT(M.nrow > 0, "M.nrow btf");
      BASKER_ASSERT(n > 0, "M.nrow btf");
      //MALLOC_INT_1DARRAY(perm,M.nrow);
      // MALLOC_INT_1DARRAY(perm, n);
      for(l_Int i = 0; i < n; i++)
        {
          perm[p[i]] = i;
        }
      BASKER_ASSERT((nblks+1) > 0, "nblks+1 btf");
      //MALLOC_INT_1DARRAY(CC, nblks+1);
      for(l_Int i = 0; i < nblks+1; i++)
        {
          CC[i] = r[i];
        }
      
      return 0;
    }
  }; //end BaskerSSWraper template <int>

  template <>
  class BaskerSSWrapper <long>
  {
  public:
    static 
    inline
    int my_strong_component 
    (   
     long           &n,
     long           *col_ptr,
     long          *row_idx,
     long           &nblks,
     long           *perm,
     long           *perm_in,
     long          *CC
        )
    {
      typedef long  l_Int;
      
      l_Int p[n]; //output row_per
      l_Int r[n+1]; //comp_tabs
      //We will want to add option to use q in the future
      
      l_Int work[n*4];
      
      //printf("before amesos call \n");
      /*
        nblks = amesos_btf_l_strongcomp(M.ncol,&(M.col_ptr[0]),
        &(M.row_idx[0]), 
        &(perm_in[0]), p, r, work);
      */
      nblks = amesos_btf_l_strongcomp(n,
                                      col_ptr,
                                      row_idx, 
                                      perm_in, p, r, work);
      //printf("after amesos call \n");
      

      
#ifdef BASKER_DEBUG_ORDER_BTF
      
      printf("\nBTF perm: \n");
    for(Int i=0; i <n; i++)
      {
	printf("%d, ", p[i]);
      }
    
    printf("\n\nBTF tabs: <right> \n");
    for(l_Int i=0; i < nblks+1; i++)
      {
	printf("%d, ", r[i]);
      }
    printf("\n");
#endif

    BASKER_ASSERT(n > 0, "M.nrow btf");
    //MALLOC_INT_1DARRAY(perm,M.nrow);
    for(l_Int i = 0; i < n; i++)
      {
	perm[p[i]] = i;
      }
    BASKER_ASSERT((nblks+1) > 0, "nblks+1 btf");
    //MALLOC_INT_1DARRAY(CC, nblks+1);
    for(l_Int i = 0; i < nblks+1; i++)
      {
	CC[i] = r[i];
      }

    return 0;
  }//strong_component<long int, Entry, Exe_Space>

  }; //end BaskerSSWrapper <long>


}//end BaskerNS
#endif
