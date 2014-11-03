#ifndef BASKER_DEF_HPP
#define BASKER_DEF_HPP

#include "basker_decl.hpp"
#include "basker.hpp"

//#include <assert.h>
#include <iostream>
#include <stdio.h>

using namespace std;

//#define BASKER_DEBUG 1
//#undef UDEBUG

namespace basker{

  template<class Int>
  inline
 bool basker_ordinal_type_valid
  (
   Int a
   )
  {
    if(a == ((Int)-1)) return false;
    else return true;
  }



  template <class Int, class Entry>
  void basker_dfs
  ( 
   Int n,
   Int j, 
   Int Li [], 
   Int Lp [], 
   Int color [], 
   Int pattern [], /* o/p */
   Int *top,       /* o/p */ 
   Int k,
   Int pinv [],
   Int stack []
  )
  {
    //int head;
    int have_elements = 1;
    Int i, t, i1, head ; 
   Int start, end, done, *store ;

   store = stack + n ;
   head = 0;
   stack[head] = j;

   //while (head >= 0)
        
   while(basker_ordinal_type_valid<Int>(head))
   {
       j = stack[head] ;
#ifdef BASKER_DEBUG
       std::cout << "DFS: " << j << "COLOR: " << color[j] << std::endl;
#endif
       t = pinv [j] ;
       if (color[j] == 0)
       {
           /* Seeing this column for first time */
           color[j] = 1 ;
           start = Lp[t] ;
       }
       else
       {
	 ASSERT (color[j] == 1) ; /* color cannot be 2 when we are here */
           start = store[j];
       }
       done = 1;

       if ( t != -1 )
       {
           end = Lp[t+1] ;
           for ( i1 = start ; i1 < end ; i1++ )
           {
               i = Li[i1] ;
               if ( color[i] == 0 )
               {
                   stack[++head] = i;
                   store[j] = i1+1;
                   done = 0;
                   break;
               }
           }
           /*if (i1 == end)
               done = 1;*/
       }
       if (done)
       {
	 // std::cout << "done called " << std::endl;
           pattern[--*top] = j ;
           color[j] = 2 ;
	   head--;
       }
   }
#ifdef BASKER_DEBUG
   std::cout << "Out of DFS: " << j << std::endl;
#endif
}

/* ==================== basker function ============== */
  template <class Int, class Entry>
  Int basker
  (
   Int Ap [],
   Int Ai [],
   Entry Ax [],
   Int anrow,
   Int ancol,
   Int ws [],
   Entry X [],
   Int *Lp, 
   Int **Li_p,
   Entry **Lx_p,
   Int *Up,
   Int **Ui_p,
   Entry **Ux_p,
   Int *llnnz_p,
   Int *uunnz_p,
   Int *pinv
   )   
  {

    Int i, j, k;
    Int *tptr, *color, *pattern, *stack ;
    Int *Li, *Ui ;
    Int top, top1, maxindex, t, j1, j2, llnnz, uunnz ;
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop ;
    Int pp, p2, p ;
    Int newsize;
    Entry pivot, value, xj ;
    Entry absv, maxv ;
    Entry *Lx, *Ux ;

    llnnz = *llnnz_p;
    uunnz = *uunnz_p;
    Li = *Li_p ;
    Lx = *Lx_p ;
    Ui = *Ui_p ;
    Ux = *Ux_p ;

    tptr = ws ;

    color = ws ;
    tptr += ancol ;

    
    pattern = tptr ;
    tptr += anrow ;

    stack = tptr ;
    tptr += 2*anrow ;

    cu_ltop = 0 ;
    cu_utop = 0 ;
    top = ancol ;
    top1 = ancol ;
    lnnz = 0 ;
    unnz = 0 ;
    

    for (k = 0; k<ancol ; k++)
    {
        pinv[k] = -1 ;
    }

    
    /* for all the columns in A ... */
    for (k = 0; k<ancol ; k++)
    {
      //std::cout << "K: " << k << std::endl;
#ifdef BASKER_DEBUG
      std::cout << "k = " << k << std::endl;
#endif

        value = 0.0 ;
        pivot = 0.0 ; 
        maxindex = -1 ;
        j1 = 0 ;
        j2 = 0 ;
        lcnt = 0;
        ucnt = 0;

#ifdef BASKER_DEBUG
        ASSERT( top == ancol ) ;

        for ( i = 0 ; i < anrow ;i++)
        {
          ASSERT( X[i] == 0 ) ;
        }

        for ( i = 0 ; i < ancol ;i++)
        {
          ASSERT ( color[i] == 0 ) ;
        }
#endif
	//printf("here \n");
	//return 1;

        /* reachability for every non zero in Ak */
        for (i = Ap[k] ;  i < Ap[k+1] ; i++ )
        {
            j = Ai[i] ;
            X [j] = Ax[i] ;

            if ( color[j] == 0 )
            {
                /*BASKER(dfs) (j, Li, Lp, color, pattern, &top, k, pinv) ;*/
	      //std::cout << "Before dfs \n";
	      basker_dfs<Int, Entry> (anrow, j, Li, Lp, color, pattern, &top, k,
                        pinv, stack) ;
	      //std::cout << "After dfs \n";

	      
            }    
        }

       
	

        xnnz = ancol - top ;
        //assert ( xnnz <= anrow ) ;

        /* Lx = b where x will be the column k in in L and U */
        top1 = top ;
        for ( pp = 0 ; pp < xnnz ; pp++ )
        {
            j = pattern[top1++] ;
            color[j] = 0 ;
            t = pinv[j] ;

            if (basker_ordinal_type_valid<Int>(t))
            {
                xj = X [j] ;
                p2 = Lp [t+1] ; /* TBV */
                for (p = Lp [t]+1 ;  p < p2 ; p++)
                {
                    X [Li [p]] -= Lx [p] * xj ;
                }
            }
        }

        /* get the pivot */
        maxv = 0.0 ; 
        for ( i = top ; i < anrow ;i++)
        {
            j = pattern[i] ;
            t = pinv[j] ;
            value =  X[j] ;
            absv = ( value < 0.0 ? -value : value ) ;

            if (!basker_ordinal_type_valid<Int>(t))
            {
                lcnt++;
                if( absv > maxv)
                {
                    maxv = absv ;
                    pivot = value ;
                    maxindex = j ;
                }
            }
        }
        ucnt = anrow - top - lcnt  + 1;
	
	//cout << "Max Index" << maxindex << std::endl;
        if (!basker_ordinal_type_valid<Int>(maxindex) || pivot == 0 )
        {
	  std::cout << "Matrix is singular at index: " << maxindex << " pivot: " << pivot << std::endl;       
            return 1;
        }

        pinv[maxindex] = k ;
        if (maxindex != k)
        {
#ifdef BASKER_DEBUG
	  std::cout << "Permuting pivot: " << k << " for row: " << maxindex << std::endl;
	  //std::cout
#endif
     
        }

       if (lnnz + lcnt >= llnnz)
        {
            /* reallocate space for L */
            /* Note that there can be atmost anrow - top - 1 entries in L
             * from the for loop below as P[maxindex] != -1.
             * The index and value for maxindex entry is assigned outside
             * the loop. */
             newsize = llnnz * 1.1 + 2 * anrow + 1;
#ifdef BASKER_DEBUG
	     std::cout << "Out of memory -- Reallocating: OldSize: " << llnnz << "NewSize: " << newsize << std::endl;
	     std::cout << "MY NOTES:  Check reallocs in C++" << std::endl;
#endif
	     Li = (Int *)REALLOC(Li, newsize*sizeof(Int));
             if (!Li)
             {
	       std::cout << "Cannot Realloc Memory" << std::endl;
               return 1;
             }
             //Lx = REALLOC(Lx, newsize * sizeof(double));
	     Lx = (Entry *)REALLOC(Lx, newsize*sizeof(Entry));
             if (!Lx)
             {
	       std::cout << "Cannot Realloc Memory" << std::endl;
               return 1;
             }
             llnnz = newsize;
        }

        if (unnz + ucnt >= uunnz)
        {
            /* reallocate space for U */
             newsize = uunnz * 1.1 + 2 * anrow + 1;
#ifdef BASKER_DEBUG
	     std::cout << "Out of memory -- Reallocating U:  OldSize: " << uunnz << " NewSize " << newsize << std::endl;
#endif
             //PRINT(("Reallocating L oldsize=%d, newsize=%d\n", uunnz, newsize));
             Ui = (Int *) REALLOC(Ui, newsize * sizeof(Int));
             
	     //Ui = (Int *)realloc(Ui, newsize*sizeof(Int));
	     if (!Ui)
             {
	       std::cout << "Cannot allocate memory\n";
	       //printf("Cannot allocate memory\n");
                 return 1;
             }
	     
             Ux = (Entry *) REALLOC(Ux, newsize * sizeof(Entry));
	     //Ux = (Entry *)realloc(Ux, newsize*sizeof(Entry));
             if (!Ux)
             {
	       std::cout << "Cannot allocate memory\n";
	       //printf("Cannot allocate memory\n");
                 return 1;
             }
             uunnz = newsize;
        }

        /* L(k,k) = 1 */
        //assert(lnnz < llnnz);
        Li[lnnz] = maxindex ;
        Lx[lnnz] = 1.0 ;
        lnnz++;

        for ( i = top ; i < anrow ;i++ )
        {
            j = pattern[i] ;
            t = pinv[j] ;

            /*  chk for numerical cancellation */
            if ( X[j] != 0 )
            {
	      if ( basker_ordinal_type_valid<Int>(t) )
                {
                    if ( unnz >= uunnz )
                    {
		      std::cout << "basker: Insufficient memory for U" << std::endl;
		      //printf ("basker : Insufficient memory for U %d %d \n", unnz, uunnz); 
                        return 1;
                    }
                    /* ASSERT(unnz < uunnz ) ; */
                    Ui[unnz] = pinv[j] ;
                    Ux[unnz] = X[j] ;
                    unnz++ ;
                }
	      else if (!basker_ordinal_type_valid<Int>(t))
                {
                    if ( lnnz >= llnnz )
                    {
		      std::cout << "basker: Insufficient memroy for L" << std::endl;
		      //printf ("basker : Insufficient memory for L \n"); 
                        return 1;
                    }
		    // assert(lnnz < llnnz ) ;
            /*printf("I am assigning Li[%d]=%d  Lx[%d]=%g t=%d, j =%d\n", lnnz, j, lnnz, X[j]/pivot, t, j) ;*/
                    Li[lnnz] = j ;
                    Lx[lnnz] = X[j]/pivot ;
                    lnnz++;
                }
            }
            X[j] = 0 ;
        }

        xnnz = 0;
        top = ancol ;

        Lp[k] = cu_ltop ;
        Lp[k+1] = lnnz ;
        cu_ltop = lnnz ;

        Up[k] = cu_utop ;
        Up[k+1] = unnz ;
        cu_utop = unnz ;

    }

#ifdef BASKER_DEBUG
    for (k = 0; k<lnnz ; k++)
    {
        printf("L[%d]=%g", k, Lx[k]);
    }
    printf("\n");

    for (k = 0; k<lnnz ; k++)
    {
        printf("Li[%d]=%d", k, Li[k]);
    }
    printf("\n");
#endif
    /*repermute L*/
    for(i = 0; i < ancol; i++)
      {
	for(k = Lp[i]; k < Lp[i+1]; k++)
	  {
	    Li[k] =pinv[Li[k]]; 	    
	  }
      }





    *llnnz_p = llnnz;
    *uunnz_p = uunnz;
    *Li_p = Li ;
    *Lx_p = Lx ;
    *Ui_p = Ui ;
    *Ux_p = Ux ;

    return 0 ;
}
  }//end namespace
#endif
