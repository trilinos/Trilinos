// @HEADER
// *****************************************************************************
//                   Basker: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the Basker contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/*======================== basker ==============================================*/
/* Finds the LU factorization of a sparse A such that A=L*U.
 * Requires 2 workspaces. The integer workspace of size ancol+2*anrow and a
 * double workspace of size 2*anrow. A is expected in compressed column form.
 *
 * The output L and U are also in compressed column form. The pointers for both
 * L and U should be preallocated. The size of expected L and U should be passed
 * to the method as lnnz and unnz. basker will return an error as the size of
 * computed L and U exceeds the expected L and U.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "basker.h"



/* ==================== dfs function ============== */
/* Does a depth first search of a trapezoidal L that is partly computed and
 * partly equal to the identity. The result of the dfs is in the pattern.
 */

void BASKER(dfs)
(
   Int j,
   Int Li [],
   Int Lp [],
   Int color [],
   Int pattern [],
   Int *top,
   Int k,
   Int pinv []
)
{
   Int i, t, i1 ;
   Int start, end ;

   PRINT(("DFS : %d ***********************\n", j));
   color[j] = 1 ;
   t = pinv [j] ;

   if ( t != -1 )
   {
       start = Lp[t] ;
       end = Lp[t+1] ;
       for ( i1 = start ; i1 < end ; i1++ )
       {
           i = Li[i1] ;
           if ( color[i] == 0 )
           {
               BASKER(dfs)(i, Li, Lp, color, pattern, top, k, pinv) ;
           }
       }
   }
   pattern[--*top] = j ;
   color[j] = 2 ;
}

void BASKER(dfs_iter)
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
   Int i, t, head, i1 ;
   Int start, end, done, *store ;

   store = stack + n ;
   head = 0;
   stack[head] = j;

   while (head >= 0)
   {
       j = stack[head] ;
       PRINT(("DFS : %d ***********************\n", j));
       t = pinv [j] ;
       if (color[j] == 0)
       {
           /* Seeing this column for first time */
           color[j] = 1 ;
           start = Lp[t] ;
       }
       else
       {
           BASKERASSERT (color[j] == 1) ; /* color cannot be 2 when we are here */
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
           pattern[--*top] = j ;
           color[j] = 2 ;
           head--; /* pop j */
       }
   }
   PRINT(("Out of DFS : %d ***********************\n", j));

}

/* ==================== basker function ============== */
Int BASKER(basker)
(
   Int Ap [],
   Int Ai [],
   double Ax [],
   Int anrow,
   Int ancol,
   Int ws [],
   double X [],
   Int *Lp,
   Int **Li_p,
   double **Lx_p,
   Int *Up,
   Int **Ui_p,
   double **Ux_p,
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
    double pivot, value, xj ;
    double absv, maxv ;
    double *Lx, *Ux ;

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
        PRINT(("k = %d ****************** \n", k));
        value = 0.0 ;
        pivot = 0.0 ;
        maxindex = -1 ;
        j1 = 0 ;
        j2 = 0 ;
        lcnt = 0;
        ucnt = 0;

#ifdef DEBUG
        BASKERASSERT( top == ancol ) ;

        for ( i = 0 ; i < anrow ;i++)
        {
          BASKERASSERT( X[i] == 0 ) ;
        }

        for ( i = 0 ; i < ancol ;i++)
        {
          BASKERASSERT ( color[i] == 0 ) ;
        }
#endif

        /* reachability for every non zero in Ak */
        for (i = Ap[k] ;  i < Ap[k+1] ; i++ )
        {
            j = Ai[i] ;
            X [j] = Ax[i] ;

            if ( color[j] == 0 )
            {
                /*BASKER(dfs) (j, Li, Lp, color, pattern, &top, k, pinv) ;*/
                BASKER(dfs_iter) (anrow, j, Li, Lp, color, pattern, &top, k,
                        pinv, stack) ;
            }
        }


        xnnz = ancol - top ;
        BASKERASSERT ( xnnz <= anrow ) ;

        /* Lx = b where x will be the column k in in L and U */
        top1 = top ;
        for ( pp = 0 ; pp < xnnz ; pp++ )
        {
            j = pattern[top1++] ;
            color[j] = 0 ;
            t = pinv[j] ;

            if ( t != -1 )
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

            if ( t == -1 )
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

        if ( maxindex == -1 || pivot == 0 )
        {
            printf(" matrix is singular maxindex=%d pivot=%g\n", maxindex, pivot ) ;
            return 1;
        }

        pinv[maxindex] = k ;
        if (maxindex != k)
        {
            PRINT(("********* Permuting pivot %d as row %d\n", k, maxindex));
        }


        if (lnnz + lcnt >= llnnz)
        {
            /* reallocate space for L */
            /* Note that there can be atmost anrow - top - 1 entries in L
             * from the for loop below as P[maxindex] != -1.
             * The index and value for maxindex entry is assigned outside
             * the loop. */
             newsize = llnnz * 1.1 + 2 * anrow + 1;
             PRINT(("Reallocating L oldsize=%d, newsize=%d\n", llnnz, newsize));
             Li = BASKERREALLOC(Li, newsize * sizeof(Int));
             if (!Li)
             {
                 printf("Cannot allocate memory\n");
                 return 1;
             }
             Lx = BASKERREALLOC(Lx, newsize * sizeof(double));
             if (!Lx)
             {
                 printf("Cannot allocate memory\n");
                 return 1;
             }
             llnnz = newsize;
        }

        if (unnz + ucnt >= uunnz)
        {
            /* reallocate space for U */
             newsize = uunnz * 1.1 + 2 * anrow + 1;
             PRINT(("Reallocating L oldsize=%d, newsize=%d\n", uunnz, newsize));
             Ui = BASKERREALLOC(Ui, newsize * sizeof(Int));
             if (!Ui)
             {
                 printf("Cannot allocate memory\n");
                 return 1;
             }
             Ux = BASKERREALLOC(Ux, newsize * sizeof(double));
             if (!Ux)
             {
                 printf("Cannot allocate memory\n");
                 return 1;
             }
             uunnz = newsize;
        }

        /* L(k,k) = 1 */
        BASKERASSERT(lnnz < llnnz);
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
                if ( t  != -1 )
                {
                    if ( unnz >= uunnz )
                    {
                        printf ("basker : Insufficient memory for U %d %d \n", unnz, uunnz);
                        return 1;
                    }
                    /* BASKERASSERT(unnz < uunnz ) ; */
                    Ui[unnz] = pinv[j] ;
                    Ux[unnz] = X[j] ;
                    unnz++ ;
                }
                else if ( t == -1 )
                {
                    if ( lnnz >= llnnz )
                    {
                        printf ("basker : Insufficient memory for L \n");
                        return 1;
                    }
                    BASKERASSERT(lnnz < llnnz ) ;
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

#ifdef DEBUG
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

    *llnnz_p = llnnz;
    *uunnz_p = uunnz;
    *Li_p = Li ;
    *Lx_p = Lx ;
    *Ui_p = Ui ;
    *Ux_p = Ux ;

    return 0 ;
}
