// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_IKLU_Utils.hpp"

//-----------------------------------------------------------------
// Allocation and Destruction
//-----------------------------------------------------------------

/* allocate a sparse matrix (triplet form or compressed-column form) */
csr *csr_spalloc (int m, int n, int nzmax, int values, int triplet)
{
  csr *A = (csr*)calloc (1, sizeof (csr)) ;    /* allocate the csr struct */
  if (!A) return (NULL) ;                 /* out of memory */
  A->m = m ;                              /* define dimensions and nzmax */
  A->n = n ;
  A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
  A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.row */
  A->p = (int*)malloc (triplet ? CS_MAX(nzmax,1) * sizeof (int) : CS_MAX(m+1,1) * sizeof (int)) ;
  A->j = (int*)malloc (CS_MAX(nzmax,1) * sizeof (int)) ;
  A->x = values ? (double*)malloc (CS_MAX(nzmax,1) * sizeof (double)) : NULL ;
  return ((!A->p || !A->j || (values && !A->x)) ? csr_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
int csr_sprealloc (csr *A, int nzmax)
{
    int ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    nzmax = (nzmax <= 0) ? (A->p [A->m]) : nzmax ;
    A->j = (int*)csr_realloc (A->j, nzmax, sizeof (int), &oki) ;
    if (CS_TRIPLET (A)) A->p = (int*)csr_realloc (A->p, nzmax, sizeof (int), &okj) ;
    if (A->x) A->x = (double*)csr_realloc (A->x, nzmax, sizeof (double), &okx) ;
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* wrapper for realloc */
void *csr_realloc (void *p, int n, size_t size, int *ok)
{
    void *pnew ;
    pnew = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

/* free a sparse matrix */
csr *csr_spfree (csr *A)
{
  if (!A) return (NULL);     /* do nothing if A already NULL */
  if (A->p) free(A->p);
  if (A->j) free(A->j);
  if (A->x) free(A->x);
  if (A) free(A);
  return (NULL) ;      /* free the csr struct and return NULL */
}

/* free a symbolic factorization */
css *csr_sfree (css *S)
{
  if (!S) return (NULL) ;     /* do nothing if S already NULL */
  if (S->pinv) free(S->pinv);
  if (S->q) free(S->q);
  if (S->parent) free(S->parent);
  if (S->cp) free(S->cp);
  if (S->leftmost) free(S->leftmost);
  if (S) free(S);
  return (NULL) ;      /* free the css struct and return NULL */
}

csrn *csr_nfree (csrn *N)
{
  if (!N) return (NULL) ;     /* do nothing if N already NULL */
  csr_spfree (N->L) ;
  csr_spfree (N->U) ;
  if (N->pinv) free(N->pinv);
  if (N->perm) free(N->perm);
  if (N->B) free(N->B);
  if (N) free(N);
  return (NULL) ;      /* free the csn struct and return NULL */
}

/* free workspace and return a sparse matrix result */
csr *csr_done (csr *C, void *w, void *x, int ok)
{
  if (w) free(w);                       /* free workspace */
  if (x) free(x);
  return (ok ? C : csr_spfree (C)) ;  /* return result if OK, else free it */
}

/* free workspace and return a numeric factorization (Cholesky, LU, or QR) */
csrn *csr_ndone (csrn *N, csr *C, void *w, void *x, int ok)
{
  csr_spfree (C) ;                    /* free temporary matrix */
  if (w) free(w);                     /* free workspace */
  if (x) free(x);
  return (ok ? N : csr_nfree (N)) ;   /* return result if OK, else free it */
}

/* free workspace and return int array result */
int *csr_idone (int *p, csr *C, void *w, int ok)
{
  csr_spfree (C) ;                    /* free temporary matrix */
  if (w) free(w);                     /* free workspace */
  /* return result if OK, else free it */
  if (ok)
    return p;
  else {
    if (p) free(p);
    return NULL;
  }
}

//-----------------------------------------------------------------
// Basic CSR math routines
//-----------------------------------------------------------------

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
double csr_cumsum (int *p, int *c, int n)
{
  int i, nz = 0 ;
  double nz2 = 0 ;
  if (!p || !c) return (-1) ;     /* check inputs */
  for (i = 0 ; i < n ; i++)
    {
      p [i] = nz ;
      nz += c [i] ;
      nz2 += c [i] ;              /* also in double to avoid int overflow */
      c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
  p [n] = nz ;
  return (nz2) ;                  /* return sum (c [0..n-1]) */
}

/* x = x + alpha * B(i,:), where x is a dense vector and B(i,:) is sparse */
int csr_scatter (const csr *B, int i, double alpha, int *w, double *x, int mark,
		 csr *C, int nz)
{
  int j, p, *Bp, *Bj, *Cj ;
  double *Bx ;
  if (!CS_CSC (B) || !w || !CS_CSC (C)) return (-1) ;     /* check inputs */
  Bp = B->p ; Bj = B->j ; Bx = B->x ; Cj = C->j ;
  for (p = Bp [i] ; p < Bp [i+1] ; p++)
    {
      j = Bj [p] ;                            /* B(i,j) is nonzero */
      if (w [j] < mark)
        {
	  w [j] = mark ;                      /* j is new entry in row i */
	  Cj [nz++] = j ;                     /* add j to pattern of C(i,:) */
	  if (x) x [j] = alpha * Bx [p] ;     /* x(j) = alpha*B(i,j) */
        }
      else if (x) x [j] += alpha * Bx [p] ;   /* j exists in C(i,:) already */
    }
  return (nz) ;
}

/* C = alpha*A + beta*B */
csr *csr_add (const csr *A, const csr *B, double alpha, double beta)
{
  int p, j, nz = 0, anz, *Cp, *Cj, *Bp, m, n, bnz, *w, values ;
  double *x, *Bx, *Cx ;
  csr *C ;
  if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;         /* check inputs */
  if ( (A->m != B->m) || (A->n != B->n) ) return (NULL);
  m = A->m ; anz = A->p [m] ;
  n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [m] ;
  w = (int*)calloc (CS_MAX(n,1), sizeof (int)) ;                       /* get workspace */
  values = (A->x != NULL) && (Bx != NULL) ;
  x = values ? (double*)malloc (n * sizeof (double)) : NULL ;    /* get workspace */
  C = csr_spalloc (m, n, anz + bnz, values, 0) ;          /* allocate result*/
  if (!C || !w || (values && !x)) return (csr_done (C, w, x, 0)) ;
  Cp = C->p ; Cj = C->j ; Cx = C->x ;
  for (j = 0 ; j < n ; j++)
    {
      Cp [j] = nz ;                   /* row j of C starts here */
      nz = csr_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(j,:)*/
      nz = csr_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(j,:) */
      if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Cj [p]] ;
    }
  Cp [m] = nz ;                       /* finalize the last row of C */
  csr_sprealloc (C, 0) ;              /* remove extra space from C */
  return (csr_done (C, w, x, 1)) ;    /* success; free workspace, return C */
}

/* C = A' */
csr *csr_transpose (const csr *A, int values)
{
  int p, q, i, *Cp, *Cj, n, m, *Ap, *Aj, *w ;
  double *Cx, *Ax ;
  csr *C ;
  if (!CS_CSC (A)) return (NULL) ;    /* check inputs */
  m = A->m ; n = A->n ; Ap = A->p ; Aj = A->j ; Ax = A->x ;
  C = csr_spalloc (n, m, Ap [m], values && Ax, 0) ;      /* allocate result */
  w = (int*)calloc (CS_MAX(n,1), sizeof (int)) ;         /* get workspace */
  if (!C || !w) return (csr_done (C, w, NULL, 0)) ;      /* out of memory */
  Cp = C->p ; Cj = C->j ; Cx = C->x ;
  for (p = 0 ; p < Ap [m] ; p++) w [Aj [p]]++ ;          /* col counts */
  csr_cumsum (Cp, w, n) ;                                 /* col pointers */
  for (i = 0 ; i < m ; i++)
    {
      for (p = Ap [i] ; p < Ap [i+1] ; p++)
        {
	  Cj [q = w [Aj [p]]++] = i ; /* place A(i,j) as entry C(j,i) */
	  if (Cx) Cx [q] = Ax [p] ;
        }
    }
  return (csr_done (C, w, NULL, 1)) ; /* success; free w and return C */
}

/* C = A*B */
csr *csr_multiply (const csr *A, const csr *B)
{
  int p, j, nz = 0, anz, *Cp, *Cj, *Ap, m, k, n, bnz, *w, values, *Aj;
  double *x, *Ax, *Cx ;
  csr *C ;
  if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;      /* check inputs */
  k = A->n;
  if (k != B->m ) return(NULL);
  m = A->m ; anz = A->p [A->m] ;
  n = B->n ; Ap = A->p ; Aj = A->j ; Ax = A->x ; bnz = B->p [k] ;
  w = (int*)calloc (CS_MAX(n,1), sizeof (int)) ;            /* get workspace */
  values = (Ax != NULL) && (B->x != NULL) ;
  x = values ? (double*)malloc (n * sizeof (double)) : NULL ; /* get workspace */
  C = csr_spalloc (m, n, anz + bnz, values, 0) ;       /* allocate result */
  if (!C || !w || (values && !x)) return (csr_done (C, w, x, 0)) ;
  Cp = C->p ;
  for (j = 0 ; j < m ; j++)
    {
      if (nz + n > C->nzmax && !csr_sprealloc (C, 2*(C->nzmax)+m))
        {
	  return (csr_done (C, w, x, 0)) ;            /* out of memory */
        }
      Cj = C->j ; Cx = C->x ;         /* C->j and C->x may be reallocated */
      Cp [j] = nz ;                   /* row j of C starts here */
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
	  nz = csr_scatter (B, Aj [p], Ax ? Ax [p] : 1, w, x, j+1, C, nz) ;
        }
      if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Cj [p]] ;
    }
  Cp [m] = nz ;                       /* finalize the last column of C */
  csr_sprealloc (C, 0) ;              /* remove extra space from C */
  return (csr_done (C, w, x, 1)) ;    /* success; free workspace, return C */
}


//-----------------------------------------------------------------
// Symbolic analysis
//-----------------------------------------------------------------

/* symbolic ordering and analysis for LU */
css *csr_sqr (int order, const csr *A )
{
    int n, ok = 1;
    css *S ;
    if (!CS_CSC (A)) return (NULL) ;        /* check inputs */
    n = A->n ;
    S = (css*)calloc(1, sizeof (css)) ;       /* allocate result S */
    if (!S) return (NULL) ;                 /* out of memory */
    S->q = csr_amd (order, A) ;             /* fill-reducing ordering */
    if (!S->q)
    {
        printf(" csr_sqr error no permutation\n");
    }
    if (order && !S->q) return (csr_sfree (S)) ;

    /* LU factorization */
    S->unz = (double) CS_MIN(4*(A->p [n]) + n, n * n );
    S->lnz = S->unz ;                       /* guess nnz(L) and nnz(U) */
    return (ok ? S : csr_sfree (S)) ;        /* return result S */
}


/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
int csr_reach (csr *G, const csr *B, int k, int *xi, const int *pinv)
{
    int p, m, top, *Bp, *Bj, *Gp ;
    if (!CS_CSC (G) || !CS_CSC (B) || !xi) return (-1) ;    /* check inputs */
    m = G->m ; Bp = B->p ; Bj = B->j ; Gp = G->p ;
    top = m ;
    for (p = Bp [k] ; p < Bp [k+1] ; p++)
    {
        if (!CS_MARKED (Gp, Bj [p]))    /* start a dfs at unmarked node i */
        {
            top = csr_dfs (Bj [p], G, top, xi, xi+m, pinv) ;
        }
    }
    for (p = top ; p < m ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */
    return (top) ;
}

/* depth-first-search of the graph of a csr matrix, starting at node j */
int csr_dfs (int j, csr *G, int top, int *xi, int *pstack, const int *pinv)
{
    int i, p, p2, done, jnew, head = 0, *Gp, *Gj ;
    if (!CS_CSC (G) || !xi || !pstack) return (-1) ;    /* check inputs */
    Gp = G->p ; Gj = G->j ;
    xi [0] = j ;                /* initialize the recursion stack */
    while (head >= 0)
    {
        j = xi [head] ;         /* get j from the top of the recursion stack */
        jnew = pinv ? (pinv [j]) : j ;
        if (!CS_MARKED (Gp, j))
        {
            CS_MARK (Gp, j) ;       /* mark node j as visited */
            pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ;
        }
        done = 1 ;                  /* node j done if no unvisited neighbors */
        p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;
        for (p = pstack [head] ; p < p2 ; p++)  /* examine all neighbors of j */
        {
            i = Gj [p] ;            /* consider neighbor node i */
            if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */
            pstack [head] = p ;     /* pause depth-first search of node j */
            xi [++head] = i ;       /* start dfs at node i */
            done = 0 ;              /* node j is not done */
            break ;                 /* break, to start dfs (i) */
        }
        if (done)               /* depth-first search at node j is done */
        {
            head-- ;            /* remove j from the recursion stack */
            xi [--top] = j ;    /* and place in the output stack */
        }
    }
    return (top) ;
}

/* depth-first search and postorder of a tree rooted at node j */
int csr_tdfs (int j, int k, int *head, const int *next, int *post, int *stack)
{
  int i, p, top = 0 ;
  if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
  stack [0] = j ;                 /* place j on the stack */
  while (top >= 0)                /* while (stack is not empty) */
    {
      p = stack [top] ;           /* p = top of stack */
      i = head [p] ;              /* i = youngest child of p */
      if (i == -1)
        {
	  top-- ;                 /* p has no unordered children left */
	  post [k++] = p ;        /* node p is the kth postordered node */
        }
      else
        {
	  head [p] = next [i] ;   /* remove i from children of p */
	  stack [++top] = i ;     /* start dfs on child node i */
        }
    }
  return (k) ;
}

//-----------------------------------------------------------------
// LU factorization
//-----------------------------------------------------------------

/*
 * Given sparse A,  
 * [L,U,pinv]=lu(A, [q lnz unz]). lnz and unz can be guesses 
 * Hypotheses: m=n, Lj[Lp[i+1]-1]=i and Uj[Lp[i]]=i
 */
csrn *csr_lu (const csr *A, const css *S, double tol)
{
    csr *L, *U ;
    csrn *N ;
    double pivot, *Lx, *Ux, *x,  a, t;
    int *Lp, *Lj, *Up, *Uj, *pinv, *xi, *q, n, ipiv, k, top, p, i, row, lnz, unz;
    int debug = 0;

    if (!CS_CSC (A) ) printf(" error csrlu: A not csc\n");
    if (!CS_CSC (A) || !S) return (NULL) ;	    /* check inputs */
    n = A->n ;
    if (n != A->m) return (NULL) ;	    /* check inputs */
    q = S->q ; lnz = (int)S->lnz ; unz = (int)S->unz ;
    x = (double*)malloc(n * sizeof (double)) ;	    /* get double workspace */
    xi = (int*)malloc (2*n * sizeof (int)) ;	    /* get int workspace */
    N = (csrn*)calloc (1, sizeof (csrn)) ;		    /* allocate result */
    if (!(A) ) printf(" error csrlu: allocation of N failed\n");
    if (!x || !xi || !N) return (csr_ndone (N, NULL, xi, x, 0)) ;

    N->L = L = csr_spalloc (n, n, lnz, 1, 0) ;	    /* allocate result L */
    N->U = U = csr_spalloc (n, n, unz, 1, 0) ;	    /* allocate result U */
    N->pinv = pinv = (int*)malloc (n * sizeof (int)) ;  /* allocate result pinv */
    N->perm = (int*)malloc (n * sizeof (int)) ;         /* allocate result perm */
    if (!L || !U || !pinv) return (csr_ndone (N, NULL, xi, x, 0)) ;
    Lp = L->p ; Up = U->p ;
    for (i = 0 ; i < n ; i++) x [i] = 0 ;	    /* clear workspace */
    for (i = 0 ; i < n ; i++) pinv [i] = -1 ;	    /* no rows pivotal yet */
    for (k = 1 ; k <= n ; k++) Up [k] = 0 ;	    /* no rows of U yet */
    for (k = 1 ; k <= n ; k++) Lp [k] = 0 ;	    /* no rows of L yet either */
    lnz = unz = 0 ;
    if( debug )
    {
       printf ("A:\n") ; csr_print (A, 0) ;  
    }
    for (k = 0 ; k < n ; k++)	    /* compute L(:,k) and U(:,k) */
    {
	/* --- Triangular solve --------------------------------------------- */
	Lp [k] = lnz ;		    /* L(:,k) starts here */
	Up [k] = unz ;		    /* U(:,k) starts here */
	if ((lnz + n > L->nzmax && !csr_sprealloc (L, 2*L->nzmax + n)) ||
	    (unz + n > U->nzmax && !csr_sprealloc (U, 2*U->nzmax + n)))
	{
	    return (csr_ndone (N, NULL, xi, x, 0)) ;
	}
        Lj = L->j ; Lx = L->x ; Uj = U->j ; Ux = U->x ;
        row = q ? (q [k]) : k ;
        if( debug > 1 )
        {
            printf("--------------------------------\n");
            printf(" %d spsolve row=%d \n",k, row);
            printf(" pinv = %d %d %d %d \n", pinv[0], pinv[1], pinv[2], pinv[3]);
        }
        top = csr_spsolve (U, A, row, xi, x, pinv, 1) ; /* x = U\A(row,:) */
        if( debug > 1 ) printf("top=%d  x = %g %g %g %g \n", top,x[0],x[1],x[2],x[3]);
        /* --- Find pivot --------------------------------------------------- */
        ipiv = -1 ;
        a = -1 ;
        for (p = top ; p < n ; p++)
        {
            i = xi [p] ;            /* x(i) is nonzero */
            if (pinv [i] < 0)       /* row i is not yet pivotal */
            {
                if ((t = fabs (x [i])) > a)
                {
                    a = t ;         /* largest pivot candidate so far */
                    ipiv = i ;
                }
            }
            else                    /* x(i) is the entry L(pinv[i],k) */
            {
                Lj [lnz] = pinv [i] ;
                Lx [lnz++] = x [i] ;
            }
        }
        if (ipiv == -1 || a <= 0) return (csr_ndone (N, NULL, xi, x, 0)) ;
        if (pinv [row] < 0 && fabs (x [row]) >= a*tol) ipiv = row;
        pivot = x [ipiv] ;          /* the chosen pivot */
        Lj [lnz] = k ;              /* last entry in L(:,k) is L(k,k) */
        Lx [lnz++] = pivot ;
        if( debug > 1 ) { printf ("L:") ; csr_print (L, 0) ;  }

        /* --- Divide by pivot ---------------------------------------------- */
        pinv [ipiv] = k ;           /* ipiv is the kth pivot row */
        Uj [unz] = ipiv ;           /* first entry in U(:,k) is U(k,k) = 1 */
        Ux [unz++] = 1 ;
        for (p = top ; p < n ; p++) /* U(k+1:n,k) = x / pivot */
        {
            i = xi [p] ;
            if (pinv [i] < 0)       /* x(i) is an entry in U(:,k) */
            {
                Uj [unz] = i ;      /* save unpermuted row in U */
                Ux [unz++] = x [i] / pivot ;    /* scale pivot row */
            }
            x [i] = 0 ;             /* x [0..n-1] = 0 for next k */
        }
        if( debug > 1 )
        {
            printf ("U:") ; csr_print (U, 0) ;  
            printf("------------------------------------\n");
        }
    }
    /* --- Finalize U and L ------------------------------------------------- */
    Lp [n] = lnz ;
    if( debug ) { printf ("L:") ; csr_print (L, 0) ;  }
    Up [n] = unz ;
    Uj = U->j ;                     /* fix column indices of U for final pinv */
    for (p = 0 ; p < unz ; p++) Uj [p] = pinv [Uj [p]] ;

    csr_sprealloc (L, 0) ;	    /* remove extra space from L and U */
    csr_sprealloc (U, 0) ;
    if( debug ) { printf ("U:") ; csr_print (U, 0) ;  }
    return (csr_ndone (N, NULL, xi, x, 1)) ;	/* success */
}

//-----------------------------------------------------------------
// Triangular solves
//-----------------------------------------------------------------

/* solve xG=b(k,:), where G is either upper (up=1) or lower (up=0) triangular */
int csr_spsolve (csr *G, const csr *B, int k, int *xi, double *x, const int *pinv, int up)
{
  int i, I, p, q, px, top, n, *Gp, *Gj, *Bp, *Bj ;
  int debug = 0;
  double *Gx, *Bx ;
  if (!CS_CSC (G) || !CS_CSC (B) || !xi || !x) return (-1) ;
  Gp = G->p ; Gj = G->j ; Gx = G->x ; n = G->n ;
  Bp = B->p ; Bj = B->j ; Bx = B->x ;
  top = csr_reach (G, B, k, xi, pinv) ;       /* xi[top..n-1]=Reach(B(:,k)) */
  for (p = top ; p < n ; p++) x [xi [p]] = 0 ;    /* clear x */
  for (p = Bp [k] ; p < Bp [k+1] ; p++) x [Bj [p]] = Bx [p] ; /* scatter B */
  if( debug )
    printf("solve k=%d   x= %g %g %g %g \n", k, x[0], x[1], x[2], x[3]);
  if( debug )
    printf("top=%d   xi= %d %d %d %d \n", top , xi[0], xi[1], xi[2], xi[3]);
  for (px = top ; px < n ; px++)
    {
      i = xi [px] ;                               /* x(i) is nonzero */
      I = pinv ? (pinv [i]) : i ;                 /* i maps to col I of G */
      if (I < 0) continue ;                       /* row I is empty */
      /* dmd */
      x [i] /= Gx [up ? (Gp [I]) : (Gp [I+1]-1)] ;/* x(i) /= G(i,i) */
      p = up ? (Gp [I]+1) : (Gp [I]) ;            /* up: L(i,i) 1st entry */
      q = up ? (Gp [I+1]) : (Gp [I+1]-1) ;        /* up: U(i,i) last entry */
      for ( ; p < q ; p++)
        {
	  if( debug )
	    printf("%d %d solve %d %g %g \n", px, i ,p,  Gx [p] , x [Gj [p]]  );
	  
	  x [Gj[p]] -= Gx [p] * x [i] ;           /* x(i) -= G(i,j) * x(j) */
        }
      if( debug )
	printf(" x= %g %g %g %g \n", x[0], x[1], x[2], x[3]);
    }
  return (top) ;                                  /* return top of stack */
}

//-----------------------------------------------------------------
// AMD routine
//-----------------------------------------------------------------

/*
 * amd for csr matrices , and  order =1
 * Hypothesis:  m = n
 */
/* clear w */
static int csr_wclear (int mark, int lemax, int *w, int n)
{
    int k ;
    if (mark < 2 || (mark + lemax < 0))
    {
        for (k = 0 ; k < n ; k++) if (w [k] != 0) w [k] = 1 ;
        mark = 2 ;
    }
    return (mark) ;     /* at this point, w [0..n-1] < mark holds */
}

/* keep off-diagonal entries; drop diagonal entries */
static int csr_diag (int i, int j, double aij, void *other) { return (i != j) ; }

/* p = amd(A+A') if symmetric is true, or amd(A'A) otherwise */
int *csr_amd (int order, const csr *A)  /* order 0:natural, 1:Chol, 2:LU, 3:QR */
{
  csr *C, *A2, *AT ;
  int *Cp, *Cj, *last, *W, *len, *nv, *next, *P, *head, *elen, *degree, *w,
    *hhead, *ATp, *ATj, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
    k2, k3, jlast, ln, dense, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
    ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, n, m, t ;
  unsigned int h ;
  /* --- Construct matrix C ----------------------------------------------- */
  if (!CS_CSC (A) || order <= 0 || order > 3) return (NULL) ; /* check */
  AT = csr_transpose (A, 0) ;             /* compute A' */
  if (!AT) return (NULL) ;
  m = A->m ; n = A->n ;
  if ( n != m) return(NULL); /* For rectangular matrices, use csr_amd */
  dense = (int)CS_MAX (16, 10 * sqrt ((double) n)) ;   /* find dense threshold */
  dense = CS_MIN (n-2, dense) ;
  if (order == 1 && n == m)
    {
      C = csr_add (A, AT, 0, 0) ;         /* C = A+A' */
    }
  else if (order == 2)
    {
      ATp = AT->p ;                       /* drop dense columns from AT */
      ATj = AT->j ;
      for (p2 = 0, j = 0 ; j < m ; j++)
        {
	  p = ATp [j] ;                   /* column j of AT starts here */
	  ATp [j] = p2 ;                  /* new column j starts here */
	  if (ATp [j+1] - p > dense) continue ;   /* skip dense col j */
	  for ( ; p < ATp [j+1] ; p++) ATj [p2++] = ATj [p] ;
        }
      ATp [m] = p2 ;                      /* finalize AT */
      A2 = csr_transpose (AT, 0) ;        /* A2 = AT' */
      C = A2 ? csr_multiply (AT, A2) : NULL ; /* C=A'*A with no dense rows */
      csr_spfree (A2) ;
    }
  else
    {
      C = csr_multiply (AT, A) ;          /* C=A'*A */
    }
  csr_spfree (AT) ;
  if (!C) return (NULL) ;
  csr_fkeep (C, &csr_diag, NULL) ;         /* drop diagonal entries */
  Cp = C->p ;
  cnz = Cp [n] ;
  P = (int*)malloc (CS_MAX(n+1,1) * sizeof (int)) ;     /* allocate result */
  W = (int*)malloc (CS_MAX(8*(n+1),1) * sizeof (int)) ; /* get workspace */
  t = cnz + cnz/5 + 2*n ;                 /* add elbow room to C */
  if (!P || !W || !csr_sprealloc (C, t)) return (csr_idone (P, C, W, 0)) ;
  len  = W           ; nv     = W +   (n+1) ; next   = W + 2*(n+1) ;
  head = W + 3*(n+1) ; elen   = W + 4*(n+1) ; degree = W + 5*(n+1) ;
  w    = W + 6*(n+1) ; hhead  = W + 7*(n+1) ;
  last = P ;                              /* use P as workspace for last */
  /* --- Initialize quotient graph ---------------------------------------- */
  for (k = 0 ; k < n ; k++) len [k] = Cp [k+1] - Cp [k] ;
  len [n] = 0 ;
  nzmax = C->nzmax ;
  Cj = C->j ;
  for (i = 0 ; i <= n ; i++)
    {
      head [i] = -1 ;                     /* degree list i is empty */
      last [i] = -1 ;
      next [i] = -1 ;
      hhead [i] = -1 ;                    /* hash list i is empty */
      nv [i] = 1 ;                        /* node i is just one node */
      w [i] = 1 ;                         /* node i is alive */
      elen [i] = 0 ;                      /* Ek of node i is empty */
      degree [i] = len [i] ;              /* degree of node i */
    }
  mark = csr_wclear (0, 0, w, n) ;         /* clear w */
  elen [n] = -2 ;                         /* n is a dead element */
  Cp [n] = -1 ;                           /* n is a root of assembly tree */
  w [n] = 0 ;                             /* n is a dead element */
  /* --- Initialize degree lists ------------------------------------------ */
  for (i = 0 ; i < n ; i++)
    {
      d = degree [i] ;
      if (d == 0)                         /* node i is empty */
        {
	  elen [i] = -2 ;                 /* element i is dead */
	  nel++ ;
	  Cp [i] = -1 ;                   /* i is a root of assemby tree */
	  w [i] = 0 ;
        }
      else if (d > dense)                 /* node i is dense */
        {
	  nv [i] = 0 ;                    /* absorb i into element n */
	  elen [i] = -1 ;                 /* node i is dead */
	  nel++ ;
	  Cp [i] = CS_FLIP (n) ;
	  nv [n]++ ;
        }
      else
        {
	  if (head [d] != -1) last [head [d]] = i ;
	  next [i] = head [d] ;           /* put node i in degree list d */
	  head [d] = i ;
        }
    }
  while (nel < n)                         /* while (selecting pivots) do */
    {
      /* --- Select node of minimum approximate degree -------------------- */
      for (k = -1 ; mindeg < n && (k = head [mindeg]) == -1 ; mindeg++) ;
      if (next [k] != -1) last [next [k]] = -1 ;
      head [mindeg] = next [k] ;          /* remove k from degree list */
      elenk = elen [k] ;                  /* elenk = |Ek| */
      nvk = nv [k] ;                      /* # of nodes k represents */
      nel += nvk ;                        /* nv[k] nodes of A eliminated */
      /* --- Garbage collection ------------------------------------------- */
      if (elenk > 0 && cnz + mindeg >= nzmax)
        {
	  for (j = 0 ; j < n ; j++)
            {
	      if ((p = Cp [j]) >= 0)      /* j is a live node or element */
                {
		  Cp [j] = Cj [p] ;       /* save first entry of object */
		  Cj [p] = CS_FLIP (j) ;  /* first entry is now CS_FLIP(j) */
                }
            }
	  for (q = 0, p = 0 ; p < cnz ; ) /* scan all of memory */
            {
	      if ((j = CS_FLIP (Cj [p++])) >= 0)  /* found object j */
                {
		  Cj [q] = Cp [j] ;       /* restore first entry of object */
		  Cp [j] = q++ ;          /* new pointer to object j */
		  for (k3 = 0 ; k3 < len [j]-1 ; k3++) Cj [q++] = Cj [p++] ;
                }
            }
	  cnz = q ;                       /* Cj [cnz...nzmax-1] now free */
        }
      /* --- Construct new element ---------------------------------------- */
      dk = 0 ;
      nv [k] = -nvk ;                     /* flag k as in Lk */
      p = Cp [k] ;
      pk1 = (elenk == 0) ? p : cnz ;      /* do in place if elen[k] == 0 */
      pk2 = pk1 ;
      for (k1 = 1 ; k1 <= elenk + 1 ; k1++)
        {
	  if (k1 > elenk)
            {
	      e = k ;                     /* search the nodes in k */
	      pj = p ;                    /* list of nodes starts at Cj[pj]*/
	      ln = len [k] - elenk ;      /* length of list of nodes in k */
            }
	  else
            {
	      e = Cj [p++] ;              /* search the nodes in e */
	      pj = Cp [e] ;
	      ln = len [e] ;              /* length of list of nodes in e */
            }
	  for (k2 = 1 ; k2 <= ln ; k2++)
            {
	      i = Cj [pj++] ;
	      if ((nvi = nv [i]) <= 0) continue ; /* node i dead, or seen */
	      dk += nvi ;                 /* degree[Lk] += size of node i */
	      nv [i] = -nvi ;             /* negate nv[i] to denote i in Lk*/
	      Cj [pk2++] = i ;            /* place i in Lk */
	      if (next [i] != -1) last [next [i]] = last [i] ;
	      if (last [i] != -1)         /* remove i from degree list */
                {
		  next [last [i]] = next [i] ;
                }
	      else
                {
		  head [degree [i]] = next [i] ;
                }
            }
	  if (e != k)
            {
	      Cp [e] = CS_FLIP (k) ;      /* absorb e into k */
	      w [e] = 0 ;                 /* e is now a dead element */
            }
        }
      if (elenk != 0) cnz = pk2 ;         /* Cj [cnz...nzmax] is free */
      degree [k] = dk ;                   /* external degree of k - |Lk\i| */
      Cp [k] = pk1 ;                      /* element k is in Cj[pk1..pk2-1] */
      len [k] = pk2 - pk1 ;
      elen [k] = -2 ;                     /* k is now an element */
      /* --- Find set differences ----------------------------------------- */
      mark = csr_wclear (mark, lemax, w, n) ;  /* clear w if necessary */
      for (pk = pk1 ; pk < pk2 ; pk++)    /* scan 1: find |Le\Lk| */
        {
	  i = Cj [pk] ;
	  if ((eln = elen [i]) <= 0) continue ;/* skip if elen[i] empty */
	  nvi = -nv [i] ;                      /* nv [i] was negated */
	  wnvi = mark - nvi ;
	  for (p = Cp [i] ; p <= Cp [i] + eln - 1 ; p++)  /* scan Ei */
            {
	      e = Cj [p] ;
	      if (w [e] >= mark)
                {
		  w [e] -= nvi ;          /* decrement |Le\Lk| */
                }
	      else if (w [e] != 0)        /* ensure e is a live element */
                {
		  w [e] = degree [e] + wnvi ; /* 1st time e seen in scan 1 */
                }
            }
        }
      /* --- Degree update ------------------------------------------------ */
      for (pk = pk1 ; pk < pk2 ; pk++)    /* scan2: degree update */
        {
	  i = Cj [pk] ;                   /* consider node i in Lk */
	  p1 = Cp [i] ;
	  p2 = p1 + elen [i] - 1 ;
	  pn = p1 ;
	  for (h = 0, d = 0, p = p1 ; p <= p2 ; p++)    /* scan Ei */
            {
	      e = Cj [p] ;
	      if (w [e] != 0)             /* e is an unabsorbed element */
                {
		  dext = w [e] - mark ;   /* dext = |Le\Lk| */
		  if (dext > 0)
                    {
		      d += dext ;         /* sum up the set differences */
		      Cj [pn++] = e ;     /* keep e in Ei */
		      h += e ;            /* compute the hash of node i */
                    }
		  else
                    {
		      Cp [e] = CS_FLIP (k) ;  /* aggressive absorb. e->k */
		      w [e] = 0 ;             /* e is a dead element */
                    }
                }
            }
	  elen [i] = pn - p1 + 1 ;        /* elen[i] = |Ei| */
	  p3 = pn ;
	  p4 = p1 + len [i] ;
	  for (p = p2 + 1 ; p < p4 ; p++) /* prune edges in Ai */
            {
	      j = Cj [p] ;
	      if ((nvj = nv [j]) <= 0) continue ; /* node j dead or in Lk */
	      d += nvj ;                  /* degree(i) += |j| */
	      Cj [pn++] = j ;             /* place j in node list of i */
	      h += j ;                    /* compute hash for node i */
            }
	  if (d == 0)                     /* check for mass elimination */
            {
	      Cp [i] = CS_FLIP (k) ;      /* absorb i into k */
	      nvi = -nv [i] ;
	      dk -= nvi ;                 /* |Lk| -= |i| */
	      nvk += nvi ;                /* |k| += nv[i] */
	      nel += nvi ;
	      nv [i] = 0 ;
	      elen [i] = -1 ;             /* node i is dead */
            }
	  else
            {
	      degree [i] = CS_MIN (degree [i], d) ;   /* update degree(i) */
	      Cj [pn] = Cj [p3] ;         /* move first node to end */
	      Cj [p3] = Cj [p1] ;         /* move 1st el. to end of Ei */
	      Cj [p1] = k ;               /* add k as 1st element in of Ei */
	      
	      len [i] = pn - p1 + 1 ;     /* new len of adj. list of node i */
	      h %= n ;                    /* finalize hash of i */
	      next [i] = hhead [h] ;      /* place i in hash bucket */
	      hhead [h] = i ;
	      last [i] = h ;              /* save hash of i in last[i] */
            }
        }                                   /* scan2 is done */
      degree [k] = dk ;                   /* finalize |Lk| */
      lemax = CS_MAX (lemax, dk) ;
      mark = csr_wclear (mark+lemax, lemax, w, n) ;    /* clear w */
      /* --- Supernode detection ------------------------------------------ */
      for (pk = pk1 ; pk < pk2 ; pk++)
        {
	  i = Cj [pk] ;
	  if (nv [i] >= 0) continue ;         /* skip if i is dead */
	  h = last [i] ;                      /* scan hash bucket of node i */
	  i = hhead [h] ;
	  hhead [h] = -1 ;                    /* hash bucket will be empty */
	  for ( ; i != -1 && next [i] != -1 ; i = next [i], mark++)
            {
	      ln = len [i] ;
	      eln = elen [i] ;
	      for (p = Cp [i]+1 ; p <= Cp [i] + ln-1 ; p++) w [Cj [p]] = mark;
	      jlast = i ;
	      for (j = next [i] ; j != -1 ; ) /* compare i with all j */
                {
		  ok = (len [j] == ln) && (elen [j] == eln) ;
		  for (p = Cp [j] + 1 ; ok && p <= Cp [j] + ln - 1 ; p++)
		    {
		      if (w [Cj [p]] != mark) ok = 0 ;    /* compare i and j*/
                    }
		  if (ok)                     /* i and j are identical */
                    {
		      Cp [j] = CS_FLIP (i) ;  /* absorb j into i */
		      nv [i] += nv [j] ;
		      nv [j] = 0 ;
		      elen [j] = -1 ;         /* node j is dead */
		      j = next [j] ;          /* delete j from hash bucket */
		      next [jlast] = j ;
                    }
		  else
                    {
		      jlast = j ;             /* j and i are different */
		      j = next [j] ;
                    }
                }
            }
        }
      /* --- Finalize new element------------------------------------------ */
      for (p = pk1, pk = pk1 ; pk < pk2 ; pk++)   /* finalize Lk */
        {
	  i = Cj [pk] ;
	  if ((nvi = -nv [i]) <= 0) continue ;/* skip if i is dead */
	  nv [i] = nvi ;                      /* restore nv[i] */
	  d = degree [i] + dk - nvi ;         /* compute external degree(i) */
	  d = CS_MIN (d, n - nel - nvi) ;
	  if (head [d] != -1) last [head [d]] = i ;
	  next [i] = head [d] ;               /* put i back in degree list */
	  last [i] = -1 ;
	  head [d] = i ;
	  mindeg = CS_MIN (mindeg, d) ;       /* find new minimum degree */
	  degree [i] = d ;
	  Cj [p++] = i ;                      /* place i in Lk */
        }
      nv [k] = nvk ;                      /* # nodes absorbed into k */
      if ((len [k] = p-pk1) == 0)         /* length of adj list of element k*/
        {
	  Cp [k] = -1 ;                   /* k is a root of the tree */
	  w [k] = 0 ;                     /* k is now a dead element */
        }
      if (elenk != 0) cnz = p ;           /* free unused space in Lk */
    }
  /* --- Postordering ----------------------------------------------------- */
  for (i = 0 ; i < n ; i++) Cp [i] = CS_FLIP (Cp [i]) ;/* fix assembly tree */
  for (j = 0 ; j <= n ; j++) head [j] = -1 ;
  for (j = n ; j >= 0 ; j--)              /* place unordered nodes in lists */
    {
      if (nv [j] > 0) continue ;          /* skip if j is an element */
      next [j] = head [Cp [j]] ;          /* place j in list of its parent */
      head [Cp [j]] = j ;
    }
  for (e = n ; e >= 0 ; e--)              /* place elements in lists */
    {
      if (nv [e] <= 0) continue ;         /* skip unless e is an element */
      if (Cp [e] != -1)
        {
	  next [e] = head [Cp [e]] ;      /* place e in list of its parent */
	  head [Cp [e]] = e ;
        }
    }
  for (k = 0, i = 0 ; i <= n ; i++)       /* postorder the assembly tree */
    {
      if (Cp [i] == -1) k = csr_tdfs (i, k, head, next, P, w) ;
    }
  return (csr_idone (P, C, W, 1)) ;
}


//-----------------------------------------------------------------
// Misc utilities
//-----------------------------------------------------------------

/* print a sparse matrix */
int csr_print (const csr *A, int brief)
{
  int p, j, m, n, nzmax, nz, nnz, *Ap, *Aj ;
  double *Ax ;
  if (!A) { printf ("(null)\n") ; return (0) ; }
  m = A->m ; n = A->n ; Ap = A->p ; Aj = A->j ; Ax = A->x ;
  nzmax = A->nzmax ; nz = A->nz ; nnz = Ap [m];
  if (nz < 0)
    {
      if( nnz == 0)
        {  /* print intermeidate matrices from csr_lu */
	  while ((Ap[m] == 0) && (m > 0))
	    {
              --m;
	    }
	  nnz = Ap [m];
        }
      if( nnz > 0)
        {
	  printf ("%d-by-%d, nzmax: %d nnz: %d, mxnorm: %g\n", m, n, nzmax,
		  Ap [m], csr_norm (A)) ;
	  for (j = 0 ; j < m ; j++)
            {
	      printf ("    row %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
	      for (p = Ap [j] ; p < Ap [j+1] ; p++)
                {
		  printf ("      %d : %g\n", Aj [p], Ax ? Ax [p] : 1) ;
		  if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
                }
            }
        }else{
	printf ("%d-by-%d, zero matrix with nzmax: %d\n", m, n, nzmax);
      }
    }
  else
    {
      printf ("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
      for (p = 0 ; p < nz ; p++)
        {
	  printf ("    %d %d : %g\n", Aj [p], Ap [p], Ax ? Ax [p] : 1) ;
	  if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
        }
    }
  return (1) ;
}

/* infinity-norm of a sparse matrix = norm(A,inf), max row sum */
double csr_norm (const csr *A)
{
    int p, j, m, *Ap ;
    double *Ax,  norm = 0, s ;
    if (!CS_CSC (A) || !A->x) return (-1) ;        /* check inputs */
    m = A->m ; Ap = A->p ; Ax = A->x ;
    for (j = 0 ; j < m ; j++)
    {
        for (s = 0, p = Ap [j] ; p < Ap [j+1] ; p++) s += fabs (Ax [p]);
        norm = CS_MAX (norm, s) ;
    }
    return (norm) ;
}

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
int csr_fkeep (csr *A, int (*fkeep) (int, int, double, void *), void *other)
{
    int j, p, nz = 0, m, *Ap, *Aj ;
    double *Ax ;
    if (!CS_CSC (A) || !fkeep) return (-1) ;    /* check inputs */
    m = A->m ; Ap = A->p ; Aj = A->j ; Ax = A->x ;
    for (j = 0 ; j < m ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        Ap [j] = nz ;                       /* record new location of col j */
        for ( ; p < Ap [j+1] ; p++)
        {
            if (fkeep (Aj [p], j, Ax ? Ax [p] : 1, other))
            {
                if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
                Aj [nz++] = Aj [p] ;
            }
        }
    }
    Ap [m] = nz ;                           /* finalize A */
    csr_sprealloc (A, 0) ;                  /* remove extra space from A */
    return (nz) ;
}
