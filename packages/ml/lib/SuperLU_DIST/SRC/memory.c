/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include <malloc.h>
#include "superlu_ddefs.h"


/*
 * Global variables
 */
ExpHeader *expanders; /* Array of pointers to 4 types of memory */
LU_stack_t stack;
int_t no_expand;


/*
 * Prototype
 */
static int_t memory_usage(const int_t, const int_t, const int_t);
static void *expand(int_t *, MemType, int_t, int_t,
		    Glu_freeable_t *);

/*
 * Internal prototypes
 */
void  SetupSpace (void *, int_t, LU_space_t *);


void
superlu_abort_and_exit(char *msg)
{
    /*fprintf(stderr, msg);
    fflush(stderr);*/
    printf(msg);
    exit (-1);
}

int_t superlu_malloc_total = 0;

#if ( DEBUGlevel>=1 )           /* Debug malloc/free. */

#define PAD_FACTOR  2
#define DWORD  (sizeof(double)) /* Be sure it's no smaller than double. */

void *superlu_malloc(int_t size)
{
    char *buf;
    int iam;

    MPI_Comm_rank(MPI_COMM_WORLD, &iam);
    buf = (char *) malloc(size + DWORD);
    if ( !buf ) {
	printf("(%d) superlu_malloc fails: malloc_total %.0f MB, size %d\n",
	       iam, superlu_malloc_total*1e-6, size);
	ABORT("superlu_malloc: out of memory");
    }

    ((int_t *) buf)[0] = size;
    superlu_malloc_total += size + DWORD;

    return (void *) (buf + DWORD);
}

void superlu_free(void *addr)
{
    char *p = ((char *) addr) - DWORD;

    if ( !addr )
	ABORT("superlu_free: tried to free NULL pointer");

    if ( !p )
	ABORT("superlu_free: tried to free NULL+DWORD pointer");

    { 
	int_t n = ((int_t *) p)[0];
	
	if ( !n )
	    ABORT("superlu_free: tried to free a freed pointer");
	*((int_t *) p) = 0; /* Set to zero to detect duplicate free's. */
	
	superlu_malloc_total -= (n + DWORD);

	if ( superlu_malloc_total < 0 )
	    ABORT("superlu_malloc_total went negative!");
	
	/*free (addr);*/
	free (p);
    }

}

#else  /* The production mode. */

void *superlu_malloc(int_t size)
{
    void *buf;
    buf = (void *) malloc(size);
    return (buf);
}

void superlu_free(void *addr)
{
    free (addr);
}

#endif  /* End debug malloc/free. */



void
copy_mem_int(int_t howmany, void *old, void *new)
{
    register int_t i;
    int_t *iold = old;
    int_t *inew = new;
    for (i = 0; i < howmany; i++) inew[i] = iold[i];
}


void
user_bcopy(char *src, char *dest, int_t bytes)
{
    char *s_ptr, *d_ptr;

    s_ptr = src + bytes - 1;
    d_ptr = dest + bytes - 1;
    for (; d_ptr >= dest; --s_ptr, --d_ptr ) *d_ptr = *s_ptr;
}



int_t *intMalloc(int_t n)
{
    int_t *buf;
    buf = (int_t *) SUPERLU_MALLOC(n * sizeof(int_t));
    return (buf);
}

int_t *intCalloc(int_t n)
{
    int_t *buf;
    register int_t i;
    buf = (int_t *) SUPERLU_MALLOC(n * sizeof(int_t));
    if ( buf )
	for (i = 0; i < n; ++i) buf[i] = 0;
    return (buf);
}


void *user_malloc(int_t bytes, int_t which_end)
{
    void *buf;
    
    if ( StackFull(bytes) ) return (NULL);

    if ( which_end == HEAD ) {
	buf = (char*) stack.array + stack.top1;
	stack.top1 += bytes;
    } else {
	stack.top2 -= bytes;
	buf = (char*) stack.array + stack.top2;
    }
    
    stack.used += bytes;
    return buf;
}

void user_free(int_t bytes, int_t which_end)
{
    if ( which_end == HEAD ) {
	stack.top1 -= bytes;
    } else {
	stack.top2 += bytes;
    }
    stack.used -= bytes;
}


/*
 * Setup the memory model to be used for factorization.
 *    lwork = 0: use system malloc;
 *    lwork > 0: use user-supplied work[] space.
 */
void SetupSpace(void *work, int_t lwork, LU_space_t *MemModel)
{
    if ( lwork == 0 ) {
	*MemModel = SYSTEM; /* malloc/free */
    } else if ( lwork > 0 ) {
	*MemModel = USER;   /* user provided space */
	stack.used = 0;
	stack.top1 = 0;
	stack.top2 = (lwork/4)*4; /* must be word addressable */
	stack.size = stack.top2;
	stack.array = (void *) work;
    }
}


/*
 * Set up pointers for integer working arrays.
 */
/************************************************************************/
void SetIWork
/************************************************************************/
(
 int_t m, int_t n, int_t panel_size, int_t *iworkptr, int_t **segrep,
 int_t **parent, int_t **xplore, int_t **repfnz, int_t **panel_lsub,
 int_t **xprune, int_t **marker
 )
{
    *segrep = iworkptr;
    *parent = iworkptr + m;
    *xplore = *parent + m;
    *repfnz = *xplore + m;
    *panel_lsub = *repfnz + panel_size * m;
    *xprune = *panel_lsub + panel_size * m;
    *marker = *xprune + n;
    ifill (*repfnz, n * panel_size, EMPTY);
    ifill (*panel_lsub, m * panel_size, EMPTY);
}




/*
 * Allocate storage for the data structures common to symbolic factorization
 * routines. For those unpredictable size, make a guess as FILL * nnz(A).
 * Return value:
 *     If lwork = -1, return the estimated amount of space required, plus n;
 *     otherwise, return the amount of space actually allocated when
 *     memory allocation failure occurred.
 */
/************************************************************************/
int_t symbfact_SubInit
/************************************************************************/
(
 fact_t fact, void *work, int_t lwork, int_t m, int_t n, int_t annz,
 Glu_persist_t *Glu_persist, Glu_freeable_t *Glu_freeable
 )
{
    int_t  iword;
    int_t  *xsup, *supno;
    int_t  *lsub, *xlsub;
    int_t  *usub, *xusub;
    int_t  nzlmax, nzumax;
    int_t  FILL = sp_ienv(6);

#if ( DEBUGlevel>=1 )
    int iam;
    MPI_Comm_rank( MPI_COMM_WORLD, &iam );
    CHECK_MALLOC(iam, "Enter symbfact_SubInit()");
#endif

    no_expand = 0;
    iword     = sizeof(int_t);

    expanders = (ExpHeader *) SUPERLU_MALLOC( NO_MEMTYPE*sizeof(ExpHeader) );
    if ( !expanders ) ABORT("SUPERLU_MALLOC fails for expanders");
    
    if ( fact == DOFACT || fact == SamePattern ) {
	/* Guess for L\U factors */
	nzlmax = nzumax = FILL * annz;

	if ( lwork == -1 ) {
	    return ( GluIntArray(n) * iword + TempSpace(m,1)
		    + (nzlmax+nzumax)*iword + n );
        } else {
	    SetupSpace(work, lwork, &Glu_freeable->MemModel);
	}
	
#if ( DEBUGlevel>=1 )
	printf(".. symbfact_SubInit(): annz %d, MemModel %d, nzumax %d, nzlmax %d\n", 
		annz, Glu_freeable->MemModel, nzumax, nzlmax);
#endif	
	
	/* Integer pointers for L\U factors */
	if ( Glu_freeable->MemModel == SYSTEM ) {
	    xsup   = intMalloc(n+1);
	    supno  = intMalloc(n+1);
	    xlsub  = intMalloc(n+1);
	    xusub  = intMalloc(n+1);
	} else {
	    xsup   = (int_t *)user_malloc((n+1) * iword, HEAD);
	    supno  = (int_t *)user_malloc((n+1) * iword, HEAD);
	    xlsub  = (int_t *)user_malloc((n+1) * iword, HEAD);
	    xusub  = (int_t *)user_malloc((n+1) * iword, HEAD);
	}

	lsub  = (int_t *) expand(&nzlmax, LSUB, 0, 0, Glu_freeable);
	usub  = (int_t *) expand(&nzumax, USUB, 0, 0, Glu_freeable);

	while ( !lsub || !usub ) {
	    if ( Glu_freeable->MemModel == SYSTEM ) {
		SUPERLU_FREE(lsub); 
		SUPERLU_FREE(usub);
	    } else {
		user_free((nzlmax+nzumax)*iword, HEAD);
	    }
	    nzlmax /= 2;
	    nzumax /= 2;
	    if ( nzumax < annz/2 ) {
		printf("Not enough memory to perform factorization.\n");
		return (memory_usage(nzlmax, nzumax, n) + n);
	    }
	    lsub  = (int_t *) expand( &nzlmax, LSUB, 0, 0, Glu_freeable );
	    usub  = (int_t *) expand( &nzumax, USUB, 0, 1, Glu_freeable );
	}

	Glu_persist->xsup    = xsup;
	Glu_persist->supno   = supno;
	Glu_freeable->lsub   = lsub;
	Glu_freeable->xlsub  = xlsub;
	Glu_freeable->usub   = usub;
	Glu_freeable->xusub  = xusub;
	Glu_freeable->nzlmax = nzlmax;
	Glu_freeable->nzumax = nzumax;
    } else {
	/* fact == SamePattern_SameRowPerm */
	if ( lwork == -1 ) {
	    return ( GluIntArray(n) * iword + TempSpace(m, 1)
		    + (nzlmax+nzumax)*iword + n );
        } else if ( lwork == 0 ) {
	    Glu_freeable->MemModel = SYSTEM;
	} else {
	    Glu_freeable->MemModel = USER;
	    stack.top2 = (lwork/4)*4; /* must be word-addressable */
	    stack.size = stack.top2;
	}
	
	expanders[USUB].mem = Glu_freeable->usub;
	expanders[LSUB].mem = Glu_freeable->lsub;
	expanders[USUB].size = nzumax;
	expanders[LSUB].size = nzlmax;
    }

    ++no_expand;

#if ( DEBUGlevel>=1 )
    /* Memory allocated but not freed: xsup, supno */
    CHECK_MALLOC(iam, "Exit symbfact_SubInit()");
#endif

    return 0;
    
} /* SYMBFACT_SUBINIT */

/*
 * Expand the data structures for L and U during the factorization.
 * Return value:   0 - successful return
 *               > 0 - number of bytes allocated when run out of space
 */
/************************************************************************/
int_t symbfact_SubXpand
/************************************************************************/
(
 int_t n,           /* total number of columns */
 int_t jcol,        /* current column */
 int_t next,        /* number of elements currently in the factors */
 MemType mem_type,  /* which type of memory to expand  */
 int_t *maxlen,     /* modified - maximum length of a data structure */
 Glu_freeable_t *Glu_freeable  /* modified - global LU data structures */
 )
{
    void   *new_mem;
    
#if ( DEBUGlevel>=1 )
    printf("symbfact_SubXpand(): jcol %d, next %d, maxlen %d, MemType %d\n",
	   jcol, next, *maxlen, mem_type);
#endif    

    new_mem = expand(maxlen, mem_type, next, 0, Glu_freeable);
    
    if ( !new_mem ) {
	int_t    nzlmax  = Glu_freeable->nzlmax;
	int_t    nzumax  = Glu_freeable->nzumax;
    	fprintf(stderr, "Can't expand MemType %d: jcol %d\n", mem_type, jcol);
    	return (memory_usage(nzlmax, nzumax, n) + n);
    }

    if ( mem_type == LSUB ) {
	Glu_freeable->lsub   = (int_t *) new_mem;
	Glu_freeable->nzlmax = *maxlen;
    } else if ( mem_type == USUB ) {
	Glu_freeable->usub   = (int_t *) new_mem;
	Glu_freeable->nzumax = *maxlen;
    } else ABORT("Tries to expand nonexisting memory type.\n");
    
    return 0;
    
} /* LUSUB_XPAND */

/*
 * Deallocate storage of the data structures common to symbolic
 * factorization routines.
 */
/************************************************************************/
int_t symbfact_SubFree(Glu_freeable_t *Glu_freeable)
/************************************************************************/
{
#if ( DEBUGlevel>=1 )
    int iam;
    MPI_Comm_rank( MPI_COMM_WORLD, &iam );
    CHECK_MALLOC(iam, "Enter symbfact_SubFree()");
#endif
    
    SUPERLU_FREE(expanders);
    SUPERLU_FREE(Glu_freeable->lsub);
    SUPERLU_FREE(Glu_freeable->xlsub);
    SUPERLU_FREE(Glu_freeable->usub);
    SUPERLU_FREE(Glu_freeable->xusub);

#if ( DEBUGlevel>=1 )    
    CHECK_MALLOC(iam, "Exit symbfact_SubFree()");
#endif
    return 0;
} /* SYMBFACT_SUBFREE */


/*
 * Expand the existing storage to accommodate more fill-ins.
 */
/************************************************************************/
static void *expand
/************************************************************************/
(
 int_t *prev_len,   /* length used from previous call */
 MemType type,    /* which part of the memory to expand */
 int_t len_to_copy, /* size of the memory to be copied to new store */
 int_t keep_prev,   /* = 1: use prev_len;
		     = 0: compute new_len to expand */
 Glu_freeable_t *Glu_freeable  /* modified - global LU data structures */
 )
{
    float    EXPAND = 1.5;
    float    alpha;
    void     *new_mem;
    int_t    new_len, tries, lword, extra, bytes_to_copy;

    alpha = EXPAND;
    lword = sizeof(int_t);

    if ( no_expand == 0 || keep_prev ) /* First time allocate requested */
        new_len = *prev_len;
    else {
	new_len = alpha * *prev_len;
    }

    if ( Glu_freeable->MemModel == SYSTEM ) {
	new_mem = (void *) SUPERLU_MALLOC(new_len * lword);
	/*new_mem = (void *) calloc(new_len, lword); */
	if ( no_expand != 0 ) {
	    tries = 0;
	    if ( keep_prev ) {
		if ( !new_mem ) return (NULL);
	    } else {
		while ( !new_mem ) {
		    if ( ++tries > 10 ) return (NULL);
		    alpha = Reduce(alpha);
		    new_len = alpha * *prev_len;
		    new_mem = (void *) SUPERLU_MALLOC(new_len * lword); 
		    /* new_mem = (void *) calloc(new_len, lword); */
		}
	    }
	    copy_mem_int(len_to_copy, expanders[type].mem, new_mem);
	    SUPERLU_FREE (expanders[type].mem);
	}
	expanders[type].mem = (void *) new_mem;
	
    } else { /* MemModel == USER */
	if ( no_expand == 0 ) {
	    new_mem = user_malloc(new_len * lword, HEAD);
	    expanders[type].mem = (void *) new_mem;
	}
	else {
	    tries = 0;
	    extra = (new_len - *prev_len) * lword;
	    if ( keep_prev ) {
		if ( StackFull(extra) ) return (NULL);
	    } else {
		while ( StackFull(extra) ) {
		    if ( ++tries > 10 ) return (NULL);
		    alpha = Reduce(alpha);
		    new_len = alpha * *prev_len;
		    extra = (new_len - *prev_len) * lword;	    
		}
	    }

	    if ( type != USUB ) {
		new_mem = (void*)((char*)expanders[type + 1].mem + extra);
		bytes_to_copy = (char*)stack.array + stack.top1
		    - (char*)expanders[type + 1].mem;
		user_bcopy(expanders[type+1].mem, new_mem, bytes_to_copy);

		if ( type < USUB ) {
		    Glu_freeable->usub = expanders[USUB].mem =
			(void*)((char*)expanders[USUB].mem + extra);
		}
		if ( type < LSUB ) {
		    Glu_freeable->lsub = expanders[LSUB].mem =
			(void*)((char*)expanders[LSUB].mem + extra);
		}
		stack.top1 += extra;
		stack.used += extra;
		
	    } /* if ... */

	} /* else ... */
    }

    expanders[type].size = new_len;
    *prev_len = new_len;
    if ( no_expand ) ++no_expand;
    
    return (void *) expanders[type].mem;
    
} /* EXPAND */

/*
 * mem_usage consists of the following fields:
 *    - for_lu (float)
 *      The amount of space used in bytes for the L\U data structures.
 *    - total (float)
 *      The amount of space needed in bytes to perform factorization.
 *    - expansions (int)
 *      Number of memory expansions during the LU factorization.
 */
/************************************************************************/
int_t QuerySpace(int_t n, int_t lsub_size, Glu_freeable_t *Glu_freeable,
		 mem_usage_t *mem_usage)
/************************************************************************/
{
    register int_t iword = sizeof(int_t);
    extern int_t no_expand;

    /* For the adjacency graphs of L and U. */
    /*mem_usage->for_lu = (float)( (4*n + 3) * iword +
				Glu_freeable->xlsub[n]*iword );*/
    mem_usage->for_lu = (float)( (4*n + 3) * iword +
				lsub_size * iword );
    mem_usage->for_lu += (float)( (n + 1) * iword +
				 Glu_freeable->xusub[n]*iword );

    /* Working storage to support factorization */
    mem_usage->total = mem_usage->for_lu + 9*n*iword;

    mem_usage->expansions = --no_expand;
    return 0;
} /* QUERYSPACE */

static int_t
memory_usage(const int_t nzlmax, const int_t nzumax, const int_t n)
{
    register int_t iword = sizeof(int_t);
    return (10*n*iword + (nzlmax+nzumax)*iword);
}


#if 0

LBlock_t *New_Lblock()
{
    LBlock_t *buf;
    if ( !(buf = (LBlock_t *) SUPERLU_MALLOC(sizeof(LBlock_t))) )
        ABORT("New_Lblock fails.");
    buf->next = NULL;
    return buf;
}

void Free_Lblock(LBlock_t *b)
{
    SUPERLU_FREE(b);
}
    
UBlock_t *New_Ublock()
{
    UBlock_t *buf;
    if ( !(buf = (UBlock_t *) SUPERLU_MALLOC(sizeof(UBlock_t))) )
        ABORT("New_Ublock fails.");
    buf->next = NULL;
    return buf;
}

void Free_Ublock(UBlock_t *b)
{
    SUPERLU_FREE(b);
}

#endif
