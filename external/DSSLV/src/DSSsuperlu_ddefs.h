/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#ifndef __SUPERLU_dDEFS /* allow multiple inclusions */
#define __SUPERLU_dDEFS

/*
 * File name:	superlu_ddefs.h
 * Purpose:     Distributed SuperLU data types and function prototypes
 * History:
 */
#ifdef _CRAY
#include <fortran.h>
#include <string.h>
#endif

#include <mpi.h>
// #include "DSS_Cnames.h"
// #include "supermatrix.h"
#include "DSSutil.h"

#ifdef _CRAY
typedef short int_t;
#undef int       /* Revert back to int of default size. */
#define mpi_int_t   MPI_SHORT
#else
typedef int int_t;
#define mpi_int_t   MPI_INT
#endif


/***********************************************************************
 * Enumerated types
 ***********************************************************************/
/*typedef enum {FALSE, TRUE}                                      boolean_t;*/
typedef enum {NO, YES}                                          yes_no_t;
typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;
typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD, MY_PERMC}colperm_t;
typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
typedef enum {NOREFINE, DOUBLE = 1, EXTRA}                      IterRefine_t;
typedef enum {LUSUP, UCOL, LSUB, USUB}                          MemType;
typedef enum {HEAD, TAIL}                                       stack_end_t;
typedef enum {SYSTEM, USER}                                     LU_space_t;


/***********************************************************************
 * Constants
 ***********************************************************************/
/* 
 * For each block column of L, the index[] array contains both the row 
 * subscripts and the integers describing the size of the blocks.
 * The organization of index[] looks like:
 *
 *     [ BLOCK COLUMN HEADER (size BC_HEADER)
 *           number of blocks 
 *           number of row subscripts, i.e., LDA of nzval[]
 *       BLOCK 0                                        <----
 *           BLOCK DESCRIPTOR (of size LB_DESCRIPTOR)  |
 *               block number (global)                      |
 *               number of full rows in the block           |
 *           actual row subscripts                          |
 *       BLOCK 1                                            | Repeat ...
 *           BLOCK DESCRIPTOR                               | number of blocks
 *               block number (global)                      | 
 *               number of full rows in the block           |
 *           actual row subscripts                          |
 *       .                                                  |
 *       .                                                  |
 *       .                                              <----
 *     ]
 *
 * For each block row of U, the organization of index[] looks like:
 *
 *     [ BLOCK ROW HEADER (of size BR_HEADER)
 *           number of blocks 
 *           number of entries in nzval[]
 *           number of entries in index[]
 *       BLOCK 0                                        <----
 *           BLOCK DESCRIPTOR (of size UB_DESCRIPTOR)  |
 *               block number (global)                      |
 *               number of nonzeros in the block            |
 *           actual fstnz subscripts                        |
 *       BLOCK 1                                            | Repeat ...
 *           BLOCK DESCRIPTOR                               | number of blocks
 *               block number (global)                      |
 *               number of nonzeros in the block            |
 *           actual fstnz subscripts                        |
 *       .                                                  |
 *       .                                                  |
 *       .                                              <----
 *     ]
 *
 */
#define BC_HEADER      2
#define LB_DESCRIPTOR  2
#define BR_HEADER      3
#define UB_DESCRIPTOR  2
#define NBUFFERS       5

/*
 * Communication tags
 */
    /* For numeric factorization. */
#define NTAGS    10000
#define UjROW    10
#define UkSUB    11
#define UkVAL    12
#define LkSUB    13
#define LkVAL    14
#define LkkDIAG  15
    /* For triangular solves. */
#define XK_H     1  /* The header preceeding each X block. */
#define LSUM_H   1  /* The header preceeding each MOD block. */
#define GSUM     20 
#define Xk       21
#define Yk       22
#define LSUM     23

#define YES      1

/* 
 * Communication scopes
 */
#define COMM_ALL      100
#define COMM_COLUMN   101
#define COMM_ROW      102

/*
 * Matrix distribution for sparse matrix-vector multiplication
 */
#define SUPER_LINEAR     11
#define SUPER_BLOCK      12

/*
 * No of marker arrays used in the symbolic factorization, each of size n
 */
#define NO_MARKER     3



/***********************************************************************
 * Macros
 ***********************************************************************/
#define IAM(comm)    { int rank; MPI_Comm_rank ( comm, &rank ); rank};
#define MYROW(iam,grid) ( (iam) / grid->npcol )
#define MYCOL(iam,grid) ( (iam) % grid->npcol )
#define BlockNum(i)     ( supno[i] )
#define FstBlockC(bnum) ( xsup[bnum] )
#define SuperSize(bnum) ( xsup[bnum+1]-xsup[bnum] )
#define LBi(bnum,grid)  ( (bnum)/grid->nprow )/* Global to local block rowwise */
#define LBj(bnum,grid)  ( (bnum)/grid->npcol )/* Global to local block columnwise*/
#define PROW(bnum,grid) ( (bnum) % grid->nprow )
#define PCOL(bnum,grid) ( (bnum) % grid->npcol )
#define PNUM(i,j,grid)  ( (i)*grid->npcol + j ) /* Process number at coord(i,j) */
#define CEILING(a,b)    ( ((a)%(b)) ? ((a)/(b) + 1) : ((a)/(b)) )
    /* For triangular solves */
#define RHS_ITERATE(i)                    \
        for (i = 0; i < nrhs; ++i)
#define X_BLK(i)                          \
        ilsum[i] * nrhs + (i+1) * XK_H
#define LSUM_BLK(i)                       \
        ilsum[i] * nrhs + (i+1) * LSUM_H


#if ( VAMPIR>=1 ) 
#define VT_TRACEON    VT_traceon()
#define VT_TRACEOFF   VT_traceoff()
#else
#define VT_TRACEON 
#define VT_TRACEOFF
#endif


/***********************************************************************
 * New data types
 ***********************************************************************/

/* 
 *   Define the 2D mapping of matrix blocks to process grid.
 *
 *   Process grid:
 *     Processes are numbered (0 : P-1).
 *     P = Pr x Pc, where Pr, Pc are the number of process rows and columns.
 *     (pr,pc) is the coordinate of IAM; 0 <= pr < Pr, 0 <= pc < Pc.
 *
 *   Matrix blocks:
 *     Matrix is partitioned according to supernode partitions, both
 *     column and row-wise. 
 *     The k-th block columns (rows) contains columns (rows) (s:t), where
 *             s=xsup[k], t=xsup[k+1]-1.
 *     Block A(I,J) contains
 *             rows from (xsup[I]:xsup[I+1]-1) and
 *             columns from (xsup[J]:xsup[J+1]-1)
 *
 *  Mapping of matrix entry (i,j) to matrix block (I,J):
 *     (I,J) = ( supno[i], supno[j] )
 *
 *  Mapping of matrix block (I,J) to process grid (pr,pc):
 *     (pr,pc) = ( MOD(I,NPROW), MOD(J,NPCOL) )
 *  
 *  (xsup[nsupers],supno[n]) are replicated on all processors.
 *
 */

/*-- Communication subgroup */
typedef struct {
    MPI_Comm comm;        /* MPI communicator */
    int Np;               /* number of processes */
    int Iam;              /* my process number */
} superlu_scope_t;

/*-- Process grid definition */
typedef struct {
    MPI_Comm comm;        /* MPI communicator */
    superlu_scope_t rscp; /* row scope */
    superlu_scope_t cscp; /* column scope */
    int iam;              /* my process number in this scope */
    int_t nprow;          /* number of process rows */
    int_t npcol;          /* number of process columns */
} gridinfo_t;


/*
 *-- The structures are determined by SYMBFACT and used thereafter.
 *
 * (xsup,supno) describes mapping between supernode and column:
 *	xsup[s] is the leading column of the s-th supernode.
 *      supno[i] is the supernode no to which column i belongs;
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 * This is allocated during symbolic factorization SYMBFACT.
 */
typedef struct {
    int_t     *xsup;
    int_t     *supno;
} Glu_persist_t;

/*
 *-- The structures are determined by SYMBFACT and used by DDISTRIBUTE.
 * 
 * (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 * (xusub,usub): xusub[i] points to the starting location of column i
 *      in usub[]. For each U-segment, only the row index of first nonzero
 *      is stored in usub[].
 *
 *      Each U column consists of a number of full segments. Each full segment
 *      starts from a leading nonzero, running up to the supernode (block)
 *      boundary. (Recall that the column-wise supernode partition is also
 *      imposed on the rows.) Because the segment is full, we don't store all
 *      the row indices. Instead, only the leading nonzero index is stored.
 *      The rest can be found together with xsup/supno pair.
 *      For example, 
 *          usub[xsub[j+1]] - usub[xsub[j]] = number of segments in column j.
 *          for any i in usub[], 
 *              supno[i]         = block number in which i belongs to
 *  	        xsup[supno[i]+1] = first row of the next block
 *              The nonzeros of this segment are: 
 *                  i, i+1 ... xsup[supno[i]+1]-1 (only i is stored in usub[])
 *
 */
typedef struct {
    int_t     *lsub;     /* compressed L subscripts */
    int_t     *xlsub;
    int_t     *usub;     /* compressed U subscripts */
    int_t     *xusub;
    int_t     nzlmax;    /* current max size of lsub */
    int_t     nzumax;    /*    "    "    "      usub */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
} Glu_freeable_t;


/* 
 *-- The structure used to store matrix A of the linear system and
 *   several vectors describing the transformations done to matrix A.
 *
 * A      (SuperMatrix*)
 *        Matrix A in A*X=B, of dimension (A->nrow, A->ncol).
 *        The number of linear equations is A->nrow. The type of A can be:
 *        Stype = NC; Dtype = _D; Mtype = GE.
 *         
 * DiagScale  (DiagScale_t)
 *        Specifies the form of equilibration that was done.
 *        = NOEQUIL: No equilibration.
 *        = ROW:  Row equilibration, i.e., A was premultiplied by diag(R).
 *        = COL:  Column equilibration, i.e., A was postmultiplied by diag(C).
 *        = BOTH: Both row and column equilibration, i.e., A was replaced 
 *                 by diag(R)*A*diag(C).
 *
 * R      double*, dimension (A->nrow)
 *        The row scale factors for A.
 *        If DiagScale = ROW or BOTH, A is multiplied on the left by diag(R).
 *        If DiagScale = NOEQUIL or COL, R is not defined.
 *
 * C      double*, dimension (A->ncol)
 *        The column scale factors for A.
 *        If DiagScale = COL or BOTH, A is multiplied on the right by diag(C).
 *        If DiagScale = NOEQUIL or ROW, C is not defined.
 *         
 * perm_r (int*) dimension (A->nrow)
 *        Row permutation vector which defines the permutation matrix Pr,
 *        perm_r[i] = j means row i of A is in position j in Pr*A.
 *
 * perm_c (int*) dimension (A->ncol)
 *	  Column permutation vector, which defines the 
 *        permutation matrix Pc; perm_c[i] = j means column i of A is 
 *        in position j in A*Pc.
 *
 */
typedef struct {
    DiagScale_t DiagScale;
    double *R;
    double *C; 
    int_t  *perm_r;
    int_t  *perm_c;
} ScalePermstruct_t;

/* 
 * On each processor, the blocks in L are stored in compressed block
 * column format, the blocks in U are stored in compressed block row format.
 */
typedef struct {
    int_t   **Lrowind_bc_ptr; /* size ceil(NSUPERS/Pc)                 */
    double  **Lnzval_bc_ptr;  /* size ceil(NSUPERS/Pc)                 */
    int_t   **Ufstnz_br_ptr;  /* size ceil(NSUPERS/Pr)                 */
    double  **Unzval_br_ptr;  /* size ceil(NSUPERS/Pr)                 */
#if 0
    int_t   *Lsub_buf;        /* Buffer for the remote subscripts of L */
    double  *Lval_buf;        /* Buffer for the remote nonzeros of L   */
#endif
    int_t   *Lsub_buf_2[2];   /* Buffers for the remote subscripts of L*/
    double  *Lval_buf_2[2];   /* Buffers for the remote nonzeros of L  */
    int_t   *Usub_buf;        /* Buffer for the remote subscripts of U */
    double  *Uval_buf;        /* Buffer for the remote nonzeros of U   */
    double  *ujrow;           /* used in panel factorization.          */
    int_t   bufmax[NBUFFERS]; /* Buffer size; 5 entries
			       *     0 : size of Lsub_buf[]
			       *     1 : size of Lval_buf[]
			       *     2 : size of Usub_buf[] 
			       *     3 : size of Uval_buf[]
			       *     4 : size of tempv[LDA]
			       */

    /*-- Record communication schedule for factorization. --*/
    int_t   *ToRecv;          /* Recv from no one (0), left (1), and up (2).*/
    int_t   *ToSendD;         /* Whether need to send down block row.       */
    int_t   **ToSendR;        /* List of processes to send right block col. */

    /*-- Record communication schedule for solves. --*/
    int_t   *fmod;            /* Modification count for L-solve.            */
    int_t   **fsendx_plist;   /* Column process list to send down Xk.       */
    int_t   *frecv;           /* Modifications to be recv'd in proc row.    */
    int_t   nfrecvx;          /* Number of Xk I will receive in L-solve.    */
    int_t   *bmod;            /* Modification count for U-solve.            */
    int_t   **bsendx_plist;   /* Column process list to send down Xk.       */
    int_t   *brecv;           /* Modifications to be recv'd in proc row.    */
    int_t   nbrecvx;          /* Number of Xk I will receive in U-solve.    */

    /*-- Auxiliary arrays used for solves. --*/
    int_t   *ilsum;           /* Starting position of each supernode in lsum
				 (local)  */
    int_t   ldalsum;          /* LDA of lsum (local) */
} LocalLU_t;

typedef struct {
    int_t *etree;
    Glu_persist_t *Glu_persist;
    LocalLU_t *Llu;
} LUstruct_t;

/*-- Auxiliary data type used in PxGSTRS/PxGSTRS1. */
typedef struct {
    int_t lbnum;  /* Row block number (local).      */
    int_t indpos; /* Starting position in Uindex[]. */
} Ucb_indptr_t;


/* 
 *-- This contains the options used to control the solve process.
 *
 * Fact   (fact_t)
 *        Specifies whether or not the factored form of the matrix
 *        A is supplied on entry, and if not, how the matrix A should
 *        be factorizaed.
 *        = FACTORED: On entry, L, U, perm_r and perm_c contain the 
 *              factored form of A. If DiagScale is not NOEQUIL, the matrix
 *              A has been equilibrated with scaling factors R and C.
 *              A, L, U, perm_r are not modified.
 *        = DOFACT: The matrix A will be factored, and the factors will be
 *              stored in L and U.
 *        = EQUILIBRATE: The matrix A will be equilibrated if necessary, then
 *              factored into L and U.
 *
 * Trans  (trans_t)
 *        Specifies the form of the system of equations:
 *        = NOTRANS: A * X = B        (No transpose)
 *        = TRANS:   A**T * X = B     (Transpose)
 *        = CONJ:    A**H * X = B     (Transpose)
 *
 * refact (refact_t)
 *        Specifies whether this is first time or subsequent factorization.
 *        = NO: this factorization is treated as the first one;
 *        = SamePattern: it means that a factorization of a matrix with
 *              the same sparsity pattern was performed prior to
 *              this one. Therefore, this factorization will re-use
 *              column permutation LUstruct->perm_c;
 *        = SamePattern_SameRowPerm: it means that a factorization of
 *              a matrix with the same sparsity pattern and similar
 *              numerical values was performed prior to this one.
 *              Therefore, this factorization will reuse the following
 *              data structures:
 *                  ScalePermstruct : DiagScale, R, C, perm_r, perm_c.
 *                  LUstruct : Glu_persist, Llu.
 *        This is useful only when the option 'fact' is not FACTORED.
 *
 * Equil  (yes_no_t)
 *        Specifies whether to equilibrate the system.
 *
 * RowPerm (rowperm_t)
 *        Specifies whether to permute rows of the original matrix.
 *        = NO: not to permute the rows
 *        = LargeDiag: make the diagonal large relative to the off-diagonal
 *        = MY_PERMR: use the permutation specified in ScalePermstruct->perm_r[]
 *           
 * ColPerm (colperm_t)
 *        Specifies what type of column permutation to use to reduce fill.
 *        = NATURAL: use the natural ordering 
 *        = MMD_ATA: use minimum degree ordering on structure of A'*A
 *        = MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
 *        = COLAMD: use approximate minimum degree column ordering
 *        = MY_PERMC: use the ordering specified in ScalePermstruct->perm_c[]
 *         
 * ReplaceTinyPivot (yes_no_t)
 *        Specifies whether to replace the tiny diagonals by
 *        sqrt(epsilon)*||A|| during LU factorization.
 *
 * IterRefine (IterRefine_t)
 *        Specifies whether to perform iterative refinement.
 *        = NO: no iterative refinement
 *        = WorkingPrec: perform iterative refinement in working precision
 *        = ExtraPrec: perform iterative refinement in extra precision
 *
 */
typedef struct {
    fact_t Fact;
    trans_t Trans;
    yes_no_t Equil;
    rowperm_t RowPerm;
    colperm_t ColPerm;
    yes_no_t ReplaceTinyPivot;
    IterRefine_t IterRefine;
} superlu_options_t;


typedef struct {
    float for_lu;
    float total;
    int_t   expansions;
} mem_usage_t;


/***********************************************************************
 * Function prototypes
 ***********************************************************************/

extern "C" {


/* Supernodal LU factor related */
extern void
dCreate_CompCol_Matrix_dist(SuperMatrix *, int_t, int_t, int_t, double *,
			    int_t *, int_t *, SluDist_Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_CompCol_Matrix_dist(SuperMatrix *, SuperMatrix *);
extern void
dCreate_Dense_Matrix_dist(SuperMatrix *, int_t, int_t, double *, int_t,
			  SluDist_Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix_dist(SuperMatrix *, int_t, int_t, int_t, double *, 
			      int_t *, int_t *, int_t *, int_t *, int_t *,
			      SluDist_Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_Dense_Matrix_dist(int_t, int_t, double *, int_t, double *, int_t);

extern void    Destroy_CompCol_Matrix_dist(SuperMatrix *);
extern void    Destroy_SuperNode_Matrix_dist(SuperMatrix *);
extern void    Destroy_CompCol_Permuted_dist(SuperMatrix *);

extern void    dallocateA_dist (int_t, int_t, double **, int_t **, int_t **);
extern void    sp_colorder (superlu_options_t*, SuperMatrix*, int_t*, int_t*,
			    SuperMatrix*);
extern int_t   sp_coletree_dist (int_t *, int_t *, int_t *, int_t, int_t,
				 int_t *);
extern void    countnz_dist (const int_t, int_t *, int_t *, int_t *,
			     Glu_persist_t *, Glu_freeable_t *);
extern int_t   fixupL_dist (const int_t, const int_t *, Glu_persist_t *,
			    Glu_freeable_t *);
extern int_t   *TreePostorder_dist (int_t, int_t *);
extern void    dGenXtrue_dist (int_t, int_t, double *, int_t);
extern void    dFillRHS_dist (char *, int_t, double *, int_t, SuperMatrix *,
			      double *, int_t);


/* Driver related */

extern void    dgsequ_dist (SuperMatrix *, double *, double *, double *,
			    double *, double *, int_t *);
extern double  dlangs_dist (char *, SuperMatrix *);
extern void    dlaqgs_dist (SuperMatrix *, double *, double *, double,
			    double, double, char *);
extern int     sp_dtrsv_dist (char *, char *, char *, SuperMatrix *,
			      SuperMatrix *, double *, int *);
extern int     sp_dgemv_dist (char *, double, SuperMatrix *, double *,
			      int, double, double *, int);
extern int     sp_dgemm (char *, char *, int, int, int, double,
			SuperMatrix *, double *, int, double, 
			double *, int);
extern double  dlamch_(char *);
extern double  slamch_(char *);

/* Memory-related */

extern void    *superlu_malloc_dist (int_t);
extern void    superlu_free_dist (void*);
extern int_t   *intMalloc_dist (int_t);
extern int_t   *intCalloc_dist (int_t);
extern double  *doubleMalloc_dist(int_t);
extern double  *doubleCalloc_dist(int_t);

/* Auxiliary routines */
extern double  SuperLU_timer_ ();
extern void    superlu_abort_and_exit_dist(char *);
extern int_t   sp_ienv_dist (int_t);
extern int     lsame_ (char *, char *);
extern int     xerbla_ (char *, int *);
extern void    ifill_dist (int_t *, int_t, int_t);
extern void    dfill_dist (double *, int_t, double);
extern void    dinf_norm_error_dist (int_t, int_t, double*, int_t, double*,
				     int_t);
extern void    super_stats_dist (int_t, int_t *);
extern void    PrintPerf (SuperMatrix *, SuperMatrix *, mem_usage_t *,
			 double, double, double *, double *, char *);


extern void    Destroy_LU(int_t, gridinfo_t *, LUstruct_t *);
extern void    ScalePermstructInit(const int_t, const int_t, 
				   ScalePermstruct_t *);
extern void    ScalePermstructFree(ScalePermstruct_t *);
extern void    LUstructInit(const int_t, const int_t, LUstruct_t *);
extern void    LUstructFree(LUstruct_t *);
extern void  dreadhb_dist (int, FILE *, int_t *, int_t *, int_t *, 
			   double **, int_t **, int_t **);
extern void  superlu_gridinit(MPI_Comm, int_t, int_t, gridinfo_t *);
extern void  superlu_gridmap(MPI_Comm, int_t, int_t, int_t [], int_t,
			     gridinfo_t *);
extern void  superlu_gridexit(gridinfo_t *);
extern void  get_perm_c_dist(int_t, int_t, SuperMatrix *, int_t *);
extern void  bcast_tree(void *, int, MPI_Datatype, int, int,
			gridinfo_t *, int, int *);
extern int_t symbfact(int, SuperMatrix *, int_t *, int_t *,
		      Glu_persist_t *, Glu_freeable_t *);
extern int_t symbfact_SubInit(fact_t, void *, int_t, int_t, int_t, int_t,
			      Glu_persist_t *, Glu_freeable_t *);
extern int_t symbfact_SubXpand(int_t, int_t, int_t, MemType, int_t *,
			       Glu_freeable_t *);
extern int_t symbfact_SubFree(Glu_freeable_t *);
extern int_t ddistribute(fact_t, int_t, SuperMatrix *, Glu_freeable_t *, 
			 LUstruct_t *, gridinfo_t *);
extern void  pdgssvx_ABglobal(superlu_options_t *, SuperMatrix *, 
			      ScalePermstruct_t *, double *,
			      int, int, gridinfo_t *, LUstruct_t *, double *,
			      SuperLUStat_t *, int *);
extern void  set_default_options(superlu_options_t *);
extern void  dldperm(int_t, int_t, int_t, int_t [], int_t [],
		     double [], int_t *, double [], double []);
extern void  pdgstrf(superlu_options_t *, int, int, double,
		     LUstruct_t*, gridinfo_t*, SuperLUStat_t*, int*);
extern void  pdgstrs_Bglobal(int_t, LUstruct_t *, gridinfo_t *,
			     double *, int_t, int, SuperLUStat_t *, int *);
extern void dlsum_fmod(double *, double *, double *, double *,
		       int, int, int_t , int_t *, int_t, int_t, int_t,
		       int_t *, gridinfo_t *, LocalLU_t *, 
		       MPI_Request [], SuperLUStat_t *);
extern void dlsum_bmod(double *, double *, double *, int, int_t, int_t *,
		       int_t *Urbs, Ucb_indptr_t **, int_t **, int_t *,
		       gridinfo_t *, LocalLU_t *,
		       MPI_Request [], SuperLUStat_t *);
extern void  pdgsrfs_ABXglobal(int_t, SuperMatrix *, double, LUstruct_t *,
			       gridinfo_t *, double [], int_t, double [],
			       int_t, int, double *, SuperLUStat_t *, int *);
/*extern void  az_gsmv_setup(SuperMatrix *, Glu_persist_t *, gridinfo_t *,
			   int [], int *[], int *[], int *[],
			   double *[], int *[], int *[],
			   int *[], int *[], int *[], int_t []);*/
extern void  get_diag_procs(int_t, Glu_persist_t *, gridinfo_t *, int_t *,
			    int_t **, int_t **);
extern void  pdgsmv_setup(SuperMatrix *, Glu_persist_t *, gridinfo_t *,
			  int_t *, int_t *[], double *[], int *[], int_t []);
extern void  pdgsmv(int_t, int_t [], double [], int [], double [], double []);
extern void  pdgsmv_abs(int_t, int_t [], double [], int [], double [], double []);
extern void  GenXtrueRHS(int, SuperMatrix *, Glu_persist_t *, gridinfo_t *,
			 double **, int *, double **, int *);
extern void  Destroy_LU(int_t, gridinfo_t *, LUstruct_t *);
extern void  *duser_malloc_dist (int_t, int_t);
extern void  duser_free_dist (int_t, int_t);
extern int_t QuerySpace_dist(int_t, int_t, Glu_freeable_t *, mem_usage_t *);
extern int_t dQuerySpace_dist(int_t, LUstruct_t *, gridinfo_t *, mem_usage_t *);
extern int   xerbla_ (char *, int *);
extern void  pxerbla (char *, gridinfo_t *, int_t);

/* Routines for debugging */
extern void  print_panel_seg_dist(int_t, int_t, int_t, int_t, int_t *, int_t *);
extern void  check_tempv(int_t, double *);
extern void  check_repfnz_dist(int_t, int_t, int_t, int_t *);
extern int_t CheckZeroDiagonal(int_t, int_t *, int_t *, int_t *);

extern void  PrintLblocks(int_t, int_t, gridinfo_t *, Glu_persist_t *,
			  LocalLU_t *);
extern void  PrintUblocks(int_t, int_t, gridinfo_t *, Glu_persist_t *,
			  LocalLU_t *);
extern void  PStatInit(SuperLUStat_t *);
extern void  PStatFree(SuperLUStat_t *);
extern void  PStatPrint(SuperLUStat_t *, gridinfo_t *);
extern void  dPrint_CompCol_Matrix_dist(SuperMatrix *);
extern void  dPrint_Dense_Matrix_dist(SuperMatrix *);
extern void  PrintDouble5(char *, int_t, double *);
extern void  PrintInt10(char *, int_t, int_t *);


  }

#endif /* __SUPERLU_dDEFS */

