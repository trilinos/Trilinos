

/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 */

#ifndef __SUPERLU_dDEFS /* allow multiple inclusions */
#define __SUPERLU_dDEFS

/*
 * File name:	superlu_ddefs.h
 * Purpose:     Distributed SuperLU data types and function prototypes
 * History:
 */

#include "superlu_defs.h"

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

/*-- Data structure for communication during matrix-vector multiplication. */
typedef struct {
    int_t *extern_start;
    int_t *ind_tosend;    /* X indeices to be sent to other processes */
    int_t *ind_torecv;    /* X indeices to be received from other processes */
    int_t *ptr_ind_tosend;/* Printers to ind_tosend[] (Size procs)
			     (also point to val_torecv) */
    int_t *ptr_ind_torecv;/* Printers to ind_torecv[] (Size procs)
			     (also point to val_tosend) */
    int   *SendCounts;    /* Numbers of X indices to be sent
			     (also numbers of X values to be received) */
    int   *RecvCounts;    /* Numbers of X indices to be received
			     (also numbers of X values to be sent) */
    double *val_tosend;   /* X values to be sent to other processes */
    double *val_torecv;   /* X values to be received from other processes */
    int_t TotalIndSend;   /* Total number of indices to be sent
			     (also total number of values to be received) */
    int_t TotalValSend;   /* Total number of values to be sent.
			     (also total number of indices to be received) */
} pdgsmv_comm_t;

/*-- Data structure for redistribution of B and X --*/
typedef struct {
    int  *B_to_X_SendCnt;
    int  *X_to_B_SendCnt;
    int  *ptr_to_ibuf, *ptr_to_dbuf;
} pxgstrs_comm_t;

/*-- Data structure holding the information for the solution phase --*/
typedef struct {
    int_t *row_to_proc;
    int_t *inv_perm_c;
    int_t num_diag_procs, *diag_procs, *diag_len;
    pdgsmv_comm_t *gsmv_comm;
    pxgstrs_comm_t *gstrs_comm;
    int_t *A_colind_gsmv; /* After p@(pre)gsmv_init(), the global column
			     indices of A are translated into the relative
			     positions in the gathered x-vector.
			     This is re-used in repeated calls to p@(pre)gsmv() */
} SOLVEstruct_t;


/***********************************************************************
 * Function prototypes
 ***********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif


/* Supernodal LU factor related */
extern void
dCreate_CompCol_Matrix_dist(SuperMatrix *, int_t, int_t, int_t, double *,
			    int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompRowLoc_Matrix_dist(SuperMatrix *, int_t, int_t, int_t, int_t,
			       int_t, double *, int_t *, int_t *,
			       Stype_t, Dtype_t, Mtype_t);
extern void
dCompRow_to_CompCol_dist(int_t, int_t, int_t, double *, int_t *, int_t *,
                         double **, int_t **, int_t **);
extern int
dCompRow_loc_to_CompCol_global(int_t, SuperMatrix *, gridinfo_t *,
			       SuperMatrix *);
extern void
dCopy_CompCol_Matrix_dist(SuperMatrix *, SuperMatrix *);
extern void
dCreate_Dense_Matrix_dist(SuperMatrix *, int_t, int_t, double *, int_t,
			  Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix_dist(SuperMatrix *, int_t, int_t, int_t, double *, 
			      int_t *, int_t *, int_t *, int_t *, int_t *,
			      Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_Dense_Matrix_dist(int_t, int_t, double *, int_t,
                        double *, int_t);

extern void    dallocateA_dist (int_t, int_t, double **, int_t **, int_t **);
extern void    dGenXtrue_dist (int_t, int_t, double *, int_t);
extern void    dFillRHS_dist (char *, int_t, double *, int_t,
                              SuperMatrix *, double *, int_t);
extern int     dcreate_matrix(SuperMatrix *, int, double **, int *, 
			      double **, int *, FILE *, gridinfo_t *);

/* Driver related */
extern void    dgsequ_dist (SuperMatrix *, double *, double *, double *,
			    double *, double *, int_t *);
extern double  dlangs_dist (char *, SuperMatrix *);
extern void    dlaqgs_dist (SuperMatrix *, double *, double *, double,
			    double, double, char *);
extern void    pdgsequ (SuperMatrix *, double *, double *, double *,
			double *, double *, int_t *, gridinfo_t *);
extern double  pdlangs (char *, SuperMatrix *, gridinfo_t *);
extern void    pdlaqgs (SuperMatrix *, double *, double *, double,
			double, double, char *);
extern int     pdPermute_Dense_Matrix(int_t, int_t, int_t [], int_t[],
				      double [], int, double [], int, int,
				      gridinfo_t *);

extern int     sp_dtrsv_dist (char *, char *, char *, SuperMatrix *,
			      SuperMatrix *, double *, int *);
extern int     sp_dgemv_dist (char *, double, SuperMatrix *, double *,
			      int, double, double *, int);
extern int     sp_dgemm (char *, char *, int, int, int, double,
			SuperMatrix *, double *, int, double, 
			double *, int);

extern int_t ddistribute(fact_t, int_t, SuperMatrix *, Glu_freeable_t *, 
			 LUstruct_t *, gridinfo_t *);
extern void  pdgssvx_ABglobal(superlu_options_t *, SuperMatrix *, 
			      ScalePermstruct_t *, double *,
			      int, int, gridinfo_t *, LUstruct_t *, double *,
			      SuperLUStat_t *, int *);
extern int_t pddistribute(fact_t, int_t, SuperMatrix *, 
			 ScalePermstruct_t *, Glu_freeable_t *, 
			 LUstruct_t *, gridinfo_t *);
extern void  pdgssvx(superlu_options_t *, SuperMatrix *, 
		     ScalePermstruct_t *, double *,
		     int, int, gridinfo_t *, LUstruct_t *,
		     SOLVEstruct_t *, double *, SuperLUStat_t *, int *);
extern int  dSolveInit(superlu_options_t *, SuperMatrix *, int_t [], int_t [],
		       int_t, LUstruct_t *, gridinfo_t *, SOLVEstruct_t *);
extern int_t pxgstrs_init(int_t, int_t, int_t, int_t,
	                  int_t [], int_t [], gridinfo_t *grid,
	                  Glu_persist_t *, SOLVEstruct_t *);
extern void pxgstrs_finalize(pxgstrs_comm_t *);
extern void dSolveFinalize(superlu_options_t *, SOLVEstruct_t *);
extern void dldperm(int_t, int_t, int_t, int_t [], int_t [],
		    double [], int_t *, double [], double []);
extern void pdgstrf(superlu_options_t *, int, int, double,
		    LUstruct_t*, gridinfo_t*, SuperLUStat_t*, int*);
extern void pdgstrs_Bglobal(int_t, LUstruct_t *, gridinfo_t *,
			     double *, int_t, int, SuperLUStat_t *, int *);
extern void pdgstrs(int_t, LUstruct_t *, ScalePermstruct_t *, gridinfo_t *,
		    double *, int_t, int_t, int_t, int, SOLVEstruct_t *,
		    SuperLUStat_t *, int *);
extern void dlsum_fmod(double *, double *, double *, double *,
		       int, int, int_t , int_t *, int_t, int_t, int_t,
		       int_t *, gridinfo_t *, LocalLU_t *, 
		       MPI_Request [], SuperLUStat_t *);
extern void dlsum_bmod(double *, double *, double *,
                       int, int_t, int_t *, int_t *, Ucb_indptr_t **,
                       int_t **, int_t *, gridinfo_t *, LocalLU_t *,
		       MPI_Request [], SuperLUStat_t *);
extern void pdgsrfs(int_t, SuperMatrix *, double, LUstruct_t *,
		    ScalePermstruct_t *, gridinfo_t *,
		    double [], int_t, double [], int_t, int,
		    SOLVEstruct_t *, double *, SuperLUStat_t *, int *);
extern int   pdgsmv_AXglobal_setup(SuperMatrix *, Glu_persist_t *,
				   gridinfo_t *, int_t *, int_t *[],
				   double *[], int_t *[], int_t []);
extern int  pdgsmv_AXglobal(int_t, int_t [], double [], int_t [],
	                       double [], double []);
extern int  pdgsmv_AXglobal_abs(int_t, int_t [], double [], int_t [],
				 double [], double []);
extern void pdgsmv_init(SuperMatrix *, int_t *, gridinfo_t *,
			pdgsmv_comm_t *);
extern void pdgsmv(int_t, SuperMatrix *, gridinfo_t *, pdgsmv_comm_t *,
		   double x[], double ax[]);

/* Memory-related */
extern double  *doubleMalloc_dist(int_t);
extern double  *doubleCalloc_dist(int_t);
extern void  *duser_malloc_dist (int_t, int_t);
extern void  duser_free_dist (int_t, int_t);
extern int_t dQuerySpace_dist(int_t, LUstruct_t *, gridinfo_t *, mem_usage_t *);
extern void    Destroy_LU(int_t, gridinfo_t *, LUstruct_t *);
extern void    LUstructInit(const int_t, const int_t, LUstruct_t *);
extern void    LUstructFree(LUstruct_t *);

/* Auxiliary routines */
extern void    dfill_dist (double *, int_t, double);
extern void    dinf_norm_error_dist (int_t, int_t, double*, int_t,
                                     double*, int_t, gridinfo_t*);
extern void    pdinf_norm_error(int, int_t, int_t, double [], int_t,
				double [], int_t , gridinfo_t *);
extern void  dreadhb_dist (int, FILE *, int_t *, int_t *, int_t *, 
			   double **, int_t **, int_t **);

/* Routines for debugging */
extern void  dPrintLblocks(int_t, int_t, gridinfo_t *, Glu_persist_t *,
		 	   LocalLU_t *);
extern void  dPrintUblocks(int_t, int_t, gridinfo_t *, Glu_persist_t *,
			   LocalLU_t *);
extern void  dPrint_CompCol_Matrix_dist(SuperMatrix *);
extern void  dPrint_Dense_Matrix_dist(SuperMatrix *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_dDEFS */

