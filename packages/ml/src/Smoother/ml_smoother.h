/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* and disclaimer.                                                      */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_Smoother structure                             */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLSMOOTHER__
#define __MLSMOOTHER__

/* ******************************************************************** */
/* data structure type definition                                       */
/* ******************************************************************** */
typedef enum {ML_GS_standard=0,ML_GS_symmetric,ML_GS_efficient_symmetric} ML_GS_SWEEP_TYPE;
/* 0 - Standard G-S or Block G-S
   1 - Symmetric G-S or Block G-S
   2 - Efficient Symmetric G-S or Block G-S
   (i.e. pre=forward, post=backward) */


typedef struct ML_SmootherFunc_Struct ML_SmootherFunc;
typedef struct ML_Smoother_Struct ML_Smoother;
typedef struct ML_Sm_BGS_Data_Struct ML_Sm_BGS_Data;
typedef struct ML_Sm_ILUT_Data_Struct ML_Sm_ILUT_Data;
typedef struct ML_Sm_Hiptmair_Data_Struct ML_Sm_Hiptmair_Data;
typedef struct ML_Sm_BlockHiptmair_Data_Struct ML_Sm_BlockHiptmair_Data;

/* ******************************************************************** */
/* local include files                                                  */
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_memory.h"
#include "ml_1level.h"
#include "ml_operator.h"
#include "ml_comminfoop.h"
#include "ml_csolve.h"
#include "ml_struct.h"
#include <math.h>

#ifdef WKC
#include <Epetra_MultiVector.h>
#endif

/* ******************************************************************** */
/* data definition for the ML_Smoother Class                            */
/* ******************************************************************** */
/* -------------------------------------------------------------------- */
/* These data structures define the smoother object.                    */
/* -------------------------------------------------------------------- */

struct ML_SmootherFunc_Struct
{
   int ML_id;
   int (*func_ptr)(ML_Smoother *, int, double *, int, double *);
   void *data;
};

/*******************************************************************************
   pre_or_post      flag that the smoother toggles to determine whether it is
                    in pre or post smoothing phase
   envelope         message-passing information that is used in
                    ML_exchange_bdry
*******************************************************************************/

struct ML_Smoother_Struct
{
   int                     ML_id;
   struct ML_1Level_Struct *my_level;
   int                     ntimes;
   int                     init_guess;
   double                  omega;
   double                  tol;
   ML_SmootherFunc         *smoother;
   void                    (*data_destroy)(void *);
   double                  build_time, apply_time;
   int                     times_applied;
   char                    *label;
   int                     pre_or_post;
   ML_Comm_Envelope        *envelope;
   int                     output_level;
   ML_GS_SWEEP_TYPE        gs_sweep_type;

};

struct ML_Sm_BGS_Data_Struct
{
   double ** blockfacts;
   int    ** perms;
   int    blocksize;
   int    *blocklengths;
   int    *blockmap;
   int    *blockOffset;
   int    Nblocks;
  int    optimized;
   double **trid_dl;
   double **trid_d;
   double **trid_du;
   double **trid_du2;
   int    **trid_ipiv;
};

struct ML_Sm_ILUT_Data_Struct
{
   int           Nrows;
   int           *mat_ia;
   int           *mat_ja;
   double        *mat_aa;
   ML_CommInfoOP *getrow_comm;
   int           fillin;
   double        threshold;
};
struct DinvA_widget {
  int ML_id;
  int (*func_ptr)(ML_Operator *, int, double *, int, double *);
  void *data;
  ML_Operator *Amat;
};

#ifdef out
#if defined(SUPERLU)
#include "dsp_defs.h"
#include "util.h"
#endif
#ifdef DSUPERLU
#include <mpi.h>
#include "superlu_ddefs.h"
#endif
typedef struct ML_Sm_Schwarz_Data_Struct ML_Sm_Schwarz_Data;

struct ML_Sm_Schwarz_Data_Struct
{
   int           Nrows;
   int           **bmat_ia;
   int           **bmat_ja;
   double        **bmat_aa;
   int           **aux_bmat_ia;
   int           **aux_bmat_ja;
   double        **aux_bmat_aa;
   ML_CommInfoOP *getrow_comm;
   int           nblocks;
   int           *blk_info;
   int           *blk_size;
   int           **blk_indices;
   int           **perm_r;
   int           **perm_c;
#if defined(SUPERLU)
   SuperMatrix   **slu_Amat;
   SuperMatrix   **slu_Lmat;
   SuperMatrix   **slu_Umat;
#endif
};
#endif

/*******************************************************************************
Hiptmair Smoother data structure
    sm_nodal    pointer to nodal smoother structure
    max_eig     eigenvalue calculated for damping parameter
    omega       damping parameter for edge and/or nodal smoother inside
                Hiptmair smoother
*******************************************************************************/

struct ML_Sm_Hiptmair_Data_Struct
{
   ML_Operator *Tmat;
   ML_Operator *Tmat_trans;
   ML_Operator *ATmat_trans;
   double      *TtAT_diag;
   ML_Operator *TtATmat;
   ML_Smoother *sm_nodal;
   double max_eig;
   double omega;
   double output_level;
   ML    *ml_nodal;
   ML    *ml_edge;
   int   reduced_smoother;
   int   external_TtATmat; /* Set if T^TAT is generated elsewhere */
};

#define FULL_HIPTMAIR 0
             /* smoothes on edges, nodes and then edges */

#define HALF_HIPTMAIR 1
             /* smoothes on edges then node for pre-smoothing */
             /* smoothes on nodes then edges fr post-smoothing */

/* for block hiptmair */

struct ML_Sm_BlockHiptmair_Data_Struct
{
   ML_Operator *Tmat;
   ML_Operator *Tmat_trans;
   ML_Operator *ATmat_trans;
   double      *TtAT_diag;
   ML_Operator *TtATmat;
   ML_Smoother *sm_nodal;
   double *res_edge;
   double *res_edge1;
   double *res_edge2;
   double *rhs_nodal1;
   double *rhs_nodal2;
   double *x_nodal1;
   double *x_nodal2;
   double *edge_update1;
   double *edge_update2;
   double max_eig;
   double omega;
   double output_level;
   ML    *ml_nodal;
   ML    *ml_edge;
   int   reduced_smoother;
};

/* ******************************************************************** */
/* ******************************************************************** */
/*      User Interface Proto-types                                      */
/* ******************************************************************** */
/* ******************************************************************** */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern  int ML_Smoother_Create(ML_Smoother **, ML_1Level *level_ptr );
extern  int ML_Smoother_Init(ML_Smoother *, ML_1Level *level_ptr);
extern  int ML_Smoother_Destroy(ML_Smoother **);
extern  int ML_Smoother_Clean(ML_Smoother *);
extern  int ML_Smoother_Set_Label( ML_Smoother *smoo, char *label);
extern  int ML_Smoother_Apply(ML_Smoother *,int,double *,int,double*,int);

#ifdef WKC
/* WKC Added proto for Epetra stuff! */
extern  int ML_Smoother_Apply(ML_Smoother *,int,Epetra_MultiVector &,
                              int,Epetra_MultiVector &,int);
#endif


extern  int ML_Smoother_Set(ML_Smoother *, void *,
                 int (*func_ptr)(ML_Smoother*,int,double*,int,double *),
                 int, double, char *);
extern  int ML_Smoother_Jacobi(ML_Smoother *, int, double *x, int, double *);
extern  int ML_Smoother_GaussSeidel(ML_Smoother *, int, double *, int, double *);
extern  int ML_Smoother_SGSSequential(ML_Smoother *, int,double *, int, double *);
extern  int ML_Smoother_SGS(ML_Smoother *, int, double *, int, double *);
extern  int ML_Smoother_BlockGS(ML_Smoother *, int, double *, int, double *);
extern int ML_Smoother_NewGS(ML_Smoother *sm,int inlen,double x[],int outlen,
                        double rhs[]);

extern int ML_BlockScaledApply(ML_Operator *Amat, int inlen, double in[],
			       int outlen, double out[]);

extern  int ML_Smoother_Cheby_Apply(ML_Smoother *, int, double *, int, double *);
extern int ML_Cheby(ML_Smoother *sm, int inlen, double x[], int outlen, double rhs[]);

#ifdef WKC
/* WKC -- double * are actually Epetra_MultiVectors */
/*extern int ML_Cheby_WKC(ML_Smoother *sm, int inlen, double *x, int outlen, double *rhs);*/
extern int ML_Cheby_WKC(void *sm, int inlen, double *x, int outlen, double *rhs);
#endif

extern int ML_complex_Cheby(ML_Smoother *sm, int inlen, double x[], int outlen,
			    double rhs[]);
extern int ML_DiagScaled_1stepKrylov(ML_Smoother *sm, int inlen, double x[],
				     int outlen, double rhs[]);
extern  int ML_Smoother_ParaSails(ML_Smoother *, int, double *, int, double *);
extern  int ML_Smoother_ParaSailsSym(ML_Smoother *, int, double *, int, double *);
extern  int ML_Smoother_ParaSailsTrans(ML_Smoother *, int, double *, int, double *);
extern  int ML_Smoother_VBlockJacobi(ML_Smoother *,int,double *x,int, double *);
extern  int ML_Smoother_LineJacobi(ML_Smoother *,int,double *x,int, double *);
extern  int ML_Smoother_LineGS(ML_Smoother *,int,double *x,int, double *);
extern  int ML_Smoother_VBlockKrylovJacobi(ML_Smoother *,int,double*,int,double*);
extern  int ML_Smoother_VBlockSGS(ML_Smoother *, int, double *x, int, double *);
extern  int ML_Smoother_VBlockSGSSequential(ML_Smoother*,int,double*,int,double*);
extern  int ML_Smoother_OverlappedILUT(ML_Smoother *,int,double *x,int,double *);
extern  int ML_Smoother_VBlockAdditiveSchwarz(ML_Smoother *,int,double*,int,double*);
extern  int ML_Smoother_VBlockMultiplicativeSchwarz(ML_Smoother *,int,double*,int,double*);
extern  int ML_Smoother_Hiptmair(ML_Smoother *, int, double *, int, double *);
extern  int ML_Smoother_BlockHiptmair(ML_Smoother *, int, double *, int, double *);
extern int ML_Smoother_ApplySubdomainOverlap(ML_Smoother *sm, int inlen,
					    double x[],int outlen, double b[]);
#ifdef HAVE_ML_PETSC
extern int ML_Smoother_Petsc(ML_Smoother *sm, int inlen, double x[], int outlen,
                      double rhs[]);
#endif

extern void ML_Smoother_DestroySubdomainOverlap(void *data);

extern int ML_EyeMinusIterationOperator_Matvec(ML_Operator *Amat, int ilen,
		       double p[], int olen, double ap[]);

extern int ML_Smoother_ComputeOmegaViaSpectralradius(ML_Operator *Amat,
    int (*smoothing_function)(ML_Smoother *, int, double *, int, double *),
    void *data, double *spectral_radius, double *omega);


/* ******************************************************************** */
/* ******************************************************************** */
/* private functions                                                   */
/* ******************************************************************** */
/* ******************************************************************** */

extern  int ML_Smoother_Create_Hiptmair_Data(ML_Sm_Hiptmair_Data **data);
extern  int ML_Smoother_Create_BlockHiptmair_Data(ML_Sm_BlockHiptmair_Data **data);
extern  int ML_Smoother_Gen_Hiptmair_Data(ML_Sm_Hiptmair_Data**,
                         ML_Operator*, ML_Operator*, ML_Operator*,
                         ML_Operator*, ML_Operator*, ML_Operator *,
                         int, int*, void *, void **, void *, void **);
extern  int ML_Smoother_Gen_BlockHiptmair_Data(ML_Sm_BlockHiptmair_Data**,
                         ML_Operator*, ML_Operator*, ML_Operator*,
                         ML_Operator*, int, int*, void *, void **,
					  void *, void **);
extern int ML_Smoother_HiptmairSubsmoother_Create(ML **ml_subproblem,
					   ML_Operator *Amat, void *smoother,
						  void **args, double default_omega);

extern void ML_Smoother_Destroy_Hiptmair_Data(void *data);
extern void ML_Smoother_Destroy_BlockHiptmair_Data(void *data);
extern  int ML_Smoother_Create_BGS_Data(ML_Sm_BGS_Data **data);
extern void ML_Smoother_Destroy_BGS_Data(void *data);
extern void ML_Smoother_Clean_BGS_Data(void *data);
extern  int ML_Smoother_Create_ILUT_Data(ML_Sm_ILUT_Data **data);
extern void ML_Smoother_Destroy_ILUT_Data(void *data);
extern  int ML_Smoother_Gen_BGSFacts(ML_Sm_BGS_Data **, ML_Operator *,int);
extern  int ML_Smoother_Gen_VBGSFacts(ML_Sm_BGS_Data**,ML_Operator*,int,int*);
extern  int ML_Smoother_Gen_LineSmootherFacts(ML_Sm_BGS_Data**, ML_Operator*, int, int*, int*);
extern void ML_Smoother_Destroy_Schwarz_Data(void *data);
extern void ML_Smoother_Clean_ParaSails(void *data);
extern struct MLSthing *ML_Smoother_Create_MLS(void);
extern int ML_BlockDinv(ML_Sm_BGS_Data *BGS_Data, int outlen, double out[]);

extern void ML_Smoother_Destroy_MLS(void *data);
extern void **ML_Smoother_Arglist_Create(int nargs);
extern int ML_Smoother_Arglist_Set(void **arglist, int which_arg, void *ptr);
extern void *ML_Smoother_Arglist_Get(void **arglist, int which_arg);
extern int ML_Smoother_Arglist_Delete(void ***arglist);
extern int ML_Smoother_Arglist_Nargs(void **arglist);



extern  int ML_Smoother_ILUTDecomposition(ML_Sm_ILUT_Data *, ML_Operator *,
                    ML_Comm *, int, int *,int*,double *,int *, int *,int);
#ifdef out
extern  int ML_Smoother_Create_Schwarz_Data(ML_Sm_Schwarz_Data **data);
extern  int ML_Smoother_VBlockSchwarzDecomposition(ML_Sm_Schwarz_Data *,
                    ML_Operator *, ML_Comm *, int, int *,int*,double *,int *,
                    int *,int);
#endif

extern  int ML_Smoother_GetOffProcRows(ML_CommInfoOP *, ML_Comm *,
                  ML_Operator *,int,int *,int,int *,int **,double **);
extern  int ML_Smoother_GetRowLengths(ML_CommInfoOP *, ML_Comm *,
                  ML_Operator *, int *, int **);
extern  int ML_Smoother_ComposeOverlappedMatrix(ML_Operator *, ML_Comm *,
                  int *, int **, int **, double **, int **, int **, int *);
extern  ML *ML_Smoother_Get_Hiptmair_nodal(ML *ml, int level, int);
extern  int ML_dgetrs_special(int blocksize, double *ablock, int *ipiv,
			      double *correc);
extern  int ML_dgetrs_trans_special(int blocksize, double *ablock, int *ipiv,
			      double *correc);
extern  int ML_permute_for_dgetrs_special(double *Z[], int Nblocks, int blocksize,
                               ML_Sm_BGS_Data *block_data_widget);


/* -------------------------------------------------------------------- */
/* Ray's functions                                                      */
/* -------------------------------------------------------------------- */

extern  int ML_MSR_SGSextra(ML_Smoother *, int , double *, int , double *);
extern void ML_MSR_GSextra_Clean(void *data);
extern  int ML_Smoother_BackGS(void *, int, double *, int, double *);
extern void ML_Smoother_Clean_OrderedSGS(void *data);
extern  int ML_Smoother_Gen_Ordering(ML_Operator *Amat, int **data_ptr);
extern  int ML_Smoother_OrderedSGS(ML_Smoother *sm,int inlen,double x[],int outlen,
                                   double rhs[]);
extern  int ML_Smoother_MSR_SGS(ML_Smoother *, int, double *, int, double *);
extern int ML_Smoother_MSR_SGSnodamping(ML_Smoother *,int ,double *,int , double *);
extern int ML_Smoother_MSR_GSforwardnodamping(void *sm,int inlen,double x[],
					      int outlen, double rhs[]);
extern int ML_Smoother_MSR_GSbackwardnodamping(void *sm,int inlen,double x[],
					       int outlen, double rhs[]);
extern int ML_Smoother_MSR_SGSdamping(void *,int ,double *,int , double *);
extern void ML_Smoother_Clean_MSR_GS(void *data);

extern int DinvA(ML_Operator *data,  int in, double p[], int out, double ap[]);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#endif

