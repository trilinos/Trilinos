/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/*
 * Include file for inclusion in any routine which will call the solver
 * library. Contains necessary constants and prototypes.
 *
 * Author:  Scott A. Hutchinson, SNL
 *          John  N. Shadid,     SNL
 *          Ray   S. Tuminaro,   SNL
 */
#ifndef __AZTECH__

/* Set variable to indicate that this file has already been included */

#define __AZTECH__


/* Some older codes use AZ_MPI to set MPI mode for AztecOO/Aztec.
 * Check to see if AZ_MPI is defined, and define AZTEC_MPI if
 * it is not already defined.
 */

#if defined(AZ_MPI) && !defined(AZTEC_MPI) 
#define AZTEC_MPI
#endif

/* Force AZTEC_MPI to be defined if ML_MPI is defined */

#if defined(ML_MPI) && !defined(AZTEC_MPI) 
#define AZTEC_MPI
#endif

/* The definition of MPI_AZRequest and MPI_AZComm depend on
 * whether or not we are using MPI.
 * NOTE:  This technique can cause problems if az_aztec.h (this file)
 *        is included with AZTEC_MPI undefined in some files and with
 *        AZTEC_MPI defined in other files.  Therefore, users must
 *        make sure that AZTEC_MPI is either defined or undefined for
 *        all files that are compiled and including az_aztec.h.
 */

#ifdef AZTEC_MPI
#include <mpi.h>
#define MPI_AZRequest MPI_Request
#define MPI_AZComm    MPI_Comm
#else
#define MPI_AZRequest int
#define MPI_AZComm    int 
#endif


/*structure definitions*/


struct AZ_MATRIX_STRUCT {

             /* Used to represent matrices. In particular, two structures are  */
             /* passed into AZ_iterate:                                      */
             /*    AZ_iterate(..., AZ_MATRIX *Amat, AZ_MATRIX *Precond, ...) */
             /* corresponding to matrix-vector products and preconditioners. */
             /*                                                              */
             /* For matrix-vector products, a subroutine Amat.'user_function'*/
             /* can be supplied. 'Amat' is be passed to this routine and thus*/
             /* relevant data can be placed in this structure. If a matrix-  */
             /* vector product is not supplied, either an MSR or VBR matrix  */
             /* must be used by specifying the arrays bindx,indx,rpntr,cpntr,*/
             /* bpntr, and val. In this case, Aztec supplies the matrix-     */
             /* vector product.                                              */
             /*                                                              */
             /* NOTE: Fortran users never explicitly see this structure but  */
             /* instead pass matrix-vector product and preconditioning       */
             /* information through parameters which Aztec copies to this    */
             /* structure.                                                   */

   int              matrix_type;  /* Indicates whether the matrix is MSR,   */
                                  /* VBR, or user-supplied.                 */
                                  /*                                        */
   int              N_local,      /* Number of local and ghost unknowns     */
                    N_ghost;      /*                                        */
                                  /*                                        */
  int           mat_create_called;/* =1 indicates that AZ_matrix_create()   */
                                  /* was invoked.                           */
  int          must_free_data_org;/* =1 indicates that data_org was created */
                                  /* via the matrix_free set functions and  */
                                  /* needs to be freed during destroy oper. */
                                  /*                                        */
   int              *rpntr,*cpntr,/* arrays to support MSR & VBR formats    */
                    *bpntr,*bindx;/*                                        */ 
   int              *indx;        /*                                        */
   double *val;                   /*                                        */
   int              *data_org;    /* array to support matvec communication  */
                                  /*                                        */
/* Begin Aztec 2.1 mheroux mod */
   int              N_update;     /* Number of nodes updated on this proc   */
                                  /*                                        */
   int              *update;      /* array containing global indices map    */
                                  /*                                        */
   int       has_global_indices;  /* true/false for say bindx has global    */
                                  /*                                        */
/* End Aztec 2.1 mheroux mod */
   void (*matvec)(double *,       /* function ptr to user-defined routine   */
                  double *, struct 
                  AZ_MATRIX_STRUCT *, 
                  int *); 
/*********************************************************************/
/*********************************************************************/
   int (*getrow)(int columns[], double values[], int row_lengths[],
		 struct AZ_MATRIX_STRUCT *Amat, int N_requested,
                 int requested_rows[], int allocated_space);
                                  /* function ptr to user-defined routine   */
/* Get some matrix rows ( requested_rows[0 ... N_requested_rows-1] ) */
/* from the user's matrix and return this information  in            */
/* 'row_lengths, columns, values'.  If there is not enough space to  */
/* complete this operation, return 0.  Otherwise, return 1.          */
/*                                                                   */
/* Parameters                                                        */
/* ==========                                                        */
/* data             On input, points to user's data containing       */
/*                  matrix values.                                   */
/* N_requested_rows On input, number of rows for which nonzero are   */
/*                  to be returned.                                  */
/* requested_rows   On input, requested_rows[0...N_requested_rows-1] */
/*                  give the row indices of the rows for which       */
/*                  nonzero values are returned.                     */
/* row_lengths      On output, row_lengths[i] is the number of       */
/*                  nonzeros in the row 'requested_rows[i]'          */
/*                  ( 0 <= i < N_requested_rows). NOTE: this         */
/*                  array is of size 'N_requested_rows'.             */
/* columns,values   On output, columns[k] and values[k] contains the */
/*                  column number and value of a matrix nonzero where*/
/*                  all the nonzeros for requested_rows[0] appear    */
/*                  first followed by the nonzeros for               */
/*                  requested_rows[1], etc. NOTE: these arrays are   */
/*                  of size 'allocated_space'.                       */
/* allocated_space  On input, indicates the space available in       */
/*                  'columns' and 'values' for storing nonzeros. If  */
/*                  more space is needed, return 0.                  */
/*********************************************************************/
/*********************************************************************/
   int  (*user_comm)(double *, struct AZ_MATRIX_STRUCT *);
                                  /* user communication routine before */
                                  /* doing matvecs. Only used when doing    */
                                  /* matrix-free.                           */ 

   double matrix_norm;            /* norm of the matrix A used in the case  */
                                  /* of least square preconditioning if the */
                                  /* matrix A is of type AZ_USER_MATRIX */
                                  /*                                        */
                                  /*                                        */
   int              **aux_ival;   /* integer, double precision, function,   */
   double           **aux_dval;   /* generic, and matrix pointers at the    */
   void              *aux_ptr;    /* product routine: 'matvec()'.           */
   void              *matvec_data;
   void              *getrow_data;
   struct AZ_MATRIX_STRUCT        
                    **aux_matrix; 
   int              N_nz, max_per_row, /* Total number of nonzeros, maximum */
		    largest_band;      /* nonzeros per row, and bandwidth.  */
                                       /* ONLY used for matrix-free         */
};

struct grid_level {
          int                     N;
          struct AZ_MATRIX_STRUCT *transfer_to_prev_grid;
          struct AZ_MATRIX_STRUCT *transfer_to_next_grid;
          struct AZ_MATRIX_STRUCT *discretization_op;
          struct AZ_PREC_STRUCT *smoother1;
          struct AZ_PREC_STRUCT *smoother2;
          void                    *mesh;
       };

struct AZ_PREC_STRUCT {

             /* Used to represent preconditioners. In particular,            */
             /* two structures  are                                          */
             /* passed into AZ_iterate:                                      */
             /* AZ_iterate(..., AZ_MATRIX *Amat, AZ_PRECOND *Precond, ...) */
             /* corresponding to matrix and preconditioner descriptions. */
             /*                                                              */
             /* For matrix-vector products, a subroutine Amat.'matvec'*/
             /* can be supplied. 'Amat' is be passed to this routine and thus*/
             /* relevant data can be placed in this structure. If a matrix-  */
             /* vector product is not supplied, either an MSR or VBR matrix  */
             /* must be used by specifying the arrays bindx,indx,rpntr,cpntr,*/
             /* bpntr, and val. In this case, Aztec supplies the matrix-     */
             /* vector product as well as a number of preconditioners.       */
             /*                                                              */
             /* Likewise, a preconditioner can be supplied via the routine   */
             /* 'Precond.prec_function'. In this case options[AZ_precond]    */
             /* must be set to "AZ_user_precond". Otherwise                  */
             /* options[AZ_precond] must be set to one of the preconditioners*/
             /* supplied by Aztec and the matrix must be a MSR or VBR format */
             /* The matrix used as preconditionner is descibed in a AZ_MATRIX*/
             /* structure which could be either the same as Amat             */
             /* (precond.Pmat = Amat) or  a different matrix described       */
             /* by the arrays bindx,indx,rpntr,cpntr, bpntr, and val.        */
             /*                                                              */
             /* NOTE: Fortran users never explicitly see these structures but*/
             /* instead pass matrix and preconditioning  information through */
             /* parameters which Aztec copies to this  structure.            */
             /*                                                              */

   struct AZ_MATRIX_STRUCT *Pmat;     /* matrix used by the precondtioner  */
                                      /* when not using multilevel stuff   */
                                      /*                                   */
  int           prec_create_called;/* =1 indicates that AZ_precond_create() */
                                  /* was invoked.                           */

   void    (*prec_function)(double *, /* function ptr to user-defined      */
               int *, int *, double *,/* preconditioning routine           */ 
               struct AZ_MATRIX_STRUCT  *, 
               struct AZ_PREC_STRUCT *);

   int                   *options;    /* used to determine preconditioner  */
   double                *params;     /* when options[AZ_precond] is set   */
   struct AZ_PREC_STRUCT *next_prec;  /* to AZ_multilevel. The series of   */
                                      /* preconditioners is done in a      */
                                      /* multiplicative fashion.           */
   struct context        *context;
   struct grid_level     grid_levels[10]; /* multilevel stuff                 */
   void *ml_ptr;         /* MLDIFF */
   double timing[2];     /* preconditioner timing array */
   void *precond_data;
   char *print_string;
};

   
typedef struct AZ_MATRIX_STRUCT AZ_MATRIX;
typedef struct AZ_PREC_STRUCT   AZ_PRECOND;

struct AZ_CONVERGE_STRUCT {
   double r0_norm, A_norm, b_norm;
   int    total_N;
   int    not_initialized;
   struct AZ_SCALING *scaling;
};
   





struct aztec_choices {
   int *options;
   double *params;
};
struct context {                       /* This structure is used to  */
   int      *iu, *iflag, *ha, *ipvt;   /* hold variables specific to */
   int      *dblock, *space_holder;    /* the preconditioner */
   int      extra_fact_nz_per_row;  
   int      N_large_int_arrays, N_large_dbl_arrays;
   int      N_nz_factors,N_nz_matrix, N_blk_rows, max_row;
   double   *pivot;
   struct   AZ_MATRIX_STRUCT     *A_overlapped;
   struct   aztec_choices        *aztec_choices;
   double   *x_pad, *ext_vals, *x_reord;
   int      *padded_data_org, *map, *ordering, *inv_ordering;
   int      N, N_unpadded, N_nz, N_nz_allocated;
   char     *tag;
   int      *proc_config;
   int      Pmat_computed;                 /* indicates that the has    */
                                           /* been called at least once */
                                           /* before with this context. */
/* Begin Aztec 2.1 mheroux mod */
   void    *precon;
/* End Aztec 2.1 mheroux mod */
};

 			    
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

struct AZ_SCALING {  /* Left and right matrices to to scale    */
                     /* the problem                            */
   int    action;
   double A_norm;
   int    mat_name;
   int    scaling_opt;
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

struct grid {                      /* used to define a grid. Still under */
                                   /* construction                       */
   int    *element_vertex_lists;
   int    *Nvertices_per_element;
   int    Nelements;
   int    Nvertices;
   double *vertices;
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


/* Aztec's previous AZ_solve() is renamed to AZ_oldsolve() with 3 new  */
/* parameters appended to it: Amat, precond, scaling.                  */
/* This routine is never called directly by an application. It is only */
/* used internally by Aztec.                                           */

#ifdef __cplusplus
extern "C" {
#endif
extern void AZ_oldsolve(double x[], double b[], int options[], double params[],
	double status[], int proc_config[], AZ_MATRIX *Amat, 
	AZ_PRECOND *precond, struct AZ_SCALING *scaling);
#ifdef __cplusplus
}
#endif




/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


/* This is the new Aztec interface. This routine calls AZ_oldsolve() passing */
/* in for example Amat->indx for the indx[] parameter in AZ_oldsolve().      */
/*                                                                           */
/* NOTE: User's can still invoke AZ_solve() in the old Aztec way. AZ_solve   */
/*       also calls AZ_oldsolve(). However, matrix-free and coarse grid      */
/*       capabilities are not available via AZ_solve().                      */

#ifdef __cplusplus
extern "C" {
#endif
extern void AZ_iterate(double x[], double b[], int options[], double params[],
                 double status[], int proc_config[],
                 AZ_MATRIX *Amat, AZ_PRECOND *precond, struct AZ_SCALING *scaling);
#ifdef __cplusplus
}
#endif

#ifdef next_release
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* This is the new fortran interface to Aztec. There is a wrapper so that  */
/* fortran user's can invoke AZ_iterate (like 'C' users), however the      */
/* parameters are different as the information that is captured in the     */
/* structures Amat, precond, and scaling must now be passed as parameters. */
/*                                                                         */
/* Note: Multilevel stuff is currently not provided to Fortran users.      */

extern void AZ_fortransolve(double x[], double b[], int options[], 
   double params[], int data_org[], double status[], int proc_config[],

   int indx[],  int bindx[],       /* VBR & MSR arrays. Note: all of these   */
   int rpntr[], int cpntr[],       /* arrays are passed to 'user_Avec' and   */
   int bpntr[], double val[],      /* to 'user_precond' if supplied.         */
                                   /*                                        */
   AZ_FUNCTION_PTR user_Avec,      /* user's matrix-free matvec              */
                                   /* If doing matrix-free, the following    */
                                   /* arrays and functions are passed to the */
				   /* user's matvec:                         */
   int    A_ival0[],       int    A_ival1[],
   int    A_ival2[],       int    A_ival3[],
   double A_dval0[],       double A_dval1[],
   double A_dval2[],       double A_dval3[],
   AZ_FUNCTION_PTR A_fun0, AZ_FUNCTION_PTR A_fun1, 
   AZ_FUNCTION_PTR A_fun2, AZ_FUNCTION_PTR A_fun3, 

   AZ_FUNCTION_PTR user_precond,   /* user's preconditioning routine         */
                                   /*                                        */
                                   /* The following arrays and functions are */
                                   /* passed to the user's preconditioner:   */
   int    M_ival0[],       int    M_ival1[],
   int    M_ival2[],       int    M_ival3[],
   double M_dval0[],       double M_dval1[],
   double M_dval2[],       double M_dval3[],
   AZ_FUNCTION_PTR M_fun0, AZ_FUNCTION_PTR M_fun1, 
   AZ_FUNCTION_PTR M_fun2, AZ_FUNCTION_PTR M_fun3, 

                                   /* Vectors used to scale the problem.     */

   double left_scale[], double right_scale[]);
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************

 The C user's matrix-vector product must look like this:

   'Amat.user_function'(double x[], double b[], int options[], double params[], 
                        AZ_MATRIX *Amat, int proc_config[])

 where on output b = A * x.  The user can put what he wants in 'Amat'
 when he calls AZ_iterate() so that he can use it inside this function.

 The C user's preconditioner must look like this:

   'Prec.user_precond'(double x[], int options[], double params[], 
                       AZ_MATRIX *Prec, int proc_config[])

 where on output x = M * x.  The user can put what he wants in 'Prec'
 when he calls AZ_iterate() so that he can use it inside this function.



 The Fortran user's matrix-vector product must look like this:

     'user_matvec'(x, b, options, params, data_org, proc_config,
       indx, bindx, rpntr, cpntr, bpntr,  val,                              
       A_ival0, A_ivals1, A_ival2, A_ival3,       
       A_dval0, A_dvals1, A_dval2, A_dval3,      
       A_fun0,  A_fun1,   A_fun2,  A_fun3)     

  The user can put what he wants in the integer arrays (A_ival*),
  double precision arrays (A_dval*) and functions (A_fun*) when invoking
  Aztec so that he can use them inside this function. 

  NOTE: Additionally, if not using MSR or VBR matrices, the arrays
  indx,bindx,rpntr,cpntr,bpntr, and val can also be used to pass information.

  
 The Fortran user's preconditioner must look like this:

     'user_precond'(x, options, params, data_org, proc_config,
       indx, bindx, rpntr, cpntr, bpntr,  val,                              
       M_ival0, M_ivals1, M_ival2, M_ival3,       
       M_dval0, M_dvals1, M_dval2, M_dval3,      
       M_fun0,  M_fun1,   M_fun2,  M_fun3)     

  The user can put what he wants in the integer arrays (M_ival*),
  double precision arrays (M_dval*) and functions (M_fun*) when invoking
  Aztec so that he can use them inside this function. 

  NOTE: Additionally, if not using MSR or VBR matrices, the arrays
  indx,bindx,rpntr,cpntr,bpntr, and val can also be used to pass information.

*/


/* Finally, some crude code that defines a function mymatvec() and 
   myprecond() each of which use the function myfun() as well as the
   vbr arrays to generate a result. 

   C:
     Amat.bindx         = bindx;
     Amat.indx          = indx;
     Amat.rpntr         = rpntr;
     Amat.cpntr         = cpntr;
     Amat.bpntr         = bpntr;
     Amat.val           = val;
     Amat.user_function = mymatvec;
     Amat.aux_funs      = (AZ_FUNCTION_PTR *) calloc(1,sizeof(AZ_FUNCTION_PTR));
     Amat.aux_funs[0]   = myfun;
 
     Prec.bindx         = bindx;
     Prec.indx          = indx;
     Prec.rpntr         = rpntr;
     Prec.cpntr         = cpntr;
     Prec.bpntr         = bpntr;
     Prec.val           = val;
     Prec.user_function = myprecond;
     Prec.aux_funs      = (AZ_FUNCTION_PTR *) calloc(1,sizeof(AZ_FUNCTION_PTR));
     Prec.aux_funs[0]   = myfun;
 
     AZ_iterate(x, ax, options, params, status, proc_config,
                &Amat, &Prec, &Scaling);

   Fortran: 

     call AZ_iterate(x,b, options, params, data_org, status,
    $           proc_config, NULL,bindx,NULL,NULL,NULL, val,
    $           mymatvec,  NULL, NULL, NULL, NULL, NULL, NULL,
    $           NULL, NULL, myfun, NULL, NULL, NULL,
    $           myprecond, NULL, NULL, NULL, NULL, NULL, NULL,
    $           NULL, NULL, myfun, NULL, NULL, NULL,
    $           NULL, NULL)

*/
/* constants */


/* function definitions */

#ifndef __cplusplus
#ifndef max
#define max(x,y) (( (x) > (y) ) ?  (x) : (y))     /* max function  */
#endif
#ifndef min
#define min(x,y) (( (x) < (y) ) ?  (x) : (y))     /* min function  */
#endif
#ifndef sgn
#define sgn(x)   (( (x) < 0.0 ) ? -1.0 : 1.0)  /* sign function */
#endif
#endif

/*
 * There are different conventions for external names for fortran subroutines.
 * In addition, different compilers return differing caluse for a fortran
 * subroutine call. In this section we take these into account.
 */

#if   defined(caps)
#ifdef T3E
#include <fortran.h>
#endif
#ifdef T3E
#   define dgemv_(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12) \
       DGEMV(_cptofcd(t1,1), t2, t3, t4, t5, t6, t7, t8, t9, t10, t11) 
#   define dgemm_(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15) \
       DGEMM(_cptofcd(t1,1),_cptofcd(t2,1),t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13)
#   define dgetrs_(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10) \
           DGETRS(_cptofcd(t1,1),t2,t3,t4,t5,t6,t7,t8,t9)
#   define dtrsm_(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15) \
           DTRSM(_cptofcd(t1,1),_cptofcd(t2,1),_cptofcd(t3,1),_cptofcd(t4,1), \
                 t5,t6,t7,t8,t9,t10,t11)
#   define dpotrf_(t1,t2,t3,t4,t5,t6) \
           DPOTRF(_cptofcd(t1,1),t2,t3,t4,t5)
#   define dtrmm_(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15) \
	   DTRMM(_cptofcd(t1,1),_cptofcd(t2,1),_cptofcd(t3,1),_cptofcd(t4,1), \
                 t5,t6,t7,t8,t9,t10,t11)
#else
#   define dgemv_                     DGEMV
#   define dgemm_                     DGEMM
#   define dgetrs_                    DGETRS
#   define dtrsm_                     DTRSM
#   define dpotrf_                    DPOTRF
#   define dtrmm_                     DTRMM
#endif
#   define dgetrf_                    DGETRF
#   define dgetri_                    DGETRI
#   define dgeco_                     DGECO
#   define dgedi_                     DGEDI
#   define y12mbf_                    Y12MBF
#   define y12mcf_                    Y12MCF
#   define y12mdf_                    Y12MDF
#   define idamax_                    IDAMAX
#   define dasum_                     DASUM
#   define daxpy_                     DAXPY
#   define dlaswp_                    DLASWP
#   define ddot_                      DDOT
#   define dcopy_                     DCOPY
#   define dscal_                     DSCAL
#   define dlaic1_                    DLAIC1
#   define fnroot_                    FNROOT
#   define mc64ad_                    MC64AD
#   define rcm_                       RCM
#   define az_broadcast_              AZ_BROADCAST
#   define az_check_input_            AZ_CHECK_INPUT
#   define az_check_msr_              AZ_CHECK_MSR
#   define az_check_vbr_              AZ_CHECK_VBR
#   define az_defaults_               AZ_DEFAULTS
#   define az_exchange_bdry_          AZ_EXCHANGE_BDRY
#   define az_find_index_             AZ_FIND_INDEX
#   define az_find_local_indices_     AZ_FIND_LOCAL_INDICES
#   define az_find_procs_for_externs_ AZ_FIND_PROCS_FOR_EXTERNS
#   define az_free_memory_            AZ_FREE_MEMORY
#   define az_gavg_double_            AZ_GAVG_DOUBLE
#   define az_gdot_                   AZ_GDOT
#   define az_gmax_double_            AZ_GMAX_DOUBLE
#   define az_gmax_int_               AZ_GMAX_INT
#   define az_gmax_matrix_norm_       AZ_GMAX_MATRIX_NORM
#   define az_gmax_vec_               AZ_GMAX_VEC
#   define az_gmin_double_            AZ_GMIN_DOUBLE
#   define az_gmin_int_               AZ_GMIN_INT
#   define az_gsum_double_            AZ_GSUM_DOUBLE
#   define az_gsum_int_               AZ_GSUM_INT
#   define az_gsum_vec_               AZ_GSUM_VEC
#   define az_gvector_norm_           AZ_GVECTOR_NORM
#   define az_init_quick_find_        AZ_INIT_QUICK_FIND
#   define az_invorder_vec_           AZ_INVORDER_VEC
#   define az_matvec_mult_            AZ_MATVEC_MULT
#   define az_VBR_matvec_mult_        AZ_VBR_MATVEC_MULT
#   define az_MSR_matvec_mult_        AZ_MSR_MATVEC_MULT
#   define az_msr2vbr_                AZ_MSR2VBR
#   define az_order_ele_              AZ_ORDER_ELE
#   define az_pr_error_               AZ_PR_ERROR
#   define az_print_out_              AZ_PRINT_OUT
#   define az_processor_info_         AZ_PROCESSOR_INFO
#   define az_quick_find_             AZ_QUICK_FIND
#   define az_read_msr_matrix_        AZ_READ_MSR_MATRIX
#   define az_read_update_            AZ_READ_UPDATE
#   define az_reorder_matrix_         AZ_REORDER_MATRIX
#   define az_reorder_vec_            AZ_REORDER_VEC
#   define az_get_comm_               AZ_GET_COMM
#   define az_set_comm_               AZ_SET_COMM
#   define az_set_message_info_       AZ_SET_MESSAGE_INFO
#   define az_set_proc_config_        AZ_SET_PROC_CONFIG
#   define az_sort_                   AZ_SORT
#   define az_solve_                  AZ_SOLVE
#   define az_transform_              AZ_TRANSFORM
#elif defined(matched)
#   define dgemm_                     dgemm
#   define dgemv_                     dgemv
#   define dgetrf_                    dgetrf
#   define dgetri_                    dgetri
#   define dgetrs_                    dgetrs
#   define dgeco_                     dgeco
#   define dgedi_                     dgedi
#   define y12mbf_                    y12mbf
#   define y12mcf_                    y12mcf
#   define y12mdf_                    y12mdf
#   define idamax_                    idamax
#   define dasum_                     dasum
#   define daxpy_                     daxpy
#   define dtrmm_                     dtrmm
#   define dlaswp_                    dlaswp
#   define ddot_                      ddot
#   define dcopy_                     dcopy
#   define dscal_                     dscal
#   define dtrsm_                     dtrsm
#   define dpotrf_                    dpotrf
#   define dlaic1_                    dlaic1
#   define fnroot_                    fnroot
#   define mc64ad_                    mc64ad
#   define rcm_                       rcm
#   define az_broadcast_              az_broadcast
#   define az_check_input_            az_check_input
#   define az_check_msr_              az_check_msr
#   define az_check_vbr_              az_check_vbr
#   define az_defaults_               az_defaults
#   define az_exchange_bdry_          az_exchange_bdry
#   define az_find_index_             az_find_index
#   define az_find_local_indices_     az_find_local_indices
#   define az_find_procs_for_externs_ az_find_procs_for_externs
#   define az_free_memory_            az_free_memory
#   define az_gavg_double_            az_gavg_double
#   define az_gdot_                   az_gdot
#   define az_gmax_double_            az_gmax_double
#   define az_gmax_int_               az_gmax_int
#   define az_gmax_matrix_norm_       az_gmax_matrix_norm
#   define az_gmax_vec_               az_gmax_vec
#   define az_gmin_double_            az_gmin_double
#   define az_gmin_int_               az_gmin_int
#   define az_gsum_double_            az_gsum_double
#   define az_gsum_int_               az_gsum_int
#   define az_gsum_vec_               az_gsum_vec
#   define az_gvector_norm_           az_gvector_norm
#   define az_init_quick_find_        az_init_quick_find
#   define az_invorder_vec_           az_invorder_vec
#   define az_matvec_mult_            az_matvec_mult
#   define az_VBR_matvec_mult_        az_VBR_matvec_mult
#   define az_MSR_matvec_mult_        az_MSR_matvec_mult
#   define az_msr2vbr_                az_msr2vbr
#   define az_order_ele_              az_order_ele
#   define az_pr_error_               az_pr_error
#   define az_print_out_              az_print_out
#   define az_processor_info_         az_processor_info
#   define az_quick_find_             az_quick_find
#   define az_read_msr_matrix_        az_read_msr_matrix
#   define az_read_update_            az_read_update
#   define az_reorder_matrix_         az_reorder_matrix
#   define az_reorder_vec_            az_reorder_vec
#   define az_get_comm_               az_get_comm
#   define az_set_comm_               az_set_comm
#   define az_set_message_info_       az_set_message_info
#   define az_set_proc_config_        az_set_proc_config
#   define az_sort_                   az_sort
#   define az_solve_                  az_solve
#   define az_transform_              az_transform
#endif

#ifndef FSUB_TYPE
#  if defined(ncube)
#     define  FSUB_TYPE void
#  elif defined(paragon)
#     define  FSUB_TYPE void
#  elif defined(hp)
#     define  FSUB_TYPE void
#  else
#     define  FSUB_TYPE void
/*#     define  FSUB_TYPE int*/
#  endif
#endif

#ifdef __cplusplus
#include <stdio.h>
extern "C" {
#endif

/* BLAS, LAPACK stuff */

extern double dasum_(int *, double *, int *);

extern FSUB_TYPE  daxpy_(int *n1, double *a1, double *x, int *skipx1, double *y,
                         int *skipy1);

extern FSUB_TYPE  dcopy_(int *n1, double *src, int *src1, double *dst,
                         int *dst1);

extern double     ddot_(int *n1, double *v1, int *dum11, double *v2,
                        int *dum21);

extern FSUB_TYPE  dgeco_(double *, int *, int *, int *, double *, double *);

extern FSUB_TYPE  dgedi_(double *, int *, int * , int *, double *, double *,
                         int *);

extern FSUB_TYPE  dlaswp_(int *, double *, int *, int *, int *, int *, int *);

#ifdef T3E
extern FSUB_TYPE  DGEMM(_fcd, _fcd, int *, int *, int *, double *, double *, 
                        int *, double *, int *, double *, double *, int *);

extern void       DGEMV(_fcd, int *, int *, double *, double *, int *, double *,
                        int *, double *, double *, int *);

extern FSUB_TYPE  DGETRS(_fcd, int *, int *, double *, int *, int *,
                          double *, int *, int *);

extern FSUB_TYPE  DTRMM(_fcd, _fcd, _fcd, _fcd, int *, int *,
			 double *, double *, int *, double *, int *);

extern FSUB_TYPE  DTRSM(_fcd, _fcd, _fcd, _fcd, int *, int *,
                         double *, double *, int *, double *, int *);

extern FSUB_TYPE  DPOTRF(_fcd, int *, double *,int *, int *);

#else
extern FSUB_TYPE  dgemm_(char *, char *, int *, int *, int *, double *, 
         double *, int *, double *, int *, double *, double *, int *
#ifdef SPARSEBLAS
);
#else
, unsigned int, unsigned int);
#endif

extern void   dgemv_(char *, int *, int *, double *, double *, int *, double *,
                     int *, double *, double *, int *, unsigned int);

extern FSUB_TYPE  dgetrs_(char *, int *, int *, double *, int *, int *,
                          double *, int *, int *, unsigned int);

extern FSUB_TYPE  dtrmm_(char *, char *, char *, char *, int *, int *,
			 double *, double *, int *, double *, int *, 
			 unsigned int, unsigned int, unsigned int,
			 unsigned int);

extern FSUB_TYPE  dtrsm_(char *, char *, char *, char *, int *, int *,
                        double *, double *, int *, double *, int *,
                        unsigned int, unsigned int, unsigned int, unsigned int);

extern FSUB_TYPE  dpotrf_(char *, int *, double *,int *, int *, unsigned int);

#endif

extern FSUB_TYPE  dgetrf_(int *, int *, double *, int *, int *, int *);

extern FSUB_TYPE  dgetri_(int *, double *, int *, int *, double *, int *,
                          int *);

extern FSUB_TYPE  dlaic1_(int * , int *, double *, double *, double *, double *,
			  double *, double *, double *);

extern FSUB_TYPE  dscal_(int *n1, double *a1, double *r, int *stride1);

extern int    idamax_(int *, double *, int *);


/* Aztec function prototypes that can be called by the user */

extern void AZ_solve(
        double x[],     /* On input 'x' contains the initial guess. On output*/
                        /* 'x' contains the solution to our linear system.   */
                        /* NOTE: THis vector must be of size >= N + NExt     */
        double b[],     /* right hand side of linear system.                 */
                        /* NOTE: This vector must be of size >= N            */
        int options[],
        double params[],
        int indx[],     /* The ith element of indx points to the location in */
                        /* val of the (0,0) entry of the ith block entry. The*/
                        /* last element is the number of nonzero entries of  */
                        /* matrix A plus one.                                */
        int bindx[],    /* Contains the block column indices of the non-zero */
                        /* block entries.                                    */
        int rpntr[],    /* The ith element of rpntr indicates the first point*/
                        /* row in the ith block row. The last element is the */
                        /* number of block rows plus one.                    */
        int cpntr[],    /* The jth element of cpntr indicates the first point*/
                        /* column in the jth block column. The last element  */
                        /* is the number of block columns plus one.          */
        int bpntr[],    /* The ith element of bpntr points to the first block*/
                        /* entry of the ith row in bindx. The last element is*/
                        /* the number of nonzero blocks of matrix A plus one.*/
        double val[],   /* matrix A in sparse format (VBR)  .                */
                        /* Indicates current level of factorization          */
                        /* factor_flag =                                     */
                        /*      1: indicates first call to precond. routine  */
                        /*      that performs some type of factorization     */
                        /*      preprocessing such as an incomplete LU.      */
                        /*                                                   */
                        /*      2: use preprocessing info. from a previous   */
                        /*      call. Implies some further change in the     */
                        /*      the numerical entries rather than the sparse */
                        /*      pattern.                                     */
                        /*                                                   */
                        /*      3: use precondtioner from last level 1 or 2  */
                        /*      call to precond. (see specific precondioner  */
                        /*      for more info)                               */
        int data_org[], double status[], int proc_config[]);

extern int AZ_initialize(double x[], double b[], int options[],
        double params[], double status[], int proc_config[], AZ_MATRIX *Amat,
        AZ_PRECOND *precond, int save_old_values[], struct AZ_SCALING *);

extern void AZ_finalize(double x[], double b[], int options[], int
   proc_config[], AZ_MATRIX *Amat, AZ_PRECOND *precond, int save_old_values[],
   struct AZ_SCALING *scaling);

extern void AZ_iterate_setup(int options[], double params[], int proc_config[],
        AZ_MATRIX *Amat, AZ_PRECOND *precond);

extern void AZ_iterate_finish(int options[], AZ_MATRIX *Amat,
        AZ_PRECOND *precond);

extern int AZ_oldsolve_setup(double x[], double b[], int options[],
        double params[], double status[], int proc_config[], AZ_MATRIX *Amat,
        AZ_PRECOND *precond, int save_old_values[], struct AZ_SCALING *);

extern void AZ_oldsolve_finish(double x[], double b[], int options[],
        int proc_config[], AZ_MATRIX *Amat, int save_old_values[],
	struct AZ_SCALING *);

extern void AZ_abs_matvec_mult (double *b, double *c,AZ_MATRIX *Amat,int proc_config[]);

extern void   AZ_add_new_ele(int cnptr[], int col, int blk_row, int bindx[],
                             int bnptr[], int indx[], double val[], int therow,
                             double new_ele, int maxcols, int blk_space,
                             int nz_space, int blk_type);

extern void   AZ_add_new_row(int therow, int *nz_ptr, int *current, double
                             **val, int **bindx, char *input,FILE *dfp,
                             int *msr_len, int *column0);

extern int AZ_adjust_N_nz_to_fit_memory(int N, int , int);

extern char *AZ_allocate(unsigned int iii);

extern char *AZ_allocate_or_free(void *ptr, unsigned int size, int action);

extern void   AZ_backsolve(double newa[], double pivot[], double x[], int snr[],
                  int ha[], int iflag[], int *ifail, int *nn, int *n, int *iha);

extern void AZ_block_diagonal_scaling(int action, AZ_MATRIX *Amat, double val[],
        int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
        int data_org[], double b[], int options[], int proc_config[],
        struct AZ_SCALING *scaling);


extern int    AZ_breakdown_f(int N, double v[], double w[], double inner,
                             int proc_config[]);

extern void AZ_broadcast(char *ptr, int length, int proc_config[], int action);

extern unsigned int AZ_broadcast_info(char buffer[], int proc_config[], 
			        unsigned int length);

extern void   AZ_calc_blk_diag_inv(double *val, int *indx, int *bindx,
                                  int *rpntr, int *cpntr, int *bpntr,
                                  double *d_inv, int *d_indx, int *d_bindx,
                                  int *d_rpntr, int *d_bpntr, int data_org[]);

extern void AZ_calc_blk_diag_LU(double *val, int *indx, int *bindx, int *rpntr,
                          int *cpntr, int *bpntr, double *d_inv, int *d_indx,
                          int *d_bindx, int *d_rpntr, int *d_bpntr,
                          int *data_org, int *ipvt);

extern double AZ_calc_iter_flops(int solver_flag, double inner_flops,
                                 double daxpy_flops, double matvec_flops,
                                 int total_its, double gnnz, double K);

extern double AZ_calc_precond_flops(int solver_flag, int options[],
                              double daxpy_flops, double matvec_flops,
                              int total_its, int gn, double gnnz, 
                              int data_org[], int proc_config[]);

extern double AZ_calc_solve_flops(int options[], int, double , int , double,
                                  int data_org[], int proc_config[]);

extern void   AZ_change_it(int indx[], int length, int *first, int *total,
                           int b[]);

extern void   AZ_change_sign(double *lambda_max, double val[], int indx[],
                             int bindx[], int rpntr[], int cpntr[], int bpntr[],
                             int data_org[]);
extern void AZ_check_block_sizes(int bindx[], int cpntr[], int Nrows,
                          int *new_block);

extern int  AZ_check_input(int data_org[], int options[], double params[],
                           int proc_config[]);

extern void AZ_check_msr(int *bindx, int N_update, int N_external,
                         int option, int *proc_config);

extern int    AZ_check_options(int * , int ,int data_org[], int,double *,
                                AZ_MATRIX *, AZ_PRECOND *);

extern void AZ_check_vbr(int N_update, int N_external, int option,
                         int bindx[], int bnptr[], int cnptr[], int rnptr[],
                         int proc_config[]);

extern void AZ_combine_overlapped_values(int sym_flag, int data_org[],
	int options[], double x[], int map[], double ext_vals[], int name,
	int proc_config[]);

extern int AZ_compare_update_vs_soln(int N, double, double alpha, double p[], 
        double x[], 
        double update_reduction, int ouput_flag, int proc_config[], int *first_time);

extern int AZ_compress_msr(int *ibindx[], double *ival[], int allocated, 
			int needed, int name, struct context *context);

extern void   AZ_compute_global_scalars(AZ_MATRIX *Amat, 
                                        double x[], double b[], double r[],
                                        double w[], double *r_norm,
                                        double *scaled_r_norm, int option_i[],
                                        int data_org[], int proc_config[],
                                        int *use_r, double v1[], double v2[],
                                        double *value, 
					struct AZ_CONVERGE_STRUCT *);

extern void AZ_compute_matrix_size(AZ_MATRIX *Amat, int options[], 
	int N_nz_unpadded, int N_unpadded, int *N_nz_padded, int N_external, 
	int *max_row, int *N, int *N_nz, double fill,int *extra_fact_nz_per_row,
	int Nb_unpadded, int *bandwidth);

extern int AZ_compute_max_nz_per_row(AZ_MATRIX *Amat, int N, int Nb, 
	int *largest_band);


extern void   AZ_compute_residual( double b[], double u[], double r[],
                                  int proc_config[], AZ_MATRIX *);

extern void   AZ_convert_ptrs_to_values(int array[], int length);

extern void   AZ_convert_values_to_ptrs(int array[], int length, int start);

extern struct AZ_CONVERGE_STRUCT *AZ_converge_create(void);
extern void AZ_converge_destroy(struct AZ_CONVERGE_STRUCT **temp);

extern AZ_MATRIX *AZ_create_matrix(int local, int additional, int matrix_type, 
                                   int local_blks, int *not_using);

extern void AZ_defaults(int options[], double params[]);

extern void AZ_delete_matrix(AZ_MATRIX *ptr);

extern void   AZ_dgemv2(int m, int n, double *a, double *x, double *y);

extern void   AZ_dgemv3(int m, int n, double *a, double *x, double *y);

extern void   AZ_direct_sort(int b[], int indx[], char buffer[], char a[],
                             int *start, int buf_len, int *ind_index,
                             int *the_first, int *real_lists, int *pre_mid);

extern void   AZ_divide_block(int i, int j, double val[], int indx[],
                              int bindx[], int cpntr[], double *z,
                              double *blockj, double *blocki, int *ipvt);

extern void   AZ_divide_block0(int i, int j, double val[], int indx[],
                                  int bindx[], int cpntr[], int *ipvt);

extern void AZ_domain_decomp(double x[], AZ_MATRIX *Amat, int options[],
                   int proc_config[], double params[],
                   struct context *context);

extern void   AZ_dtrans(int *, int *, double *);

extern void AZ_equil_scaling(int action, AZ_MATRIX *Amat, 
                              double b[],
                              double x[], int options[],
                              int proc_config[], struct AZ_SCALING *scaling);

extern void AZ_exchange_bdry(double x[], int data_org[], int proc_config[]);

extern void   AZ_exchange_local_info(int N_neighbors, int proc_num_neighbor[],
                                  char *message_send_add[],
                                  unsigned int message_send_length[],
                                  char *message_recv_add[],
                                  unsigned int message_recv_length[], 
				  int type, int proc_config[]);

extern int AZ_exit(int input);

extern int AZ_extract_comm_info(int **idata_org, int (*user_comm)(double *,
	AZ_MATRIX *), AZ_MATRIX *, 
	int proc_config[], int N_cols, int Nghost);

extern void AZ_fact_bilu(int new_blks, AZ_MATRIX *A_overlapped, 
				 int *diag_block, int *pivot);

extern void AZ_fact_chol(int bindx[], double val[], int N,
                           double rthresh, double athresh);

extern void  AZ_fact_ilut( int *, AZ_MATRIX *, double *a, int *ja,
		  double drop, int extra_fact_nz_per_row, int shift,
		  int *iu, double *cr, double *unorm, int *ind, 
		  int *nz_used, int *jnz,
            double rthresh, double athresh);
     
extern void   AZ_fact_lu(double x[], AZ_MATRIX *A_overlapped, double *aflag, double *pivot,
                         int *rnr, int *ha, int *iflag, int *, int*, int *, int *, int *);

extern void AZ_fact_rilu(int N, int *nz_used, int *iu, int *iw,
                         AZ_MATRIX *A_overlapped, double omega,
                         double rthresh, double athresh);

extern void AZ_factor_subdomain(struct context *context, int N, 
	int N_nz, int *nz_used);

extern int    AZ_fill_sparsity_pattern(struct context *context, int ifill, 
			     int bindx[], double val[], int N);

extern int    AZ_find_block_col(int cnptr[], int column, int maxcols,
                                int blk_type);

extern int    AZ_find_block_in_row(int bindx[], int bnptr[], int i, int blk_col,
                                   int indx[], int, double val[], int blk_space,
                                   int nz_space);

extern void AZ_find_MSR_ordering(int bindx2[],int **ordering,int N,
        int **inv_ordering, int name, struct context *);

extern int    AZ_find_closest_not_larger(int key, int list[], int length);

extern int  AZ_find_index(int key, int list[], int length);

extern void AZ_find_local_indices(int N_update, int bindx[], int update[],
                                  int **external, int *N_external, int mat_type,
                                  int bpntr[]);

extern void AZ_find_procs_for_externs(int N_update, int update[],
                                      int external[], int N_external,
                                      int proc_config[], int **extern_proc);

extern int AZ_find_simple(int, int *, int, int *, int, int *, int *);

extern void AZ_find_global_ordering(int proc_config[], AZ_MATRIX *Amat,
                             int **global_bindx, int **update);

extern void AZ_revert_to_global(int proc_config[], AZ_MATRIX *Amat,
                             int **global_bindx, int **update);


extern void AZ_fix_pt(double *, double *, double *, int *, double * , int * ,
	double *, AZ_MATRIX *, AZ_PRECOND *, struct AZ_CONVERGE_STRUCT *);

extern void AZ_flop_rates(int data_org[],int indx[],int bpntr[], int bindx[],
              int options[], double status[], double total_time,
              int proc_config[]);

extern void AZ_free(void *ptr);

extern void AZ_free_memory(int label);

extern void AZ_free_space_holder(struct context *variables);

extern void   AZ_gappend(int vals[], int *cur_length, int total_length,
                         int proc_config[]);

extern double AZ_gavg_double(double var, int proc_config[]);

extern double AZ_gdot(int N, double r[], double z[], int proc_config[]);

extern void   AZ_gdot_vec(int N, double dots[], double dots2[],
                          int proc_config[]);

extern int    AZ_get_block(int j, int k, int bindx[], int bpntr[], int *ptr_j);

extern unsigned int AZ_get_sol_param_size(void);

extern int    AZ_get_new_eps(double *epsilon, double, double,
                             int proc_config[]);

extern void   AZ_get_poly_coefficients(int power, double b, double c[],
                                       int param_flag);

extern int    AZ_get_sym_indx(int, int, int *, int *, int *);

extern void   AZ_get_x_incr(int options[], int data_org[], int proc_config[], 
			    double params[], int i, double **hh, double *rs, 
                            double *trash, double **ss, AZ_MATRIX *, 
                            AZ_PRECOND *, double *, int *, int *, int);

extern double AZ_gmax_double(double, int proc_config[]);

extern int    AZ_gmax_int(int val, int proc_config[]);

extern double AZ_gmax_matrix_norm(double val[], int indx[], int bindx[],
                                  int rpntr[], int cpntr[], int bpntr[],
                                  int proc_config[], int data_org[]);

extern double AZ_gmax_vec(int N, double vec[], int proc_config[]);

extern double AZ_gmin_double(double var, int proc_config[]);

extern int AZ_gmin_int(int val, int proc_config[]);

extern double AZ_gsum_double(double , int proc_config[]);

extern int    AZ_gsum_int(int totals, int proc_config[]);

extern void   AZ_gsum_vec_int(int vals[], int vals2[], int length,
                              int proc_config[]);

extern double AZ_gvector_norm(int n, int p, double *x, int *);

extern void AZ_hold_space(struct context *context, int N);

extern void   AZ_init_quick_find(int list[], int length, int *shift, int *bins);

extern void   AZ_init_subdomain_solver(struct context *context);

extern void   AZ_invorder_vec(double vec[], int data_org[], int update_index[],
			     int rpntr[],double newarray[]);

extern void AZ_list_print(int ivec[] , int length, double dvec[], int length2);

extern void AZ_loc_avg(AZ_MATRIX *Amat, double r[], double newr[], int N_fixed,
	   int fixed_pts[], int proc_config[]);

extern void AZ_lower_triang_vbr_solve(int Nrows, int cpntr[], int bpntr[],
           int indx[], int bindx[], double val[], double b[]);

extern void  AZ_lower_icc(int bindx[],double val[],int N, double rhs[]);

extern void  AZ_lower_tsolve(double x[],int  , double l[], int il[], 
                            int jl[],  double y[] );

extern double *AZ_manage_memory(unsigned int size, int action, int type, 
				char *name, int *status);

extern struct AZ_MATRIX_STRUCT *AZ_matrix_create(int local);

extern struct AZ_MATRIX_STRUCT *AZ_submatrix_create(AZ_MATRIX *Amat, int Nsub_rows, 
				int sub_rows[], int Nsub_cols, int sub_cols[], int proc_config[]);

void AZ_submatrix_destroy(AZ_MATRIX **submat);

extern struct AZ_MATRIX_STRUCT *AZ_blockmatrix_create(AZ_MATRIX **submat_list, int Nsub_mats, 
        int **submat_locs, int Nblock_rows, int Nblock_cols, int Nsub_rows[], int **sub_rows, 
        int Nsub_cols[], int **sub_cols, int proc_config[]);

void AZ_blockmatrix_destroy(AZ_MATRIX **blockmat);

extern void AZ_matrix_init(AZ_MATRIX *Amat, int local);

extern struct AZ_PREC_STRUCT   *AZ_precond_create(struct AZ_MATRIX_STRUCT *Pmat,
		void (*prec_fun)( double *, int *, int *, double *, 
		      struct AZ_MATRIX_STRUCT  *, struct AZ_PREC_STRUCT *),
		void *data);

extern void AZ_matrix_destroy( struct AZ_MATRIX_STRUCT **Amat);
extern void AZ_precond_destroy(struct AZ_PREC_STRUCT **precond);



extern void AZ_matfree_Nnzs(AZ_MATRIX *Amat);

extern void AZ_matfree_2_msr(AZ_MATRIX *Amat,double *val, int *bindx, int N_nz);

extern void AZ_mat_colperm(int N, int bindx2[], double val2[],
               int **inv_ordering, int name, struct context *);

extern void AZ_mat_reorder(int n, int bindx[], double val[], int perm[],
               int invp[]);

extern void   AZ_matvec_mult(double *val, int *indx, int *bindx, int *rpntr,
                             int *cpntr, int *bpntr, double *b, double *c,
                             int exchange_flag, int *data_org);

extern void AZ_mk_context(int options[], double params[], int data_org[], 
			  AZ_PRECOND *precond, int proc_config[]);

extern void AZ_mk_identifier(double *params, int *options,
        int *data_org, char *tag);

extern void AZ_MSR_mult_patterns(int *bindx, int N, int *work1, int length,
	int *work2);

extern void AZ__MPI_comm_space_ok(void);

extern int  AZ_MSR_getrow(int columns[], double values[], int row_lengths[],
        struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space);

extern int  AZ_VBR_getrow(int columns[], double values[], int row_lengths[],
        struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space);




extern void   AZ_msr2lu(int oldN, AZ_MATRIX *A_overlapped, int *rnr);

extern void   AZ_msr2vbr(double val[], int indx[], int rnptr[], int cnptr[],
                         int bnptr[], int bindx[], int msr_bindx[],
                         double msr_val[], int total_blk_rows,
                         int total_blk_cols, int blk_space, int nz_space,
                         int blk_type);

extern void AZ_msr2vbr_mem_efficient(int N, int **ibindx, double **ival,
			  int **icpntr, int **ibpntr, int **iindx, 
			  int *N_blk_rows, int name, char *label, int);

extern void   AZ_order(int M, double *val_old, double *val_new, int *bindx,
                       int *indx_old, int *indx_new, int *bpntr,
                       int *diag_bloc);

extern void   AZ_order_ele(int update_index[], int extern_index[],
                           int *internal, int *border, int N_update,
                           int msr_bindx[], int bindx[], int extern_proc[],
                           int N_external, int option, int m_type);

extern void   AZ_p_error(char *str, int proc);

extern void AZ_pad_matrix(struct context *context, int proc_config[],
   int N_unpadded, int *N, int **map, int **padded_data_org,
   int *N_nz, int estimated_requirements);

extern void   AZ_pbicgstab(double *, double *, double *, int *, double *,
	int *, double *, AZ_MATRIX *, AZ_PRECOND *, struct AZ_CONVERGE_STRUCT *);

extern void   AZ_pcg_f(double *, double *, double *, int *, double * , int * ,
	double *, AZ_MATRIX *, AZ_PRECOND *, struct AZ_CONVERGE_STRUCT *);

extern void   AZ_pcgs(double *, double *, double *, int *, double * , int * ,
	double *, AZ_MATRIX *, AZ_PRECOND *, struct AZ_CONVERGE_STRUCT *);

extern void AZ_perror(char *string);

extern void   AZ_pgmresr(double *, double *, double *, int *, double * , int * ,
	double *, AZ_MATRIX *, AZ_PRECOND *, struct AZ_CONVERGE_STRUCT *);

extern void   AZ_pgmres(double *, double *, double *, int *, double * , int * ,
	double *, AZ_MATRIX *, AZ_PRECOND *, struct AZ_CONVERGE_STRUCT *);

extern void   AZ_polynomial_expansion( double z[], int options[],
                                      int proc_config[], AZ_PRECOND *);

extern int AZ_pos( int , int bindx[] , int position[], int inv_ordering[],
        double , int );


extern void   AZ_precondition(double x[], int options[], int proc_config[],
                              double params[], AZ_MATRIX *, AZ_PRECOND *);

extern void   AZ_pqmrs(double *, double *, double *, int *, double * , int * ,
	double *, AZ_MATRIX *, AZ_PRECOND *, struct AZ_CONVERGE_STRUCT *);

extern void   AZ_print_call_iter_solve(int * , double *, int , int, 
			AZ_PRECOND *);

extern void   AZ_print_error(int error_code);

extern void AZ_print_header(int options[], int mem_overlapped,
                          int mem_orig, int mem_factor);

extern void   AZ_capture_matrix(double val[], int indx[], int bindx[],
                               int rpntr[], int cpntr[], int bpntr[],
                               int proc_config[], int data_org[], double b[]);

extern void AZ_print_out(int update_index[], int extern_index[], int update[],
                        int external[],
                        double val[], int indx[],  int
                        bindx[], int rpntr[], int cpntr[], int bpntr[], int
                        proc_config[], int choice, int matrix, int N_update,
                        int N_external, int off_set );

extern void   AZ_print_sync_start(int proc,int do_print_line,int proc_config[]);

extern void   AZ_print_sync_end(int proc_config[], int do_print_line);

extern void   AZ_processor_info(int proc_config[]);

extern void AZ_put_in_dbl_heap(int *row, double vals[], int heap[],
			       int *length);
extern void AZ_put_in_heap(int heap[], int *val,int *length);

extern int    AZ_quick_find(int key, int list[],int length, int shift,
                            int bins[]);

extern void   AZ_random_vector(double u[], int data_org[], int proc_config[]);

extern void   AZ_read_msr_matrix(int update[], double **val, int **bindx,
                                 int N_update, int proc_config[]);

extern void   AZ_read_update(int *N_update_blks, int *update_blks[],
                                 int proc_config[], int bigN, int chunk,
                                 int input_option);

extern void   AZ_input_msr_matrix(char datafile[], int update[], double **val, int **bindx,
                                 int N_update, int proc_config[]);

extern void   AZ_input_update(char datafile[], int *N_update_blks, int *update_blks[],
                                 int proc_config[], int bigN, int chunk,
                                 int input_option);

extern char *AZ_realloc(void *ptr, unsigned int size);

extern void AZ_recover_sol_params(int instance, int **sub_options, 
	double **sub_params, double **sub_status, AZ_MATRIX **sub_matrix,
        AZ_PRECOND **sub_precond, struct AZ_SCALING **);

extern void   AZ_reorder_matrix(int N_update, int bindx[], double val[],
                                int update_index[], int extern_index[],
                                int indx[], int rnptr[], int bnptr[],
                                int N_external, int cnptr[], int option,
                                int);

extern void   AZ_reorder_vec(double vec[], int data_org[], int update_index[],
			     int rpntr[]);

extern void AZ_global2local(int data_org[], int bindx[], int update[], 
			    int update_index[], int externs[], int extern_index[]);

extern void AZ_restore_unreordered_bindx(int bindx[], double val[], int update[],
					 int update_index[], int external[],
					 int extern_index[], int data_org[]);

extern void   AZ_reverse_it(int indx[], int length, int first, int total,
                            int b[]);

extern void AZ_rm_context(int options[], double params[], int data_org[]);

extern void AZ_rm_dbl_heap_root(int heap[], double vals[], int *length);
 
extern void AZ_rm_heap_root(int heap[], int *length);

extern void AZ_rm_duplicates(int array[], int *N);

extern void AZ_row_sum_scaling(int action, AZ_MATRIX *Amat, 
	double b[], int options[], 
	struct AZ_SCALING *scaling);

extern void AZ_scale_f(int action, AZ_MATRIX *Amat, int options[], double b[], 
                     double x[], int proc_config[], struct AZ_SCALING *scaling);

extern struct AZ_SCALING *AZ_scale_matrix_only(AZ_MATRIX *Amat, int options[],
                        int proc_config[]);

extern void AZ_scale_rhs_only(double b[], AZ_MATRIX *Amat, int options[], 
                       int proc_config[], struct AZ_SCALING *scaling);

extern void AZ_scale_sol_only(double x[], AZ_MATRIX *Amat, int options[], 
                       int proc_config[], struct AZ_SCALING *scaling);

extern void AZ_scale_rhs_sol_before_iterate(double x[], double b[], 
 AZ_MATRIX *Amat, int options[], int proc_config[], struct AZ_SCALING *scaling);

extern void AZ_unscale_after_iterate(double x[], double b[], AZ_MATRIX *Amat,
                              int options[], int proc_config[],
                              struct AZ_SCALING *scaling);

extern void AZ_clean_scaling(struct AZ_SCALING **scaling);


extern void   AZ_scale_true_residual(double x[], double b[], double v[],
                                     double w[], double *actual_residual,
                                     double *scaled_r_norm, int options[],
                                     int data_org[], int proc_config[], 
                                     AZ_MATRIX *Amat, 
				     struct AZ_CONVERGE_STRUCT *);

extern struct AZ_SCALING *AZ_scaling_create(void);

extern void AZ_scaling_destroy(struct AZ_SCALING **temp);

extern double AZ_second(void);

extern MPI_AZComm *AZ_get_comm(int proc_config[]);

extern void AZ_set_comm(int proc_config[], MPI_AZComm );

extern void AZ_set_MATFREE_name(AZ_MATRIX *Amat, int name);

extern void AZ_set_MATFREE_matrix_norm(AZ_MATRIX *Amat, double mat_norm);

extern void AZ_set_MATFREE(AZ_MATRIX *Amat, void *data,
    void (*matvec)(double *, double *, struct AZ_MATRIX_STRUCT *, int *));

extern void AZ_set_MATFREE_getrow(AZ_MATRIX *Amat, void *data, 
   int  (*getrow)(int *, double *, int *, struct AZ_MATRIX_STRUCT *, int ,
	int *, int),
   int  (*user_comm)(double *, AZ_MATRIX *), int N_ghost, int proc_config[]);

extern void AZ_set_MSR(AZ_MATRIX *Amat, int bindx[], double val[], 
	int data_org[], int N_update, int update[], int option);

extern void AZ_set_VBR(AZ_MATRIX *Amat, int rpntr[], int cpntr[], int bpntr[],
        int indx[], int bindx[], double val[], int data_org[],
        int N_update, int update[], int option);


extern void   AZ_set_message_info(int N_external, int extern_index[],
                                  int N_update, int external[],
                                  int extern_proc[], int update[],
                                  int update_index[], int proc_config[],
                                  int cnptr[], int *data_org[], int);

extern void AZ_set_precond_print_string(struct AZ_PREC_STRUCT *precond, 
				char str[]);

extern void AZ_set_proc_config(int proc_config[], MPI_AZComm );

extern int AZ_set_solver_parameters(double *params, int *options, AZ_MATRIX *Amat,
				AZ_PRECOND *Pmat, struct AZ_SCALING *S);

extern void AZ_setup_dd_olap_msr(int N_rows, int *New_N_rows, int *bindx,
   double *val, int olap_size, int *proc_config, int *data_org[], int **map3,
   int bindx_length,int name, int *prev_data_org,int estimated_requirements,
   struct context *context);
 
extern double AZ_condest(int N, struct context *context);

extern void AZ_solve_subdomain(double x[],int N, struct context *context);

extern void   AZ_sort(int list[], int N, int list2[], double list3[]);

extern void   AZ_sort_dble(char a[], int indx[], int start, int end, int b[],
                           int *mid, int real_lists, char buffer[], int buf_len,
                           int afirst, int );

extern void   AZ_sort_ints(char a[], int indx[], int start, int end, int b[],
                           int *mid, int real_lists, char buffer[], int buf_len,
                           int afirst, int );

extern void AZ_sort_msr(int bindx[], double val[], int N);

extern void   AZ_sortqlists(char a[], int b[], int lists[], int length,
                            int type_length, int ind_length);

extern void AZ_space_for_factors(double input_fill, int N_nz, int N,
        int *extra_factor_nonzeros, int options[],int bandwidth, int );

extern void AZ_space_for_kvecs(int request, int **kvec_sizes, double ***saveme, 
	double **ptap, int *options, int *data_org, char *suffix, int proc, double **);

extern void AZ_space_for_padded_matrix(int overlap, int N_nonzeros, int N,
    int *extra_rows, int *extra_nonzeros, int N_external, int *largest);

extern void   AZ_splitup_big_msg(int num_neighbors, char *buffer, char *buf2,
			         unsigned int element_size, 
				 int *start_send_proc, 
				 int *actual_send_length,int *num_nonzeros_recv,
				 int *proc_num_neighbor, int type, 
				 int *total_num_recv, int *proc_config);

extern double AZ_srandom1(int *seed);

extern void AZ_sum_bdry(double x[], int data_org[], int proc_config[]);


extern void   AZ_sym_block_diagonal_scaling(double val[], int indx[],
                                            int bindx[], int rpntr[],
                                            int cpntr[], int bpntr[],
                                            double b[], int options[],
                                            int data_org[],
					    int proc_config[]
				            /* struct AZ_SCALING * */);

extern void AZ_sym_diagonal_scaling(int action, AZ_MATRIX *Amat, 
                double b[], double x[], int options[],
                int proc_config[], struct AZ_SCALING *scaling);


extern void   AZ_sym_gauss_seidel(void);

extern void   AZ_sym_gauss_seidel_sl(double val[], int bindx[], double x[],
                 int data_org[], int options[], struct context *,
		 int proc_config[]);

extern void AZ_sym_reinvscale_sl(double x[], int data_org[], int options[],
        int proc_config[], struct AZ_SCALING *scaling);

extern void   AZ_sym_rescale_sl(double x[], int data_org[], int options[],
				int proc_config[],struct AZ_SCALING * );

extern void AZ_sym_row_sum_scaling(int action, AZ_MATRIX *Amat, 
                              double b[],
                              double x[], int options[],
                              int proc_config[], struct AZ_SCALING *scaling);

extern void   AZ_sync(int proc_config[]);

extern void   AZ_terminate_status_print(int situation, int iter,
                                        double status[], double rec_residual,
                                        double params[], double scaled_r_norm,
                                        double actual_residual, int options[],
                                        int proc_config[]);

extern void   AZ_transform(int proc_config[], int *external[], int bindx[],
                           double val[], int update[], int *update_index[],
                           int *extern_index[], int *data_org[], int N_update,
                           int indx[], int bnptr[], int rnptr[], int *cnptr[],
                           int mat_type);

extern void   AZ_update_block(int i, int k, int j, double val[], int indx[],
                              int bindx[], int cpntr[]);

extern void  AZ_upper_icc( int bindx[],double val[],int N, double rhs[]);

extern void AZ_upper_triang_vbr_solve(int Nrows, int cpntr[], int bpntr[], 
    int indx[], int bindx[], double val[], double b[], int piv[], int dblock[]);

extern void  AZ_upper_tsolve( double x[],int ,double u[],int iui[], 
                             int ju[]);

extern void   AZ_vb2msr(int m, double val[], int indx[], int bindx[],
                        int rpntr[], int cpntr[], int bpntr[], double msr_val[],
                        int msr_bindx[]);

void AZ_zero_out_context(struct context *);

void AZ_version(char string[]);

extern void   AZ_MSR_matvec_mult(double x[], double b[], AZ_MATRIX *Amat, 
				 int proc_config[]); 

extern void   AZ_VBR_matvec_mult(double x[], double b[], AZ_MATRIX *Amat, 
			         int proc_config[]); 

extern void PAZ_compose_external(int, int*, int *, int *, int **);

extern void PAZ_find_local_indices(int,int*,int*,int*,int,int*);

extern void PAZ_order_ele(int*,int,int*, int, int*, int*, int);

extern void PAZ_set_message_info(int, int, int*, int*, int*, int*, 
                                 int **, int ,int,int ,struct context*);

extern int  PAZ_sorted_search(int, int, int*);

/*****************************************************************************/
/*                    IFPACK interface routine
*/
/*****************************************************************************/
#ifdef IFPACK
extern void az2ifp_blockmatrix (void **bmat, AZ_MATRIX *Amat);
#endif

/*****************************************************************************/
/*                    Machine Dependent communication routines               */
/*****************************************************************************/
extern unsigned int md_wrap_iread(void *, unsigned int, int *, int *, MPI_AZRequest *);

extern unsigned int md_wrap_iwrite(void *,unsigned int, int , int ,int *, MPI_AZRequest *);

extern unsigned int md_wrap_wait(void *, unsigned int, int *, int *,int *,MPI_AZRequest *);

extern unsigned int md_wrap_write(void *, unsigned int , int , int , int *);

extern unsigned int md_wrap_request_free(MPI_AZRequest *);

#define mdwrap_request_free(a)   md_wrap_request_free(a)
#ifdef AZTEC_MPI
#define mdwrap_wait(a,b,c,x,y,z)   md_mpi_wait(a,b,c,(x),(y),(z),proc_config)
#define mdwrap_iwrite(a,b,c,x,y,z) md_mpi_iwrite(a,b,c,(x),(y),(z),proc_config)
#define mdwrap_iread(a,b,c,x,y)   md_mpi_iread((a),(b),(c),(x),(y),proc_config)
#define mdwrap_write(a,b,c,x,y)   md_mpi_write((a),(b),(c),(x),(y),proc_config)

extern unsigned int md_mpi_iread(void *, unsigned int, int *, int *, 
		MPI_AZRequest *, int *);

extern unsigned int md_mpi_iwrite(void *,unsigned int, int , int ,int *, 
		MPI_AZRequest *, int *);

extern unsigned int md_mpi_wait(void *, unsigned int, int *, int *,int *,
		MPI_AZRequest *, int *);

extern unsigned int md_mpi_write(void *, unsigned int ,int , int , int *,int *);
#else
#define mdwrap_wait(a,b,c,x,y,z)   md_wrap_wait(a,b,c,(x),(y),(z))
#define mdwrap_iwrite(a,b,c,x,y,z) md_wrap_iwrite(a,b,c,(x),(y),(z))
#define mdwrap_iread(a,b,c,x,y)   md_wrap_iread((a),(b),(c),(x),(y))
#define mdwrap_write(a,b,c,x,y)   md_wrap_write((a),(b),(c),(x),(y))
#endif
/*****************************************************************************/
/*                    Auxilliary fortran rroutines needed by Aztec           */
/*****************************************************************************/

extern void fnroot_(int *,int *,int *,int *, int *, int *, int *);

extern void mc64ad_(int *, int *, int *, int *, int *, double*,
                    int *, int *, int *, int *, int *, double*,
                    int *, int *);

extern void rcm_(int *, int *,int *, int *,int *, int *, int *);

extern void   y12mbf_(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, int ha[], int *iha, double aflag[],
                      int iflag[], int *ifail);

extern void   y12mcf_(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, double pivot[], double b[], int ha[],
                      int *iha, double aflag[], int iflag[], int *ifail);

extern void   y12mdf_(int *n, double val[], int *nn, double b[], double pivot[],
                      int snr[], int ha[], int *iha, int iflag[], int *ifail);


/*****************************************************************************/
/*                    Auxilliary routines available to users                 */
/*****************************************************************************/

extern void AZ_check_update(int update[], int N_update, int proc_config[]);

extern void AZ_clear_solver_parameters(int handle);

extern void   AZ_mysleep(int i);

extern void   AZ_output_matrix(double val[], int indx[], int bindx[],
                               int rpntr[], int cpntr[], int bpntr[],
                               int proc_config[], int data_org[]);

extern void AZ_print_vbr_matrix(
        int matrix_flag, /* = 0 no matrix output, = 1 output matrix */
        int Proc,        /* Processor number                  */
        int itotal_nodes,/* Number of internal + border nodes */
        int ext_nodes,   /* Number of external nodes          */
        double  val[],   /* matrix A in sparse format (VBR)   */
        int  indx[],     /* The ith element of indx points to the location in */
                         /* val of the (0,0) entry of the ith block entry. The*/
                         /* last element is the number of nonzero entries of  */
                         /* matrix A plus one.                                */
        int bindx[],     /* Contains the block column indices of the non-zero */
                         /* block entries.                                    */
        int rpntr[],     /* The ith element of rpntr indicates the first point*/
                         /* row in the ith block row. The last element is the */
                         /* number of block rows plus one.                    */
        int bpntr[]      /* The ith element of bpntr points to the first block*/
                         /* entry of the ith row in bindx. The last element is*/
                         /* the number of nonzero blocks of matrix A plus one.*/
        );


extern double AZ_sync_timer(int proc_config[]);







/*****************************************************************************/
/*                    Routines just used locally at Sandia                   */
/*****************************************************************************/

#ifdef Sandia
extern void   AZ_dvbr_diag_sparax(int m, double *val, int *rpntr, int *bpntr,
				  double *b, double *c);

extern void   AZ_transpose(int N, double l[], int ijl[], double lt[],
                           int ijlt[], int row_counter[]);

extern void   AZ_psymmlq(double *, double *, double *, int *, double *, int * ,
                       double *, AZ_MATRIX *, AZ_PRECOND *);



extern void   AZ_gather_mesg_info(double x[],int data_org[],char **, char **,
                                  int *, int *);

extern void   AZ_read_local_info(int data_org[], char *message_recv_add[],
                                 int message_recv_length[]);
extern void   AZ_write_local_info(int data_org[], char *message_recv_add[],
                                  char *message_send_add[],
                                  int message_recv_length[],
                                  int message_send_length[]);

#endif

/*****************************************************************************/
/*                    Timing Routine                                         */
/*****************************************************************************/

#ifdef TIME_VB
extern void   AZ_time_kernals(int , int , double , double *, int *, int *,
                              int *, int *, int *, double *, double *, int, 
                              double *,AZ_MATRIX *);
#endif

#ifdef next_version
extern void   AZ_sym_rescale_vbr(double x[], int data_org[], int options[]);
#endif

/* When calling this fortran routine from C we need to include an extra     */
/* parameter on the end indicating the string length of the first parameter */
#include "az_aztec_defs.h"


#if defined (hp)
extern void   dgemvnsqr_(int *, double *, double *, double *);
extern void   vec_$dcopy(double *, double *, int *);
extern void   blas_$dgemm(char *, char *, int *, int *, int *, double *,
                          double *, int *, double *, int *, double *, double *,
                          int *, int, int);
#endif



#ifdef __cplusplus
}
#endif

#endif
