/*
#define AZ_CONVR0
#define AZ_CONVRHS
#define ML_partition
#define ReuseOps
*/
/*****************************************************************************/
/* Copyright 1999, Sandia Corporation. The United States Government retains  */
/* a nonexclusive license in this software as prescribed in AL 88-1 and AL   */
/* 91-7. Export of this program may require a license from the United States */
/* Government.                                                               */
/*****************************************************************************/

/*****************************************************************************/
/* Sample driver for AZTEC/ML package. The software is tested by reading     */ 
/* input to drive ML from a file (ml_inputfile). Matrices can be created     */
/* by using Aztec's AZ_matrix_capture capability and then running            */
/* ml/util/az_capt2read.c (see comments for ways using matlab to create      */
/* matrices).                                                                */ 
/*****************************************************************************/

/*******************************************************************************
Here is a sample input file:
#
#  Test input file ML 
#       
# Notes: 
#   1) Capitalization should not matter
#   2) Lines starting with # are comments
#   3) Parallel Partitioning File is not used
#      when only 1 processor is present.
#   4) comments in [] indicate available options.
#   5) The matrix must be stored in a file '.data' according
#      to Aztec's AZ_read_msr() format.
#   6) Including the file 'rhsfile' will cause this 
#      data (stored in Aztec's AZ_read_msr() format)
#      to be used as righthand side.
#   7) Including the file 'initguessfile' will cause this 
#      data (stored in Aztec's AZ_read_msr() format)
#      to be used as righthand side.
#   8) rigid body mode information can be input by
#      keeping files 'rigid_body_mode%d' (stored
#      in Aztec's AZ_read_msr() format) where %d
#      starts from 0 and increases. Each file
#      should contain 1 rigid body mode.
#   9) The utility ml/util/az_capt2read.c (see comments)
#      can be used to convert matlab/Aztec type data into 
#      AZ_read_msr() format.
#
-----------------------------------------------
      General Problem Specifications
-----------------------------------------------
Number of DOF per node       = 1
Parallel Partitioning File   = myfile  
Output Frequency             = 2       
Tolerance                    = 1.0e-11

-----------------------------------------------
      Solution Specifications
-----------------------------------------------
Max Number of Levels         = 4
Type of Smoother             = SymGaussSeidel 
#                              [Parasails, GaussSeidel, SymGaussSeidel, 
#                               BlockGaussSeidel, Aggregate, Jacobi, Metis]
Smoother steps per level     = 7
Coarse grid solver           = SuperLU
#                              [Parasails, GaussSeidel, SymGaussSeidel, 
#                               BlockGaussSeidel, Aggregate, Jacobi, Metis,
#                               SuperLU]
Coarse Grid iterations       = 1
Outer Iteration              = Cg
#                              [Cg, Bicgstab, Tfqmr, Gmres] 

-----------------------------------------------
      Aggregation Specifications
-----------------------------------------------
Type of Aggregation          = Mis
#                              [Mis, Uncoupled, Coupled]
Aggregate threshold          = 0.0
Max coarse size              = 30
Smoothed aggregation damping = 1.5 
Spectral norm calculation    = Anorm
#                              [Anorm, Calc]
# end of sample inputfile
*******************************************************************************/

/*******************************************************************************
Output files:
    PPn_1 (n = 1,...,numproc)   Matrix-vector product.  Matrix is the
                                prolongator P, and vector is the vector whose
                                i_th entry is the aggregate that i belongs to. 

    aggn_0 (n=1,...,numproc)    (Node,aggregate) pairs.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "az_aztec.h"
#include "ml_include.h"
#include "ml_read_utils.h"
#include <math.h>
#include "ml_mat_formats.h"
#include "ml_aztec_utils.h"

extern int AZ_using_fortran;
int    parasails_factorized = 0;
int    parasails_sym        = 1;
double parasails_thresh     = 0.01;
int    parasails_nlevels    = 0;
double parasails_filter     = 0.;
double parasails_loadbal    = 0.;
int my_proc_id;

#define HORIZONTAL 1
#define VERTICAL   2
#define OUT        3
extern int inv2dindex(int index, int *i, int *j, int n, int *hor_or_vert);
extern int inv2dnodeindex(int index, int *i, int *j, int n);
extern int northeast2d(int i, int j, int n);
extern int northwest2d(int i, int j, int n);
extern int southeast2d(int i, int j, int n);
extern int southwest2d(int i, int j, int n);

extern int north2d(int i, int j, int n);
extern int south2d(int i, int j, int n);
extern int east2d(int i, int j, int n);
extern int west2d(int i, int j, int n);
extern int southwestfront3d(int i, int j, int k, int n);
extern int northwestback3d(int i, int j, int k, int n);
extern int southwestback3d(int i, int j, int k, int n);
extern int northwestfront3d(int i, int j, int k, int n);
extern int southeastback3d(int i, int j, int k, int n);
extern int southeastfront3d(int i, int j, int k, int n);
extern int northeastback3d(int i, int j, int k, int n);
extern int northeastfront3d(int i, int j, int k, int n);
extern void compress_matrix(double val[], int bindx[], int N_update);




int *Ke_data_org = NULL, *Kn_data_org = NULL, *global_edge_inds = NULL,
  *global_node_inds = NULL, *global_edge_externs = NULL,
  *global_node_externs = NULL;
int    *reordered_glob_edges = NULL, *reordered_edge_externs = NULL;
int    *reordered_glob_nodes = NULL, *reordered_node_externs = NULL;
int    *cpntr = NULL, *Ke_bindx = NULL, *Kn_bindx = NULL, Nlocal_edges, iii, *Tmat_bindx;
int    *update, *update_index;
int *external, *extern_index;
struct reader_context *context;
int mls_poly_degree, eig_ratio_tol = 27.;
int reduced_smoother_flag = 0;

#ifdef debugSmoother
double *xxx;
#endif
/******************************************************************************/
extern int  ML_make_block_matrix(ML_Operator *blockmat, ML_Operator *original1,
                          ML_Operator *original2);

/*
extern int ML_Gen_Hierarchy_ComplexMaxwell(ML *ml_edges,ML_Operator **Tmat_array, 
					   ML_Operator **Tmat_trans_array,
					   ML **newml, ML_Operator *M);
*/
struct ml_operator_wrapper {
  int (*diag_matvec)(void *, int, double *, int, double *);
  int (*diag_getrow)(void *,int,int*,int,int*,double*,int*);
  void *diag_matvec_data, *diag_getrow_data;
  int (*offdiag_matvec)(void *, int, double *, int, double *);
  int (*offdiag_getrow)(void *,int,int*,int,int*,double*,int*);
  void *offdiag_matvec_data, *offdiag_getrow_data;
  /* work vectors for block matrix getrow */
  int *cols;
  double *vals;
};
void aztec_block_matvec(double *x, double *y, AZ_MATRIX *Amat,
                             int proc_config[]);
struct aztec_block_data {
  int N;
  AZ_MATRIX *Ke;
  AZ_MATRIX *M;
};
ML_Operator *globbie;

/***************************************************************************/

int main(int argc, char *argv[])
{

  int num_PDE_eqns=1, N_levels=3, nsmooth=2;
  int Nglobal_edges, Nglobal_nodes, Nlocal_edges, Nlocal_nodes;
  int level, coarsest_level;

  /* See Aztec User's Guide for more information on the */
  /* variables that follow.                             */

  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE], status[AZ_STATUS_SIZE];

  /* data structure for matrix corresponding to the fine grid */

  double *Ke_val = NULL, *Kn_val, *Tmat_val = NULL, *rhs, solve_time,
    setup_time, start_time, *yyy, *vvv, *zzz;
  /* A Josh hack - double *blockrhs */
  double *blockrhs;
  
  double *M_val = NULL;
  int *M_bindx = NULL, *M_data_org = NULL;
  AZ_MATRIX *M_mat, *blockmat;
  ML *ml_M;
  ML_Operator *Mmat_ml;
  struct aztec_block_data aztec_block_data;
  int *junk1 = NULL, *junk2 = NULL, *junk3 = NULL, *junk4 = NULL;
  AZ_MATRIX *Ke_mat, *Kn_mat;
  ML_Operator *Tmat, *Tmat_trans, *Tmat_transbc=NULL, *Amat;

  /* Operator for use in Hiptmair smoother only */
  ML_Operator *Tmatbc = NULL;
  double *Tmatbc_val, sum;
  int *Tmatbc_bindx;

  /* Used in zeroing out rows of Tmat. */
  int allocated_space = 0, *cols = NULL, length = 0;
  double *vals = NULL;
  struct ML_CSR_MSRdata *data;
  int *row_ptr;
  double *val_ptr;
  int *BCindices, BCcount;



  ML_Operator **Tmat_array, **Tmat_trans_array;
  AZ_PRECOND *Pmat = NULL;
  ML *ml_edges, *ml_nodes;
  ML *ml_block;
  FILE *fp;
  int i, j, Nrigid, *garbage= NULL, nblocks, *blocks = NULL;
  struct AZ_SCALING *scaling;
  ML_Aggregate *ag;
  double *mode, *rigid, alpha;
  char filename[80];
  int    one = 1;
  int lower, Ncols;
  struct ML_CSR_MSRdata *csr_data;
  int *bindx = NULL;
  int row_length;
  double dtemp, dtemp2;
  int Nexterns;
  int Nnz, nx, offset, ii, jj, kk, index, horv, stencil[40],off2set, nz_ptr;
  double sigma, dcenter, doffdiag, sten_val[40], c1, c2, c3, c4, c5;
  int mg_cycle_type;
  double omega, nodal_omega, edge_omega;
  void **edge_args, **nodal_args, *edge_smoother, *nodal_smoother;
  int  edge_its, nodal_its, *temp1, *temp2;
  double edge_eig_ratio[25], nodal_eig_ratio[25];
#ifndef debugSmoother
  double *xxx;
#endif
  int Ndirichlet, *dirichlet;
  struct ML_Operator_blockmat_data *block_data;
   ML_Krylov   *kdata;
   double *M_diag, *diag_of_M, *diag_of_S, lambda, tdiag, ttdiag;


#ifdef ML_partition
  FILE *fp2;
  int *block_list=NULL;
  int count, *pcounts, *proc_id;
  int Tmat_size, *proc_assignment, p1, p2, col1, col2, row_start;
  int itemp;
  int evenodd[2]; int Iwin[4];
  int tiebreaker;

  nblocks = -1;
  if (argc < 2) {
    printf("Usage: ml_readfile num_processors\n");
    exit(1);
  }
  else sscanf(argv[1],"%d",&nblocks);
  if (nblocks == -1) {
    printf("Usage: ml_readfile num_processors\n");
    exit(1);
  }
#endif /* ifdef ML_partition */
  temp1 = (int *) malloc(sizeof(int)*4);
  temp2 = (int *) malloc(sizeof(int)*4);

#ifdef ML_MPI
  MPI_Init(&argc,&argv);

  /* get number of processors and the name of this processor */

  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif
  my_proc_id = proc_config[AZ_node];

#ifndef ML_partition
  if (proc_config[AZ_node] == 0)
  {
    printf("Reading settings from 'ml_inputfile'\n");
    ML_Reader_ReadInput("ml_inputfile", &context);
  }
  else context = (struct reader_context *) malloc(sizeof(struct reader_context));
  AZ_broadcast((char *) context,  sizeof(struct reader_context), proc_config,
	       AZ_PACK);
  AZ_broadcast((char *) NULL        ,   0          , proc_config, AZ_SEND);

  N_levels = context->N_levels;
  nsmooth   = context->nsmooth;
  num_PDE_eqns = context->N_dofPerNode;
  ML_Set_PrintLevel(context->output_level);
#else
  context = (struct reader_context *) malloc(sizeof(struct reader_context));
  ML_Reader_InitContext(context);
#endif /* ifndef ML_partition */

  /* read in the number of edge matrix equations */
  Nglobal_edges = 0;
  Nglobal_nodes = 0;

#if defined(HARDWIRE3D) || defined(HARDWIRE2D)
  if (proc_config[AZ_node] == 0)
  {
    printf("Enter the total number of elements on a side\n");
    scanf("%d",&Nglobal_edges);
#ifdef HARDWIRE3D
    Nglobal_nodes = Nglobal_edges*Nglobal_edges*Nglobal_edges;
Nglobal_nodes = (Nglobal_edges+1)*(Nglobal_edges+1)*(Nglobal_edges+1); /* rst dirichlet */
    Nglobal_edges = Nglobal_edges*Nglobal_edges*Nglobal_edges*3;

#else
    Nglobal_nodes = Nglobal_edges*Nglobal_edges;
    Nglobal_edges = Nglobal_edges*Nglobal_edges*2;
#endif
  }
#else
  if (proc_config[AZ_node] == 0) {
#    ifdef binary
    fp=fopen("Ke_mat.az","rb");
#    else
    fp=fopen("Ke_mat.az","r");
#    endif
    if (fp==NULL) {
      printf("couldn't open file Ke_mat.az\n");
      exit(1);
    }
#    ifdef binary
    fread(&Nglobal_edges, sizeof(int), 1, fp);
#    else
    fscanf(fp,"%d",&Nglobal_edges);
#    endif
    fclose(fp);
  }

  /* read in the number of node matrix equations */

  if (proc_config[AZ_node] == 0) {
#    ifdef binary
    fp=fopen("Kn_mat.az","rb");
#    else
    fp=fopen("Kn_mat.az","r");
#    endif
    if (fp==NULL) {
      printf("couldn't open file Kn_mat.az\n");
      exit(1);
    }
#    ifdef binary
    fread(&Nglobal_nodes, sizeof(int), 1, fp);
#    else
    fscanf(fp,"%d",&Nglobal_nodes);
#    endif
    fclose(fp);
  }
#endif
  Nglobal_edges = AZ_gsum_int(Nglobal_edges, proc_config);
  Nglobal_nodes = AZ_gsum_int(Nglobal_nodes, proc_config);

  /*******************************************************************/
  /* initialize the list of global indices indicating which rows are */
  /* stored on which processor. NOTE: the list of global             */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly.                         */
  /*-----------------------------------------------------------------*/

  printf("N procs = %d\n",proc_config[AZ_N_procs]);
  if (proc_config[AZ_N_procs] == 1) i = AZ_linear;
  else i = AZ_file;
  AZ_input_update("edge_partition",&Nlocal_edges, &global_edge_inds, 
		  proc_config, Nglobal_edges, num_PDE_eqns,i);

  /********************************************************************/
  /*                      Set up Ke_mat                               */
  /*  1) read it in using Aztec.                                      */
  /*  2) call AZ_transform to convert global column indices to local  */
  /*     indices and to set up Aztec's communication structure.       */
  /*  3) Stuff the arrays into an Aztec matrix.                       */
  /*------------------------------------------------------------------*/

  if (proc_config[AZ_node] == 0)
  {
     printf("Reading edge matrix.\n"); fflush(stdout);
  }
#if defined(HARDWIRE3D) || defined(HARDWIRE2D)
  sigma = 1.0e-14;
  sigma = 0.000000001;
  sigma = 0.001;
#ifdef HARDWIRE2D
  nx = (int) sqrt( ((double) Nglobal_nodes) + .00001);
  Ke_bindx = (int    *) malloc((8*Nlocal_edges+5)*sizeof(int));
  Ke_val   = (double *) malloc((8*Nlocal_edges+5)*sizeof(double));
  Ke_bindx[0] = Nlocal_edges+1;

  dcenter  = 2 + 2.*sigma/((double) ( 3 * nx * nx));
  doffdiag = -1 + sigma/((double) ( 6 * nx * nx));

  for (i = 0; i < Nlocal_edges; i++) {
    Ke_val[global_edge_inds[i]] = dcenter;
    Ke_bindx[i+1] = Ke_bindx[i] + 6;
    inv2dindex(global_edge_inds[i], &ii, &jj, nx, &horv);
    stencil[0] = north2d(ii,jj,nx);
    stencil[1] = west2d(ii,jj,nx);
    stencil[2] = east2d(ii,jj,nx);
    stencil[3] = south2d(ii,jj,nx);
    if (horv == HORIZONTAL) {
      sten_val[0] = doffdiag;
      sten_val[1] =  1.;
      sten_val[2] = -1.;
      sten_val[3] = dcenter;
    }
    else {
      sten_val[0] = -1.;
      sten_val[1] =  dcenter;
      sten_val[2] = doffdiag;
      sten_val[3] = 1.;
    }
    if (horv == HORIZONTAL) {
      if (jj == 0) jj = nx-1;
      else jj--;
    }
    else {
      if (ii == 0) ii = nx-1;
      else ii--;
    }
    stencil[4] = west2d(ii,jj,nx);
    stencil[5] = south2d(ii,jj,nx);
    if (horv == HORIZONTAL) stencil[6] = east2d(ii,jj,nx);
    else stencil[6] = north2d(ii,jj,nx);
    if (horv == HORIZONTAL) {
      sten_val[4] = -1.;
      sten_val[5] = doffdiag;
      sten_val[6] =  1.;
    }
    else {
      sten_val[4] = doffdiag;
      sten_val[5] = -1.;
      sten_val[6] =  1.;
    }
    AZ_sort(stencil,7,NULL,sten_val);

    off2set = 0;
    for (ii = 0; ii < 7; ii++) {
      if (stencil[ii] != global_edge_inds[i]) {
	Ke_val[Ke_bindx[global_edge_inds[i]]+off2set] = sten_val[ii];
	Ke_bindx[Ke_bindx[global_edge_inds[i]]+off2set] = stencil[ii];
	off2set++;
      }
    }
  }
#else
  nx = (int) pow( ((double) Nglobal_edges/3) + .00001,.333333333334);
  Ke_bindx = (int    *) malloc((35*Nlocal_edges+5)*sizeof(int));
  Ke_val   = (double *) malloc((35*Nlocal_edges+5)*sizeof(double));
  Ke_bindx[0] = Nlocal_edges+1;
  nz_ptr = Ke_bindx[0];

  dcenter  = 2 + 2.*sigma/((double) ( 3 * nx * nx * nx));
  doffdiag = -1 + sigma/((double) ( 6 * nx * nx * nx));


  c1 = 1./6.;
  c2 = 2./3.;
  c3 = 8./3. + 16.*sigma/((double) (36 * nx *nx));
  c4 = -1./3. + 4.*sigma/((double) (36 * nx *nx));
  c5 = -1./3. + sigma/((double) (36 * nx * nx));

  for (i = 0; i < Nlocal_edges; i++) {
    Ke_val[global_edge_inds[i]] = c3;
    if (Ke_bindx[i] != nz_ptr) {
      printf("problems %d %d\n",i,nz_ptr); exit(1);
    }
    Ke_bindx[i+1] = Ke_bindx[i] + 32;
    inv3dindex(global_edge_inds[i], &ii, &jj, &kk, nx, &horv);

    if (horv == HORIZONTAL) {
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] = northback3d(ii,jj,kk,nx); /* H -,1,0 */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] = westback3d(ii,jj,kk,nx);  /* V 0,0,0 */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx);  /* V 1,0,0 */
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  northfront3d(ii,jj,kk,nx);/* H -,1,1 */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  westfront3d(ii,jj,kk,nx); /* V 0,0,1 */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  eastfront3d(ii,jj,kk,nx); /* V 1,0,1 */
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx); /* H -,0,1 */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx);  /* I 0,0,0 */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx);  /* I 1,0,0 */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx); /* I 0,1,0 */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northeast3d(ii,jj,kk,nx); /* I 1,1,0 */

      if (jj == 0) jj = nx - 1;
      else jj--;
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx);  /* H -,-1,0 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx);   /* V 0,-1,0 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx);   /* V 1,-1,0 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  westfront3d(ii,jj,kk,nx); /* V 0,-1,1 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  eastfront3d(ii,jj,kk,nx); /* V 1,-1,1 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx); /* H -1,1 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx); /* I 0,-1,0 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx); /* I 1,-1,0 */
if (jj == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */

      if (kk == 0) kk = nx - 1;
      else kk--;
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx); /* H -1,-1 */
if ((jj == nx-1) || (kk == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx); /* V 0,-1,-1 */
if ((jj == nx-1) || (kk == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx); /* V 1,-1,-1 */
if ((jj == nx-1) || (kk == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  northback3d(ii,jj,kk,nx); /* H 0,-1 */
if ((jj == nx-1) || (kk == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx); /* I 0,-1,-1 */
if ((jj == nx-1) || (kk == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx); /* I 1,-1,-1 */
if ((jj == nx-1) || (kk == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */

      if (jj == nx-1) jj = 0;
      else jj++;
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  northback3d(ii,jj,kk,nx); /* H 1,-1 */
if (kk == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx); /* V 0,0,-1 */
if (kk == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx); /* V 1,0,-1 */
if (kk == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx); /* I 0,0,-1 */
if (kk == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx); /* I 1,0,-1 */
if (kk == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx); /* I 0,1,-1 */
if (kk == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northeast3d(ii,jj,kk,nx); /* I 1,1,-1 */
if (kk == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
 if (kk == nx-1) kk = 0;
 else kk++;

if ((jj==0) || (kk==0)) { /* rst dirichlet */
  for (ii = Ke_bindx[i] ; ii < Ke_bindx[i+1]; ii++) Ke_bindx[ii]=-1;
}
    }

    if (horv == VERTICAL) {
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] = eastback3d(ii,jj,kk,nx); /* H -,1,0 */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] = northback3d(ii,jj,kk,nx);  /* V 0,0,0 */
      Ke_val[nz_ptr] =  c2;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx);  /* V 1,0,0 */

      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  eastfront3d(ii,jj,kk,nx);/* H -,1,1 */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northfront3d(ii,jj,kk,nx); /* V 0,0,1 */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx); /* V 1,0,1 */
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  westfront3d(ii,jj,kk,nx); /* H -,0,1 */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx);  /* I 0,0,0 */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx);  /* I 1,0,0 */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northeast3d(ii,jj,kk,nx); /* I 0,1,0 */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx); /* I 1,1,0 */

      if (ii == 0) ii = nx - 1;
      else ii--;
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx);  /* H -,-1,0 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  northback3d(ii,jj,kk,nx);   /* V 0,-1,0 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx);   /* V 1,-1,0 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northfront3d(ii,jj,kk,nx); /* V 0,-1,1 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx); /* V 1,-1,1 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  westfront3d(ii,jj,kk,nx); /* H -1,1 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx); /* I 0,-1,0 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx); /* I 1,-1,0 */
if (ii == nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      if (kk == 0) kk = nx - 1;
      else kk--;
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx); /* H -1,-1 */
if ((kk==nx-1) || (ii == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northback3d(ii,jj,kk,nx); /* V 0,-1,-1 */
if ((kk==nx-1) || (ii == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx); /* V 1,-1,-1 */
if ((kk==nx-1) || (ii == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx); /* H 0,-1 */
if ((kk==nx-1) || (ii == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx); /* I 0,-1,-1 */
if ((kk==nx-1) || (ii == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx); /* I 1,-1,-1 */
if ((kk==nx-1) || (ii == nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */

      if (ii == nx-1) ii = 0;
      else ii++;
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx); /* H 1,-1 */
if (kk==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northback3d(ii,jj,kk,nx); /* V 0,0,-1 */
if (kk==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx); /* V 1,0,-1 */
if (kk==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */

      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx); /* I 0,0,-1 */
if (kk==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx); /* I 1,0,-1 */
if (kk==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northeast3d(ii,jj,kk,nx); /* I 0,1,-1 */
if (kk==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx); /* I 1,1,-1 */
if (kk==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
if (kk == nx-1) kk = 0;
else kk++;
if ((ii==0) || (kk==0)) { /* rst dirichlet */
  for (jj = Ke_bindx[i] ; jj < Ke_bindx[i+1]; jj++) Ke_bindx[jj]=-1;
}
    }

    if (horv == OUT) {
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] = northwest3d(ii,jj,kk,nx); /* H -,1,0 */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] = westfront3d(ii,jj,kk,nx);  /* V 0,0,0 */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx);  /* V 1,0,0 */
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  northeast3d(ii,jj,kk,nx);/* H -,1,1 */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  eastfront3d(ii,jj,kk,nx); /* V 0,0,1 */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx); /* V 1,0,1 */
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx); /* H -,0,1 */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx);  /* I 0,0,0 */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx);  /* I 1,0,0 */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northfront3d(ii,jj,kk,nx); /* I 0,1,0 */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northback3d(ii,jj,kk,nx); /* I 1,1,0 */

      if (jj == 0) jj = nx - 1;
      else jj--;
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx);  /* H -,-1,0 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  westfront3d(ii,jj,kk,nx);   /* V 0,-1,0 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx);   /* V 1,-1,0 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  eastfront3d(ii,jj,kk,nx); /* V 0,-1,1 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  eastback3d(ii,jj,kk,nx); /* V 1,-1,1 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  southeast3d(ii,jj,kk,nx); /* H -1,1 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx); /* I 0,-1,0 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx); /* I 1,-1,0 */
if (jj==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */

      if (ii == 0) ii = nx - 1;
      else ii--;
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  southwest3d(ii,jj,kk,nx); /* H -1,-1 */
if ((ii==nx-1)||(jj==nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  westfront3d(ii,jj,kk,nx); /* V 0,-1,-1 */
if ((ii==nx-1)||(jj==nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx); /* V 1,-1,-1 */
if ((ii==nx-1)||(jj==nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c4;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx); /* H 0,-1 */
if ((ii==nx-1)||(jj==nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx); /* I 0,-1,-1 */
if ((ii==nx-1)||(jj==nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx); /* I 1,-1,-1 */
if ((ii==nx-1)||(jj==nx-1)) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */

      if (jj == nx-1) jj = 0;
      else jj++;
      Ke_val[nz_ptr] = c5;
      Ke_bindx[nz_ptr++] =  northwest3d(ii,jj,kk,nx); /* H 1,-1 */
if (ii==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  westfront3d(ii,jj,kk,nx); /* V 0,0,-1 */
if (ii==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  westback3d(ii,jj,kk,nx); /* V 1,0,-1 */
if (ii==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c2;
      Ke_bindx[nz_ptr++] =  southfront3d(ii,jj,kk,nx); /* I 0,0,-1 */
if (ii==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c2;
      Ke_bindx[nz_ptr++] =  southback3d(ii,jj,kk,nx); /* I 1,0,-1 */
if (ii==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = c1;
      Ke_bindx[nz_ptr++] =  northfront3d(ii,jj,kk,nx); /* I 0,1,-1 */
if (ii==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
      Ke_val[nz_ptr] = -c1;
      Ke_bindx[nz_ptr++] =  northback3d(ii,jj,kk,nx); /* I 1,1,-1 */
if (ii==nx-1) Ke_bindx[nz_ptr-1] = -1; /* rst dirichlet */
if (ii == nx-1) ii = 0;
else ii++;
if ((jj==0) || (ii==0)) { /* rst dirichlet */
  for (kk = Ke_bindx[i] ; kk < Ke_bindx[i+1]; kk++) Ke_bindx[kk]=-1;
}
    }
  }
  compress_matrix(Ke_val, Ke_bindx, Nlocal_edges);

#endif
#else
  AZ_input_msr_matrix("Ke_mat.az", global_edge_inds, &Ke_val, &Ke_bindx, 
                      Nlocal_edges, proc_config);
#endif

  if (proc_config[AZ_node] == 0)
  {
     printf("Done reading edge matrix.\n"); fflush(stdout);
  }

  AZ_transform_norowreordering(proc_config, &global_edge_externs, Ke_bindx,
			       Ke_val, global_edge_inds, &reordered_glob_edges, 
			       &reordered_edge_externs, &Ke_data_org, Nlocal_edges, 0, 0, 0, 
			       &cpntr,	       AZ_MSR_MATRIX);


  update = global_edge_inds;
  update_index = reordered_glob_edges;

  Ke_mat = AZ_matrix_create( Nlocal_edges );
  AZ_set_MSR(Ke_mat, Ke_bindx, Ke_val, Ke_data_org, 0, NULL, AZ_LOCAL);
  Ke_mat->matrix_type    = Ke_data_org[AZ_matrix_type];
  Ke_data_org[AZ_N_rows] = Ke_data_org[AZ_N_internal] + 
    Ke_data_org[AZ_N_border];

  /* rst: Read in the mass matrix */

  AZ_input_msr_matrix("M_mat.az", global_edge_inds, &M_val, &M_bindx, 
                      Nlocal_edges, proc_config);
  AZ_transform_norowreordering(proc_config, &junk1, M_bindx,
			       M_val, global_edge_inds, &junk2,
			       &junk3, &M_data_org, Nlocal_edges, 0, 0, 0, 
			       &junk4,	       AZ_MSR_MATRIX);

  /* Create an Aztec matrix wrapper */
  M_mat = AZ_matrix_create( Nlocal_edges );
  AZ_set_MSR(M_mat, M_bindx, M_val, M_data_org, 0, NULL, AZ_LOCAL);

  M_mat->matrix_type    = M_data_org[AZ_matrix_type];
  M_data_org[AZ_N_rows] = M_data_org[AZ_N_internal] + M_data_org[AZ_N_border];

  /* rst: Make an ML matrix. Unfortunately, we have to create a bogus */
  /* hierarchy to make this work. After these two lines,         */
  /* &(ml_M->Amat[0]) is an ML matrix corresponding to M.        */
  //leak!!!
  ML_Create(&ml_M,1);
  AZ_ML_Set_Amat(ml_M, 0, Nlocal_edges, Nlocal_edges, M_mat, proc_config);
  Mmat_ml = &(ml_M->Amat[0]);

  /* Set rhs */

  /*rhs=(double *) ML_allocate(Nlocal_edges*sizeof(double));*/
 
  if(0)
  {
    printf("%d: rhsfile file pointer is NULL\n",proc_config[AZ_node]); fflush(stdout);
    if (proc_config[AZ_node] == 0 && 0.5 < ML_Get_PrintLevel())
       printf("taking zero vector for rhs\n");
    fflush(stdout);
    rhs = (double *)
	      ML_allocate(2*(Nlocal_edges + Ke_mat->data_org[AZ_N_external])
                      *sizeof(double)); 
    for (i = 0; i < 2*Nlocal_edges; i++) rhs[i] = 0.0;
    rhs[Nlocal_edges/2] = 1.;
    rhs[Nlocal_edges/2+1] = -1.;
    xxx = (double *) ML_allocate(2*(3 + Nlocal_edges + Ke_mat->data_org[AZ_N_external])*sizeof(double));

#ifdef garbage
for (i = 0; i < -Nlocal_edges; i++) {
  inv3dindex(global_edge_inds[i], &ii, &jj, &kk, nx, &horv);
  rhs[i] =
    sin(( (double) ii*1)/( (double) nx) ) *
    sin(( (double) jj*1)/( (double) nx) ) *
    sin(( (double) kk*1)/( (double) nx) );
    if (horv == HORIZONTAL) {
      if ((jj==0) || (kk==0)) rhs[i] = 0.;/* rst dirichlet */
    }
    if (horv == VERTICAL) {
      if ((ii==0) || (kk==0)) rhs[i] = 0.;/* rst dirichlet */
    }
    if (horv == OUT) {
      if ((ii==0) || (jj==0)) rhs[i] = 0.;/* rst dirichlet */
    }
}
#endif
printf("should not be doing this\n");

     AZ_random_vector(rhs, Ke_data_org, proc_config);
for (i = 0; i < Nlocal_edges; i++) rhs[i+Nlocal_edges] = rhs[i];

/*     AZ_random_vector(&(rhs[Nlocal_edges]), Ke_data_org, proc_config); */

Ke_mat->matvec(rhs, xxx, Ke_mat, proc_config);
/* for (i = 0; i < Nlocal_edges; i++) rhs[i] = xxx[i]; */
free(xxx);
  }

  else
  {
    rhs = (double *)
	  ML_allocate(2*(Nlocal_edges + Ke_mat->data_org[AZ_N_external])
          *sizeof(double)); 
/* read in real rhs here */
    fp = fopen("rhsfile_re","r");
    if (proc_config[AZ_node] == 0)
    {
       printf("%d: reading real part of rhs from a file\n",proc_config[AZ_node]);
       fflush(stdout);
    }
    AZ_input_msr_matrix("rhsfile_re", global_edge_inds, &rhs, &garbage,
                                    Nlocal_edges, proc_config);
    
    dtemp = sqrt(abs(ML_gdot(Nlocal_edges, rhs, rhs, ml_M->comm)));
    printf("\n\n\n%20.15e\n\n\n",dtemp);


    fclose(fp);
    if (proc_config[AZ_node] == 0)
    {
       printf("%d: Done reading real part of rhs from a file\n",proc_config[AZ_node]);
       fflush(stdout);
    }

/* create block rhs 
   double *blockrhs; */
   blockrhs = (double *)
	      ML_allocate(2*(Nlocal_edges + Ke_mat->data_org[AZ_N_external])
                      *sizeof(double)); 
   for (i=0; i < Nlocal_edges; i++) blockrhs[i] = rhs[i];

/* read in imaginary rhs here */
    fp = fopen("rhsfile_im","r");
    if (proc_config[AZ_node] == 0)
    {
       printf("%d: reading imaginary part of rhs from a file\n",proc_config[AZ_node]);
       fflush(stdout);
    }
    AZ_input_msr_matrix("rhsfile_im", global_edge_inds, &rhs, &garbage,
                                    Nlocal_edges, proc_config);

    fclose(fp);
    if (proc_config[AZ_node] == 0)
    {
       printf("%d: Done reading imaginary part of rhs from a file\n",proc_config[AZ_node]);
       fflush(stdout);
    }
  }
//  }
/* and shove into blockrhs */

   for (i=0; i < Nlocal_edges; i++) blockrhs[i+Nlocal_edges+Ke_mat->data_org[AZ_N_external]] = rhs[i];

   ML_free(rhs);
   rhs = blockrhs;

/*
#define ZEROOUTDIRICHLET
*/
#ifdef ZEROOUTDIRICHLET
  if (proc_config[AZ_node] == 0 && 0.5 < ML_Get_PrintLevel() )
     printf("Zeroing out Dirichlet columns\n");
  AZ_zeroDirichletcolumns(Ke_mat, rhs, proc_config);
#else
  if (proc_config[AZ_node] == 0 && 0.5 < ML_Get_PrintLevel() )
     printf("Not zeroing out Dirichlet columns\n");
#endif

  /*******************************************************************/
  /* initialize the list of global indices indicating which rows are */
  /* stored on which processor. NOTE: the list of global             */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly.                         */
  /*-----------------------------------------------------------------*/

  if (proc_config[AZ_N_procs] == 1) i = AZ_linear;
  else i = AZ_file;
  AZ_input_update("node_partition",&Nlocal_nodes, &global_node_inds, 
		  proc_config, Nglobal_nodes, num_PDE_eqns,i);

  /********************************************************************/
  /*                      Set up Kn_mat                               */
  /*  1) read it in using Aztec.                                      */
  /*  2) call AZ_transform to convert global column indices to local  */
  /*     indices and to set up Aztec's communication structure.       */
  /*  3) Stuff the arrays into an Aztec matrix.                       */
  /*------------------------------------------------------------------*/

  if (proc_config[AZ_node] == 0)
  {
     printf("Reading node matrix.\n"); fflush(stdout);
  }
#if defined(HARDWIRE3D) || defined(HARDWIRE2D)
  Kn_bindx = (int    *) malloc((27*Nlocal_nodes+5)*sizeof(int));
  Kn_val   = (double *) malloc((27*Nlocal_nodes+5)*sizeof(double));
  Kn_bindx[0] = Nlocal_nodes+1;
#ifdef HARDWIRE2D
  nx = (int) sqrt( ((double) Nglobal_nodes) + .00001);
  for (i = 0; i < Nlocal_nodes; i++) {
    Kn_bindx[i+1] = Kn_bindx[i] + 8;
    Kn_val[i] = 8.;
    inv2dnodeindex(global_node_inds[i], &ii, &jj, nx);
    Kn_bindx[Kn_bindx[i]+0] =  southeast2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+0] = -1.;
    Kn_bindx[Kn_bindx[i]+1] =  northwest2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+1] = -1.;
    Kn_bindx[Kn_bindx[i]+4] =  northeast2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+4] = -.00000001;
    if (ii == 0) ii = nx-1;
    else ii--;
    Kn_bindx[Kn_bindx[i]+5] =  northwest2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+5] = -.00000001;
    if (jj == 0) jj = nx-1;
    else jj--;
    Kn_bindx[Kn_bindx[i]+2] =  southeast2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+2] = -1.;
    Kn_bindx[Kn_bindx[i]+3] =  northwest2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+3] = -1.;
    Kn_bindx[Kn_bindx[i]+6] =  southwest2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+6] = -.00000001;
    if (ii == nx-1) ii = 0;
    else ii++;
    Kn_bindx[Kn_bindx[i]+7] =  southeast2d(ii,jj,nx);
    Kn_val[  Kn_bindx[i]+7] = -.00000001;
    AZ_sort( &(Kn_bindx[Kn_bindx[i]]),8,NULL,&(Kn_val[Kn_bindx[i]]));
  }
#else
  nx = (int) pow( ((double) Nglobal_nodes) + .00001,.3333333333334);
  for (i = 0; i < Nlocal_nodes; i++) {
    Kn_bindx[i+1] = Kn_bindx[i] + 26;
    Kn_val[i] = 26.;
    for (ii = 0; ii < 26; ii++)  Kn_val[Kn_bindx[i]+ii] = -1.;
    for (ii = 0; ii < 26; ii++)  Kn_bindx[Kn_bindx[i]+ii] = i;
    inv3dnodeindex(global_node_inds[i], &ii, &jj, &kk, nx);

    Kn_bindx[Kn_bindx[i]+0] = northwestback3d(ii,jj,kk,nx);
    Kn_bindx[Kn_bindx[i]+1] = southeastback3d(ii,jj,kk,nx);
    Kn_bindx[Kn_bindx[i]+2] = southwestfront3d(ii,jj,kk,nx);  
    Kn_bindx[Kn_bindx[i]+6] = northeastback3d(ii,jj,kk,nx);  
    Kn_bindx[Kn_bindx[i]+7] = southeastfront3d(ii,jj,kk,nx);  
    Kn_bindx[Kn_bindx[i]+8] = northeastfront3d(ii,jj,kk,nx);  
    Kn_bindx[Kn_bindx[i]+9] = northwestfront3d(ii,jj,kk,nx);  
    if (ii == 0) ii = nx -1;
    else ii--;
    Kn_bindx[Kn_bindx[i]+10] = northwestfront3d(ii,jj,kk,nx);  
if (ii==nx-1) Kn_bindx[Kn_bindx[i]+10] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+11] = northwestback3d(ii,jj,kk,nx);  
if (ii==nx-1) Kn_bindx[Kn_bindx[i]+11] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+12] = southwestfront3d(ii,jj,kk,nx);  
if (ii==nx-1) Kn_bindx[Kn_bindx[i]+12] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+13] = southwestback3d(ii,jj,kk,nx);  
if (ii==nx-1) Kn_bindx[Kn_bindx[i]+13] = -1; /* rst dirichlet */

    if (jj == 0) jj = nx -1;
    else jj--;
    Kn_bindx[Kn_bindx[i]+14] = southwestback3d(ii,jj,kk,nx);  
if ((ii==nx-1) || (jj==nx-1)) Kn_bindx[Kn_bindx[i]+14] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+15] = southwestfront3d(ii,jj,kk,nx);  
if ((ii==nx-1) || (jj==nx-1)) Kn_bindx[Kn_bindx[i]+15] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+16] = southeastfront3d(ii,jj,kk,nx);  
if (ii==nx-1) Kn_bindx[Kn_bindx[i]+16] = southwestfront3d(0,jj,kk,nx);  
if (jj==nx-1) Kn_bindx[Kn_bindx[i]+16] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+17] = southeastback3d(ii,jj,kk,nx);  
if (ii==nx-1) Kn_bindx[Kn_bindx[i]+17] = southwestback3d(0,jj,kk,nx);  
if (jj==nx-1) Kn_bindx[Kn_bindx[i]+17] = -1; /* rst dirichlet */

    if (ii == nx-1) ii = 0;
    else ii++;
    Kn_bindx[Kn_bindx[i]+18] = southeastback3d(ii,jj,kk,nx);  
if (jj==nx-1) Kn_bindx[Kn_bindx[i]+18] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+19] = southeastfront3d(ii,jj,kk,nx);  
if (jj==nx-1) Kn_bindx[Kn_bindx[i]+19] = -1; /* rst dirichlet */

    if (kk == 0) kk = nx -1;
    else kk--;
    Kn_bindx[Kn_bindx[i]+20] = northeastback3d(ii,jj,kk,nx);
if (jj==nx-1) Kn_bindx[Kn_bindx[i]+20] = southeastback3d(ii,0,kk,nx);
if (kk==nx-1) Kn_bindx[Kn_bindx[i]+20] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+21] = northwestback3d(ii,jj,kk,nx);
if (jj==nx-1) Kn_bindx[Kn_bindx[i]+21] = southwestback3d(ii,0,kk,nx);
if (kk==nx-1) Kn_bindx[Kn_bindx[i]+21] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+22] = southeastback3d(ii,jj,kk,nx);
if ((jj==nx-1)||(kk==nx-1)) Kn_bindx[Kn_bindx[i]+22] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+23] = southwestback3d(ii,jj,kk,nx);
if ((jj==nx-1)||(kk==nx-1)) Kn_bindx[Kn_bindx[i]+23] = -1; /* rst dirichlet */

    if (ii == 0) ii = nx -1;
    else ii--;
    Kn_bindx[Kn_bindx[i]+24] = southwestback3d(ii,jj,kk,nx);
if ((ii==nx-1)||(jj==nx-1)||(kk==nx-1)) Kn_bindx[Kn_bindx[i]+24] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+25] = northwestback3d(ii,jj,kk,nx);
if (jj==nx-1)  Kn_bindx[Kn_bindx[i]+25] = southwestback3d(ii,0,kk,nx);
if ((ii==nx-1)||(kk==nx-1)) Kn_bindx[Kn_bindx[i]+25] = -1; /* rst dirichlet */
    if (jj == nx-1) jj = 0;
    else jj++;
    Kn_bindx[Kn_bindx[i]+3] = northwestback3d(ii,jj,kk,nx);
if ((ii==nx-1)||(kk==nx-1)) Kn_bindx[Kn_bindx[i]+3] = -1; /* rst dirichlet */
    Kn_bindx[Kn_bindx[i]+4] = northeastback3d(ii,jj,kk,nx);
if(ii==nx-1) Kn_bindx[Kn_bindx[i]+4] = northwestback3d(0,jj,kk,nx);
if (kk==nx-1) Kn_bindx[Kn_bindx[i]+4] = -1; /* rst dirichlet */
    if (ii == nx-1) ii = 0;
    else ii++;
    Kn_bindx[Kn_bindx[i]+5] = northeastback3d(ii,jj,kk,nx);
if (kk==nx-1) Kn_bindx[Kn_bindx[i]+5] = -1; /* rst dirichlet */
  }
#endif
  compress_matrix(Kn_val, Kn_bindx, Nlocal_nodes);

  for (i = 0; i < Nlocal_nodes; i++) {
    sum = 0.;
    for (j = Kn_bindx[i]; j < Kn_bindx[i+1]; j++) {
      sum += Kn_val[j];
    }
    Kn_val[i] = -sum;
  }

#else
  AZ_input_msr_matrix("Kn_mat.az", global_node_inds, &Kn_val, &Kn_bindx, 
		      Nlocal_nodes, proc_config);
  if (proc_config[AZ_node] == 0)
  {
     printf("Done reading node matrix.\n"); fflush(stdout);
  }
#endif

  AZ_transform_norowreordering(proc_config, &global_node_externs, Kn_bindx, 
			       Kn_val, global_node_inds, &reordered_glob_nodes, 
			       &reordered_node_externs, &Kn_data_org, Nlocal_nodes, 0, 0, 0, 
			       &cpntr, AZ_MSR_MATRIX);
  Kn_mat = AZ_matrix_create( Nlocal_nodes );
  AZ_set_MSR(Kn_mat, Kn_bindx, Kn_val, Kn_data_org, 0, NULL, AZ_LOCAL);
  Kn_mat->matrix_type  = Kn_data_org[AZ_matrix_type];
  Kn_data_org[AZ_N_rows]  = Kn_data_org[AZ_N_internal] + Kn_data_org[AZ_N_border];

  /********************************************************************/
  /*                      Set up T_mat                                */
  /*  1) read it in using Aztec.                                      */
  /*  2) convert it to CSR                                            */
  /*      NOTE: we store the diagonal with the offdiagonals (wasting  */
  /*      the diagonal storage). This is so the manual mapping from   */
  /*      global to local indices does not get confusing.             */
  /*  3) Make it into an ML Operator                                  */
  /*     note: since ML_Operator_Create needs a communicator, I call  */
  /*     ML_Create() now to get one. I could have just created a      */
  /*     a communicator and done the ML_Create() later.               */
  /*------------------------------------------------------------------*/

  /* This copy of Tmat does not contain boundary conditions. */

  if (proc_config[AZ_node] == 0)
  {
     printf("Reading T matrix\n"); fflush(stdout);
  }
#if defined(HARDWIRE3D) || defined(HARDWIRE2D)
  nx = (int) pow( ((double) Nglobal_edges/3) + .00001,.333333333334);
#ifdef HARDWIRE2D 
  nx = (int) sqrt( ((double) Nglobal_nodes) + .00001);
#endif
  Tmat_bindx = (int    *) malloc((3*Nlocal_edges+5)*sizeof(int));
  Tmat_val   = (double *) malloc((3*Nlocal_edges+5)*sizeof(double));
  Tmat_bindx[0] = Nlocal_edges + 1;
  for (i = 0; i < Nlocal_edges; i++) {
    Tmat_bindx[i+1] = Tmat_bindx[i] + 2;
    Tmat_val[i] = 0.0;
#ifdef HARDWIRE2D
    inv2dindex(global_edge_inds[i], &ii, &jj, nx, &horv);
    if (horv == HORIZONTAL) {
      Tmat_bindx[Tmat_bindx[i]] = southwest2d(ii,jj,nx);
      Tmat_val[Tmat_bindx[i]] = -1.;
      Tmat_bindx[Tmat_bindx[i]+1] = southeast2d(ii,jj,nx);
      Tmat_val[Tmat_bindx[i]+1] = 1.;
    }
    else {
      Tmat_bindx[Tmat_bindx[i]] = northwest2d(ii,jj,nx);
      Tmat_val[Tmat_bindx[i]] = -1.;
      Tmat_bindx[Tmat_bindx[i]+1] = southwest2d(ii,jj,nx);
      Tmat_val[Tmat_bindx[i]+1] = 1.;
    }
#else
    inv3dindex(global_edge_inds[i], &ii, &jj, &kk, nx, &horv);
nx++; /* rst dirichlet */
    if (horv == HORIZONTAL) {
      Tmat_bindx[Tmat_bindx[i]] = southwestback3d(ii,jj,kk,nx); 
      Tmat_val[Tmat_bindx[i]] = -1.;
if ((ii==0)||(jj==0)||(kk==0)) Tmat_val[Tmat_bindx[i]]=0.; /* rst dirichlet */
      Tmat_bindx[Tmat_bindx[i]+1] = southeastback3d(ii,jj,kk,nx);
      Tmat_val[Tmat_bindx[i]+1] = 1.;
if ((ii==nx-2)||(jj==0)||(kk==0)) Tmat_val[Tmat_bindx[i]+1]=0.; /* rst dirichlet */
      if ((jj==0) || (kk==0)) Tmat_bindx[i+1]=Tmat_bindx[i];/* rst dirichlet */
    }
    else if (horv == VERTICAL) {
      Tmat_bindx[Tmat_bindx[i]] = northwestback3d(ii,jj,kk,nx);
      Tmat_val[Tmat_bindx[i]] = -1.;
if ((ii==0)||(jj==nx-2)||(kk==0)) Tmat_val[Tmat_bindx[i]]=0.; /* rst dirichlet */
      Tmat_bindx[Tmat_bindx[i]+1] = southwestback3d(ii,jj,kk,nx);
      Tmat_val[Tmat_bindx[i]+1] = 1.;

if ((ii==0)||(jj==0)||(kk==0)) Tmat_val[Tmat_bindx[i]+1]=0.; /* rst dirichlet */
      if ((ii==0) || (kk==0)) Tmat_bindx[i+1]=Tmat_bindx[i];/* rst dirichlet */
    }
    else {
      Tmat_bindx[Tmat_bindx[i]] = southwestback3d(ii,jj,kk,nx);
      Tmat_val[Tmat_bindx[i]] = 1.;
if ((ii==0)||(jj==0)||(kk==0)) Tmat_val[Tmat_bindx[i]]=0.; /* rst dirichlet */
      Tmat_bindx[Tmat_bindx[i]+1] = southwestfront3d(ii,jj,kk,nx);
      Tmat_val[Tmat_bindx[i]+1] = -1.;
if ((ii==0)||(jj==0)||(kk==nx-2)) Tmat_val[Tmat_bindx[i]+1]=0.; /* rst dirichlet */
      if ((ii==0) || (jj==0)) Tmat_bindx[i+1]=Tmat_bindx[i];/* rst dirichlet */
    }
nx = nx--; /* rst dirichlet */
#endif
    AZ_sort(&(Tmat_bindx[Tmat_bindx[i]]),Tmat_bindx[i+1]-Tmat_bindx[i],NULL,&(Tmat_val[Tmat_bindx[i]]));
  }
#else
  AZ_input_msr_matrix_nodiag("Tmat.az", global_edge_inds, &Tmat_val, 
			     &Tmat_bindx,  Nlocal_edges, proc_config);
  if (proc_config[AZ_node] == 0)
  {
     printf("Done reading T matrix\n"); fflush(stdout);
  }
#endif

  /* compress out any zeros which might occur due to empty rows  */

  lower = Tmat_bindx[0];
  Nnz = Tmat_bindx[Nlocal_edges];
  for (i = 0; i < Nlocal_edges; i++) {
    row_length = 0;
    for (j = lower; j < Tmat_bindx[i+1]; j++) {
      if (Tmat_val[j] != 0.0) row_length++;
    }
    lower = Tmat_bindx[i+1];
    Tmat_bindx[i+1] = Tmat_bindx[i] + row_length;
  }
  lower = Tmat_bindx[0];
  for (i = Tmat_bindx[0]; i < Nnz; i++ ) {
    if (Tmat_val[i] != 0.0) {
      Tmat_val[lower] = Tmat_val[i];
      Tmat_bindx[lower] = Tmat_bindx[i];
      lower++;
    }
  }

  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
							  ML_CSR_MSRdata));
  csr_data->columns = Tmat_bindx;
  csr_data->values  = Tmat_val;
  ML_MSR2CSR(csr_data, Nlocal_edges, &Ncols);
  ML_Create(&ml_edges, N_levels);
  ML_Create(&ml_nodes, N_levels);

  Nexterns = Kn_data_org[AZ_N_external];

  AZ_Tmat_transform2ml(Nexterns, global_node_externs, reordered_node_externs,
		       Tmat_bindx, Tmat_val, csr_data->rowptr, Nlocal_nodes,
		       global_node_inds, ml_edges->comm, Nlocal_edges, &Tmat);
  ML_free(csr_data);
  Tmat->data_destroy = ML_CSR_MSRdata_Destroy;

  start_time = AZ_second();
  options[AZ_scaling] = AZ_none;

  /********************************************************************/
  /* Convert the two discretization matrices (Ke_mat and Kn_mat) to   */
  /* ML matrices and stuff them into the ml grid hierarchies.         */
  /*------------------------------------------------------------------*/

  
  AZ_ML_Set_Amat(ml_edges, N_levels-1, Nlocal_edges, Nlocal_edges, Ke_mat, 
		 proc_config);
  AZ_ML_Set_Amat(ml_nodes, N_levels-1, Nlocal_nodes, Nlocal_nodes, Kn_mat, 
		 proc_config);
  ml_edges->Amat[N_levels-1].N_nonzeros = Ke_bindx[Nlocal_edges];

/* Check symmetry of Ke. */
  Amat = &(ml_edges->Amat[N_levels-1]);
  yyy = (double *) ML_allocate( Amat->invec_leng * sizeof(double) );
  vvv = (double *) ML_allocate( Amat->invec_leng * sizeof(double) );
  zzz = (double *) ML_allocate( Amat->invec_leng * sizeof(double) );

  AZ_random_vector(yyy, Ke_data_org, proc_config);
  AZ_random_vector(vvv, Ke_data_org, proc_config);

  ML_Operator_Apply(Amat, Amat->invec_leng, yyy,Amat->outvec_leng,zzz);
  dtemp = sqrt(abs(ML_gdot(Amat->outvec_leng, vvv, zzz, ml_edges->comm)));
  ML_Operator_Apply(Amat, Amat->invec_leng, vvv,Amat->outvec_leng,zzz);
  dtemp2 =  sqrt(abs(ML_gdot(Amat->outvec_leng, yyy, zzz, ml_edges->comm)));

  if (abs(dtemp-dtemp2) > 1e-15)
  {
     if (proc_config[AZ_node]== 0 && 0.5 < ML_Get_PrintLevel())
     {
        printf("\n\n*****************\n"
                       "WARNING: Edge matrix may not be symmetric.\n");
        printf("\n              ||vvv^{t} * Ke_mat * yyy|| = %20.15e\n",
                dtemp);
        printf("\n              ||yyy^{t} * Ke_mat * vvv|| = %20.15e\n",
                dtemp2);
        printf("               (vvv and yyy are random)\n");
        printf("*****************\n\n");
        fflush(stdout);
     }

     /* If Amat is not symmetric (i.e., Dirichlet rows have been zeroed with
        a one on diagonal, but corresponding columns have not been zeroed)
        zero out the Dirichlet rows in Tmat. This is necessary to have a
        potential matrix (T^{*}AT) that is equivalent to the potential matrix
        when Amat is symmetric.  See the function
        ML_Smoother_Gen_Hiptmair_Data to see how an equivalent potential
        matrix is created. */

     /* Read in a copy of Tmat that will contain boundary conditions. */
   
     if (proc_config[AZ_node] == 0)
     {
        printf("Reading T matrix\n"); fflush(stdout);
     }
     AZ_input_msr_matrix_nodiag("Tmat.az", global_edge_inds, &Tmatbc_val, 
   			     &Tmatbc_bindx,  Nlocal_edges, proc_config);
     if (proc_config[AZ_node] == 0)
     {
        printf("Done reading T matrix\n"); fflush(stdout);
     }
   
     /* compress out any zeros which might occur due to empty rows  */
   
     lower = Tmatbc_bindx[0];
     Nnz = Tmatbc_bindx[Nlocal_edges];
     for (i = 0; i < Nlocal_edges; i++) {
       row_length = 0;
       for (j = lower; j < Tmatbc_bindx[i+1]; j++) {
         if (Tmatbc_val[j] != 0.0) row_length++;
       }
       lower = Tmatbc_bindx[i+1];
       Tmatbc_bindx[i+1] = Tmatbc_bindx[i] + row_length;
     }
     lower = Tmatbc_bindx[0];
     for (i = Tmatbc_bindx[0]; i < Nnz; i++ ) {
       if (Tmatbc_val[i] != 0.0) {
         Tmatbc_val[lower] = Tmatbc_val[i];
         Tmatbc_bindx[lower] = Tmatbc_bindx[i];
         lower++;
       }
     }
   
     csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
   							  ML_CSR_MSRdata));
     csr_data->columns = Tmatbc_bindx;
     csr_data->values  = Tmatbc_val;
     ML_MSR2CSR(csr_data, Nlocal_edges, &Ncols);
     Nexterns = Kn_data_org[AZ_N_external];
   
     AZ_Tmat_transform2ml(Nexterns, global_node_externs, reordered_node_externs,
   		       Tmatbc_bindx, Tmatbc_val, csr_data->rowptr, Nlocal_nodes,
   		       global_node_inds, ml_edges->comm, Nlocal_edges, &Tmatbc);
     ML_free(csr_data);
     Tmatbc->data_destroy = ML_CSR_MSRdata_Destroy;

     if (proc_config[AZ_node]== 0) printf("Zeroing out rows of Tmat that "
                                       "correspond to Dirichlet points.\n");
     printf("\nProcessor %d owns %d rows of Ke.\n",proc_config[AZ_node],
             Amat->outvec_leng);
   
     data = (struct ML_CSR_MSRdata *) (Tmatbc->data);
     row_ptr = data->rowptr;
     bindx = data->columns;
     val_ptr = data->values;
     BCindices = (int *) ML_allocate( Amat->outvec_leng * sizeof(int) );
     BCcount = 0;

     /* Step through edge matrix Amat to find Dirichlet rows. */
     for (i=0; i < Amat->outvec_leng; i++)
     {
        ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                          &length, 0);
        Nnz = 0;
        for (j=0; j<length; j++)
           /*if ( abs(vals[j]) != 0.0 ) Nnz++;*/
           if ( abs(vals[j]) > 1e-12 ) Nnz++;
        if (Nnz <=1)   /* Dirichlet row */
        {
           /* Zero out corresponding row in Tmat. */
           /*printf("%d: %d is a Dirichlet row\n",proc_config[AZ_node],i+1);*/
           for (j = row_ptr[i]; j < row_ptr[i+1]; j++)
              val_ptr[j] = 0.0;
           BCindices[BCcount] = i;
           BCcount++;
        }
     }
     ML_Set_BoundaryTypes(ml_edges,N_levels-1,ML_BDRY_DIRICHLET,
                          BCcount,BCindices);
     ML_free(cols);
     ML_free(vals);
     ML_free(BCindices);

     ML_CommInfoOP_Clone(&(Tmatbc->getrow->pre_comm),
                         ml_nodes->Amat[N_levels-1].getrow->pre_comm);
     Tmat_transbc = ML_Operator_Create(ml_edges->comm);
     ML_Operator_Transpose_byrow(Tmatbc, Tmat_transbc);
  }
  else
  {
     if (proc_config[AZ_node]== 0 && 0.5 < ML_Get_PrintLevel())
     {
        printf("Edge matrix passed symmetry check.\n");
        fflush(stdout);
     }
  }
  ML_free(yyy); ML_free(zzz); ML_free(vvv);
  /* end of Ke symmetry check. */

  ML_CommInfoOP_Clone(&(Tmat->getrow->pre_comm),
                      ml_nodes->Amat[N_levels-1].getrow->pre_comm);


  /********************************************************************/
  /*                      Set up Tmat_trans                           */
  /*------------------------------------------------------------------*/

  Tmat_trans = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Transpose_byrow(Tmat, Tmat_trans);
  globbie = Tmat_trans;

#ifdef ML_partition

  evenodd[0] = 0; evenodd[1] = 0;
  for (i=0; i<nblocks; i++) Iwin[i] = 0;

  /* this code is meant to partition the matrices so that things can be */
  /* run in parallel later.                                             */
  /* It is meant to be run on only one processor.                       */

  ML_Operator_AmalgamateAndDropWeak(&(ml_nodes->Amat[N_levels-1]),
                                    num_PDE_eqns, 0.0);
  ML_Gen_Blocks_Metis(ml_nodes, N_levels-1, &nblocks, &block_list);

  fp2 = fopen("node_partition","w");
  for (i = 0; i < nblocks; i++)
  {
    count = 0;
    for (j = 0; j < ml_nodes->Amat[N_levels-1].outvec_leng; j++)
    {
      if (block_list[j] == i) count++;
    }
    fprintf(fp2,"   %d\n",count*num_PDE_eqns);
    for (j = 0; j < ml_nodes->Amat[N_levels-1].outvec_leng; j++)
    {
      if (block_list[j] == i)
      {
	     for (k = 0; k < num_PDE_eqns; k++)
            fprintf(fp2,"%d\n",j*num_PDE_eqns+k);
      }
    }
  }
  fclose(fp2);
  ML_Operator_UnAmalgamateAndDropWeak(&(ml_nodes->Amat[N_levels-1]),
                                      num_PDE_eqns,0.0);

  /* Need to partition Tmat given a partitioning of Ke_mat */
  /* Here's how it should work:                          */
  /*    1) if both column points are within processor,   */
  /*       assign row to this processor.                 */
  /*    2) if columns are from different procs,          */
  /*        a) take lower if product of columns is even  */
  /*        b) take upper if product of columns is odd.  */

  csr_data = (struct ML_CSR_MSRdata *) Tmat->data;
  proc_assignment = (int *) AZ_allocate( Tmat->outvec_leng*sizeof(int));

  j = 0;
  for (i = 0; i < Tmat->outvec_leng; i++)
  {
    /* Calculate the actual number of nonzeros (<=2) in the row. */
    itemp = (csr_data->rowptr)[i+1] - (csr_data->rowptr)[i];
    row_start = (csr_data->rowptr)[i];
    if (itemp == 2)
    {
      /* Entry could be zero if node is a Dirichlet point. */
      if ( (csr_data->values)[row_start+1] == 0.0) itemp--;
    }
    if (itemp > 0)
    {
      /* Entry could be zero if node is a Dirichlet point. */
       if ( (csr_data->values)[row_start] == 0.0)
       {
          itemp--;
          row_start++;
       }
    }
    if ( itemp > 2) 
      pr_error("Too many nonzeros per row in Tmat   %d\n", itemp);

    if (itemp == 1)
    {
      col1 = (csr_data->columns)[row_start];
      proc_assignment[i] = block_list[col1];
    }
    else if (itemp == 2)
    {
      col1 = (csr_data->columns)[row_start];
      col2 = (csr_data->columns)[row_start+1];
      p1   = block_list[col1];
      p2   = block_list[col2];
      /* 1+(int) (  2.0*rand()/(RAND_MAX+1.0)); */
      /* Randomly pick an integer between 1 and 10. */
      tiebreaker =  1+(int) (  10.0*rand()/(RAND_MAX+1.0));
      if ( tiebreaker > 5  )
      /*if ( (col1*col2)%2 == 0 ) */
      {
         evenodd[0]++;
         if (p1 < p2)
         {
            proc_assignment[i] = p1;
            Iwin[p1]++;
         }
         else
         {
            proc_assignment[i] = p2;
            Iwin[p2]++;
         }
      }
      else
      {
         evenodd[1]++;
         if (p2 < p1)
         {
            proc_assignment[i] = p1;
            Iwin[p1]++;
         }
         else
         {
            proc_assignment[i] = p2;
            Iwin[p2]++;
         }
      }
    }
    else
    {
       proc_assignment[i] = -1;
       j = 1;
    }
  }
  if (j==0) printf("\aWarning: no Dirichlet edges found in Tmat.\n");
  printf("even tie breakers = %d\n",evenodd[0]);
  printf("odd tie breakers = %d\n",evenodd[1]);
  for (i=0; i<nblocks; i++)
     printf("proc. %d won %d times\n",i,Iwin[i]);
  pcounts = (int *) ML_allocate( sizeof(int)*nblocks);
  proc_id = (int *) ML_allocate( sizeof(int)*nblocks);
  for (i = 0; i < nblocks; i++) pcounts[i] = 0;
  for (i = 0; i < nblocks; i++) proc_id[i] = i;

  count = 0;
  for (i = 0; i < Tmat->outvec_leng ; i++)
  {
    if (proc_assignment[i] != -1) pcounts[proc_assignment[i]]++;
    else count++;
  }

  printf("%d:Tmat->invec_leng = %d\n",Tmat->comm->ML_mypid,Tmat->invec_leng);
  printf("%d:Tmat->outvec_leng = %d\n",Tmat->comm->ML_mypid,Tmat->outvec_leng);
  fflush(stdout);
  
  ML_az_sort(pcounts, nblocks, proc_id, NULL);

  i = 0; j = 0;
  while ( (i < nblocks) && (pcounts[i] < pcounts[nblocks-1]) &&
          (j < Tmat->outvec_leng) )
  {
    if ( proc_assignment[j] == -1)
    {
      proc_assignment[j] = proc_id[i];
      pcounts[i]++;
      if (pcounts[i] == pcounts[nblocks-1]) i++;
    }
    j++;
  }
  for (i = j; i < Tmat->outvec_leng; i++) {
    if (proc_assignment[i] == -1) proc_assignment[i] = i%nblocks;
  }

  fp2 = fopen("edge_partition","w");

  for (i = 0; i < nblocks; i++) {
    count = 0;
    for (j = 0; j < Tmat->outvec_leng; j++) {
      if (proc_assignment[j] == i) count++;
    }
    fprintf(fp2,"   %d\n",count*num_PDE_eqns);
    for (j = 0; j < Tmat->outvec_leng; j++) {
      if (proc_assignment[j] == i) {
	for (k = 0; k < num_PDE_eqns; k++)  fprintf(fp2,"%d\n",j*num_PDE_eqns+k);
      }
    }
  }
  fclose(fp2);


  exit(1);
#endif /*ifdef ML_partition*/

  /********************************************************************/
  /* Set some ML parameters.                                          */
  /*------------------------------------------------------------------*/
	
  ML_Set_ResidualOutputFrequency(ml_edges, context->output);
  ML_Set_Tolerance(ml_edges, context->tol);
  ML_Aggregate_Create( &ag );
  if (ML_strcmp(context->agg_coarsen_scheme,"Mis") == 0) {
    ML_Aggregate_Set_CoarsenScheme_MIS(ag);
  }
  else if (ML_strcmp(context->agg_coarsen_scheme,"Uncoupled") == 0) {
    ML_Aggregate_Set_CoarsenScheme_Uncoupled(ag);
  }
  else if (ML_strcmp(context->agg_coarsen_scheme,"Coupled") == 0) {
    ML_Aggregate_Set_CoarsenScheme_Coupled(ag);
  }
  else {
    printf("ML: Unknown aggregation scheme %s\n",context->agg_coarsen_scheme);
  }
  ML_Aggregate_Set_DampingFactor(ag, context->agg_damping);
  ML_Aggregate_Set_MaxCoarseSize( ag, context->maxcoarsesize);
  ML_Aggregate_Set_Threshold(ag, context->agg_thresh);

  if (ML_strcmp(context->agg_spectral_norm,"Calc") == 0) {
    ML_Aggregate_Set_SpectralNormScheme_Calc(ag);
  }
  else if (ML_strcmp(context->agg_spectral_norm,"Anorm") == 0) {
    ML_Aggregate_Set_SpectralNormScheme_Anorm(ag);
  }
  else {
    printf("ML: Unknown spectral norm scheme %s\n",context->agg_spectral_norm);
  }

  /* read in the rigid body modes */

  Nrigid = 0;
  if (proc_config[AZ_node] == 0) {
    sprintf(filename,"rigid_body_mode%d",Nrigid+1);
    while( (fp = fopen(filename,"r")) != NULL) {
      fclose(fp);
      Nrigid++;
      sprintf(filename,"rigid_body_mode%d",Nrigid+1);
    }
  }
  Nrigid = AZ_gsum_int(Nrigid,proc_config);

  if (Nrigid != 0) {
    rigid = (double *) ML_allocate( sizeof(double)*Nrigid*(Nlocal_edges+1) );
    if (rigid == NULL) {
      printf("Error: Not enough space for rigid body modes\n");
    }
  }

  //blockrhs ?
  AZ_reorder_vec(rhs, Ke_data_org, reordered_glob_edges, NULL);

  for (i = 0; i < Nrigid; i++) {
    sprintf(filename,"rigid_body_mode%d",i+1);
    AZ_input_msr_matrix(filename, global_edge_inds, &mode, &garbage, Nlocal_edges, 
			proc_config);

    /*
     *  Rescale matrix/rigid body modes and checking 
     *
     AZ_sym_rescale_sl(mode, Ke_mat->Ke_data_org, options, proc_config, scaling);
     Ke_mat->matvec(mode, rigid, Ke_mat, proc_config);
     for (j = 0; j < Nlocal_edges; j++) printf("this is %d %e\n",j,rigid[j]);
     */

    for (j = 0; j < i; j++) {
      alpha = -AZ_gdot(Nlocal_edges, mode, &(rigid[j*Nlocal_edges]), proc_config)/
	AZ_gdot(Nlocal_edges, &(rigid[j*Nlocal_edges]), &(rigid[j*Nlocal_edges]), 
		proc_config);
      daxpy_(&Nlocal_edges, &alpha,  &(rigid[j*Nlocal_edges]),  &one, mode, &one);
    }
   
    /* rhs orthogonalization */

    alpha = -AZ_gdot(Nlocal_edges, mode, rhs, proc_config)/
      AZ_gdot(Nlocal_edges, mode, mode, proc_config);
    daxpy_(&Nlocal_edges, &alpha,  mode,  &one, rhs, &one);

    for (j = 0; j < Nlocal_edges; j++) rigid[i*Nlocal_edges+j] = mode[j];
    free(mode);
    if (garbage != NULL) free(garbage);
    garbage = NULL;
  }

  for (j = 0; j < Nrigid; j++) {
    alpha = -AZ_gdot(Nlocal_edges, rhs, &(rigid[j*Nlocal_edges]), proc_config)/
      AZ_gdot(Nlocal_edges, &(rigid[j*Nlocal_edges]), &(rigid[j*Nlocal_edges]), 
	      proc_config);
    daxpy_(&Nlocal_edges, &alpha,  &(rigid[j*Nlocal_edges]),  &one, rhs, &one);
  }

  if (Nrigid != 0) {
    ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, Nrigid, rigid, Nlocal_edges);
  }

#ifdef HierarchyCheck
  xxx = (double *) ML_allocate(2*(3+Nlocal_edges + Ke_mat->data_org[AZ_N_external])
                               *sizeof(double)); 


  for (iii = 0; iii < 2*Nlocal_edges; iii++) xxx[iii] = 0.0; 

  /* Set xxx */

  if (proc_config[AZ_node]== 0)
     printf("putting in an edge based initial guess\n");
  fp = fopen("initguessfile","r");
  if (fp != NULL)
  {
    fclose(fp);
    free(xxx);
    if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
    AZ_input_msr_matrix("initguessfile", global_edge_inds, &xxx, &garbage,
			Nlocal_edges, proc_config);
    printf("this is not done yet\n");
    exit(1);
    options[AZ_conv] = AZ_expected_values;
    printf("done reading initial guess\n");
  }
#endif /* ifdef HierarchyCheck */
  if (Tmat_transbc != NULL)
     coarsest_level = ML_Gen_MGHierarchy_UsingReitzinger(ml_edges, ml_nodes,
						         N_levels-1, ML_DECREASING, ag, Tmatbc,
                                 Tmat_transbc, &Tmat_array, &Tmat_trans_array,
							 ML_NO, 1.5);
  else
     coarsest_level = ML_Gen_MGHierarchy_UsingReitzinger(ml_edges, ml_nodes,
						         N_levels-1, ML_DECREASING, ag, Tmat,
                                 Tmat_trans, &Tmat_array, &Tmat_trans_array,
							 ML_NO, 1.5);

  ml_edges->ML_finest_level = N_levels-1;

  printf("before complex\n");
  ML_Gen_Hierarchy_ComplexMaxwell(ml_edges, /* Tmat_array, Tmat_trans_array, */ &ml_block, Mmat_ml);
  printf("after Complex\n");
   // rst: comment this out to use the original multigrid stuff
      ml_edges = ml_block;

#ifdef ReuseOps
  {printf("Starting reuse\n"); fflush(stdout);}
  ML_Operator_Clean(&(ml_edges->Amat[N_levels-1]));
  ML_Operator_Init(&(ml_edges->Amat[N_levels-1]),ml_edges->comm);
  AZ_ML_Set_Amat(ml_edges, N_levels-1, Nlocal_edges, Nlocal_edges, Ke_mat, 
		 proc_config);

  ML_Gen_MGHierarchy_ReuseExistingOperators(ml_edges);
  {printf("Ending reuse\n"); fflush(stdout);}
#endif
  coarsest_level = N_levels - coarsest_level;

 /* The MLS smoother needs the largest eigenvalue of the matrix. */
  /* Normally, this happens automatically within ML by checking   */
  /* if the eigenvalue is not already defined and then calling a  */
  /* CG routine to compute it. However, CG won't work for the     */
  /* equivalent real system. Instead, we will manually compute the*/
  /* largest eigenvalue of the pieces to estimate the largest     */
  /* eigenvalue of the equivalent real system.                    */
  printf("here we are %d %d\n",
	 ml_block->ML_finest_level, ml_block->ML_coarsest_level);
  for (j = N_levels-1; j >= coarsest_level; j--) {
    printf("THE LEVEL IS %d\n",j);


     block_data = (struct ML_Operator_blockmat_data *) ml_block->Amat[j].data;

    /* compute eigenvalue of 'stiffness' matrix (1,1) block */
     kdata = ML_Krylov_Create( ml_edges->comm );
     ML_Krylov_Set_PrintFreq( kdata, 0 );
     ML_Krylov_Set_ComputeEigenvalues( kdata );
     ML_Krylov_Set_Amatrix(kdata, block_data->Ke_mat );

     ML_Krylov_Solve(kdata, block_data->Ke_mat->outvec_leng, NULL, NULL);
     ml_edges->Amat[j].lambda_max = 
                      ML_Krylov_Get_MaxEigenvalue(kdata);
     ML_Krylov_Destroy( &kdata );

     /* Get the max eigenvalue of M (a diagonal matrix) */
     ML_Operator_Get_Diag(block_data->M_mat,
                          block_data->M_mat->outvec_leng, &M_diag);

     /* diagonal of (curl,curl) matrix */
     ML_Operator_Get_Diag(block_data->Ke_mat,
                          block_data->Ke_mat->outvec_leng, &diag_of_S);

     lambda = 0.;
     tdiag = 0.0;
     for (i = 0; i < block_data->Ke_mat->outvec_leng; i++) {
       if ( fabs(M_diag[i]/diag_of_S[i]) > lambda) {
         lambda = fabs(M_diag[i]/diag_of_S[i]);
         tdiag = diag_of_S[i];
         ttdiag = M_diag[i];
       }
     }
     lambda = ML_gmax_double(lambda, ml_edges->comm);
     block_data->M_mat->lambda_max = lambda;

     /* Put in eigenvalue estimate for equivalent real form */
     ml_block->Amat[j].lambda_max     = ml_edges->Amat[j].lambda_max;
     ml_block->Amat[j].lambda_max_img = block_data->M_mat->lambda_max;
     if (ml_edges->comm->ML_mypid == 0) {
       printf("\n\nlevel %d:  (lambda_max_real,lambda_max_imag) = %e  %e\n\n",
              j,ml_edges->Amat[j].lambda_max, block_data->M_mat->lambda_max);
       printf("(level %d) numer = %e, denom = %e\n",j,ttdiag,tdiag);
       fflush(stdout);
     }
  } /* for j=ml_edges->ML_finest_level... */

  /* Here is the stuff to set the subsmoothers within the Hiptmair */
  /* smoother.                                                     */

  if (ML_strcmp(context->smoother,"Hiptmair") == 0)
    {
      nodal_its      = 1;
      edge_its       = 1;
      if (ML_strcmp(context->subsmoother,"GaussSeidel") == 0
	  || ML_strcmp(context->subsmoother,"default") == 0)
	{
	  nodal_smoother = (void *) ML_Gen_Smoother_SymGaussSeidel;
	  edge_smoother  = (void *) ML_Gen_Smoother_SymGaussSeidel;
/* edge_smoother  = (void *) ML_Gen_Smoother_VBlockSymGaussSeidel; */
	  nodal_omega    = ML_DDEFAULT;
	  edge_omega     = ML_DDEFAULT;
	  nodal_args = ML_Smoother_Arglist_Create(2);
	  ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
	  ML_Smoother_Arglist_Set(nodal_args, 1, &nodal_omega);
      /*this flag doesn't matter and is just to make the nodal and edge
        arg arrays the same length. */
//	  ML_Smoother_Arglist_Set(nodal_args, 2, &reduced_smoother_flag);

	  if (edge_smoother == (void *) ML_Gen_Smoother_VBlockSymGaussSeidel)
	    edge_args = ML_Smoother_Arglist_Create(5);
	  else
	    edge_args = ML_Smoother_Arglist_Create(2);
	  ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
	  ML_Smoother_Arglist_Set(edge_args, 1, &edge_omega);
      /* if flag is nonzero, in Hiptmair do edge/nodal combination on pre-
         smooth and nodal/edge combination on post-smooth.   This maintains
         symmetry of the preconditioner. */
//	  ML_Smoother_Arglist_Set(edge_args, 2, &reduced_smoother_flag);
	}
      else if (ML_strcmp(context->subsmoother,"MLS") == 0)
	{
	  printf("mls\n");
	  //	  edge_smoother  = (void *) ML_Gen_Smoother_MLS;
	  edge_smoother  = (void *) ML_Gen_Smoother_MLS;
	  nodal_smoother = (void *) ML_Gen_Smoother_MLS;
	  nodal_omega    = 1.0;
	  edge_omega     = 1.0;
	  nodal_args = ML_Smoother_Arglist_Create(2);
	  ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
	  //	  ML_Smoother_Arglist_Set(nodal_args, 2, &nodal_omega);
	  edge_args = ML_Smoother_Arglist_Create(2);
	  ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
	  //	  ML_Smoother_Arglist_Set(edge_args, 2, &edge_omega);
	  printf("bottom of mls\n");
	}
      else if (ML_strcmp(context->subsmoother,"ERF_1Step") == 0)
	{
	  printf("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
	  edge_smoother  = (void *) ML_Gen_Smoother_ERF_1StepKrylov;
	  nodal_smoother = (void *) ML_Gen_Smoother_MLS;
	  nodal_omega    = 1.0;
	  edge_omega     = 1.0;
	  nodal_args = ML_Smoother_Arglist_Create(2);
	  ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
	  edge_args = ML_Smoother_Arglist_Create(0);
	}
    }

  /****************************************************
  * Set up smoothers for all levels but the coarsest. *
  ****************************************************/
  for (level = N_levels-1; level > coarsest_level; level--)
  {
      num_PDE_eqns = ml_edges->Amat[level].num_PDEs;
		
      /*  Sparse approximate inverse smoother that acutally does both */
      /*  pre and post smoothing.                                     */

      if (ML_strcmp(context->smoother,"Parasails") == 0)
	  {
	  ML_Gen_Smoother_ParaSails(ml_edges , level, ML_PRESMOOTHER, nsmooth, 
				    parasails_sym, parasails_thresh, 
				    parasails_nlevels, parasails_filter,
				    parasails_loadbal, parasails_factorized);
	  }

      else if (ML_strcmp(context->smoother,"Hiptmair") == 0)
	  {
	    printf("in hip\n");
          /* Setting omega to any other value will override the automatic
             calculation in ML_Smoother_Gen_Hiptmair_Data. */
          omega = (double) ML_DEFAULT;
	  if (edge_smoother == (void *) ML_Gen_Smoother_VBlockSymGaussSeidel) {
	    nblocks = ml_edges->Amat[level].outvec_leng/2;
	    blocks = (int *) ML_allocate(sizeof(int)*
					ml_edges->Amat[level].outvec_leng);
	    for (i =0; i < nblocks; i++) blocks[i] = i;
	    for (i =0; i < nblocks; i++) blocks[nblocks+i] = i;
	    ML_Smoother_Arglist_Set(edge_args, 3, &nblocks);
	    ML_Smoother_Arglist_Set(edge_args, 4, blocks);
	  }
	  if (edge_smoother == (void *) ML_Gen_Smoother_ERF_1StepKrylov) {
	  printf("222IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");

	    nodal_eig_ratio[level] = eig_ratio_tol;
	    temp1[2] = Tmat_array[level]->invec_leng;
	    if (level != 0) temp1[3] = Tmat_array[level-1]->invec_leng;
	    else            temp1[3] = 0;

	    ML_gsum_vec_int(&temp1, &temp2, 4, ml_edges->comm);
	    if (temp1[3] != 0) {
	      nodal_eig_ratio[level] = 2.*((double) temp1[2])/ ((double) temp1[3]);
	    }
	    if (nodal_eig_ratio[level] < eig_ratio_tol) nodal_eig_ratio[level] = eig_ratio_tol;
	    ML_Smoother_Arglist_Set(nodal_args, 1, &(nodal_eig_ratio[level]));
        mls_poly_degree = 1;
        printf("\n\nChebychev polynomial degree hardwired to %d\a\a\n\n",
                mls_poly_degree);
	  }

	  if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {

	    printf("more mls %u %u\n",temp1, temp2);
	    edge_eig_ratio[level] = eig_ratio_tol;
	    nodal_eig_ratio[level] = eig_ratio_tol;
	    printf("before temp1\n");
	    temp1[0] = Tmat_array[level]->outvec_leng;
	    printf("after temp1\n");
	    temp1[2] = Tmat_array[level]->invec_leng;
	    if (level != 0) {
	      temp1[1] = Tmat_array[level-1]->outvec_leng;
	      temp1[3] = Tmat_array[level-1]->invec_leng;
	    }
	    else { 
	      temp1[1] = 0;
	      temp1[3] = 0;
	    }
	    printf("before gsum\n");
	    ML_gsum_vec_int(&temp1, &temp2, 4, ml_edges->comm);
	    printf("in hip\n");
	    if (temp1[1] != 0) {
	      edge_eig_ratio[level] = 2.*((double) temp1[0])/ ((double) temp1[1]);
	      nodal_eig_ratio[level] = 2.*((double) temp1[2])/ ((double) temp1[3]);
	    }
	    if ( edge_eig_ratio[level] < eig_ratio_tol) edge_eig_ratio[level] = eig_ratio_tol;
	    if (nodal_eig_ratio[level] < eig_ratio_tol) nodal_eig_ratio[level] = eig_ratio_tol;
	    ML_Smoother_Arglist_Set(edge_args, 1, &(edge_eig_ratio[level]));
	    ML_Smoother_Arglist_Set(nodal_args, 1, &(nodal_eig_ratio[level]));
        mls_poly_degree = 1;
        printf("\n\nChebychev polynomial degree hardwired to %d\a\a\n\n",
                mls_poly_degree);
	//	    ML_Smoother_Arglist_Set(nodal_args, 3, &mls_poly_degree);
	//	    ML_Smoother_Arglist_Set(edge_args, 3, &mls_poly_degree);
	    printf("bottom  hip\n");
	  }

      if (level == N_levels-1)
         ML_Gen_Smoother_BlockHiptmair(ml_edges, level, ML_BOTH, nsmooth,
                      Tmat_array, Tmat_trans_array, Tmatbc, edge_smoother,
                      edge_args, nodal_smoother,nodal_args, HALF_HIPTMAIR);
      else
         ML_Gen_Smoother_BlockHiptmair(ml_edges, level, ML_BOTH, nsmooth,
                      Tmat_array, Tmat_trans_array, NULL, edge_smoother,
		      edge_args, nodal_smoother,nodal_args,HALF_HIPTMAIR);
	  }
      /* This is the symmetric Gauss-Seidel smoothing that we usually use. */
      /* In parallel, it is not a true Gauss-Seidel in that each processor */
      /* does a Gauss-Seidel on its local submatrix independent of the     */
      /* other processors.                                                 */

      else if (ML_strcmp(context->smoother,"GaussSeidel") == 0)
	  {
	  ML_Gen_Smoother_GaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.);
	  }
      else if (ML_strcmp(context->smoother,"SymGaussSeidel") == 0)
	  {
	  ML_Gen_Smoother_SymGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.);
	  }
      else if (ML_strcmp(context->smoother,"BlockGaussSeidel") == 0)
  	  {
	  ML_Gen_Smoother_BlockGaussSeidel(ml_edges , level, ML_BOTH, nsmooth,1.,
					   num_PDE_eqns);
	  }
      else if (ML_strcmp(context->smoother,"Aggregate") == 0)
	  {
	  ML_Gen_Blocks_Aggregates(ag, level, &nblocks, &blocks);
	  ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , level, ML_BOTH,
					       nsmooth,1., nblocks, blocks);
	  }

      /* This is a true Gauss Seidel in parallel. This seems to work for  */
      /* elasticity problems.  However, I don't believe that this is very */
      /* efficient in parallel.                                           */       
      /*
         blocks = (int *) ML_allocate(sizeof(int)*ml_edges->Amat[level].invec_leng);
	nblocks = ml_edges->Amat[level].invec_leng;
	for (i =0; i < nblocks; i++) blocks[i] = i;
	ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml_edges , level,
	ML_PRESMOOTHER, nsmooth, 1., nblocks, blocks);
	ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml_edges, level,
	ML_POSTSMOOTHER, nsmooth, 1., nblocks, blocks);
	ML_free(blocks);
	*/

      /* Jacobi Smoothing                                                 */

      else if (ML_strcmp(context->smoother,"Jacobi") == 0)
	  {
	  /*
	  ML_Gen_Smoother_Jacobi(ml_edges , level, ML_PRESMOOTHER, nsmooth,.67);
	  ML_Gen_Smoother_Jacobi(ml_edges , level, ML_POSTSMOOTHER, nsmooth,.67);
	  */
      ML_Gen_Smoother_Jacobi(ml_edges , level, ML_BOTH, nsmooth, 0.67);
	  }

      /*  This does a block Gauss-Seidel (not true GS in parallel)        */
      /*  where each processor has 'nblocks' blocks.                      */
      /* */

      else if (ML_strcmp(context->smoother,"Metis") == 0)
	  {
	  nblocks = 250;
	  nblocks = ml_edges->Amat[level].invec_leng/25;
	  nblocks++;
	  ML_Gen_Blocks_Metis(ml_edges, level, &nblocks, &blocks);
	  ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , level, ML_BOTH,
					       nsmooth,1., nblocks, blocks);
	  ML_free(blocks);
	  }
      else
	  {
	  printf("unknown smoother %s\n",context->smoother);
	  exit(1);
	  }
  }
  printf("after fine smoothers\n");

  /*******************************************
  * Set up smoothers for the coarsest level. *
  *******************************************/
  nsmooth   = context->coarse_its;
  /*  Sparse approximate inverse smoother that actually does both */
  /*  pre and post smoothing.                                     */

  if (ML_strcmp(context->coarse_solve,"Parasails") == 0)
  {
      ML_Gen_Smoother_ParaSails(ml_edges , coarsest_level, ML_PRESMOOTHER,
				nsmooth, parasails_sym, parasails_thresh, 
                                parasails_nlevels, parasails_filter,
                                parasails_loadbal, parasails_factorized);
  }

  else if (ML_strcmp(context->coarse_solve,"Hiptmair") == 0)
  {
    /* Setting omega to any other value will override the automatic
       calculation in ML_Smoother_Gen_Hiptmair_Data. */
    edge_eig_ratio[coarsest_level] = eig_ratio_tol;
    nodal_eig_ratio[coarsest_level] = eig_ratio_tol;
    omega = (double) ML_DEFAULT;
    if (edge_smoother == (void *) ML_Gen_Smoother_ERF_1StepKrylov) {
	  printf("333IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");

       ML_Smoother_Arglist_Set(nodal_args, 1, &(nodal_eig_ratio[coarsest_level]));
       mls_poly_degree = 1;
       printf("\n\nChebychev polynomial degree hardwired to %d\a\a\n\n",
              mls_poly_degree);
    }
    if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
       ML_Smoother_Arglist_Set(edge_args, 1, &(edge_eig_ratio[coarsest_level]));
       ML_Smoother_Arglist_Set(nodal_args, 1, &(nodal_eig_ratio[coarsest_level]));
       mls_poly_degree = 1;
       printf("\n\nChebychev polynomial degree hardwired to %d\a\a\n\n",
              mls_poly_degree);
    }
    if (coarsest_level == N_levels-1)
       ML_Gen_Smoother_BlockHiptmair(ml_edges, level, ML_PRESMOOTHER, nsmooth,
				Tmat_array, Tmat_trans_array, Tmatbc, 
edge_smoother,edge_args, nodal_smoother,nodal_args,HALF_HIPTMAIR);
    else
       ML_Gen_Smoother_BlockHiptmair(ml_edges, level, ML_PRESMOOTHER, nsmooth,
				Tmat_array, Tmat_trans_array, NULL, 
				edge_smoother,edge_args, nodal_smoother,nodal_args,HALF_HIPTMAIR);
  }
  else if (ML_strcmp(context->coarse_solve,"GaussSeidel") == 0) {
    ML_Gen_Smoother_GaussSeidel(ml_edges , coarsest_level, ML_BOTH, nsmooth,1.);
  }
  else if (ML_strcmp(context->coarse_solve,"SymGaussSeidel") == 0) {
    printf("this is GS with %d\n",nsmooth);
    ML_Gen_Smoother_SymGaussSeidel(ml_edges , coarsest_level, ML_BOTH, nsmooth,1.);
  }
  else if (ML_strcmp(context->coarse_solve,"BlockGaussSeidel") == 0) {
    ML_Gen_Smoother_BlockGaussSeidel(ml_edges, coarsest_level, ML_BOTH, nsmooth,1.,
				     num_PDE_eqns);
  }
  else if (ML_strcmp(context->coarse_solve,"Aggregate") == 0) {
    ML_Gen_Blocks_Aggregates(ag, coarsest_level, &nblocks, &blocks);
    ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , coarsest_level, ML_BOTH, 
					 nsmooth,1., nblocks, blocks);
  }
  else if (ML_strcmp(context->coarse_solve,"Jacobi") == 0) {
    printf("number of steps is %d\n",nsmooth);
    ML_Gen_Smoother_Jacobi(ml_edges , coarsest_level, ML_BOTH, nsmooth,.67);
  }
  else if (ML_strcmp(context->coarse_solve,"Metis") == 0) {
    nblocks = 250;
    nblocks = ml_edges->Amat[coarsest_level].invec_leng/25;
    nblocks++;
    ML_Gen_Blocks_Metis(ml_edges, coarsest_level, &nblocks, &blocks);
    ML_Gen_Smoother_VBlockSymGaussSeidel(ml_edges , coarsest_level, ML_BOTH, 
					 nsmooth,1., nblocks, blocks);
    if (blocks != NULL) ML_free(blocks);
  }
  else if (ML_strcmp(context->coarse_solve,"SuperLU") == 0) {
    ML_Gen_CoarseSolverSuperLU( ml_edges, coarsest_level);
  }
  else {
    printf("unknown coarse grid solver %s\n",context->coarse_solve);
    exit(1);
  }


#ifdef ReuseOps

  ML_Smoother_Reinit(ml_edges);

  /***************************************************************
  *  Regenerate Hiptmair smoother data on all levels as a test.  *
  ***************************************************************/

  if (ML_strcmp(context->smoother,"Hiptmair") == 0)
  {
     if (proc_config[AZ_node] == 0)
     {
        printf("Regenerating smoother data\n");
        fflush(stdout);
     }
     nsmooth   = context->nsmooth;
     for (level = N_levels-1; level > coarsest_level; level--)
     {
         num_PDE_eqns = ml_edges->Amat[level].num_PDEs;
         /* Setting omega to any other value will override the automatic
            calculation in ML_Smoother_Gen_Hiptmair_Data. */
        omega = (double) ML_DEFAULT;
         if (level == N_levels-1)
	        ML_Gen_Smoother_Hiptmair(ml_edges, level, ML_BOTH, nsmooth,
					 Tmat_array, Tmat_trans_array, 
					 Tmat, Tmat_trans, Tmatbc,
				edge_smoother,edge_args, nodal_smoother,nodal_args);
         else
	        ML_Gen_Smoother_Hiptmair(ml_edges, level, ML_BOTH, nsmooth,
					 Tmat_array, Tmat_trans_array,
					 Tmat, Tmat_trans,NULL, 
				edge_smoother,edge_args, nodal_smoother,nodal_args);
     }
  }
  nsmooth   = context->coarse_its;
  coarsest_level = N_levels - coarsest_level;
  if (ML_strcmp(context->coarse_solve,"Hiptmair") == 0)
  {
     if (proc_config[AZ_node] == 0)
     {
        printf("Regenerating smoother data on coarsest level\n");
        fflush(stdout);
     }
     /* Setting omega to any other value will override the automatic
        calculation in ML_Smoother_Gen_Hiptmair_Data. */
     omega = (double) ML_DEFAULT;
     if (coarsest_level == N_levels-1)
	     //jes
        ML_Gen_Smoother_BlockHiptmair(ml_edges , coarsest_level, ML_BOTH,
				 nsmooth,Tmat_array, Tmat_trans_array, 
				 Tmat, Tmat_trans, Tmatbc, 
				edge_smoother,edge_args, nodal_smoother,nodal_args,HALF_HIPTMAIR);
     else
	     //jes
        ML_Gen_Smoother_BlockHiptmair(ml_edges, level, ML_BOTH, nsmooth, 
				 Tmat_array, Tmat_trans_array, Tmat, 
				 Tmat_trans, NULL, 
				edge_smoother,edge_args, nodal_smoother,nodal_args,HALF_HIPTMAIR);
  }
  else if (ML_strcmp(context->coarse_solve,"SuperLU") == 0)
  {
    ML_Gen_CoarseSolverSuperLU( ml_edges, coarsest_level);
  }
  else {
    printf("unknown coarse grid solver %s\n",context->coarse_solve);
    exit(1);
  }
  if (proc_config[AZ_node] == 0)
  {
     printf("Done regenerating data\n");
     fflush(stdout);
  }
#endif /* ifdef ReuseOps */
		
  mg_cycle_type = ML_MGV;
  ML_Gen_Solver(ml_edges, mg_cycle_type, N_levels-1, coarsest_level); 
  AZ_defaults(options, params);
	
  if (ML_strcmp(context->krylov,"Cg") == 0) {
    options[AZ_solver]   = AZ_cg;
  }
  else if (ML_strcmp(context->krylov,"Bicgstab") == 0) {
    options[AZ_solver]   = AZ_bicgstab;
  }
  else if (ML_strcmp(context->krylov,"Tfqmr") == 0) {
    options[AZ_solver]   = AZ_tfqmr;
  }
  else if (ML_strcmp(context->krylov,"Gmresr") == 0) {
    options[AZ_solver]   = AZ_GMRESR;
  }
  else if (ML_strcmp(context->krylov,"Gmres") == 0) {
    options[AZ_solver]   = AZ_gmres;
  }
  else
  {
    printf("unknown krylov method %s\n",context->krylov);
  }
  options[AZ_scaling]  = AZ_none;
  options[AZ_precond]  = AZ_user_precond;
  options[AZ_conv]     = AZ_noscaled;
#ifdef AZ_CONVR0
  options[AZ_conv]     = AZ_r0;
#endif
#ifdef AZ_CONVRHS
  options[AZ_conv]     = AZ_rhs;
#endif
  options[AZ_output]   = 1;
  options[AZ_max_iter] = 500; 
  options[AZ_poly_ord] = 5;
  options[AZ_kspace]   = 30;
  options[AZ_orthog]   = AZ_classic;
  params[AZ_tol]       = context->tol;
  options[AZ_output]   = context->output;
	
  AZ_set_ML_preconditioner(&Pmat, Ke_mat, ml_edges, options); 
  setup_time = AZ_second() - start_time;
	
  xxx = (double *) ML_allocate(2*(3+Nlocal_edges + Ke_mat->data_org[AZ_N_external])
                               *sizeof(double)); 

  for (iii = 0; iii < 2*Nlocal_edges; iii++) xxx[iii] = 0.0; 


  /* Set xxx */

  fp = fopen("initguessfile","r");
  if (fp != NULL) {
    fclose(fp);
    free(xxx);
    if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
    AZ_input_msr_matrix("initguessfile", global_edge_inds, &xxx, &garbage,
			Nlocal_edges, proc_config);
    printf("not done\n");
    exit(1);
    printf("done reading initial guess\n");
  }
  else
  {
  //  if (proc_config[AZ_node]== 0) printf("taking random initial guess \n");
  //  AZ_random_vector(xxx, Ke_data_org, proc_config);
  //  AZ_random_vector(&(xxx[Nlocal_edges]), Ke_data_org, proc_config);
  }
  //  AZ_reorder_vec(xxx, Ke_data_org, reordered_glob_edges, NULL);
  /*
  fp = fopen("randomvec","w");
  for (i=0; i<Nlocal_edges; i++)
     fprintf(fp,"%20.15f\n",xxx[i]); 
  fclose(fp);
  */

  dtemp = sqrt(ML_gdot(2*Nlocal_edges, xxx, xxx, ml_edges->comm));
  if (proc_config[AZ_node]== 0 && 5 < ML_Get_PrintLevel() )
  {
    printf("length of initial guess = %d\n",Nlocal_edges);
    printf("||xxx|| = %e\n",dtemp);
  }
  dtemp = sqrt(ML_gdot(Nlocal_edges, rhs, rhs, ml_edges->comm));
  if (proc_config[AZ_node]== 0 && 5 < ML_Get_PrintLevel() )
  {
  printf("||rhs|| = %e\n",dtemp);
  fflush(stdout);
  }

  /*
    printf("putting in an node based xxxx\n");
    fp = fopen("initguessfile","r");
    if (fp != NULL) {
    fclose(fp);
    if (proc_config[AZ_node]== 0) printf("reading initial guess from file\n");
    AZ_input_msr_matrix("initguessfile", global_node_inds, &xxx, &garbage, Nlocal_nodes, 
    proc_config);
    options[AZ_conv] = AZ_expected_values;
    }
    else if (proc_config[AZ_node]== 0) printf("taking 0 initial guess \n");

    dtemp = sqrt(ML_gdot(Nlocal_nodes, xxx, xxx, ml_edges->comm));
    printf("length of initial guess = %d\n",Nlocal_nodes);
    printf("||xxx|| = %e\n",dtemp);

    AZ_reorder_vec(xxx, Kn_data_org, reordered_glob_nodes, NULL);
    */

  fp = fopen("AZ_no_multilevel.dat","r");
  scaling = AZ_scaling_create();
  start_time = AZ_second();
  if (fp != NULL)
  {
    fclose(fp);
    options[AZ_precond] = AZ_none;
    options[AZ_scaling] = AZ_sym_diag;
    options[AZ_ignore_scaling] = AZ_TRUE;

    printf("hereeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");

    options[AZ_keep_info] = 1;

    AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, NULL, scaling); 

    /*
      options[AZ_pre_calc] = AZ_reuse;
      options[AZ_conv] = AZ_expected_values;
      if (proc_config[AZ_node] == 0) 
      printf("\n-------- Second solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, NULL, scaling); 
      if (proc_config[AZ_node] == 0) 
      printf("\n-------- Third solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, NULL, scaling); 
      */
  }
  else
  {
    options[AZ_keep_info] = 1;
    options[AZ_conv]     = AZ_noscaled;
#ifdef AZ_CONVR0
  options[AZ_conv]     = AZ_r0;
#endif
#ifdef AZ_CONVRHS
  options[AZ_conv]     = AZ_rhs;
#endif
    options[AZ_output] = 1;

    /*
      options[AZ_precond] = AZ_none;
     
      Ke_mat->matvec(xxx, rhs, Ke_mat, proc_config);
      for (i = 0; i < Nlocal_edges; i++)
         printf("%7d     %7d %20.15e %20.15e\n",i+1,i+1,xxx[i],rhs[i]);
      printf("huhhhh %e\n",Ke_mat->val[0]);
      */

  options[AZ_conv]     = AZ_noscaled;
#ifdef AZ_CONVR0
  options[AZ_conv]     = AZ_r0;
#endif
#ifdef AZ_CONVRHS
  options[AZ_conv]     = AZ_rhs;
#endif
    /*options[AZ_conv] = AZ_expected_values;*/



  //sfksdkfksdfkskdfksdfksdk
 
  /**** check various operators and vectors ****/
  if ( -1 < ML_Get_PrintLevel() )
  {

    printf("\nChecking various operators...\n\n");
    if (N_levels > 1)
    {
       Amat = &(ml_edges->Rmat[N_levels-1]);
       yyy = (double *) malloc( Amat->outvec_leng * sizeof(double) );
       ML_Operator_Apply(Amat, Amat->invec_leng, rhs,Amat->outvec_leng,yyy);
       dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
       printf("||R_e * rhs|| = %20.15e\n",dtemp);

       ML_Operator_Apply(Amat, Amat->invec_leng, xxx,Amat->outvec_leng,yyy);
       dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
       printf("||R_e * xxx|| = %20.15e\n",dtemp);
       ML_free(yyy);
    }

    dtemp = sqrt(ML_gdot(Amat->outvec_leng, rhs, rhs, ml_edges->comm));
    printf("||rhs|| = %20.15e\n",dtemp);
    
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, xxx, xxx, ml_edges->comm));
    printf("||xxx|| = %20.15e\n",dtemp);
    
    yyy = (double *) malloc( 2 * Mmat_ml->outvec_leng * sizeof(double) );
    ML_Operator_Apply(Mmat_ml, Mmat_ml->invec_leng, rhs, Mmat_ml->outvec_leng, yyy);
    dtemp = sqrt(ML_gdot(Mmat_ml->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Mmat_ml * rhs|| = %20.15e\n",dtemp);
    //ML_Operator_Print(Mmat_ml,"Mmat_ml");
    
    Amat = &(ml_edges->Amat[N_levels-1]);
    ML_Operator_Apply(Amat, Amat->invec_leng, rhs,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Ke_mat * rhs|| = %20.15e\n",dtemp);

    //ML_Operator_Print(Amat,"Amat1");
    
    ML_Operator_Apply(Amat, Amat->invec_leng, xxx,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Ke_mat * xxx|| = %20.15e\n",dtemp);
    ML_free(yyy);

    Amat = Tmat_trans;
    yyy = (double *) malloc( Amat->outvec_leng * sizeof(double) );

    ML_Operator_Apply(Amat, Amat->invec_leng, xxx,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Tmat_trans * xxx|| = %20.15e\n",dtemp);

    ML_Operator_Apply(Amat, Amat->invec_leng, rhs,Amat->outvec_leng,yyy);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, yyy, yyy, ml_edges->comm));
    printf("||Tmat_trans * rhs|| = %20.15e\n",dtemp);

    if (N_levels > 1)
    {
       Amat = &(ml_nodes->Rmat[N_levels-1]);
       vvv = (double *) malloc( Amat->outvec_leng * sizeof(double) );
       ML_Operator_Apply(Amat, Amat->invec_leng, yyy,Amat->outvec_leng,vvv);
       dtemp = sqrt(ML_gdot(Amat->outvec_leng, vvv, vvv, ml_edges->comm));
       printf("||R_n * Tmat_trans * yyy|| = %20.15e\n",dtemp);
       ML_free(vvv);
    }

    ML_free(yyy);
    printf("\nEnd of check.\n\n");
    fflush(stdout);
/**** end of check ****/

/**** check prolongation operator ****
    Amat = &(ml_edges->Pmat[N_levels-1]);
    vvv = (double *) malloc( Amat->outvec_leng * sizeof(double) );
    ML_Operator_Apply(Amat, Amat->invec_leng, yyy,Amat->outvec_leng,vvv);
    dtemp = sqrt(ML_gdot(Amat->outvec_leng, vvv, vvv, ml_edges->comm));
    printf("||P * R * xxx|| = %10.7e\n",dtemp);
    free(yyy);
    free(vvv);
**** end of check ****/

  } /*end of operator check*/

    if (proc_config[AZ_node] == 0 && 0.5 < ML_Get_PrintLevel())
    {
       if (mg_cycle_type == ML_MGV)
          printf("Cycle type = MGV\n");
       else if (mg_cycle_type == ML_MGW)
          printf("Cycle type = MGW\n");
       else
          printf("Cycle type = other\n");
    }
    fflush(stdout);
    options[AZ_scaling] = AZ_none;
    /* rst: create the aztec block 2x2 matrix. Shove in Ke and M   */
    /* into the data structure and associate the new matrix-vector */
    /* product routine with this matrix.                           */
   blockmat = AZ_matrix_create((Ke_mat->data_org[AZ_N_internal] + 
				 Ke_mat->data_org[AZ_N_border])*2);
   aztec_block_data.N = Ke_mat->data_org[AZ_N_internal] + 
     Ke_mat->data_org[AZ_N_border];
   aztec_block_data.Ke = Ke_mat;
   aztec_block_data.M = M_mat;
   AZ_set_MATFREE(blockmat, &aztec_block_data, aztec_block_matvec);

   /*
   fp = fopen("rhs.dat","a");
   xxx = ML_allocate( 2*Nlocal_edges * sizeof(double));
   for (i=0; i < 2*Nlocal_edges; i ++) xxx[i] = 0.;
   //for (i=0; i < 2*Nlocal_edges; i ++)
   i = Nlocal_edges+20;
   {
      xxx[i] = 1.;
      blockmat->matvec(xxx, rhs, blockmat, proc_config);
      for (j = 0; j < 2*Nlocal_edges; j++)
         if (rhs[j] != 0.0)  fprintf(fp,"%d %d %20.13e\n",i,j,rhs[j]);
      xxx[i] = 0.;
   }
   fclose(fp);
   */
   /*
AZ_random_vector(xxx, blockmat->data_org, proc_config);
blockmat->matvec(xxx, rhs, blockmat, proc_config);
 fp = fopen("xxx.dat","w");
 for (i = 0; i < 2*Nlocal_edges; i++) fprintf(fp,"%20.13e\n",xxx[i]);
 fclose(fp);
 fp = fopen("rhs.dat","w");
 for (i = 0; i < 2*Nlocal_edges; i++) fprintf(fp,"%20.13e\n",rhs[i]);
 fclose(fp);
   */

   /* rst: toggle these two to switch between regular and 2x2 block system */
   AZ_iterate(xxx, rhs, options, params, status, proc_config, blockmat,
   /*  AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat,*/
	      Pmat, scaling);

    options[AZ_pre_calc] = AZ_reuse;
    /*
      if (proc_config[AZ_node] == 0) 
      printf("\n-------- Second solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, scaling); 
      if (proc_config[AZ_node] == 0) 
      printf("\n-------- Third solve with improved convergence test -----\n");
      AZ_iterate(xxx, rhs, options, params, status, proc_config, Ke_mat, Pmat, scaling); 
      */
  }
  solve_time = AZ_second() - start_time;

  if (proc_config[AZ_node] == 0) 
    printf("Solve time = %e, MG Setup time = %e\n", solve_time, setup_time);

  ML_Smoother_Arglist_Delete(&nodal_args);
  ML_Smoother_Arglist_Delete(&edge_args);
  ML_Aggregate_Destroy(&ag);
  ML_Destroy(&ml_edges);
  ML_Destroy(&ml_nodes);
  AZ_free((void *) Ke_mat->data_org);
  AZ_free((void *) Ke_mat->val);
  AZ_free((void *) Ke_mat->bindx);
  AZ_free((void *) global_edge_inds);
  AZ_free((void *) global_edge_externs);
  AZ_free((void *) reordered_edge_externs);
  AZ_free((void *) reordered_glob_edges);
  AZ_scaling_destroy(&scaling);
  if (Ke_mat  != NULL) AZ_matrix_destroy(&Ke_mat);
  if (Pmat  != NULL) AZ_precond_destroy(&Pmat);
  free(xxx);
  free(rhs);
  if (Kn_mat != NULL) AZ_matrix_destroy(&Kn_mat);
  if (garbage != NULL) free(garbage);
  free(Kn_val);
  free(Kn_bindx);
  ML_Operator_Destroy(&Tmat);
  ML_Operator_Destroy(&Tmat_trans);
  ML_MGHierarchy_ReitzingerDestroy(N_levels-2, /* coarsest_level,*/ &Tmat_array, &Tmat_trans_array);


#ifdef ML_MPI
  MPI_Finalize();
#endif
		
  return 0;
		
}

int northwest2d(int i, int j, int n)
{
  return(i + (j)*n);
}
int southwest2d(int i, int j, int n)
{
  if (j == 0) j = n;
  return(i + (j-1)*n);
}
int northeast2d(int i, int j, int n)
{
  if (i == n-1) i = -1;
  return(i+1 + (j)*n);
}
int southeast2d(int i, int j, int n)
{
  if (j == 0) j = n;
  if (i == n-1) i = -1;
  return(i+1 + (j-1)*n);
}
int north2d(int i, int j, int n)
{
  return(i + j*n);
}
int south2d(int i, int j, int n)
{
  if (j == 0) j = n;
  return(i + (j-1)*n);
}
int west2d(int i, int j, int n)
{
  if (j == 0) j = n;
  return(i + (j-1)*n + n*n);
}
int east2d(int i, int j, int n)
{
  if (j == 0) j = n;
  if (i == n-1) i = -1;
  return(i+1 + (j-1)*n + n*n);
}
int inv2dindex(int index, int *i, int *j, int n, int *hor_or_vert)
{
  *hor_or_vert = HORIZONTAL;
  if (index >= n*n) {
    *hor_or_vert = VERTICAL;
    index -= n*n;
  }
  *i = (index%n);
  *j = (index - (*i))/n;
  (*j) =  ((*j)+1)%n;

  return 1;
}
int inv3dindex(int index, int *i, int *j, int *k, int n, int *hor_or_vert)
{
  *hor_or_vert = HORIZONTAL;
  if (index >= n*n*n) {
    *hor_or_vert = VERTICAL;
    index -= n*n*n;
  }
  if (index >= n*n*n) {
    *hor_or_vert = OUT;
    index -= n*n*n;
  }
  *i = (index%n);
  index = (index - (*i))/n;
  *j = (index%n);
  *k = (index - (*j))/n;
  return 1;
}
int inv2dnodeindex(int index, int *i, int *j, int n)
{
  *i = (index%n);
  *j = (index - (*i))/n;
  (*j) =  ((*j)+1)%n;

  return 1;
}
int inv3dnodeindex(int index, int *i, int *j, int *k, int n)
{
  /*n++; */ /* rst dirichlet */
  *i = (index%n);
  index = (index - (*i))/n;
  *j = (index%n);
  *k = (index - (*j))/n;

  return 1;
}

int southwestback3d(int i, int j, int k, int n)
{
  /*n++; */ /* rst dirichlet */
  return(i + j*n + k*n*n);
}
int southwestfront3d(int i, int j, int k, int n)
{
  /*n++; */ /* rst dirichlet */
if (k == n-1) return(-1); /* rst dirichlet */
  if (k == n-1) k = -1;
  return(i + j*n + (k+1)*n*n);
}
int northwestback3d(int i, int j, int k, int n)
{
  /* n++; */  /* rst dirichlet */
if (j == n-1) return(-1); /* rst dirichlet */
  if (j == n-1) j = -1;
  return(i + (j+1)*n + k*n*n);
}
int northwestfront3d(int i, int j, int k, int n)
{
  /* n++; */ /* rst dirichlet */
if (j == n-1) return(-1); /* rst dirichlet */
if (k == n-1) return(-1); /* rst dirichlet */
  if (j == n-1) j = -1;
  if (k == n-1) k = -1;
  return(i + (j+1)*n + (k+1)*n*n);
}
int southeastback3d(int i, int j, int k, int n)
{
  /* n++; */ /* rst dirichlet */
if (i == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  return(i+1 + j*n + k*n*n);
}
int southeastfront3d(int i, int j, int k, int n)
{
  /* n++; */ /* rst dirichlet */
if (i == n-1) return(-1); /* rst dirichlet */
if (k == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  if (k == n-1) k = -1;
  return(i+1 + j*n + (k+1)*n*n);
}
int northeastback3d(int i, int j, int k, int n)
{
  /* n++; */ /* rst dirichlet */
if (i == n-1) return(-1); /* rst dirichlet */
if (j == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  if (j == n-1) j = -1;
  return(i+1 + (j+1)*n + k*n*n);
}
int northeastfront3d(int i, int j, int k, int n)
{
  /* n++; */ /* rst dirichlet */
if (i == n-1) return(-1); /* rst dirichlet */
if (j == n-1) return(-1); /* rst dirichlet */
if (k == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  if (j == n-1) j = -1;
  if (k == n-1) k = -1;
  return(i+1 + (j+1)*n + (k+1)*n*n);
}

int southback3d(int i, int j, int k, int n)
{
if (k == 0) return(-1); /* rst dirichlet */
if (j == 0) return(-1); /* rst dirichlet */
  return(i + j*n + k*n*n);
}
int southfront3d(int i, int j, int k, int n)
{
if (j == 0) return(-1);   /* rst dirichlet */
if (k == n-1) return(-1); /* rst dirichlet */
  if (k == n-1) k = -1;
  return(i + j*n + (k+1)*n*n);
}
int northback3d(int i, int j, int k, int n)
{
if (k == 0) return(-1); /* rst dirichlet */
if (j == n-1) return(-1); /* rst dirichlet */
  if (j == n-1) j = -1;
  return(i + (j+1)*n + k*n*n);
}
int northfront3d(int i, int j, int k, int n)
{
if (j == n-1) return(-1); /* rst dirichlet */
if (k == n-1) return(-1); /* rst dirichlet */
  if (j == n-1) j = -1;
  if (k == n-1) k = -1;
  return(i + (j+1)*n + (k+1)*n*n);
}
int westback3d(int i, int j, int k, int n)
{
if (i == 0) return(-1); /* rst dirichlet */
if (k == 0) return(-1); /* rst dirichlet */
  return(i + j*n + k*n*n + n*n*n);
}
int westfront3d(int i, int j, int k, int n)
{
if (k == n-1) return(-1); /* rst dirichlet */
if (i == 0) return(-1); /* rst dirichlet */
  if (k == n-1) k = -1;
  return(i + j*n + (k+1)*n*n + n*n*n);
}
int eastfront3d(int i, int j, int k, int n)
{
if (k == n-1) return(-1); /* rst dirichlet */
if (i == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  if (k == n-1) k = -1;
  return(i+1 + (j)*n + (k+1)*n*n + n*n*n);
}
int eastback3d(int i, int j, int k, int n)
{
if (k == 0) return(-1); /* rst dirichlet */
if (i == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  return(i+1 + (j)*n + k*n*n + n*n*n);
}
int southwest3d(int i, int j, int k, int n)
{
if (j == 0) return(-1); /* rst dirichlet */
if (i == 0) return(-1); /* rst dirichlet */
  return(i + j*n + k*n*n + 2*n*n*n);
}
int southeast3d(int i, int j, int k, int n)
{
if (j == 0) return(-1); /* rst dirichlet */
if (i == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  return(i+1 + j*n + k*n*n + 2*n*n*n);
}
int northwest3d(int i, int j, int k, int n)
{
if (i == 0) return(-1); /* rst dirichlet */
if (j == n-1) return(-1); /* rst dirichlet */
  if (j == n-1) j = -1;
  return(i + (j+1)*n + k*n*n + 2*n*n*n);
}
int northeast3d(int i, int j, int k, int n)
{
if (i == n-1) return(-1); /* rst dirichlet */
if (j == n-1) return(-1); /* rst dirichlet */
  if (i == n-1) i = -1;
  if (j == n-1) j = -1;
  return(i+1 + (j+1)*n + k*n*n + 2*n*n*n);
}

void compress_matrix(double val[], int bindx[], int N_points)

/*
 * Take an existing MSR matrix which contains some columns labeled '-1'
 * and compress the matrix by eliminating all the columns labeled '-1'.
 * Note: this routine is usually used in conjunction with init_msr().
 */

{

  int free_ptr,start,end,i,j;
  int temp;

  free_ptr = bindx[0];
  start    = bindx[0];

  for (i = 0 ; i < N_points ; i++ ) {
    bindx[i] = bindx[i+1] - bindx[i];
  }
  for (i = 0 ; i < N_points ; i++ ) {
    end   = start + bindx[i] -1;


    for (j = start ; j <= end ; j++ ) {
      if (bindx[j] != -1) {
        val[free_ptr]   = val[j];
        bindx[free_ptr] = bindx[j];
        free_ptr++;
      }
      else {
        bindx[i]--;
      }
    }
    start = end+1;
  }

  start = N_points+1;
  for (i = 0 ; i <= N_points ; i++ ) {
    temp = bindx[i];
    bindx[i] = start;
    start += temp;
  }

} /* compress_matrix */

/****************************************************************************/
/* Generate the new grid hierarchy for 2x2 block matrices from the original */
/* hierarchy stored in 'ml_edges'.                                          */

int ML_Gxn_Hierarchy_ComplexMaxwell(ML *ml_edges, ML_Operator **Tmat_array, 
				    ML_Operator **Tmat_trans_array,
				    ML **new_ml , ML_Operator *originalM)
{

   int mesh_level, old_mesh_level, i, levels;
   ML_Operator *original, *blockmat, *mat, *newM, *lastM;
   ML  *block_ml;
   struct ml_operator_wrapper *ml_operator_wrapper;
   int scale_fact = 2;

   mesh_level = ml_edges->ML_finest_level;

   /* create a new empty hierarchy with the same number of levels */
   /* as in ml_edges.                                             */
   ML_Create(&block_ml,ml_edges->ML_num_levels);
   *new_ml = block_ml;


   /* Make A on the fine level into a 2x2 block matrix */

   levels = 1;
   original = &(ml_edges->Amat[mesh_level]);
   blockmat = &(block_ml->Amat[mesh_level]);
   ML_make_block_matrix(blockmat, original , originalM );
   printf("should not be here\n"); exit(1);
   //   blockmat->M_mat = originalM;
   /* ML_Operator_Print(blockmat,"Ablock"); */

   /* rst: I'm not sure ... but something like the following */
   /* should work for T if needed?                           */
   /*
   original = &(Tmat_array[mesh_level]);
   blockmat = &(blk_Tmat_array[mesh_level]);
   ML_make_block_matrix(blockmat, original);
   original = &(Tmat_trans_array[mesh_level]);
   blockmat = &(blk_Tmat_trans_array[mesh_level]);
   ML_make_block_matrix(blockmat, original);
   */

   lastM = originalM;
   while( ml_edges->SingleLevel[mesh_level].Rmat->to != NULL) {
     levels++;
     old_mesh_level = mesh_level;
     mesh_level = ml_edges->SingleLevel[mesh_level].Rmat->to->levelnum;

     /* Make 2x2 block diagonal P */

     original = &(ml_edges->Pmat[mesh_level]);
     blockmat = &(block_ml->Pmat[mesh_level]);
     ML_make_block_matrix(blockmat, original , NULL );
     /* This stuff sets the 'to' and 'from' field in P */
     /* which indicates from what level we interpolate */
     /* and to what level the interpolation goes.      */
     ML_Operator_Set_1Levels(blockmat, &(block_ml->SingleLevel[mesh_level]), 
			     &(block_ml->SingleLevel[old_mesh_level]));

     /* Make 2x2 block diagonal R */

     original = &(ml_edges->Rmat[old_mesh_level]);
     blockmat = &(block_ml->Rmat[old_mesh_level]);
     ML_make_block_matrix(blockmat, original, NULL );
     /* This stuff sets the 'to' and 'from' field in P */
     /* which indicates from what level we interpolate */
     /* and to what level the interpolation goes.      */
     ML_Operator_Set_1Levels(blockmat,
                             &(block_ml->SingleLevel[old_mesh_level]), 
			     &(block_ml->SingleLevel[mesh_level]));
				  
     /* Make 2x2 block diagonal A */

     original = &(ml_edges->Amat[mesh_level]);
     blockmat = &(block_ml->Amat[mesh_level]);
     /*  newM = ML_Operator_Create(ml_edges->comm);
	 ML_rap(&(ml_edges->Rmat[old_mesh_level]), original, 
	 &(ml_edges->Pmat[mesh_level]), newM, ML_CSR_MATRIX);
     //JJH
     */
     newM = ML_Operator_Create(ml_edges->comm);
     ML_rap(&(ml_edges->Rmat[old_mesh_level]), lastM, 
            &(ml_edges->Pmat[mesh_level]), newM, ML_CSR_MATRIX);
     lastM = newM;
     
     /* comment these two out if you want to do rap */
     ML_make_block_matrix(blockmat, original, newM);
     //     blockmat->M_mat = newM;
   printf("should not be here\n"); exit(1);
     blockmat->getrow->pre_comm = ML_CommInfoOP_Create(); 
              /* ugh. The superlu interface seems to not work */
              /* with an empty communication object */

     /* RAP works too */
     /* ML_Gen_AmatrixRAP(block_ml, old_mesh_level, mesh_level); */

     /* rst: I'm not sure ... but something like the following */
     /* should work for T if needed?                           */
     /*
       original = &(Tmat_array[mesh_level]);
       blockmat = &(blk_Tmat_array[mesh_level]);
       ML_make_block_matrix(blockmat, original);
       original = &(Tmat_trans_array[mesh_level]);
       blockmat = &(blk_Tmat_trans_array[mesh_level]);
       ML_make_block_matrix(blockmat, original);
     */

   }

   return levels;

}

/* ML type of matvec that corresponds to a 2x2 block matrix. */
/* We basically call the same matvec function for the (1,1)  */
/* block matrix two times                                    */
int block_matvec(void *data, int inlen, double invec[],
		      int outlen, double outvec[])
{
  ML_Operator *mat;
  struct ml_operator_wrapper *ml_operator_wrapper;
  double *z;
  int i;

  mat = (ML_Operator *) data;
  ml_operator_wrapper = (struct ml_operator_wrapper *) mat->data;
  /* workspace vector */
  z = (double *) ML_allocate (outlen * sizeof(double) );

  /* multiply by (1,1) block (stiffness)*/
  ml_operator_wrapper->diag_matvec(ml_operator_wrapper->diag_matvec_data, inlen/2,
				     invec, outlen/2, outvec);

  /* multiply by (1,2) block (-mass)*/
  if (ml_operator_wrapper->offdiag_matvec != NULL) {
    ml_operator_wrapper->offdiag_matvec(ml_operator_wrapper->offdiag_matvec_data, inlen/2,
				     &(invec[inlen/2]), outlen/2, z);
    for (i=0; i< outlen/2; i++) outvec[i] -= z[i];
  }
  /***********/
  /* multiply by (2,2) block (stiffness)*/
  ml_operator_wrapper->diag_matvec(ml_operator_wrapper->diag_matvec_data, inlen/2,
				     &(invec[inlen/2]), outlen/2, 
				     &(outvec[outlen/2]));
  /* multiply by (2,1) block (mass)*/
  if (ml_operator_wrapper->offdiag_matvec != NULL) {
    ml_operator_wrapper->offdiag_matvec(ml_operator_wrapper->offdiag_matvec_data, inlen/2,
				     invec, outlen/2, z);
    for (i=0; i < outlen/2; i++) outvec[i+outlen/2] += z[i];
  }

  ML_free(z);
}

/* ML type of getrow that corresponds to a 2x2 block matrix. */
/* We basically call the getrow function for the (1,1). If   */
/* we are getting a row in the (2,2) block, we have to first */
/* change the requested row so that it corresponds to the    */
/* (1,1) block and then we have to post-process the returned */
/* row information so that the columns are properly shifted  */
/* for the (2,2) block.                                      */
int blockdiag_getrow(void *data, int N_requested,
			  int requested_rows[], int allocated,
			  int columns[], double values[], int row_lengths[])
{
  ML_Operator *mat;
  struct ml_operator_wrapper *ml_operator_wrapper;
  int newrow, status, i;

  mat = (ML_Operator *) data;
  ml_operator_wrapper = (struct ml_operator_wrapper *) mat->data;

  if (N_requested != 1) return(1);

  if (requested_rows[0] < mat->outvec_leng/2)
    return(ml_operator_wrapper->diag_getrow(ml_operator_wrapper->diag_getrow_data,
				      N_requested, requested_rows, allocated, 
				      columns, values, row_lengths));
  else {
    /* shift requested row to (1,1) block */
    newrow = requested_rows[0] - mat->outvec_leng/2;
    status = ml_operator_wrapper->diag_getrow(ml_operator_wrapper->diag_getrow_data,
					N_requested, &newrow, allocated, 
					columns, values, row_lengths);

    /* post process data so that columns correspond to (2,2) block */
    if (status != 0) {
      for (i = 0; i < row_lengths[0]; i++) 
	columns[i] += mat->invec_leng/2;
    }
    return(status);
  }
}

/* Getrow for 2x2 block matrix */

int block_getrow(void *data, int N_requested,
			  int requested_rows[], int allocated,
			  int columns[], double values[], int row_lengths[])
{
  ML_Operator *mat;
  struct ml_operator_wrapper *ml_operator_wrapper;
  int newrow, status = 1, i;
  int *workcol;
  double *workval;
  int work_lengths[1];

  mat = (ML_Operator *) data;
  ml_operator_wrapper = (struct ml_operator_wrapper *) mat->data;

  workcol = ml_operator_wrapper->cols;
  workval = ml_operator_wrapper->vals;
  work_lengths[0] = 0;
  row_lengths[0] = 0;

  if (N_requested != 1) return(1);

  /* (1,1) and (1,2) blocks */
  if (requested_rows[0] < mat->outvec_leng/2)
  {
    /* (1,1) block */
    status = ml_operator_wrapper->diag_getrow(ml_operator_wrapper->diag_getrow_data,
				      N_requested, requested_rows, allocated, 
				      columns, values, row_lengths);
    if (status == 0) return(status);
    if (ml_operator_wrapper->offdiag_getrow != NULL) {
      /* (1,2) block */
      status = ml_operator_wrapper->offdiag_getrow(ml_operator_wrapper->offdiag_getrow_data,
				      N_requested, requested_rows, 100, 
				      workcol, workval, work_lengths);
      if (status == 0) return(status);
    }
    for (i=0; i< work_lengths[0]; i++) {
      workcol[i] += mat->invec_leng/2;
      workval[i] = -workval[i];     /* josh says -M */
    }

  }
  /* (2,1) and (2,2) blocks */
  else {
    /* shift requested row to (1,1) block */
    newrow = requested_rows[0] - mat->outvec_leng/2;
    /* (2,2) block */
    status = ml_operator_wrapper->diag_getrow(ml_operator_wrapper->diag_getrow_data,
					N_requested, &newrow, 100, 
					workcol, workval, work_lengths);
    if (status == 0) return(status);

    /* post process data so that columns correspond to (2,2) block */
    if (status != 0) {
      for (i = 0; i < work_lengths[0]; i++) 
	workcol[i] += mat->invec_leng/2;
    }

    /* (2,1) block */
    if (ml_operator_wrapper->offdiag_getrow != NULL) {
      status = ml_operator_wrapper->offdiag_getrow(
		                      ml_operator_wrapper->offdiag_getrow_data,
				      N_requested, &newrow, allocated, 
				      columns, values, row_lengths);
      if (status == 0) return(status);

    }
  }
  /*make sure columns is long enough for the concatenation*/
  if (row_lengths[0] + work_lengths[0] > allocated) {
    return(0);
  }

  for (i=0; i<work_lengths[0]; i++) {
    columns[i + row_lengths[0]] = workcol[i];
    values[i + row_lengths[0]] = workval[i];
  }
  row_lengths[0] += work_lengths[0]; 
  return(1);
}

/******************************************************************************/


/* Convert 'original' into a 2x2 block matrix of the following form: */
/*                         (original    0     )                      */
/*      blockmat =         (  0       original)                      */
/*                                                                   */
/* NOTE: This routine needs to be modified to handle matrices of the */
/* form:                                                             */
/*                         ( orig1   -orig2   )                      */
/*      blockmat =         ( orig2    orig1   )                      */

int  ML_make_block_matrix(ML_Operator *blockmat, ML_Operator *original1,
                          ML_Operator *original2)
{

  struct ml_operator_wrapper *ml_operator_wrapper;
  int scale_fact = 2; /* b/c we're building a 2x2 system */

  /* make a new data structure which `blockdiag_matvec' and `blockdiag_getrow' */
  /* will use. Shove the new data structure into blockmat along with the       */
  /* 'blockdiag_matvec' and 'blockdiag_getrow'.                                */
  ML_Operator_Init(blockmat,original1->comm);
  blockmat->max_nz_per_row = 50;

  ml_operator_wrapper = (struct ml_operator_wrapper *) AZ_allocate(
				       sizeof(struct ml_operator_wrapper));
  ml_operator_wrapper->cols = (int *) ML_allocate(100 * sizeof(int));
  ml_operator_wrapper->vals = (double *) ML_allocate(100 * sizeof(double));

  /* setup matvec for diagonal part */
  if (original1->matvec->ML_id == ML_INTERNAL) {
    ml_operator_wrapper->diag_matvec = original1->matvec->internal;
    ml_operator_wrapper->diag_matvec_data = original1;
  }
  else {
    ml_operator_wrapper->diag_matvec = original1->matvec->external;
    ml_operator_wrapper->diag_matvec_data = original1->data;
  } //if

  /* setup matvec for offdiagonal part */

  ml_operator_wrapper->offdiag_matvec = NULL;
  ml_operator_wrapper->offdiag_matvec_data = NULL;
  if (original2 != NULL) {
    if (original2->matvec->ML_id == ML_INTERNAL) {
      ml_operator_wrapper->offdiag_matvec = original2->matvec->internal;
      ml_operator_wrapper->offdiag_matvec_data = original2;
    }
    else {
      ml_operator_wrapper->offdiag_matvec = original2->matvec->external;
      ml_operator_wrapper->offdiag_matvec_data = original2->data;
    }
  } //if

  ML_Operator_Set_ApplyFuncData(blockmat, scale_fact*original1->invec_leng, 
				scale_fact*original1->outvec_leng, ML_INTERNAL,
				ml_operator_wrapper,
				scale_fact*original1->outvec_leng,
				block_matvec,0);

  /* set getrow for diagonal block */
  if (original1->getrow->ML_id == ML_INTERNAL) {
    ml_operator_wrapper->diag_getrow = original1->getrow->internal;
    ml_operator_wrapper->diag_getrow_data = original1;
  }
  else {
    ml_operator_wrapper->diag_getrow = original1->getrow->external;
    ml_operator_wrapper->diag_getrow_data = original1->data;
  }

  ml_operator_wrapper->offdiag_getrow = NULL;
  ml_operator_wrapper->offdiag_getrow_data = NULL;

  /* set getrow for offdiagonal block */
  if (original2 != NULL) {
    if (original2->getrow->ML_id == ML_INTERNAL) {
      ml_operator_wrapper->offdiag_getrow = original2->getrow->internal;
      ml_operator_wrapper->offdiag_getrow_data = original2;
    }
    else {
      ml_operator_wrapper->offdiag_getrow = original2->getrow->external;
      ml_operator_wrapper->offdiag_getrow_data = original2->data;
    }
  }

  ML_Operator_Set_Getrow(blockmat, ML_INTERNAL, scale_fact*
			 original1->outvec_leng, block_getrow /*blockdiag_getrow*/);
  
  return 1;
}

/******************************************************************************/

/* New aztec block matvec function */
void aztec_block_matvec(double *x, double *y, AZ_MATRIX *Amat,
                             int proc_config[])
{
  struct aztec_block_data *aztec_block_data;
  double *z;
  int i;

  aztec_block_data = (struct aztec_block_data *) AZ_get_matvec_data(Amat);
  z = (double *) AZ_allocate((aztec_block_data->N+1)*sizeof(double));
  aztec_block_data->Ke->matvec( x, y, aztec_block_data->Ke, proc_config);

  if (aztec_block_data->M != NULL)
  {
    aztec_block_data->M->matvec( &(x[aztec_block_data->N]), z, 
                                 aztec_block_data->M, proc_config);
    for (i = 0; i < aztec_block_data->N; i++) y[i] -= z[i];
  }

  aztec_block_data->Ke->matvec( &(x[aztec_block_data->N]), 
                                &(y[aztec_block_data->N]),
				aztec_block_data->Ke, proc_config);

  if (aztec_block_data->M != NULL)
  {
    aztec_block_data->M->matvec( x, z, aztec_block_data->M, proc_config);
    for (i = 0; i < aztec_block_data->N; i++) y[i+aztec_block_data->N] += z[i];
  }
  else
     printf("Block matrix appears to be diagonal!!\n");
  AZ_free(z);
}

/* bogus fortran functions */
/*
int f_iob() { printf("not here f_iob\n"); exit; }
int f_stop() { printf("not here f_stop\n"); exit; }
int f_powdi() { printf("not here f_powdi\n"); exit; }
int f_cpystr() { printf("not here f_cpystr\n"); exit; }
int f_concat() { printf("not here f_concat\n"); exit; }
*/

