#define AZ_MPI
#ifdef AZ_MPI
#include "mpi.h"
#endif
#include <stdio.h>
#include <math.h>
#include "az_aztec.h"
#include "ml_include.h"

#define ML_NEWNORM
#define ML_SCALEIT
/*
#define ML_JANUS_FORT(iiii) iiii##_
*/
#define ML_JANUS_FORT(iiii) iiii

/**************************************************************************/
/* structure/functions to wrap greg's routines.                           */
/**************************************************************************/
struct wrap_greg_data {
  double **ca;
  int **index, *mcol, *nrows, *ncols, *nodes, *nput, *nget, **nsol, *proc;
  int Nghost;
};
extern void aztec_matvec_wrap_of_greg(double *, double *, AZ_MATRIX *, int *);
extern int  aztec_comm_wrap_of_greg  (double *, AZ_MATRIX *);
extern int  aztec_getrow_wrap_of_greg(int *, double *, int *, AZ_MATRIX *, 
				      int , int *, int );
extern int ml_grad_matvec_wrap_of_greg(void *, int , double *, int, double *);
extern int ml_grad_comm_wrap_of_greg( double *x, void *data);
extern int ml_grad_getrow_wrap_of_greg(void *,int,int *,int,int *,double *,
				       int *);
extern int dump_greg(ML_Operator *Ke, ML_Operator *Kn, ML_Operator *M,
		     ML_Operator *T, int Nghost_nodes, int Nghost_edges,
		     double *, double *);
extern int setup_agg_MIS_dump(ML_Operator *Kn, int Nghost_nodes);
ML_JANUS_FORT(extern void ray_matvec_)(/* ML_Operator *Amat, */ double ERF_x[], double ERF_y[]);
extern void eigen_(int*, int*, int*, int*, int*, int*, char*,
                   double*, double*);

/*
clock_t gtime;
double dtime = 0.;

void time_it_(int *mypid)
{
  if (dtime == 0.) {
    gtime = clock();
    dtime = 1.;
  }
  else {
    gtime = clock();
    dtime = (double) ( gtime / 1000000);
    if (*mypid == 0) {
      printf(">>>>>>>>>>>>>>> total time = %e<<<<<<<<<<<<<<\n",dtime);
      fflush(stdout);
    }
    dtime = 0.;
  }
}
*/
ML_Operator *global_dude;
ML_Operator *globbie;     /*T transpose */
/*****************************************************************************/
/* structure/function to wrap Aztec subblocks into a large equivalent real   */
/* form Aztec matrix and then invoke Aztec using ML as a preoconditioner     */
/*****************************************************************************/
int *update_index, *update, *extern_index, *external;
ML_JANUS_FORT(void wrappers_for_mg_)(double **ca, int **index, int mcol[], int *Nedges, 
		      int nodes[], int nput[], int nget[], 
		      int **nsol, int *proc, double **ca_imag, 
		      double complex_x[], double complex_rhs[],
		      int *Nnodes, int nodespot[], int nputpot[],
		      int ngetpot[], int **nsolpot, double **grad_mat,
		      int **indexpot, int mcolpot[])
{
  struct wrap_greg_data        *wrap_real_greg_data;
  struct wrap_greg_data        *wrap_imag_greg_data;
  struct wrap_greg_data        *wrap_grad_greg_data;
  struct AZ_MAT_blockmat_data  *AZ_MAT_blockmat_data;
  AZ_MATRIX                    *Ke, *M;
  AZ_MATRIX                    *ERF_AZmat;    /* Equivalent real form matrix.*/
  AZ_PRECOND                   *Prec = NULL;  /* Aztec preconditioner. It    */
                                              /* will be properly set by     */
                                              /* AZ_set_ML_preconditioner    */
  ML_Operator                  *M_mat, *ERF_MLmat;
  ML_Operator                  *S_mat;
  ML                           *ml_edges, *ml_ERF;
  ML_Operator                  *Tmat;

  ML_Aggregate *ag;
  int Nits_per_presmooth=1;         /* # of pre & post smoothings per level */
  int smoothPe_flag = ML_YES;       /* ML_YES: smooth tentative prolongator */
                                    /* ML_NO: don't smooth prolongator      */

  double *ERF_x, *ERF_y, *grad_scale;
  /*JJH*/
  double  *vec1, *vec2, resnorm, rhsnorm, gresnorm, grhsnorm;
  int iii;
  double result, result1, xnorm;
  /*--JJH*/
  int MaxMgLevels = 10, Nlocal_edges, Nlocal_nodes, i,Nghost =0,Nghost_node =0;
  int Nlevels;
  FILE *fp,*fp1;
  double tdiag, ttdiag;

  /* Aztec stuff */

  int proc_config[AZ_PROC_SIZE];
  int options[        /* Used to pass parameters to Aztec.                   */
     AZ_OPTIONS_SIZE];/*                                                     */
double               /*                                                      */
   params[           /* Used to pass in parameters to Aztec                  */
     AZ_PARAMS_SIZE],/*                                                      */
   status[           /* Returned from AZ_iterate() indicating success or     */
     AZ_STATUS_SIZE];/* failure of the iterative solvers                     */
char str [80];
 int colInd[15], ncnt;
 double colVal[15], lambda;
   ML_Krylov   *kdata;
   int j, count;
  struct ML_CSR_MSRdata *data;
  int *rowptr;
  double *values, *diagonal;
  double *Dir_bdry, *copy_ERF_x, *diag_of_M, *diag_of_S;

  /***************************************************************************/
  /* Select Hiptmair relaxation subsmoothers for the nodal and edge problems */
  /* Choices include                                                         */
  /*   1) ML_Gen_Smoother_SymGaussSeidel: this corresponds to a processor    */
  /*      local version of symmetric Gauss-Seidel/SOR. The number of sweeps  */
  /*      can be set via either 'edge_its' or 'nodal_its'. The damping can   */
  /*      be set via 'edge_omega' or 'nodal_omega'. When set to ML_DDEFAULT, */
  /*      the damping is set to '1' on one processor. On multiple processors */
  /*      a lower damping value is set. This is needed to converge processor */
  /*      local SOR.                                                         */
  /*   2) ML_Gen_Smoother_MLS: this corresponds to polynomial relaxation.    */
  /*      The degree of the polynomial is set via 'edge_its' or 'nodal_its'. */
  /*      If the degree is '-1', Marian Brezina's MLS polynomial is chosen.  */
  /*      Otherwise, a Chebyshev polynomial is used over high frequencies    */
  /*      [ lambda_max/alpha , lambda_max]. Lambda_max is computed. 'alpha'  */
  /*      is hardwired in this example to correspond to twice the ratio of   */
  /*      unknowns in the fine and coarse meshes.                            */
  /*                                                                         */
  /* Using 'hiptmair_type' (see comments below) it is also possible to choose*/
  /* when edge and nodal problems are relaxed within the Hiptmair smoother.  */
  /***************************************************************************/

  void  *edge_smoother=(void *)     /* Edge relaxation:                     */

    ML_Gen_Smoother_MLS;
                                    /*     ML_Gen_Smoother_MLS              */
                                    /*     ML_Gen_Smoother_SymGaussSeidel   */
                                    /*ML_Gen_Smoother_SymGaussSeidelSequential*/
  void *nodal_smoother=(void *)     /* Nodal relaxation                     */
    ML_Gen_Smoother_MLS;
                                    /*     ML_Gen_Smoother_MLS              */
                                    /*     ML_Gen_Smoother_SymGaussSeidel   */
                                    /*ML_Gen_Smoother_SymGaussSeidelSequential*/

  int  edge_its = 2;                /* Iterations or polynomial degree for  */
  int  nodal_its = 2;               /* edge/nodal subsmoothers.             */
  double nodal_omega = ML_DDEFAULT,  /* SOR damping parameter for node/edge  */
         edge_omega  = ML_DDEFAULT;  /* subsmoothers (see comments above).   */
  int   hiptmair_type=HALF_HIPTMAIR;/* FULL_HIPTMAIR: each invocation       */
                                    /*     smoothes on edges, then nodes,   */
                                    /*     and then once again on edges.    */
                                    /* HALF_HIPTMAIR: each pre-invocation   */
                                    /*     smoothes on edges, then nodes.   */
                                    /*     Each post-invocation smoothes    */
                                    /*     on nodes then edges. .           */

  double       edge_coarsening_rate, node_coarsening_rate;
  void         **edge_args, **nodal_args;
  int          Nfine_edge, Ncoarse_edge, Nfine_node, Ncoarse_node;
  int          level, coarsest_level, itmp;

  /* What should be hidden from the user */

  ML *ml_nodes;
  ML_Operator  *Tmat_trans, **Tmat_array, **Tmat_trans_array;
  struct ML_Operator_blockmat_data *block_data;
  double *M_diag;

  /* variables for Greg's eigenvalue routine "eigen_".*/

  int err_flag=0;                    /* if not 0, indicates error condition */
  int io_flag=0;                   /* set to > 0 for more screen output */
  double real_eigenval;
  double imag_eigenval;
  double realsav, imagsav;
  char c[10];                    /* arpack option */
  

  /************************* start of execution *****************************/

  MPI_Barrier(MPI_COMM_WORLD);
  ML_Set_PrintLevel(5);
  Nlocal_edges = *Nedges;

  /* fill processor information for Aztec */

#ifdef AZ_MPI
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

  /* Compute the total number of ghost nodes for later use */

  for (i = 0; i < 12 ; i++) {
    if (nodes[i] != -1) Nghost += nget[i];
  }

  /* Compute the total number of nodal ghost nodes for later use */

  for (i = 0; i < 6 ; i++) {
    if (nodespot[i] != -1) Nghost_node += ngetpot[i];
  }

  /* We need to scale the gradient operator so that it contains */
  /* only 1's or -1's. We also need to scale S, M, and the      */
  /* vectors to keep things consistent. Basically, the scaled   */
  /* system corresponds to                                      */
  /*      (D^-1 S D^-1 + D^-1 M D^-1) D u = D^-1 b              */
  /* with  D^-1 S D^-1 D T = 0.                                 */

  grad_scale = (double *) malloc(sizeof(double)*(*Nedges+Nghost));
  for (i = 0; i < Nghost + *Nedges; i++) grad_scale[i] = 1.;
  i = Nghost + *Nedges;
  scale3_(grad_mat,grad_scale,indexpot,mcolpot,Nedges,nodes,
	  nput,nget,nsol,proc, &i,  ca, index, mcol, ca_imag);

  /* copy complex vectors into the right form for the equivalent real form. */
  /* Additionally, scale them to be consistent with other operators.        */

  ERF_x = (double *) malloc(sizeof(double)*2*( *Nedges + Nghost));
  ERF_y = (double *) malloc(sizeof(double)*2*( *Nedges ));
  vec1 = (double *) malloc(sizeof(double)*2*( *Nedges ));

  for (i = 0; i < *Nedges; i++) {
    ERF_x[i         ] = complex_x[2*i    ]*grad_scale[i];
    ERF_x[i+ *Nedges] = complex_x[2*i+1  ]*grad_scale[i];
    ERF_y[i         ] = complex_rhs[2*i  ]/grad_scale[i];
    ERF_y[i+ *Nedges] = complex_rhs[2*i+1]/grad_scale[i];
  }

  /* Build aztec version of the real part (curl-curl) of Greg's matrix */
  /* First bundle greg's data into one structure, then supply a        */
  /* matvec-wrapper that just calls Greg's routine.                    */

  wrap_real_greg_data = (struct wrap_greg_data *) 
                           malloc(sizeof(struct wrap_greg_data));
  wrap_real_greg_data->ca    = ca;
  wrap_real_greg_data->index = index;
  wrap_real_greg_data->mcol  = mcol;
  wrap_real_greg_data->nrows = Nedges;
  wrap_real_greg_data->ncols = Nedges;
  wrap_real_greg_data->nodes = nodes;
  wrap_real_greg_data->nput  = nput;
  wrap_real_greg_data->nget  = nget;
  wrap_real_greg_data->nsol  = nsol;
  wrap_real_greg_data->proc  = proc;
  wrap_real_greg_data->Nghost= Nghost;

  Ke = AZ_matrix_create( *Nedges );
  AZ_set_MATFREE(Ke, wrap_real_greg_data, aztec_matvec_wrap_of_greg);
  AZ_set_MATFREE_getrow(Ke, (void *) wrap_real_greg_data, 
			aztec_getrow_wrap_of_greg, aztec_comm_wrap_of_greg, 
			Nghost, proc_config);

  /* Build aztec version of the imag part (mass) of Greg's matrix      */
  /* First bundle greg's data into one structure, then supply a        */
  /* matvec-wrapper that just calls Greg's routine.                    */

  wrap_imag_greg_data = (struct wrap_greg_data *) 
                           malloc(sizeof(struct wrap_greg_data));
  wrap_imag_greg_data->ca    = ca_imag;
  wrap_imag_greg_data->index = index;
  wrap_imag_greg_data->mcol  = mcol;
  wrap_imag_greg_data->nrows = Nedges;
  wrap_imag_greg_data->ncols = Nedges;
  wrap_imag_greg_data->nodes = nodes;
  wrap_imag_greg_data->nput  = nput;
  wrap_imag_greg_data->nget  = nget;
  wrap_imag_greg_data->nsol  = nsol;
  wrap_imag_greg_data->proc  = proc;
  wrap_imag_greg_data->Nghost= Nghost;

  M = AZ_matrix_create( *Nedges );
  AZ_set_MATFREE(M, wrap_imag_greg_data, aztec_matvec_wrap_of_greg);
  AZ_set_MATFREE_getrow(M, (void *) wrap_imag_greg_data, 
			aztec_getrow_wrap_of_greg, NULL, Nghost, proc_config);

  /* Build aztec version of the equivalent real form. First build a    */
  /* structure that holds the two individual matrices, and then supply */
  /* a block matrix routine that invokes the matvecs of the subblocks  */
  /* and glues them together appropriately.                            */

  AZ_MAT_blockmat_data         = (struct AZ_MAT_blockmat_data *) 
                                   malloc(sizeof(struct AZ_MAT_blockmat_data));
  AZ_MAT_blockmat_data->N      = *Nedges;
  AZ_MAT_blockmat_data->Ke     = Ke;
  AZ_MAT_blockmat_data->M      = M;
  AZ_MAT_blockmat_data->Nghost = Nghost;
  ERF_AZmat = AZ_matrix_create( (*Nedges)*2 );
  AZ_set_MATFREE(ERF_AZmat, AZ_MAT_blockmat_data, AZ_block_matvec);
  AZ_set_MATFREE_getrow(ERF_AZmat, (void *) AZ_MAT_blockmat_data, 
			NULL, NULL, Nghost, proc_config);

  /* Start building the ML operators */

   ML_Create(&ml_edges, MaxMgLevels);

   /* Take Aztec matrix 'Ke', convert it to an ML matrix and put */
   /* it in level 'MaxMgLevels-1' in the multigrid hierarchy.    */

   AZ_ML_Set_Amat(ml_edges, MaxMgLevels-1, Nlocal_edges, Nlocal_edges, 
		  Ke, proc_config); 

   M_mat = ML_Operator_Create(ml_edges->comm);
   AZ_convert_aztec_matrix_2ml_matrix(M, M_mat, proc_config);

   S_mat = ML_Operator_Create(ml_edges->comm);
   AZ_convert_aztec_matrix_2ml_matrix(Ke, S_mat, proc_config);

   /* Build ML version of the gradient operator. First bundle */
   /* greg's data into one structure, then supply a matvec    */
   /* wrapper that calls Greg's routine.                      */

   wrap_grad_greg_data = (struct wrap_greg_data *) 
                           malloc(sizeof(struct wrap_greg_data));
   wrap_grad_greg_data->ca    = grad_mat;
   wrap_grad_greg_data->index = indexpot;
   wrap_grad_greg_data->mcol  = mcolpot;
   wrap_grad_greg_data->nrows = Nedges;
   wrap_grad_greg_data->ncols = Nnodes;
   wrap_grad_greg_data->nodes = nodespot;
   wrap_grad_greg_data->nput  = nputpot;
   wrap_grad_greg_data->nget  = ngetpot;
   wrap_grad_greg_data->nsol  = nsolpot;
   wrap_grad_greg_data->proc  = proc;
   wrap_grad_greg_data->Nghost= Nghost_node;
   
   Tmat = ML_Operator_Create(ml_edges->comm);
   
   ML_Operator_Set_ApplyFuncData(Tmat, *Nnodes, *Nedges, ML_INTERNAL,
				 (void *) wrap_grad_greg_data, *Nedges,
				 ml_grad_matvec_wrap_of_greg, 0);
   
   ML_Operator_Set_Getrow(Tmat,ML_INTERNAL,*Nedges,ml_grad_getrow_wrap_of_greg);
   ML_CommInfoOP_Generate(&(Tmat->getrow->pre_comm),ml_grad_comm_wrap_of_greg,
			  wrap_grad_greg_data, ml_edges->comm,Tmat->invec_leng,
			  Nghost_node);
   
   /* Zero out rows of T that correspond to Dirichlet points */
   /* I think we also zero out columns of Ke                 */

   Dir_bdry = (double *) malloc(sizeof(double)*(*Nedges + Nghost));
   for (i = 0; i < *Nedges; i++) {
     Dir_bdry[i] = 0.;
     ML_Operator_Getrow(&(ml_edges->Amat[MaxMgLevels-1]),1,&i,15,colInd,
			colVal,&ncnt);
     if ( ncnt == 1) { mcolpot[i] = 0; Dir_bdry[i] = 1.; }
  }
   
#ifdef out
  sprintf(str,"greg_x%d",*proc);
  fp = fopen(str,"w");
  sprintf(str,"greg_b%d",*proc);
  fp1 = fopen(str,"w");
  for (i = 0 ; i < Nlocal_edges; i++) {
    /*    if (Dir_bdry[i] == 1.) { */
      fprintf(fp,"%d %20.13e   %20.13e\n",i,complex_x[2*i],complex_x[2*i +1 ]);
      fprintf(fp1,"%d %20.13e   %20.13e\n",i,complex_rhs[2*i],complex_rhs[2*i +1 ]);
      /*    } */
  }
  fflush(fp);  fclose(fp);
  fflush(fp1);  fclose(fp1);
  printf("finished\n");
  fflush(stdout); 
  while(1 == 1);
#endif

   /* Get diagonal of M */
   for (i = 0; i < Nlocal_edges; i++) complex_x[i] = 1.;
   diag_of_M = &(complex_x[Nlocal_edges]); /* complex_x currently unused */
   ML_Operator_Apply(M_mat, M_mat->invec_leng, complex_x, 
		     M_mat->outvec_leng, diag_of_M);

   /* Need to make a copy of the first half of 'ERF_x' as there is */
   /* not enough room for ghost variables at the end of vector */

   copy_ERF_x = complex_x; /* complex_x currently unused */
   for (i = 0; i < Nlocal_edges; i++) copy_ERF_x[i] = ERF_x[i]; 

  ML_JANUS_FORT(zero_dir_)(wrap_real_greg_data->ca, Dir_bdry, wrap_real_greg_data->index,
	     wrap_real_greg_data->mcol,  wrap_real_greg_data->nrows,
	     wrap_real_greg_data->nodes, wrap_real_greg_data->nput,
	     wrap_real_greg_data->nget,  wrap_real_greg_data->nsol,
	     wrap_real_greg_data->proc,  copy_ERF_x, &(ERF_x[*Nedges]),
	     ERF_y, &(ERF_y[*Nedges]),diag_of_M);


  /* copy back into ERF_x */

  for (i = 0; i < Nlocal_edges; i++) ERF_x[i] = copy_ERF_x[i];

  ML_free(Dir_bdry);

  /********************************************************************/
  /*                      Set up Tmat_trans                           */
  /*------------------------------------------------------------------*/
  Tmat_trans = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Transpose_byrow(Tmat, Tmat_trans);

  /* More clean up for Dirichlet boundary conditions. Count   */
  /* the number of nonzeros per row in T^T and zero out any   */
  /* row that does not have 6 nonzeros.                       */

  data = (struct ML_CSR_MSRdata *) Tmat_trans->data;
  values = data->values;
  rowptr = data->rowptr;
  for (i = 0; i < Tmat_trans->outvec_leng; i++) {
    ML_Operator_Getrow(Tmat_trans,1,&i,15,colInd, colVal,&ncnt);
    count = 0;
    for (j = 0; j < ncnt; j++)
      if (colVal[j] != 0.0) count++;
    if (count != 6)
      for (j = rowptr[i]; j < rowptr[i+1]; j++) values[j] = 0.0;
  }

  /* We need to rebuild Tmat to reflect the zero'd out rows */

  ML_Operator_Destroy(&Tmat);
  Tmat = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Transpose_byrow(Tmat_trans, Tmat);

#ifdef ML_NEWNORM
  /* for applying new T norm */
  globbie = Tmat_trans;
#endif

  /* Set up nodal hierarchy */

  ML_Create(&ml_nodes, MaxMgLevels);
  ML_2matmult(Tmat_trans, Tmat, &(ml_nodes->Amat[MaxMgLevels-1]),
              ML_CSR_MATRIX);

  /********************************************************************/
  /* Set some ML parameters.                                          */
  /*------------------------------------------------------------------*/

  ML_Set_Tolerance(ml_edges, 1.0e-8);
  ML_Aggregate_Create( &ag );
   ML_Aggregate_Set_CoarsenScheme_Uncoupled(ag); 
   /*ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(ag);*/
  /* ML_Aggregate_Set_CoarsenScheme_MIS(ag); */
  ML_Aggregate_Set_DampingFactor(ag, 0.0); /* must use 0 for maxwell */
  ML_Aggregate_Set_MaxCoarseSize(ag, 30);
  ML_Aggregate_Set_Threshold(ag, 0.0);

  /********************************************************************/
  /*                      Build the MG hierarchy                      */
  /*------------------------------------------------------------------*/

  /*  setup_agg_MIS_dump(&(ml_nodes->Amat[MaxMgLevels-1]), Nghost_node);*/
  Nlevels=ML_Gen_MGHierarchy_UsingReitzinger(ml_edges, ml_nodes,MaxMgLevels-1,
                                             ML_DECREASING,ag,Tmat,Tmat_trans, 
                                             &Tmat_array,&Tmat_trans_array, 
                                             smoothPe_flag, 1.5);
  ML_Gen_Hierarchy_ComplexMaxwell(ml_edges, &ml_ERF, M_mat);
  /* ray_matvec( &(ml_ERF->Amat[MaxMgLevels-1]), ERF_x, ERF_y); */
  global_dude = &(ml_ERF->Amat[MaxMgLevels-1]);


  /*
  if (ml_edges->comm->ML_mypid == 0) printf("Skipping calls to Arpack\n");
  */
/*  --ARPACK calls follow 
  sprintf(c,"%s","LR");
  io_flag = 0;

  eigen_(&(Tmat->outvec_leng),
         &(Tmat->comm->ML_mypid),
	 &(Tmat->comm->ML_nprocs),
	 (int *) &(Tmat->comm->USR_comm),
	 &io_flag,
	 &err_flag,
	 c,
	 &real_eigenval,
         &imag_eigenval );

  realsav = real_eigenval;

  sprintf(c,"%s","LI");
  io_flag = 0;

  eigen_(&(Tmat->outvec_leng),
         &(Tmat->comm->ML_mypid),
	 &(Tmat->comm->ML_nprocs),
	 (int *) &(Tmat->comm->USR_comm),
	 &io_flag,
	 &err_flag,
	 c,
	 &real_eigenval,
         &imag_eigenval );

  real_eigenval = realsav;
*/

  /* The MLS smoother needs the largest eigenvalue of the matrix. */
  /* Normally, this happens automatically within ML by checking   */
  /* if the eigenvalue is not already defined and then calling a  */
  /* CG routine to compute it. However, CG won't work for the     */
  /* equivalent real system. Instead, we will manually compute the*/
  /* largest eigenvalue of the pieces to estimate the largest     */
  /* eigenvalue of the equivalent real system.                    */
  for (j=ml_edges->ML_finest_level; j>ml_edges->ML_coarsest_level; j--) {

     block_data = (struct ML_Operator_blockmat_data *) ml_ERF->Amat[j].data;

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
     ml_ERF->Amat[j].lambda_max     = ml_edges->Amat[j].lambda_max;
     ml_ERF->Amat[j].lambda_max_img = block_data->M_mat->lambda_max;
     if (ml_edges->comm->ML_mypid == 0) {
       printf("\n\nlevel %d:  (lambda_max_real,lambda_max_imag) = %e  %e\n\n",
              j,ml_edges->Amat[j].lambda_max, block_data->M_mat->lambda_max);
       printf("(level %d) numer = %e, denom = %e\n",j,ttdiag,tdiag);
       fflush(stdout);
     }
  } /* for j=ml_edges->ML_finest_level... */

  /* Set the Hiptmair subsmoothers */

  if (nodal_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    ML_Smoother_Arglist_Set(nodal_args, 1, &nodal_omega);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) {
    edge_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    ML_Smoother_Arglist_Set(edge_args, 1, &edge_omega);
  }
  if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    Nfine_node = Tmat_array[MaxMgLevels-1]->invec_leng;
    ML_gsum_scalar_int(&Nfine_node, &itmp, ml_ERF->comm);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
    edge_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    Nfine_edge = Tmat_array[MaxMgLevels-1]->outvec_leng;
    ML_gsum_scalar_int(&Nfine_edge, &itmp, ml_ERF->comm);
  }
  if (nodal_smoother == (void *) ML_Gen_Smoother_SymGaussSeidelSequential) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    nodal_omega = 1.0;
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    ML_Smoother_Arglist_Set(nodal_args, 1, &nodal_omega);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_SymGaussSeidelSequential) {
    edge_args = ML_Smoother_Arglist_Create(2);
    edge_omega = 1.0;
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    ML_Smoother_Arglist_Set(edge_args, 1, &edge_omega);
  }
  if (nodal_smoother == (void *) ML_Gen_Smoother_Jacobi) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    nodal_omega = 1.0;
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    ML_Smoother_Arglist_Set(nodal_args, 1, &nodal_omega);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_Jacobi) {
    edge_args = ML_Smoother_Arglist_Create(2);
    edge_omega = 1.0;
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    ML_Smoother_Arglist_Set(edge_args, 1, &edge_omega);
  }


  /***************************************************************
  * Set up smoothers for all levels. For the MLS polynomial pick 
  * parameters based on the coarsening rate. See paper 'Parallel
  * Multigrid Smoothing: Polynomial Versus Gauss-Seidel' by Adams,
  * Brezina, Hu, Tuminaro
  ****************************************************************/
  coarsest_level = MaxMgLevels - Nlevels;
  for (level = MaxMgLevels-1; level >= coarsest_level; level--)
  {
    if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
      if (level != coarsest_level) {
        Ncoarse_edge = Tmat_array[level-1]->outvec_leng;
        ML_gsum_scalar_int(&Ncoarse_edge, &itmp, ml_ERF->comm);
        edge_coarsening_rate =  2.*((double) Nfine_edge)/
                                   ((double) Ncoarse_edge);
      }
      else edge_coarsening_rate =  (double) Nfine_edge;
      if (ml_edges->comm->ML_mypid == 0) printf("hardwiring coarsening rate\n");
      edge_coarsening_rate = 27.;
      if (level == coarsest_level) edge_coarsening_rate = 300.;
      ML_Smoother_Arglist_Set(edge_args, 1, &edge_coarsening_rate);
      Nfine_edge = Ncoarse_edge;
    }
    if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
      if (level != coarsest_level) {
        Ncoarse_node = Tmat_array[level-1]->invec_leng;
        ML_gsum_scalar_int(&Ncoarse_node, &itmp, ml_ERF->comm);
        node_coarsening_rate =  2.*((double) Nfine_node)/ 
                                   ((double) Ncoarse_node);
      }
      else node_coarsening_rate = (double) Nfine_node;
      if (ml_edges->comm->ML_mypid == 0) printf("hardwiring coarsening rate\n");
      node_coarsening_rate = 27.;
      if (level == coarsest_level) node_coarsening_rate = 20.1825;

      ML_Smoother_Arglist_Set(nodal_args, 1, &node_coarsening_rate);
      Nfine_node = Ncoarse_node;
    }
    ML_Gen_Smoother_BlockHiptmair(ml_ERF, level, ML_PRESMOOTHER,
                             Nits_per_presmooth,
                             Tmat_array, Tmat_trans_array, NULL, 
                             edge_smoother, edge_args, nodal_smoother,
                             nodal_args, hiptmair_type);
  }
/*
  if (ml_ERF->comm->ML_mypid == 0)
    printf("\n\nSetting up coarse level separately\n\n");
  level = coarsest_level;
  Nits_per_presmooth = 6;
  ML_Gen_Smoother_BlockHiptmair(ml_ERF, level, ML_PRESMOOTHER,
                             Nits_per_presmooth,
                             Tmat_array, Tmat_trans_array, NULL, 
                             edge_smoother, edge_args, nodal_smoother,
                             nodal_args, hiptmair_type);
*/

  /*ML_Gen_CoarseSolverSuperLU(ml_ERF, coarsest_level);*/

  /* Must be called before invoking the preconditioner */

  ML_Gen_Solver(ml_ERF, ML_MGV, MaxMgLevels-1, coarsest_level); 


   /* initialize AZTEC options */

   AZ_defaults(options, params);
   options[AZ_precond] = AZ_none;
   options[AZ_solver] = AZ_tfqmr;
   options[AZ_kspace] = 400;
   options[AZ_max_iter] = 1000;
   options[AZ_conv] = AZ_r0;
   options[AZ_output] = 1;
   /*options[AZ_conv] = AZ_rhs;*/
   params[AZ_tol] = 1.0e-8;

/*
   if (ml_edges->comm->ML_mypid == 0) {
       for (iii=0; iii<*Nedges;iii++) {
         if (ERF_y[iii] != 0.0 || ERF_y[iii+(*Nedges)] != 0.0)
           printf("ERF_y(%d) = %e , %e\n",iii,ERF_y[iii],ERF_y[iii+(*Nedges)]);
       }
   }
*/

   
  AZ_set_ML_preconditioner(&Prec, ERF_AZmat, ml_ERF, options); 
/*    dump_greg(&(ml_edges->Amat[MaxMgLevels-1]),
  	  &(ml_nodes->Amat[MaxMgLevels-1]), M_mat,
    Tmat, Nghost_node, Nghost, ERF_x, ERF_y); */

  /* calculate an initial residual */
  ML_Operator_Apply(global_dude, global_dude->invec_leng, ERF_x,
		    global_dude->outvec_leng, vec1);
  for (i = 0; i < global_dude->outvec_leng; i++)
    vec1[i] = ERF_y[i] - vec1[i];

  fflush(stdout);
#ifdef ML_NEWNORM
  /* don't scale before applying T-norm */
  new_norm(Prec, ERF_x, &xnorm);
  /*if (ml_edges->comm->ML_mypid == 0) printf("\n\nx :  %e\n",xnorm);*/
  new_norm(Prec, vec1, &result);
  /*if (ml_edges->comm->ML_mypid == 0) printf("res :  %e\n",result);*/
  new_norm(Prec, ERF_y, &result1);
  /*if (ml_edges->comm->ML_mypid == 0) printf("rhs :  %e\n",result1);*/
  if (ml_edges->comm->ML_mypid == 0) printf("\n(|r|_t,|b|_t,x,rr) = %e %e %e %e\n\n",result,result1,xnorm,result/result1);
#endif

  xnorm = sqrt(ML_gdot(2 *(*Nedges), ERF_x, ERF_x,ml_edges->comm));
  result = sqrt(ML_gdot(2 * (*Nedges), vec1, vec1, ml_edges->comm));
  result1 = sqrt(ML_gdot(2 *(*Nedges), ERF_y, ERF_y,ml_edges->comm));
  if (ml_edges->comm->ML_mypid == 0) printf("(|r|_2,|b|_2,x,rr) = %e %e %e %e\n\n",result,result1,xnorm,result/result1);
  fflush(stdout);

#ifdef ML_SCALEIT
  /* to greg's scaling */
  for (i = 0; i < *Nedges; i++) {
    ERF_x[i         ] /= grad_scale[i];
    ERF_x[i+ *Nedges] /= grad_scale[i];
    ERF_y[i         ] *= grad_scale[i];
    ERF_y[i+ *Nedges] *= grad_scale[i];
    vec1[i         ] *= grad_scale[i];
    vec1[i+ *Nedges] *= grad_scale[i];
  }
#endif

  xnorm = sqrt(ML_gdot(2 *(*Nedges), ERF_x, ERF_x,ml_edges->comm));
  result = sqrt(ML_gdot(2 * (*Nedges), vec1, vec1, ml_edges->comm));
  result1 = sqrt(ML_gdot(2 *(*Nedges), ERF_y, ERF_y,ml_edges->comm));
  if (ml_edges->comm->ML_mypid == 0) printf("(|r|_2,|b|_2,x,rr) = %e %e %e %e\n",result,result1,xnorm,result/result1);
  fflush(stdout);

#ifdef ML_SCALEIT
  /* from greg's scaling */
  for (i = 0; i < *Nedges; i++) {
    ERF_x[i         ] *= grad_scale[i];
    ERF_x[i+ *Nedges] *= grad_scale[i];
    ERF_y[i         ] /= grad_scale[i];
    ERF_y[i+ *Nedges] /= grad_scale[i];
  }
#endif

  /* solve */
  /*time_it_(&(ml_edges->comm->ML_mypid));*/
  AZ_iterate(ERF_x, ERF_y, options, params, status,proc_config, 
  	         ERF_AZmat, Prec, NULL);
  /*time_it_(&(ml_edges->comm->ML_mypid));*/

  /* calculate a final residual */
  ML_Operator_Apply(global_dude, global_dude->invec_leng, ERF_x,
		    global_dude->outvec_leng, vec1);
  for (i = 0; i < global_dude->outvec_leng; i++)
    vec1[i] = ERF_y[i] - vec1[i];

#ifdef ML_NEWNORM
  /* don't scale before applying T-norm */
  if (ml_edges->comm->ML_mypid == 0) printf("\n\nx :");
  new_norm(Prec, ERF_x, &xnorm);
  if (ml_edges->comm->ML_mypid == 0) printf("res :");
  new_norm(Prec, vec1, &result);
  if (ml_edges->comm->ML_mypid == 0) printf("rhs :");
  new_norm(Prec, ERF_y, &result1);
  if (ml_edges->comm->ML_mypid == 0) printf("\n(|r|_t,|b|_t,x,rr) = %e %e %e %e\n\n",result,result1,xnorm,result/result1);
#endif

  xnorm = sqrt(ML_gdot(2 *(*Nedges), ERF_x, ERF_x,ml_edges->comm));
  result = sqrt(ML_gdot(2 * (*Nedges), vec1, vec1, ml_edges->comm));
  result1 = sqrt(ML_gdot(2 *(*Nedges), ERF_y, ERF_y,ml_edges->comm));
  if (ml_edges->comm->ML_mypid == 0) printf("(|r|_2,|b|_2,x,rr) = %e %e %e %e\n\n",result,result1,xnorm,result/result1);
  fflush(stdout);

#ifdef ML_SCALEIT
  /* to greg's scaling */
  for (i = 0; i < *Nedges; i++) {
    ERF_x[i         ] /= grad_scale[i];
    ERF_x[i+ *Nedges] /= grad_scale[i];
    ERF_y[i         ] *= grad_scale[i];
    ERF_y[i+ *Nedges] *= grad_scale[i];
    vec1[i         ] *= grad_scale[i];
    vec1[i+ *Nedges] *= grad_scale[i];
  }
#endif

  i = 2* (*Nedges);
  xnorm = sqrt(ML_gdot(i, ERF_x, ERF_x,ml_edges->comm));
  result = sqrt(ML_gdot(i, vec1, vec1, ml_edges->comm));
  result1 = sqrt(ML_gdot(i, ERF_y, ERF_y,ml_edges->comm));
  if (ml_edges->comm->ML_mypid == 0) printf("(|r|_2,|b|_2,x,rr) = %e %e %e %e\n\n",result,result1,xnorm,result/result1);
  fflush(stdout); fflush(stderr);

#ifdef ML_SCALEIT
  /* from greg's scaling */
  for (i = 0; i < *Nedges; i++) {
    ERF_x[i         ] *= grad_scale[i];
    ERF_x[i+ *Nedges] *= grad_scale[i];
    ERF_y[i         ] /= grad_scale[i];
    ERF_y[i+ *Nedges] /= grad_scale[i];
  }
#endif

  fflush(stdout);

  for (i = 0; i < AZ_MAT_blockmat_data->N; i++) {
    complex_x[2*i  ] = ERF_x[i          ]/grad_scale[i];
    complex_x[2*i+1] = ERF_x[i+(*Nedges)]/grad_scale[i];
  }

#ifdef out
  sprintf(str,"mlfinal%d",*proc);
  fp = fopen(str,"w");
  for (i = 0 ; i < Nlocal_edges; i++) {
    fprintf(fp,"%d %20.13e\n",i,ERF_x[i]/grad_scale[i]);
    fprintf(fp,"%d %20.13e\n",i,ERF_x[i + *Nedges]/grad_scale[i]);
  }
  fflush(fp);  fclose(fp);
  printf("finished\n");
#endif


  ML_Operator_Destroy(&Tmat_trans);
  ML_Aggregate_Destroy(&ag);
  ML_Operator_Destroy(&Tmat);
  ML_Operator_Destroy(&M_mat);
  AZ_matrix_destroy(&Ke);
  AZ_matrix_destroy(&M);
  ML_Destroy(&ml_nodes);
  ML_Destroy(&ml_edges);
  ML_MGHierarchy_ReitzingerDestroy(MaxMgLevels-2,&Tmat_array,
				  &Tmat_trans_array);
  ML_free(edge_args);
  ML_free(nodal_args);
  ML_Destroy2(&ml_ERF);
  AZ_matrix_destroy(&ERF_AZmat);
  AZ_precond_destroy(&Prec);
  ML_free(wrap_real_greg_data);
  ML_free(wrap_imag_greg_data);
  ML_free(AZ_MAT_blockmat_data);
  ML_free(ERF_x);
  ML_free(ERF_y);
  ML_free(wrap_grad_greg_data);
  ML_free(grad_scale);

  return;
}

/******************************************************************************/
/******************************************************************************/
/* Aztec functions to wrap greg's matvec, communication, and getrow.          */
/******************************************************************************/
/******************************************************************************/
void aztec_matvec_wrap_of_greg(double *x, double *y, AZ_MATRIX *Amat,
			       int proc_config[])
{
  struct wrap_greg_data *wrap_greg_data;
  int inlen, Nghost,i;
  double *temp;

  wrap_greg_data = (struct wrap_greg_data *) AZ_get_matvec_data(Amat);
  Nghost         =  wrap_greg_data->Nghost;
  inlen = Amat->data_org[AZ_N_internal] + Amat->data_org[AZ_N_border];
  temp = (double *) malloc(sizeof(double)*(inlen+Nghost));
  for (i = 0; i < inlen; i++) temp[i] = x[i];


  ML_JANUS_FORT(mat_vec_real_)(wrap_greg_data->ca, y, x, wrap_greg_data->index,
		wrap_greg_data->mcol,  wrap_greg_data->nrows,
		wrap_greg_data->nodes, wrap_greg_data->nput,
		wrap_greg_data->nget,  wrap_greg_data->nsol,
		wrap_greg_data->proc);
  ML_free(temp);
}

/******************************************************************************/

int aztec_comm_wrap_of_greg( double *x, AZ_MATRIX *Amat)
{
  struct wrap_greg_data *wrap_greg_data;
  int k;

  wrap_greg_data = (struct wrap_greg_data *) AZ_get_matvec_data(Amat);
  k = Amat->data_org[AZ_N_internal] + Amat->data_org[AZ_N_border] 
                                    + wrap_greg_data->Nghost;

  ML_JANUS_FORT(com_dbl_real_)(x, wrap_greg_data->proc, wrap_greg_data->nodes,
	       wrap_greg_data->nput, wrap_greg_data->nget, 
	       wrap_greg_data->nsol, &k);
  return 1;
}

/******************************************************************************/

int  aztec_getrow_wrap_of_greg( int columns[], double values[],
        int row_lengths[], struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
				       int requested_rows[], int allocated_space)
{
  /*
   * Supply matrix diagonals for rows requested_rows[0 ... N_requested_rows-1].
   * Return this information in 'row_lengths, columns, values'.  If there is
   * not enough space to complete this operation, return 0. Otherwise, return 1.
   *
   * Parameters
   * ==========
   * Amat             On input, points to user's data containing matrix values.
   * N_requested_rows On input, number of rows for which nonzero are to be
   *                  returned.
   * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the
   *                  row indices of the rows for which nonzero values are
   *                  returned.
   * row_lengths      On output, row_lengths[i] is the number of nonzeros in the
   *                  row 'requested_rows[i]'
   * columns,values   On output, columns[k] and values[k] contains the column
   *                  number and value of a matrix nonzero where all nonzeros for
   *                  requested_rows[i] appear before requested_rows[i+1]'s
   *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
   * allocated_space  On input, indicates the space available in 'columns' and
   *                  'values' for storing nonzeros. If more space is needed,
   *                  return 0.
   */

  struct wrap_greg_data *wrap_greg_data;

  wrap_greg_data = (struct wrap_greg_data *) AZ_get_matvec_data(Amat);

  if (allocated_space < 13) return 0;
  
  ML_JANUS_FORT(greg_getrow_)(requested_rows, columns, values, row_lengths,
	       wrap_greg_data->ca,    wrap_greg_data->index,
	       wrap_greg_data->mcol);  

  return 1;
}

/******************************************************************************/

int ml_grad_matvec_wrap_of_greg(void *data, int inlen, double x[],
				int outlen, double y[]) {

  struct wrap_greg_data *wrap_greg_data;
  int Nghost, i;
  double *temp;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;
  wrap_greg_data = (struct wrap_greg_data *) ML_Get_MyMatvecData(mat_in);
  Nghost         =  wrap_greg_data->Nghost;

  temp = (double *) malloc(sizeof(double)*(inlen+Nghost));
  for (i = 0; i < inlen; i++) temp[i] = x[i];

  ML_JANUS_FORT(mat_vec_grad_)(wrap_greg_data->ca, y, temp, wrap_greg_data->index,
		wrap_greg_data->mcol,  wrap_greg_data->nrows,
		wrap_greg_data->ncols,
		wrap_greg_data->nodes, wrap_greg_data->nput,
		wrap_greg_data->nget,  wrap_greg_data->nsol,
                &Nghost, wrap_greg_data->proc);
  ML_free(temp);

  return 0;
}

/******************************************************************************/

int ml_grad_comm_wrap_of_greg( double *x, void *data) {
  struct wrap_greg_data *wrap_greg_data;
  int k;

  wrap_greg_data = (struct wrap_greg_data *) data;
  k = *(wrap_greg_data->ncols) + wrap_greg_data->Nghost;

  ML_JANUS_FORT(com_pot_real_)(x, wrap_greg_data->proc, wrap_greg_data->nodes,
	       wrap_greg_data->nput, wrap_greg_data->nget, 
	       wrap_greg_data->nsol, &k);

  return 0;
}
int ml_grad_getrow_wrap_of_greg(void *data, int N_requested,
			  int requested_rows[], int allocated,
			  int columns[], double values[], 
			  int row_lengths[])
{
  struct wrap_greg_data *wrap_greg_data;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;

  wrap_greg_data = (struct wrap_greg_data *) ML_Get_MyGetrowData(mat_in);

  if (allocated < 13) return 0;
  
  ML_JANUS_FORT(greg_getrow_)(requested_rows, columns, values, row_lengths,
	       wrap_greg_data->ca,    wrap_greg_data->index,
	       wrap_greg_data->mcol);  


  return 1;
}

struct ml_Tmat_wrap {
  ML_Operator  *Tmat_trans, **Tmat_array, **Tmat_trans_array;
};


int  ML_Gen_MGHierarchy_viaReitzinger(ML *ml_edges, ML* ml_nodes, 
                    int fine_level, int incr_or_decrease,
                    ML_Aggregate *ag, ML_Operator *Tmat,
                    int smooth_flag, double smooth_factor)
{


struct ml_Tmat_wrap *ml_Tmat_wrap;

  
  ml_Tmat_wrap = (struct ml_Tmat_wrap *) ML_allocate(sizeof(
                                             struct ml_Tmat_wrap));
  ml_edges->void_options = (void *) ml_Tmat_wrap;

  /********************************************************************/
  /*                      Set up Tmat_trans                           */
  /*------------------------------------------------------------------*/

  ml_Tmat_wrap->Tmat_trans = ML_Operator_Create(ml_edges->comm);
  printf("3a\n"); fflush(stdout);
  ML_Operator_Transpose_byrow(Tmat, ml_Tmat_wrap->Tmat_trans);
  printf("3b\n"); fflush(stdout);

  if (ml_nodes->Amat[fine_level].getrow->ML_id == ML_EMPTY) {
    ML_2matmult(ml_Tmat_wrap->Tmat_trans, Tmat, 
		&(ml_nodes->Amat[fine_level]),
                ML_CSR_MATRIX);
  }


  /********************************************************************/
  /*                      Build the MG hierarchy                      */
  /*------------------------------------------------------------------*/

  return(ML_Gen_MGHierarchy_UsingReitzinger(ml_edges, ml_nodes,fine_level,
					    incr_or_decrease,ag,Tmat,
					    ml_Tmat_wrap->Tmat_trans, 
                                            &(ml_Tmat_wrap->Tmat_array),
		    		            &(ml_Tmat_wrap->Tmat_trans_array), 
                                            smooth_flag, smooth_factor));
}


int ML_Gen_Smoother_wrapBlockHiptmair(ML *ml_edges, int nl, int pre_or_post, 
			 int ntimes, 
			 void *edge_smoother, void **edge_args,
			 void *nodal_smoother, void **nodal_args, int type,
			 int nodal_its, double nodal_omega, int edge_its,
				      double edge_omega, int fine_level,
				      int coarsest_level,
				      int Nits_per_presmooth)
{
  int Nfine_node, Nfine_edge, Ncoarse_node, Ncoarse_edge, level;
  double edge_coarsening_rate, node_coarsening_rate;
  int itmp;
  struct ml_Tmat_wrap *ml_Tmat_wrap;
  ML_Operator   **Tmat_array, **Tmat_trans_array;


  ml_Tmat_wrap = (struct ml_Tmat_wrap *) ml_edges->void_options;
  Tmat_array = ml_Tmat_wrap->Tmat_array;
  Tmat_trans_array = ml_Tmat_wrap->Tmat_trans_array;

  /* Set the Hiptmair subsmoothers */

  if (nodal_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    ML_Smoother_Arglist_Set(nodal_args, 1, &nodal_omega);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_SymGaussSeidel) {
    edge_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    ML_Smoother_Arglist_Set(edge_args, 1, &edge_omega);
  }
  if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
    nodal_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(nodal_args, 0, &nodal_its);
    Nfine_node = Tmat_array[fine_level]->invec_leng;
    ML_gsum_scalar_int(&Nfine_node, &itmp, ml_edges->comm);
  }
  if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
    edge_args = ML_Smoother_Arglist_Create(2);
    ML_Smoother_Arglist_Set(edge_args, 0, &edge_its);
    Nfine_edge = Tmat_array[fine_level]->outvec_leng;
    ML_gsum_scalar_int(&Nfine_edge, &itmp, ml_edges->comm);
  }

  /***************************************************************
  * Set up smoothers for all levels. For the MLS polynomial pick 
  * parameters based on the coarsening rate. See paper 'Parallel
  * Multigrid Smoothing: Polynomial Versus Gauss-Seidel' by Adams,
  * Brezina, Hu, Tuminaro
  ****************************************************************/

  for (level = fine_level; level >= coarsest_level; level--) {
    if (edge_smoother == (void *) ML_Gen_Smoother_MLS) {
      if (level != coarsest_level) {
        Ncoarse_edge = Tmat_array[level-1]->outvec_leng;
        ML_gsum_scalar_int(&Ncoarse_edge, &itmp, ml_edges->comm);
        edge_coarsening_rate =  2.*((double) Nfine_edge)/
                                   ((double) Ncoarse_edge);
      }
      else edge_coarsening_rate =  (double) Nfine_edge;

      ML_Smoother_Arglist_Set(edge_args, 1, &edge_coarsening_rate);
      Nfine_edge = Ncoarse_edge;
    }
    if (nodal_smoother == (void *) ML_Gen_Smoother_MLS) {
      if (level != coarsest_level) {
        Ncoarse_node = Tmat_array[level-1]->invec_leng;
        ML_gsum_scalar_int(&Ncoarse_node, &itmp, ml_edges->comm);
        node_coarsening_rate =  2.*((double) Nfine_node)/ 
                                   ((double) Ncoarse_node);
      }
      else node_coarsening_rate = (double) Nfine_node;

      ML_Smoother_Arglist_Set(nodal_args, 1, &node_coarsening_rate);
      Nfine_node = Ncoarse_node;
    }
    ML_Gen_Smoother_Hiptmair(ml_edges, level, ML_PRESMOOTHER, Nits_per_presmooth,
                             Tmat_array, Tmat_trans_array, NULL, 
                             edge_smoother, edge_args, nodal_smoother,
                             nodal_args, type);
  }
  return 1;
}


int dump_greg(ML_Operator *Ke, ML_Operator *Kn, ML_Operator *M_mat,
	      ML_Operator *T_mat, int Nghost_nodes, int Nghost_edges,
	      double *x, double *rhs)
{
  double *global_edges, *global_Knodes, *global_Tnodes, colVal[15];
  int    N_edges, N_nodes, node_offset, edge_offset;
  int colInd[15], i, j, ncnt;
  char str[80];
  FILE *fid;
  ML_Comm *comm;
  int Nedges_global;
  
  comm = Ke->comm;

  if (comm->ML_mypid == 0) printf("DUMPING GREG's MATRICES\n");

  N_nodes = Kn->outvec_leng;
  N_edges = Ke->outvec_leng;

  node_offset = ML_gpartialsum_int(N_nodes, comm);
  edge_offset = ML_gpartialsum_int(N_edges, comm);
  Nedges_global = N_edges;
  ML_gsum_scalar_int(&Nedges_global, &i, comm);

  global_Knodes =(double *) ML_allocate(sizeof(double)*(N_nodes+Nghost_nodes));
  global_Tnodes =(double *) ML_allocate(sizeof(double)*(N_nodes+Nghost_nodes));
  global_edges  =(double *) ML_allocate(sizeof(double)*(N_edges+Nghost_edges));

  for (i = 0 ; i < N_nodes; i++) global_Knodes[i] = (double) (node_offset + i);
  for (i = 0 ; i < N_nodes; i++) global_Tnodes[i] = (double) (node_offset + i);
  for (i = 0 ; i < N_edges; i++) global_edges[i] = (double) (edge_offset + i);

  for (i = 0 ; i < Nghost_nodes; i++) global_Knodes[i+N_nodes] = -1;
  for (i = 0 ; i < Nghost_nodes; i++) global_Tnodes[i+N_nodes] = -1;
  for (i = 0 ; i < Nghost_edges; i++) global_edges[i+N_edges] = -1;

  ML_exchange_bdry(global_Tnodes,T_mat->getrow->pre_comm, 
 		 Kn->invec_leng,comm,ML_OVERWRITE,NULL);
  ML_exchange_bdry(global_Knodes,Kn->getrow->pre_comm, 
 		 Kn->invec_leng,comm,ML_OVERWRITE,NULL);
  ML_exchange_bdry(global_edges,Ke->getrow->pre_comm, 
 		 Ke->invec_leng,comm,ML_OVERWRITE,NULL);

  /* spit out Kn */

  sprintf(str,"Kn.%d",comm->ML_mypid);
  fid = fopen(str,"w");
  for (i = 0; i < Kn->outvec_leng; i++) {
    j = ML_Operator_Getrow(Kn,1,&i,15,colInd,colVal,&ncnt);
    for (j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	fprintf(fid,"%5d %5d %20.13e\n",(int) global_Knodes[i]+1,
		       (int) global_Knodes[colInd[j]]+1, colVal[j]);
      }
    }
  }
  fclose(fid);

  /* spit out Ke  */

  sprintf(str,"Ke.%d",comm->ML_mypid);
  fid = fopen(str,"w");
  for (i = 0; i < Ke->outvec_leng; i++) {
    j = ML_Operator_Getrow(Ke,1,&i,15,colInd,colVal,&ncnt);
    for (j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	fprintf(fid,"%5d %5d %20.13e\n",(int) global_edges[i]+1,
		       (int) global_edges[colInd[j]]+1, colVal[j]);
      }
    }
  }
  fclose(fid);

  /* spit out M  */

  sprintf(str,"M.%d",comm->ML_mypid);
  fid = fopen(str,"w");
  for (i = 0; i < M_mat->outvec_leng; i++) {
    j = ML_Operator_Getrow(M_mat,1,&i,15,colInd,colVal,&ncnt);
    for (j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	fprintf(fid,"%5d %5d %20.13e\n",(int) global_edges[i]+1,
		       (int) global_edges[colInd[j]]+1, colVal[j]);
      }
    }
  }
  fclose(fid);

  /* spit out T */

  sprintf(str,"T.%d",comm->ML_mypid);
  fid = fopen(str,"w");
  for (i = 0; i < T_mat->outvec_leng; i++) {
    j = ML_Operator_Getrow(T_mat,1,&i,15,colInd,colVal,&ncnt);
    for (j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	fprintf(fid,"%5d %5d %20.13e\n",(int) global_edges[i]+1,
		       (int) global_Tnodes[colInd[j]]+1, colVal[j]);
      }
    }
  }
  fclose(fid);

  /* spit out x */

  sprintf(str,"xxx.%d",comm->ML_mypid);
  fid = fopen(str,"w");
  for (i = 0; i < T_mat->outvec_leng; i++) {
    fprintf(fid,"%5d %20.13e\n",(int) global_edges[i]+1,x[i]);
    fprintf(fid,"%5d %20.13e\n",(int) global_edges[i]+1 + Nedges_global,
	    x[i+T_mat->outvec_leng]);
  }
  fclose(fid);

  /* spit out rhs */

  sprintf(str,"rhs.%d",comm->ML_mypid);
  fid = fopen(str,"w");
  for (i = 0; i < T_mat->outvec_leng; i++) {
    fprintf(fid,"%5d %20.13e\n",(int) global_edges[i]+1,rhs[i]);
    fprintf(fid,"%5d %20.13e\n",(int) global_edges[i]+1 + Nedges_global,
	    rhs[i+T_mat->outvec_leng]);
  }
  fclose(fid);


  ML_free(global_Knodes);
  ML_free(global_Tnodes);
  ML_free(global_edges);
  return 0;
}
/* Use this function to set up the global arrays that are needed */
/* to dump out aggregate information in the MIS code.            */
int setup_agg_MIS_dump(ML_Operator *Kn, int Nghost_nodes)
{
  double *dupdate;
  int     N_nodes, node_offset;
  int     i;
  

  N_nodes = Kn->outvec_leng;

  node_offset = ML_gpartialsum_int(N_nodes, Kn->comm);

  dupdate =(double *) ML_allocate(sizeof(double)*(N_nodes+Nghost_nodes));

  for (i = 0 ; i < N_nodes; i++) dupdate[i] = (double) (node_offset + i);

  for (i = 0 ; i < Nghost_nodes; i++) dupdate[i+N_nodes] = -1.;

  ML_exchange_bdry(dupdate,Kn->getrow->pre_comm, 
 		 Kn->invec_leng, Kn->comm, ML_OVERWRITE, NULL);

  update =(int *) ML_allocate(sizeof(int)*(N_nodes));
  update_index =(int *) ML_allocate(sizeof(int)*(N_nodes));
  for (i = 0 ; i < N_nodes; i++) {
     update[i]       = (int) dupdate[i];
     update_index[i] = i;
  }
  external =(int *) ML_allocate(sizeof(int)*(Nghost_nodes + 1));
  extern_index =(int *) ML_allocate(sizeof(int)*(Nghost_nodes + 1));
  for (i = 0 ; i < Nghost_nodes; i++) {
    extern_index[i] = i + N_nodes;
    external[i] =  (int) dupdate[i + N_nodes];
  }

  ML_free(dupdate);

  return 0;
}



void dumbo_(double *vec, int *ntotal)
{
  ML *ml;

  printf("in dumbo %d\n",*ntotal);
  ML_Create(&ml,1);

  ML_random_vec(vec, *ntotal, ml->comm);
  ML_random_vec(&(vec[*ntotal]), *ntotal, ml->comm);
  ML_Destroy(&ml);
}


void ML_JANUS_FORT(ray_matvec_)(/* ML_Operator *Amat,*/double complex_x[],double complex_Ax[])
{
  double *ERF_x, *ERF_y;
  int i, Nsize;
  ML_Operator *Amat;
  struct ML_Operator_blockmat_data *block_data;
  double *Ke_diag, *M_diag, denom;

  Amat = global_dude;
  Nsize = Amat->invec_leng/2;

  block_data = (struct ML_Operator_blockmat_data *) Amat->data;


  ML_Operator_Get_Diag(block_data->M_mat, Nsize, &M_diag);
  if (block_data->M_diag == NULL) {
    block_data->M_diag = ML_allocate(sizeof(double) * Nsize);
    for (i=0; i<Nsize; i++) block_data->M_diag[i] = M_diag[i];
    /*block_data->M_diag = M_diag;*/
  }

  ML_Operator_Get_Diag(block_data->Ke_mat, Nsize, &Ke_diag);
  if (block_data->Ke_diag == NULL) {
    block_data->Ke_diag = ML_allocate(sizeof(double) * Nsize);
    for (i=0; i<Nsize; i++) block_data->Ke_diag[i] = Ke_diag[i];
    /*block_data->Ke_diag = Ke_diag;*/
  }


  ERF_x = (double *) malloc(sizeof(double)*Nsize*2);
  ERF_y = (double *) malloc(sizeof(double)*Nsize*2);

  for (i = 0; i < Nsize; i++) {
    ERF_x[i         ] = complex_x[2*i    ];
    ERF_x[i+ Nsize  ] = complex_x[2*i+1  ];
  }
  ML_Operator_Apply(Amat, Amat->invec_leng, ERF_x,
		    Amat->outvec_leng, ERF_y);


  for (i = 0; i < Nsize; i++) {
    /* Jonathan ... put in a test to check that diagonal is not zero */
    denom = 1./(M_diag[i]*M_diag[i] + Ke_diag[i]*Ke_diag[i]);
    complex_Ax[2*i  ] = (Ke_diag[i]*ERF_y[i] + M_diag[i]*ERF_y[i+Nsize])*denom;
    complex_Ax[2*i+1] = (Ke_diag[i]*ERF_y[i+Nsize] - M_diag[i]*ERF_y[i])*denom;
  }

  ML_free(ERF_x);
  ML_free(ERF_y);
}
