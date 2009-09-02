/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to set up AMG multigrid structure                               */
/* ************************************************************************* */
/* ML_AMG_Create                                                             */
/* ML_AMG_Destroy                                                            */
/* ML_AMG_Set_OutputLevel                                                    */
/* ML_AMG_Set_MaxCoarseSize                                                  */
/* ML_AMG_Set_CoarsenScheme_MIS                                              */
/* ML_AMG_Set_Threshold                                                      */
/* ML_AMG_Set_MaxLevels                                                      */
/* ML_AMG_Set_CurrentLevel                                                   */
/* ML_AMG_Set_StartLevel                                                     */
/* ML_AMG_Set_Coarsen                                                        */
/* ML_AMG_Set_Smoother                                                       */
/* ML_AMG_Set_SmootherAztec                                                  */
/* ML_AMG_Set_CoarseSolve                                                    */
/* ML_AMG_Print                                                              */
/* ML_AMG_Print_Complexity                                                   */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : October, 2000                                             */
/* ************************************************************************* */
/* ************************************************************************* */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ml_amg.h"
#include "ml_operator.h"

/* ************************************************************************* */
/* Constructor                                                               */
/* ------------------------------------------------------------------------- */

int ML_AMG_Create( ML_AMG **amg )
{
   ML_memory_alloc( (void **) amg, sizeof(ML_AMG), "MG1" );
   (*amg)->ML_id                      = ML_ID_AMG;
   (*amg)->print_flag                 = 1;
   (*amg)->max_coarse_size            = 10;
   (*amg)->threshold                  = 0.5;
   (*amg)->curr_threshold             = 0.2;
   (*amg)->amg_scheme                 = ML_AMG_SCALAR;
   (*amg)->coarsen_scheme             = ML_AMG_MIS;
   (*amg)->num_PDE_eqns               = 1;
   (*amg)->blk_info                   = NULL;
   (*amg)->max_levels                 = 25;
   (*amg)->begin_level                = 25;
   (*amg)->cur_level                  = 25;
   (*amg)->operator_complexity        = 0.0;
   (*amg)->fine_complexity            = 0.0;
   (*amg)->presmoother_type           = ML_AMG_SM_SGS;
   (*amg)->postsmoother_type          = ML_AMG_SM_SGS;
   (*amg)->presmoother_ntimes         = 2;
   (*amg)->postsmoother_ntimes        = 2;
   (*amg)->presmoother_jacobiwt       = 0.0;
   (*amg)->postsmoother_jacobiwt      = 0.0;
   (*amg)->coarse_solver_type         = ML_AMG_SM_SUPERLU;
   (*amg)->coarse_solver_ntimes       = 1;
   (*amg)->coarse_solver_jacobiwt     = 0.0;
   (*amg)->fine_Amat                  = NULL;
   (*amg)->pre_aztec_options          = NULL;
   (*amg)->pre_aztec_params           = NULL;
   (*amg)->pre_aztec_proc_config      = NULL;
   (*amg)->pre_aztec_status           = NULL;
   (*amg)->pre_function               = NULL;
   (*amg)->post_aztec_options         = NULL;
   (*amg)->post_aztec_params          = NULL;
   (*amg)->post_aztec_proc_config     = NULL;
   (*amg)->post_aztec_status          = NULL;
   (*amg)->post_function              = NULL;
   return 0;
}

/* ************************************************************************* */
/* destructor                                                                */
/* ------------------------------------------------------------------------- */

int ML_AMG_Destroy( ML_AMG **amg )
{
   if ( (*amg)->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Destroy : wrong object. \n");
      exit(-1);
   }
   if ( (*amg)->blk_info != NULL ) ML_memory_free((void **) (*amg)->blk_info);
   ML_memory_free( (void **) amg );
   (*amg) = NULL;
   return 0;
}

/* ************************************************************************* */
/* set output level                                                     */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_OutputLevel( ML_AMG *amg, double level )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_OutputLevel : wrong object. \n");
      exit(-1);
   }
   amg->print_flag = level;
   return 0;
}

/* ************************************************************************* */
/* set maximum coarsest grid size                                            */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_MaxCoarseSize( ML_AMG *amg, int size  )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_MaxCoarseSize : wrong object. \n");
      exit(-1);
   }
   amg->max_coarse_size = size;
   return 0;
}

/* ************************************************************************* */
/* select AMG scheme (scalar, system (unknown), system (nodal)          */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_AMGScheme_Scalar( ML_AMG *amg  )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_AMGScheme_Scalar : wrong object. \n");
      exit(-1);
   }
   amg->amg_scheme = ML_AMG_SCALAR;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_AMG_Set_AMGScheme_SystemUnknown( ML_AMG *amg, int numPDE  )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_AMGScheme_SystemUnknown : wrong object. \n");
      exit(-1);
   }
   amg->amg_scheme   = ML_AMG_SYSTEM_UNKNOWN;
   if ( numPDE > 0 ) amg->num_PDE_eqns = numPDE;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_AMG_Set_AMGScheme_SystemNodal( ML_AMG *amg, int numPDE  )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_AMGScheme_SystemNodal : wrong object. \n");
      exit(-1);
   }
   amg->amg_scheme = ML_AMG_SYSTEM_NODAL;
   if ( numPDE > 0 ) amg->num_PDE_eqns = numPDE;
   return 0;
}

/* ************************************************************************* */
/* select coarsening scheme                                                  */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_CoarsenScheme_MIS( ML_AMG *amg  )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_CoarsenScheme_MIS : wrong object. \n");
      exit(-1);
   }
   amg->coarsen_scheme = ML_AMG_MIS;
   return 0;
}

/* ************************************************************************* */
/* set/reset strength of connection threshold                                */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_Threshold( ML_AMG *amg, double thresh )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_Threshold : wrong object. \n");
      exit(-1);
   }
   if ( thresh >= 0.0 && thresh <= 1.0 ) amg->threshold = thresh;
   else                                  amg->threshold = 0.0;
   return 0;
}

/* ************************************************************************* */
/* set max number of levels and other level information                 */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_MaxLevels( ML_AMG *amg, int level )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_MaxLevels : wrong object. \n");
      exit(-1);
   }
   amg->max_levels = level;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_AMG_Set_CurrentLevel( ML_AMG *amg, int level )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_CurrentLevel : wrong object. \n");
      exit(-1);
   }
   amg->cur_level = level;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_AMG_Set_StartLevel( ML_AMG *amg, int level )
{
   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Set_StartLevel : wrong object. \n");
      exit(-1);
   }
   amg->begin_level = level;
   return 0;
}

/* ************************************************************************* */
/* Coarsening routine                                                        */
/* ------------------------------------------------------------------------- */

int ML_AMG_Coarsen(ML_AMG *amg,ML_Operator *Amatrix,ML_Operator **Pmatrix, 
                   ML_Comm *comm)
{
   int    Ncoarse, mypid;
#ifdef ML_TIMING
   double t0;

   t0 = GetClock();
#endif

   if ( amg->ML_id != ML_ID_AMG ) 
   {
      printf("ML_AMG_Coarsen : wrong object. \n");
      exit(-1);
   }

   mypid = comm->ML_mypid;
   if (mypid == 0 && amg->print_flag < ML_Get_PrintLevel())
      printf("ML_AMG_Coarsen begins ...\n");

   Amatrix->num_PDEs = amg->num_PDE_eqns;
   switch ( amg->coarsen_scheme )
   {
      case ML_AMG_MIS :
         Ncoarse = ML_AMG_CoarsenMIS(amg,Amatrix,Pmatrix,comm);
         break;
      default :
         if ( mypid == 0 ) printf("ML_AMG_Coarsen : invalid scheme.\n");
         exit(1);
   } 

#ifdef ML_AMG_DEBUG
   i = 0;
   i = ML_gmax_int(i, comm);
   if ( mypid == 0 && amg->print_flag < ML_Get_PrintLevel()) printf("ML_AMG_Coarsen ends.\n");
#endif
#ifdef ML_TIMING
   t0 = GetClock() - t0;
   t0 = ML_gsum_double(t0, comm);
   t0 = t0/((double) comm->ML_nprocs);
   if ( mypid == 0 && ML_Get_PrintLevel() > 10)
     printf(" AMG setup time \t= %e\n",t0);
#endif

   return Ncoarse;
}

/* ************************************************************************* */
/* set smoother and smoother parameters                                      */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_Smoother(ML_AMG *amg,int smoother_type, int pre_or_post,
                        ML_Operator *Amatrix,int ntimes, double weight)
{
   amg->fine_Amat = Amatrix;
   if ( pre_or_post == ML_PRESMOOTHER )
   {
      if ( smoother_type >= ML_AMG_SM_JACOBI && 
           smoother_type <= ML_AMG_SM_MSCHWARZ)
         amg->presmoother_type = smoother_type;
      amg->presmoother_ntimes = ntimes;
      amg->presmoother_jacobiwt = weight;
   }
   else if ( pre_or_post == ML_POSTSMOOTHER )
   {
      if ( smoother_type >= ML_AMG_SM_JACOBI && 
           smoother_type <= ML_AMG_SM_MSCHWARZ)
         amg->postsmoother_type = smoother_type;
      amg->postsmoother_ntimes = ntimes;
      amg->postsmoother_jacobiwt = weight;
   }
   return 0;
}
  
/* ************************************************************************* */
/* set Aztec-related smoother parameters                                     */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_SmootherAztec(ML_AMG *amg, int pre_or_post, int *options, 
			     double *params, int *proc_config, double *status, 
			     void (*aztec_function)(int)
		     /* Trying to avoid including az_aztec.h. 
			void (*aztec_function)(double *,int *,int *,double *,
                        struct AZ_MATRIX_STRUCT*,struct AZ_PREC_STRUCT*)
		     */)
{
   if ( pre_or_post == ML_PRESMOOTHER )
   {
      amg->pre_aztec_options     = options;
      amg->pre_aztec_params      = params;
      amg->pre_aztec_proc_config = proc_config;
      amg->pre_aztec_status      = status;
      amg->pre_function          = aztec_function;
   }
   else if ( pre_or_post == ML_POSTSMOOTHER )
   {
      amg->post_aztec_options     = options;
      amg->post_aztec_params      = params;
      amg->post_aztec_proc_config = proc_config;
      amg->post_aztec_status      = status;
      amg->post_function          = aztec_function;
   }
   return 0;
}

/* ************************************************************************* */
/* set coarse grid solver                                                    */
/* ------------------------------------------------------------------------- */

int ML_AMG_Set_CoarseSolve(ML_AMG *amg, int solve_type, int ntimes, 
                           double weight)
{
   if ( solve_type >= ML_AMG_SM_JACOBI && solve_type <= ML_AMG_SM_SUPERLU)
      amg->coarse_solver_type = solve_type;
   amg->coarse_solver_ntimes = ntimes;
   if ( solve_type == ML_AMG_SM_SUPERLU) amg->coarse_solver_ntimes = 1;
   amg->coarse_solver_jacobiwt = weight;
   return 0;
}

/* ************************************************************************* */
/* Print information about current state of ML_AMG                           */
/* ------------------------------------------------------------------------- */

int ML_AMG_Print( ML_AMG *amg )
{
   char string[100];

   printf("**************************************************************\n");
   printf("* ML AMG information                                         *\n");
   printf("==============================================================\n");
   switch (amg->coarsen_scheme)
   {
      case ML_AMG_MIS :
           printf("ML_AMG : coarsen scheme     = MIS\n");
   }
   printf("ML_AMG : strong threshold   = %e\n", amg->threshold);
   printf("ML_AMG : number of PDEs     = %d\n", amg->num_PDE_eqns);
   printf("ML_AMG : max coarse size    = %d\n", amg->max_coarse_size);
   printf("ML_AMG : coarsen scheme     = MIS\n");
   printf("ML_AMG : max no. of levels  = %d\n", amg->max_levels);
   switch (amg->presmoother_type)
   {
      case ML_AMG_SM_JACOBI   : strcpy( string, "Jacobi" ); break;
      case ML_AMG_SM_GS       : strcpy( string, "Gauss Seidel" ); break;
      case ML_AMG_SM_SGS      : strcpy( string, "symmetric Gauss Seidel" ); break;
      case ML_AMG_SM_VBJACOBI : strcpy( string, "VBlock Jacobi" ); break;
      case ML_AMG_SM_VBGS     : strcpy( string, "VBlock Gauss Seidel" ); break;
      case ML_AMG_SM_VBSGS    : strcpy( string, "VBlock symmetric GS" ); break;
      case ML_AMG_SM_ASCHWARZ : strcpy( string, "additive Schwarz" ); break;
      case ML_AMG_SM_MSCHWARZ : strcpy( string, "multiplicative Schwarz" ); break;
   }
   printf("ML_AMG : presmoother type   = %s\n", string);
   printf("ML_AMG : presmoother ntimes = %d\n", amg->presmoother_ntimes);
   if ( amg->presmoother_type == ML_AMG_SM_JACOBI )
      printf("ML_AMG : damping factor     = %e\n", amg->presmoother_jacobiwt);
   switch (amg->postsmoother_type)
   {
      case ML_AMG_SM_JACOBI   : strcpy( string, "Jacobi" ); break;
      case ML_AMG_SM_GS       : strcpy( string, "Gauss Seidel" ); break;
      case ML_AMG_SM_SGS      : strcpy( string, "symmetric Gauss Seidel" ); break;
      case ML_AMG_SM_VBJACOBI : strcpy( string, "VBlock Jacobi" ); break;
      case ML_AMG_SM_VBGS     : strcpy( string, "VBlock Gauss Seidel" ); break;
      case ML_AMG_SM_VBSGS    : strcpy( string, "VBlock symmetric GS" ); break;
      case ML_AMG_SM_ASCHWARZ : strcpy( string, "additive Schwarz" ); break;
      case ML_AMG_SM_MSCHWARZ : strcpy( string, "multiplicative Schwarz" ); break;
   }
   printf("ML_AMG : postsmoother type  = %s\n", string);
   printf("ML_AMG : postsmoother ntimes= %d\n", amg->postsmoother_ntimes);
   if ( amg->postsmoother_type == ML_AMG_SM_JACOBI )
      printf("ML_AMG : damping factor     = %e\n", amg->postsmoother_jacobiwt);
   switch (amg->coarse_solver_type)
   {
      case ML_AMG_SM_JACOBI   : strcpy( string, "Jacobi" ); break;
      case ML_AMG_SM_GS       : strcpy( string, "Gauss Seidel" ); break;
      case ML_AMG_SM_SGS      : strcpy( string, "symmetric Gauss Seidel" ); break;
      case ML_AMG_SM_VBJACOBI : strcpy( string, "VBlock Jacobi" ); break;
      case ML_AMG_SM_VBGS     : strcpy( string, "VBlock Gauss Seidel" ); break;
      case ML_AMG_SM_VBSGS    : strcpy( string, "VBlock symmetric GS" ); break;
      case ML_AMG_SM_ASCHWARZ : strcpy( string, "additive Schwarz" ); break;
      case ML_AMG_SM_MSCHWARZ : strcpy( string, "multiplicative Schwarz" ); break;
      case ML_AMG_SM_SUPERLU  : strcpy( string, "SuperLU"); break;
   }
   printf("ML_AMG : coarse solver      = %s\n", string);

   printf("**************************************************************\n");
   return 1;
}

/* ************************************************************************* */
/* Print information about operator complexity                               */
/* ------------------------------------------------------------------------- */

int ML_AMG_Print_Complexity( ML_AMG *amg )
{
   if ( amg->fine_complexity != 0.0 )
      printf("AMG : operator complexity = %e.\n",
              amg->operator_complexity / amg->fine_complexity);
   else
      printf("AMG error :  fine complexity = 0.0.\n");
   return 0;
}

