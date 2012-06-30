/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Function for manipulating the ML_BdryPts structure                   */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_bdrypts.h"

/* ******************************************************************** */
/* Constructor                                                          */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Create(ML_BdryPts **bc)
{
   ML_BdryPts *ml_bc;

   ML_memory_alloc( (void**) bc, sizeof(ML_BdryPts), "BC1" );
   ml_bc = (*bc);
   ml_bc->ML_id = ML_ID_BC;
   ml_bc->Dirichlet_grid_CreateOrDup = 0;
   ml_bc->Dirichlet_grid_length = 0;
   ml_bc->Dirichlet_grid_list = NULL;
   ml_bc->Dirichlet_eqn_CreateOrDup = 0;
   ml_bc->Dirichlet_eqn_length = 0;
   ml_bc->Dirichlet_eqn_list = NULL;
   return 0;
}

/* ******************************************************************** */
/* Initialize                                                           */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Init(ML_BdryPts *ml_bc)
{
   ml_bc->ML_id = ML_ID_BC;
   ml_bc->Dirichlet_grid_CreateOrDup = 0;
   ml_bc->Dirichlet_grid_length = 0;
   ml_bc->Dirichlet_grid_list = NULL;
   ml_bc->Dirichlet_eqn_CreateOrDup = 0;
   ml_bc->Dirichlet_eqn_length = 0;
   ml_bc->Dirichlet_eqn_list = NULL;
   return 0;
}

/* ******************************************************************** */
/* Destructor                                                           */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Destroy(ML_BdryPts **bc)
{
   ML_BdryPts *ml_bc;

   ml_bc = (*bc);
   if ( ml_bc->ML_id != ML_ID_BC ) {
      printf("ML_BdryPts_Destroy : wrong object.\n");
      exit(1);
   }
   ml_bc->ML_id = -1;
   ml_bc->Dirichlet_grid_CreateOrDup = 0;
   ml_bc->Dirichlet_grid_length = 0;
   if ( ml_bc->Dirichlet_grid_CreateOrDup == 1 )
      ML_memory_free( (void**) &(ml_bc->Dirichlet_grid_list) );
   ml_bc->Dirichlet_eqn_CreateOrDup = 0;
   ml_bc->Dirichlet_eqn_length = 0;
   if ( ml_bc->Dirichlet_eqn_CreateOrDup == 1 )
      ML_memory_free( (void**) &(ml_bc->Dirichlet_eqn_list) );
   ML_memory_free( (void**) bc );
   (*bc) = NULL;
   return 0;
}

/* ******************************************************************** */
/* Destructor (only destroy internal data structures)                   */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Clean(ML_BdryPts *ml_bc)
{
   if ( ml_bc->ML_id != ML_ID_BC ) {
      printf("ML_BdryPts_Clean : wrong object.\n");
      exit(1);
   }
   ml_bc->ML_id = -1;
   ml_bc->Dirichlet_grid_length = 0;
   if ( ml_bc->Dirichlet_grid_CreateOrDup == 1 )
      ML_memory_free( (void**) &(ml_bc->Dirichlet_grid_list) );
   ml_bc->Dirichlet_grid_CreateOrDup = 0;
   ml_bc->Dirichlet_eqn_length = 0;
   if ( ml_bc->Dirichlet_eqn_CreateOrDup == 1 )
      ML_memory_free( (void**) &(ml_bc->Dirichlet_eqn_list) );
   ml_bc->Dirichlet_eqn_CreateOrDup = 0;
   return 0;
}

/* ******************************************************************** */
/* check the existence of Dirichlet grid list                           */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Check_Dirichlet_Grid(ML_BdryPts *ml_bc)
{
   if ( ml_bc->Dirichlet_grid_list != NULL ) return 1;
   else                                      return 0;
}

/* ******************************************************************** */
/* check the existence of Dirichlet eqn list                            */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Check_Dirichlet_Eqn(ML_BdryPts *ml_bc)
{
   if ( ml_bc->Dirichlet_eqn_list != NULL ) return 1;
   else                                     return 0;
}

/* ******************************************************************** */
/* Get the Dirichlet grid list pointer and length                       */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Get_Dirichlet_Grid_Info(ML_BdryPts *bc, int *n, int **list)
{
   (*n) = bc->Dirichlet_grid_length;
   (*list) = bc->Dirichlet_grid_list;
   return 0;
}

/* ******************************************************************** */
/* Get the Dirichlet equation list pointer and length                   */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Get_Dirichlet_Eqn_Info(ML_BdryPts *bc, int *n, int **list)
{
   (*n) = bc->Dirichlet_eqn_length;
   (*list) = bc->Dirichlet_eqn_list;
   return 0;
}

/* ******************************************************************** */
/* Load the Dirichlet boundary points                                   */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Load_Dirichlet_Grid(ML_BdryPts *ml_bc, int leng, int *list)
{
   int i, nbytes;

   if ( ml_bc->ML_id != ML_ID_BC ) {
      printf("ML_BdryPts_Load_Dirichlet_Grid : wrong object.\n");
      exit(1);
   }
   if ( leng < 0 ) {
      printf("ML_BdryPts_Load_Dirichlet_Grid warning : length <= 0.\n");
      exit(1);
   }
   if ( ml_bc->Dirichlet_grid_CreateOrDup == 1 )
      ML_memory_free( (void**) &(ml_bc->Dirichlet_grid_list) );

   nbytes = (leng+1)*sizeof(int);
   ML_memory_alloc((void**) &(ml_bc->Dirichlet_grid_list), (unsigned int) nbytes, "BC2");
   ml_bc->Dirichlet_grid_length = leng;
   ml_bc->Dirichlet_grid_CreateOrDup = 1;
   for ( i = 0; i < leng; i++ )
      ml_bc->Dirichlet_grid_list[i] = list[i];

   return 0;
}

/* ******************************************************************** */
/* Load the Dirichlet boundary points in equation space                 */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Load_Dirichlet_Eqn(ML_BdryPts *ml_bc, int leng, int *list)
{
   int i, nbytes;

   if ( ml_bc->ML_id != ML_ID_BC ) {
      printf("ML_BdryPts_Load_Dirichlet_Eqn : wrong object.\n");
      exit(1);
   }
   if ( leng <= 0 ) {
      printf("ML_BdryPts_Load_Dirichlet_Eqn warning : length <= 0.\n");
      exit(1);
   }
   if ( ml_bc->Dirichlet_eqn_CreateOrDup == 1 )
      ML_memory_free( (void**) &(ml_bc->Dirichlet_eqn_list) );

   nbytes = leng * sizeof(int);
   ML_memory_alloc((void**) &(ml_bc->Dirichlet_eqn_list), (unsigned int) nbytes, "BC3");
   ml_bc->Dirichlet_eqn_length = leng;
   ml_bc->Dirichlet_eqn_CreateOrDup = 1;
   for ( i = 0; i < leng; i++ )
      ml_bc->Dirichlet_eqn_list[i] = list[i];

   return 0;
}

/* ******************************************************************** */
/* Duplicate the Dirichlet boundary points from grid to equation space  */
/* -------------------------------------------------------------------- */

int ML_BdryPts_Copy_Dirichlet_GridToEqn(ML_BdryPts *ml_bc)
{

   if ( ml_bc->ML_id != ML_ID_BC ) {
      printf("ML_BdryPts_Copy_Dirichlet_GridToEqn : wrong object.\n");
      exit(1);
   }
   if ( ml_bc->Dirichlet_eqn_CreateOrDup == 1 )
      ML_memory_free( (void**) &(ml_bc->Dirichlet_eqn_list) );

   ml_bc->Dirichlet_eqn_length = ml_bc->Dirichlet_grid_length;
   ml_bc->Dirichlet_eqn_list = ml_bc->Dirichlet_grid_list; 
   ml_bc->Dirichlet_eqn_CreateOrDup = 2;
   return 0;
}

/* ******************************************************************** */
/* reset the incoming boundary points to 0 in grid space                */
/* -------------------------------------------------------------------- */

int ML_BdryPts_ApplyZero_Dirichlet_Grid(ML_BdryPts *ml_bc, double *invec)
{
   int i, length, *list;

   if ( ml_bc->ML_id != ML_ID_BC ) {
      printf("ML_BdryPts_ApplyZero_Dirichlet_Grid : wrong object.\n");
      exit(1);
   }

   length = ml_bc->Dirichlet_grid_length;
   list = ml_bc->Dirichlet_grid_list;
   for ( i = 0; i < length; i++ ) invec[list[i]] = 0.0; 
   return 0;
}

/* ******************************************************************** */
/* reset the incoming boundary points to 0 in equation space            */
/* -------------------------------------------------------------------- */

int ML_BdryPts_ApplyZero_Dirichlet_Eqn(ML_BdryPts *ml_bc, double *invec)
{
   int i, length, *list;

   if ( ml_bc->ML_id != ML_ID_BC ) {
      printf("ML_BdryPts_ApplyZero_Dirichlet_Eqn : wrong object.\n");
      exit(1);
   }

   length = ml_bc->Dirichlet_eqn_length;
   list = ml_bc->Dirichlet_eqn_list;
   for ( i = 0; i < length; i++ ) invec[list[i]] = 0.0; 
   return 0;
}


