/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for accessing a user grid                                  */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>

#include "ml_grid.h"
#include "ml_memory.h"

/* ******************************************************************** */
/* Create a ML_Grid data structure                                      */
/* -------------------------------------------------------------------- */

int ML_Grid_Create( ML_Grid ** grid ) 
{
   ML_Grid *ml_grid;
  
   ML_memory_alloc( (void**) grid, sizeof(ML_Grid), "GR1" );
   ml_grid               = (*grid);
   ml_grid->ML_id        = ML_ID_GRID;
   ml_grid->Grid         = NULL;
   ml_grid->gridfcn      = NULL;
   ml_grid->gf_SetOrLoad = 0;
   return 0;
}

/* ******************************************************************** */
/* Initialize a ML_Grid data structure                                  */
/* -------------------------------------------------------------------- */

int ML_Grid_Init( ML_Grid *ml_grid ) 
{
   ml_grid->ML_id        = ML_ID_GRID;
   ml_grid->Grid         = NULL;
   ml_grid->gridfcn      = NULL;
   ml_grid->gf_SetOrLoad = 0;
   return 0;
}

/* ******************************************************************** */
/* Destroy a ML_Grid data structure                                     */
/* -------------------------------------------------------------------- */

int ML_Grid_Destroy( ML_Grid ** grid ) 
{
   ML_Grid *ml_grid;
  
   ml_grid               = (*grid);
   ml_grid->ML_id        = -1;
   ml_grid->Grid         = NULL;
   if ( ml_grid->gridfcn != NULL && ml_grid->gf_SetOrLoad == 2 )
      ML_GridFunc_Destroy( &(ml_grid->gridfcn) );
   ml_grid->gridfcn      = NULL;
   ml_grid->gf_SetOrLoad = 0;
   ML_memory_free( (void**) grid);
   (*grid) = NULL;
   return 0;
}

/* ******************************************************************** */
/* Clean up a ML_Grid data structure                                    */
/* -------------------------------------------------------------------- */

int ML_Grid_Clean( ML_Grid *ml_grid ) 
{
   ml_grid->ML_id        = -1;
   ml_grid->Grid         = NULL;
   if ( ml_grid->gridfcn != NULL && ml_grid->gf_SetOrLoad == 2 )
      ML_GridFunc_Destroy( &(ml_grid->gridfcn) );
   ml_grid->gridfcn      = NULL;
   ml_grid->gf_SetOrLoad = 0;
   return 0;
}

/* ******************************************************************** */
/* Set the data field of the ML_Grid data structure                     */
/* -------------------------------------------------------------------- */

int ML_Grid_Set_Grid( ML_Grid *grid, void *data )
{
   if ( grid->ML_id != ML_ID_GRID ) {
      printf("ML_Grid_Set_Grid error : wrong object.\n");
      exit(-1);
   }
   grid->Grid = data;
   return 0;
}

/* ******************************************************************** */
/* Set the function field of the ML_Grid data structure                 */
/* -------------------------------------------------------------------- */

int ML_Grid_Set_GridFunc( ML_Grid *grid, ML_GridFunc *gf )
{
   if ( grid->ML_id != ML_ID_GRID ) {
      printf("ML_Grid_Set_GridFunc error : wrong object.\n");
      exit(-1);
   }
   grid->gridfcn = gf;
   return 0;
}

/* ******************************************************************** */
/* Get the function field of the ML_Grid data structure                 */
/* -------------------------------------------------------------------- */

int ML_Grid_Get_GridFunc( ML_Grid *grid, ML_GridFunc **gf )
{
   if ( grid->ML_id != ML_ID_GRID ) {
      printf("ML_Grid_Get_GridFunc error : wrong object.\n");
      exit(-1);
   }
   (*gf) = grid->gridfcn;
   return 0;
}

/* ******************************************************************** */
/* Get the function field of the ML_Grid data structure                 */
/* -------------------------------------------------------------------- */

int ML_Grid_Create_GridFunc( ML_Grid *grid )
{
   if ( grid->ML_id != ML_ID_GRID ) {
      printf("ML_Grid_Create_GridFunc error : wrong object.\n");
      exit(-1);
   }
   if ( grid->gridfcn != NULL && grid->gf_SetOrLoad == 2 )
      ML_GridFunc_Destroy( (ML_GridFunc**) &(grid->gridfcn) );
   
   ML_GridFunc_Create( &(grid->gridfcn) );
   grid->gf_SetOrLoad = 2;
   return 0;
}

