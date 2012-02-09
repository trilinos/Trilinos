/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for accessing user grid functions                          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : Sept, 1998                                           */
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>

#include "ml_gridfunc.h"
#include "ml_memory.h"

/* ******************************************************************** */
/* Create a data structure for grid access functions                    */
/* -------------------------------------------------------------------- */

int ML_GridFunc_Create( ML_GridFunc ** gf )
{
   ML_GridFunc *gf_ptr;

   ML_memory_alloc( (void**) gf, sizeof(ML_GridFunc), "GF1" );
   gf_ptr                                  = (*gf);
   gf_ptr->ML_id                           = ML_ID_GRIDFCN;
   gf_ptr->ML_MaxElmntVert                 = 8;
   gf_ptr->USR_grid_get_dimension          = NULL;
   gf_ptr->USR_grid_get_nvertices          = NULL;
   gf_ptr->USR_grid_get_nelements          = NULL;
   gf_ptr->USR_grid_get_element_global_num = NULL;
   gf_ptr->USR_grid_get_element_nvertices  = NULL;
   gf_ptr->USR_grid_get_element_vlist      = NULL;
   gf_ptr->USR_grid_get_vertex_global_num  = NULL;
   gf_ptr->USR_grid_get_vertex_coordinate  = NULL;
   gf_ptr->USR_compute_basis_coefficients  = NULL;
   gf_ptr->USR_grid_get_element_volumes    = NULL;
   gf_ptr->USR_grid_get_element_matrix     = NULL;
   return 0;
}

/* ******************************************************************** */
/* destroy a data structure for grid access functions                   */
/* -------------------------------------------------------------------- */

int ML_GridFunc_Destroy( ML_GridFunc ** gf )
{
   ML_GridFunc *gf_ptr;

   gf_ptr = (*gf);
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN ) return -1;
   gf_ptr->ML_id = -1;
   ML_memory_free( (void **) gf );
   return 0;
}

/* ******************************************************************** */
/* Check that all grid access functions have been loaded.               */
/* (return -1 if not)                                                   */
/* -------------------------------------------------------------------- */

int ML_GridFunc_Check( ML_GridFunc *gf_ptr )
{
   int ready_flag = 1;

   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Check : wrong object. \n");
      exit(-1);
   }
   if ( gf_ptr->USR_grid_get_dimension == NULL ) ready_flag = 0;
   if ( gf_ptr->USR_grid_get_nvertices == NULL ) ready_flag = 0;
   if ( gf_ptr->USR_grid_get_nelements == NULL ) ready_flag = 0;
   if ( gf_ptr->USR_grid_get_element_global_num == NULL ) ready_flag = 0;
   if ( gf_ptr->USR_grid_get_element_nvertices == NULL )  ready_flag = 0;
   if ( gf_ptr->USR_grid_get_element_vlist == NULL )      ready_flag = 0;
   if ( gf_ptr->USR_grid_get_vertex_global_num == NULL )  ready_flag = 0;
   if ( gf_ptr->USR_grid_get_vertex_coordinate == NULL )  ready_flag = 0;
   if ( gf_ptr->USR_compute_basis_coefficients == NULL )  ready_flag = 0;
   /*if ( gf_ptr->USR_grid_get_element_volumes   == NULL )  ready_flag = 0;*/
   if ( ready_flag == 1 ) return 0;
   else                   return -1;
}

/* ******************************************************************** */
/* Set grid access functions                                            */
/* -------------------------------------------------------------------- */

int ML_GridFunc_Set_MaxVertPerElmnt(ML_GridFunc *gf, int leng)
{
   if ( gf->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_MaxVertPerElmnt: wrong object. \n");
      exit(-1);
   }
   if ( leng <= 0 )
   {
      printf("ML_GridFunc_Set_MaxVertPerElmnt : value <= 0. \n");
      exit(-1);
   }
   gf->ML_MaxElmntVert = leng;
   return 0;
}

/* ******************************************************************** */
/* Set grid access functions                                            */
/* -------------------------------------------------------------------- */
#ifdef NOTSTRICT_PROTO
int ML_GridFunc_Set_Function( ML_GridFunc *gf_ptr, int ind, int (*func)())
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_Function : wrong object. \n");
      exit(-1);
   }
   switch ( ind )
   {
      case 0  : gf_ptr->USR_grid_get_dimension = func; break;
      case 1  : gf_ptr->USR_grid_get_nvertices = func; break;
      case 2  : gf_ptr->USR_grid_get_nelements = func; break;
      case 3  : gf_ptr->USR_grid_get_element_global_num = func; break;
      case 4  : gf_ptr->USR_grid_get_element_nvertices  = func; break;
      case 5  : gf_ptr->USR_grid_get_element_vlist      = func; break;
      case 6  : gf_ptr->USR_grid_get_vertex_global_num  = func; break;
      case 7  : gf_ptr->USR_grid_get_vertex_coordinate  = func; break;
      case 8  : gf_ptr->USR_compute_basis_coefficients  = func; break;
      case 9  : gf_ptr->USR_grid_get_element_volumes    = func; break;
      case 10 : gf_ptr->USR_grid_get_element_matrix     = func; break;

      /* these are the ones eventually migrated to */

      case ML_GRID_DIMENSION :
                gf_ptr->USR_grid_get_dimension = func; break;
      case ML_GRID_NVERTICES :
                gf_ptr->USR_grid_get_nvertices = func; break;
      case ML_GRID_NELEMENTS :
                gf_ptr->USR_grid_get_nelements = func; break;
      case ML_GRID_ELEM_GLOBAL :
                gf_ptr->USR_grid_get_element_global_num = func; break;
      case ML_GRID_ELEM_NVERT :
                gf_ptr->USR_grid_get_element_nvertices  = func; break;
      case ML_GRID_ELEM_VLIST :
                gf_ptr->USR_grid_get_element_vlist      = func; break;
      case ML_GRID_VERT_GLOBAL :
                gf_ptr->USR_grid_get_vertex_global_num  = func; break;
      case ML_GRID_VERT_COORD :
                gf_ptr->USR_grid_get_vertex_coordinate  = func; break;
      case ML_GRID_BASISFCNS :
                gf_ptr->USR_compute_basis_coefficients  = func; break;
      case ML_GRID_ELEM_VOLUME :
                gf_ptr->USR_grid_get_element_volumes    = func; break;
      case ML_GRID_ELEM_MATRIX :
                gf_ptr->USR_grid_get_element_matrix     = func; break;
      default :
         printf("ML_GridFunc_Set_Function : function not recognized. \n");
         exit(-1);
         break;
   }
   return 0;
}
#endif
/* ******************************************************************** */
/* Set grid access functions INDIVIDUALLY                               */
/* -------------------------------------------------------------------- */

int ML_GridFunc_Set_GetDimension( ML_GridFunc *gf_ptr, int (*func)(void *))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetDimension : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_dimension = func;
   return 0;
}

int ML_GridFunc_Set_GetNVert( ML_GridFunc *gf_ptr, int (*func)(void *))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetNVert : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_nvertices = func;
   return 0;
}

int ML_GridFunc_Set_GetNElmnts( ML_GridFunc *gf_ptr, int (*func)(void *))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetNElmnts : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_nelements = func;
   return 0;
}

int ML_GridFunc_Set_GetElmntGlobalNum(ML_GridFunc *gf_ptr, ml_big_int (*func)(void *, int))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetElmntGlobalNum : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_element_global_num = func;
   return 0;
}

int ML_GridFunc_Set_GetElmntNVert(ML_GridFunc * gf_ptr, int (*func)(void *, int))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetElmntNVert : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_element_nvertices = func;
   return 0;
}

int ML_GridFunc_Set_GetElmntVertList( ML_GridFunc *gf_ptr, int (*func)(void *, int, int *))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetElmntVertList : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_element_vlist = func;
   return 0;
}

int ML_GridFunc_Set_GetVertGlobalNum(ML_GridFunc *gf_ptr, int (*func)(void *, int))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetVertGlobalNum : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_vertex_global_num = func;
   return 0;
}

int ML_GridFunc_Set_GetVertCoordinate(ML_GridFunc *gf_ptr, int (*func)(void *, int, double *))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetVertCoordinate : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_vertex_coordinate = func;
   return 0;
}

int ML_GridFunc_Set_ComputeBasisCoef( ML_GridFunc *gf_ptr, int (*func)(void*,int,double*,int,double*,int*))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_ComputeBasisCoef : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_compute_basis_coefficients = func;
   return 0;
}

int ML_GridFunc_Set_GetElmntVolumes( ML_GridFunc *gf_ptr, int (*func)(void*,int,int*,double*))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetElmntVolumes : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_element_volumes = func;
   return 0;
}

int ML_GridFunc_Set_GetElmntMatrix( ML_GridFunc *gf_ptr, int (*func)(void*,int,double**))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetElmntMatrix : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_element_matrix = func;
   return 0;
}

int ML_GridFunc_Set_GetElmntNullSpace(ML_GridFunc *gf_ptr, int (*func)(void*,int,double*))
{
   if ( gf_ptr->ML_id != ML_ID_GRIDFCN )
   {
      printf("ML_GridFunc_Set_GetElmntNullSpace : wrong object. \n");
      exit(-1);
   }
   gf_ptr->USR_grid_get_element_nullspace = func;
   return 0;
}

