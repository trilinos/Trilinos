/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for manipulating Grid objects                              */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : Sept, 1997                                           */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_gridagx.h"

/* ******************************************************************** */
/* initialize the grid structure                                        */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Create(ML_GridAGX **ingrid)
{
   ML_GridAGX *grid;

   ML_memory_alloc( (void **) ingrid, sizeof(ML_GridAGX), "GD1" );
   grid                     = (*ingrid);
   grid->ML_id              = ML_ID_GRIDAGX;
   grid->Ndim               = 0;
   grid->Nelements          = 0;
   grid->Nvertices          = 0;
   grid->Nvertices_expanded = 0;
   grid->global_element     = 0;
   grid->global_vertex      = 0;
   grid->x                  = 0;
   grid->y                  = 0;
   grid->z                  = 0;
   grid->elmnt_proc_map     = 0;
   grid->node_proc_map      = 0;
   /* ML_IntList_Create( &(grid->ele_nodes), 0, 0 ); */
   grid->ele_nodes          = 0;
   return 0;
}

/* ******************************************************************** */
/* deallocate the grid structure                                        */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Destroy( ML_GridAGX **ingrid )
{
   ML_GridAGX *grid;

   grid = (*ingrid);
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Destroy : wrong object. \n");
      exit(1);
   }
   if ( grid->global_element != 0 )
   {
      ML_memory_free( (void **) &(grid->global_element));
   }
   if ( grid->global_vertex != 0 )
   {
      ML_memory_free( (void **) &(grid->global_vertex));
   }
   if ( grid->ele_nodes != 0 ) ML_IntList_Destroy(&(grid->ele_nodes));
   if ( grid->x != 0 ) ML_memory_free( (void**) &(grid->x) );
   if ( grid->y != 0 ) ML_memory_free( (void**) &(grid->y) );
   if ( grid->z != 0 ) ML_memory_free( (void**) &(grid->z) );
   if ( grid->elmnt_proc_map != 0 ) {
      ML_memory_free( (void **) &(grid->elmnt_proc_map) );
   }
   if ( grid->node_proc_map != 0 )
   {
      ML_memory_free( (void **) &(grid->node_proc_map) );
   }
   grid->ML_id = -1;
   ML_memory_free ( (void**) ingrid );
   return 0;
}

/* ******************************************************************** */
/* fetch all information of an element                                  */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_Element(ML_GridAGX *grid,int index,
                             ML_ElementAGX *elmnt)
{
   int    i, *nlist, ind, ncnt, nd;
   double x, y, z;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_Element : wrong object. \n");
      exit(1);
   }
   if (index < 0 || index >= grid->Nelements)
   {
      printf("ML_GridAGX_Get_Element : access error.\n");
      exit(-1);
   }

   ML_ElementAGX_Reuse( elmnt );
   ML_memory_alloc( (void **) &nlist, 30 * sizeof(int), "GD1");
   ML_IntList_Get_Sublist( grid->ele_nodes, index, &ncnt, nlist );
   if ( ncnt > 30 )
   {
      printf("Warning : Int_lists - sublist length > 30.\n");
      exit(0);
   }
   nd = ML_GridAGX_Get_Dimension( grid );
   for ( i = 0; i < ncnt; i++)
   {
      ind = nlist[i];
      x   = grid->x[ind];
      y   = grid->y[ind];
      if ( nd > 2 ) z = grid->z[ind]; else z = 0.0;
      ML_ElementAGX_Load_VertCoordinate( elmnt, ind, x, y, z );
   }
   ML_memory_free( (void **) &nlist);
   return 0;
}

/* ******************************************************************** */
/* print all information about the grid (for debug purpose)             */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Print( ML_GridAGX *grid )
{
   int i;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Print : wrong object. \n");
      exit(1);
   }
   printf("Grid : number of elements    = %d \n",grid->Nelements);
   printf("Grid : number of vertices    = %d \n",grid->Nvertices);
   for ( i = 0; i < grid->Nelements; i++ )
   {
      printf("Grid : global element %d = %d \n",i,grid->global_element[i]);
   }
   for ( i = 0; i < grid->Nvertices_expanded; i++ )
   {
      printf("Grid : global vertex %d = %d \n",i,grid->global_vertex[i]);
   }
   if ( grid->Ndim == 2 ) {
      for ( i = 0; i < grid->Nvertices_expanded; i++ )
      {
         printf("Grid : (x,y) %d = %e %e \n", i, grid->x[i], grid->y[i]);
      }
   } else {
      for ( i = 0; i < grid->Nvertices_expanded; i++ )
      {
         printf("Grid : (x,y,z) %d = %e %e %e \n", i, grid->x[i], grid->y[i],
                grid->z[i]);
      }
   }
   ML_IntList_Print( grid->ele_nodes );
   return 0;
}

/* ******************************************************************** */
/* return the grid dimension                                            */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_Dimension( ML_GridAGX *grid )
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_Dimension : wrong object. \n");
      exit(1);
   }
   return grid->Ndim;
}

/* ******************************************************************** */
/* return the number of local vertices in the grid                      */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_NVert( ML_GridAGX *grid )
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_NVert: wrong object (%d).\n", grid->ML_id);
      exit(1);
   }
   return grid->Nvertices;
}

/* ******************************************************************** */
/* return the number of local elements in the grid                      */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_NElmnts( ML_GridAGX *grid )
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_NElmnts : wrong object. \n");
      exit(1);
   }
   return grid->Nelements;
}

/* ******************************************************************** */
/* return the global number of a local element                          */
/* -------------------------------------------------------------------- */

ml_big_int ML_GridAGX_Get_ElmntGlobalNum(ML_GridAGX *grid, int index)
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_ElmntGlobalNum : wrong object. \n");
      exit(1);
   }
   return grid->global_element[index];
}

/* ******************************************************************** */
/* return the number of vertices for a given local element              */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_ElmntNVert(ML_GridAGX *grid, int index)
{
   int ndiff;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_ElmntNVert : wrong object. \n");
      exit(1);
   }
  ndiff = grid->ele_nodes->start[index+1] - grid->ele_nodes->start[index];
  return ndiff;
}

/* ******************************************************************** */
/* return the global number of a local vertex                           */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_VertGlobalNum(ML_GridAGX *grid, int index)
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_VertGlobalNum : wrong object. \n");
      exit(1);
   }
   return grid->global_vertex[index];
}

/* ******************************************************************** */
/* return the vertex list of a given local element                      */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_ElmntVertList(ML_GridAGX *grid,int index,int *vlist)
{
   int i, k, begin, end;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_ElmntVertList : wrong object. \n");
      exit(1);
   }
   begin = grid->ele_nodes->start[index];
   end   = grid->ele_nodes->start[index+1];
   k     = 0;
   for ( i = begin; i < end; i++) vlist[k++] = grid->ele_nodes->members[i];
   return (end-begin);
}

/* ******************************************************************** */
/* return the coordinate of a given local vertex                        */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Get_VertCoordinate(ML_GridAGX *grid, int index,
                                     double *coord)
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Get_VertCoordinate : wrong object. \n");
      exit(1);
   }
   coord[0] = grid->x[index];
   coord[1] = grid->y[index];
   if (grid->Ndim > 2) coord[2] = grid->z[index];
   return 0;
}

/* ******************************************************************** */
/* load the grid dimension                                              */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Set_Dimension( ML_GridAGX *grid, int ndim )
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Set_Dimension : wrong object. \n");
      exit(1);
   }
   grid->Ndim = ndim;
   return 0;
}

/* ******************************************************************** */
/* return the number of local vertices in the grid                      */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Set_NVert( ML_GridAGX *grid, int nvert )
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Set_NVert : wrong object. \n");
      exit(1);
   }
   grid->Nvertices = nvert;
   return 0;
}

/* ******************************************************************** */
/* return the number of local elements in the grid                      */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Set_NElmnts(ML_GridAGX *grid, int nelmnt, int vert_per_ele)
{
   int ncnt;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Set_NElmnts : wrong object. \n");
      exit(1);
   }
   grid->Nelements = nelmnt;
   ncnt = nelmnt * vert_per_ele;
   ML_IntList_Create( &(grid->ele_nodes), nelmnt, ncnt);
   return 0;
}

/* ******************************************************************** */
/* load the global numbers of all local elements                        */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Load_ElmntGlobalNum(ML_GridAGX *grid, int leng, ml_big_int *gnum)
{
   int i;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Load_ElmntGlobalNum : wrong object. \n");
      exit(1);
   }
   if (grid->Nelements != leng)
   {
      printf("ML_GridAGX_Load_ElmntGlobalNum : - lengths do not match. \n");
   }
   ML_memory_alloc((void**) &(grid->global_element),leng*sizeof(grid->global_element[0]),"GD2");
   for ( i = 0; i < leng; i++ ) grid->global_element[i] = gnum[i];
   return 0;
}

/* ******************************************************************** */
/* load the global numbers of all local vertices                        */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Load_VertGlobalNum(ML_GridAGX *grid, int leng, int *gnum)
{
   int i, nbytes;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Load_VertGlobalNum : wrong object. \n");
      exit(1);
   }
   grid->Nvertices_expanded = leng;
   nbytes = leng * sizeof(int);
   ML_memory_alloc((void**) &(grid->global_vertex), (unsigned int) nbytes, "GD3");
   for ( i = 0; i < leng; i++ ) grid->global_vertex[i] = gnum[i];
   return 0;
}

/* ******************************************************************** */
/* load the vertex list of a local element                              */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Load_ElmntVertList(ML_GridAGX *grid, int leng, int *vlist)
{
   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Load_ElmntVertList : wrong object. \n");
      exit(1);
   }
   ML_IntList_Load_Sublist( grid->ele_nodes, leng, vlist);
   return 0;
}

/* ******************************************************************** */
/* load the coordinates of all local vertices                           */
/* -------------------------------------------------------------------- */

int ML_GridAGX_Load_AllVertCoordinates(ML_GridAGX *grid,int leng,
                                         double *coord)
{
   int i, ndim, nbytes;

   if ( grid->ML_id != ML_ID_GRIDAGX )
   {
      printf("ML_GridAGX_Load_AllVertCoordinates : wrong object. \n");
      exit(1);
   }
   ndim = grid->Ndim;
   nbytes = leng * sizeof(double);
   ML_memory_alloc( (void**) &(grid->x), (unsigned int) nbytes, "GDX");
   ML_memory_alloc( (void**) &(grid->y), (unsigned int) nbytes, "GDY");
   if (ndim > 2)
      ML_memory_alloc( (void**) &(grid->z), (unsigned int) nbytes, "GDZ");
   for ( i = 0; i < leng; i++ )
   {
      grid->x[i] = coord[i*ndim];
      grid->y[i] = coord[i*ndim+1];
      if (ndim > 2) grid->z[i] = coord[i*ndim+2];
   }
   return 0;
}

