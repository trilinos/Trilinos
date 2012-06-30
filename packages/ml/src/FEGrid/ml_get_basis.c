/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : April, 1998                                          */
/* ******************************************************************** */

/* ******************************************************************** */
/* Usr_compute_basis_coefficients - this subroutine generates           */
/*  coefficients assuming the elements are 3-D subcubes with 8 vertices */
/*  in a unit cube.  In addition, it is assumed that the vertices for   */
/*  each element are stored in the order of x-first, y-second and       */
/*  z-last manner.                                                      */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_gridfunc.h"
#include "ml_utils.h"

extern ML_GridFunc *gridfcns_basis;


int ML_compute_basis_coefficients3D(void *grid, double *coord,
                                   int ncoord, double *coefs, int *coef_ptr)
{
   int    *vlist, the_pt, ncnt, i, j, ind, not_all_zero;
   int    max_vert_per_ele;
   double xyz[3];
   double x, y, z, xdist, ydist, zdist, xmax, ymax, zmax, xmin, ymin, zmin;
   double xwidth, ywidth, zwidth, coarse_x, coarse_y, coarse_z, local_coef[8];

   /* checking the presence of the set of grid access functions */
 
   if ( gridfcns_basis == NULL )
   {
      printf("Error in compute_basis : no grid functions available. \n");
      exit(0);
   }
 
   /* fetch the vertices (local vertex numbers) for the given coarse  */
   /* element (leng = the number of vertices for the element)         */

   max_vert_per_ele = gridfcns_basis->ML_MaxElmntVert;

   ML_memory_alloc( (void**) &vlist, max_vert_per_ele*sizeof(int), "BAS");

   /* fetch the left-bottom-front and right-top-back vertices, which  */
   /* should give information about its bounds in each dimension      */

   xmax = ymax = zmax = -1.0E10;
   xmin = ymin = zmin =  1.0E10;
   for (i = 0; i < 8; i++)
   {
      if ( vlist[i] >= 0 )
      {
         the_pt   = vlist[i];
         gridfcns_basis->USR_grid_get_vertex_coordinate(grid, the_pt, xyz);
         if (xyz[0] > xmax) xmax = xyz[0];
         if (xyz[0] < xmin) xmin = xyz[0];
         if (xyz[1] > ymax) ymax = xyz[1];
         if (xyz[1] < ymin) ymin = xyz[1];
         if (xyz[2] > zmax) zmax = xyz[2];
         if (xyz[2] < zmin) zmin = xyz[2];
      }
   }
   if (xmax == xmin || ymax == ymin || zmax == zmin)
   {
      printf("Error : get_basis - width = 0. \n");
      exit(-1);
   }
   xwidth = (double) 1.0 / (xmax - xmin);
   ywidth = (double) 1.0 / (ymax - ymin);
   zwidth = (double) 1.0 / (zmax - zmin);

   /* Now examine each incoming vertex and determine its relationship */
   /* with the coarse element vertices.                               */

   ncnt   = 0;

   for (i=0; i<ncoord; i++)
   {

      /* fetch the coordinate of the fine vertex */

      x = (double) coord[3*i]; 
      y = (double) coord[3*i+1]; 
      z = (double) coord[3*i+2];

      /* for each of the 8 coarse vertices */

      not_all_zero = 0;
      for ( j = 0; j < 8; j++ )
      {

         /* get the coarse vertex coordinate and find its */
         /* distance from the fine vertex in question     */

         ind = vlist[j];
         if ( ind >= 0 )
         {
            gridfcns_basis->USR_grid_get_vertex_coordinate(grid, ind, xyz);
            coarse_x = (double) xyz[0];
            coarse_y = (double) xyz[1];
            coarse_z = (double) xyz[2];
            xdist    = 1.0 - ML_dabs((x - coarse_x)) * xwidth;
            ydist    = 1.0 - ML_dabs((y - coarse_y)) * ywidth;
            zdist    = 1.0 - ML_dabs((z - coarse_z)) * zwidth;
            if (xdist > 0.0 && ydist > 0.0 && zdist > 0.0)
            { 
               local_coef[j] = xdist * ydist * zdist;
               if (local_coef[j] > 1.0E-6) not_all_zero++;
               else                        local_coef[j] = 0.0;
            } else local_coef[j] = 0.0;
         } else local_coef[j] = 0.0;
      }
      if (not_all_zero > 0)
      {
         for ( j = 0; j < 8; j++ ) coefs[ncnt++] = local_coef[j];
         coef_ptr[i] = 8;
      } 
      else
      {
         coefs[ncnt++] = -1.0;
         coef_ptr[i] = 1;
      }
   }
   ML_memory_free( (void **) &vlist );
   return 0;
}

/* ******************************************************************** */

int ML_compute_basis_coefficients2D(void *grid, double *coord,
                                   int ncoord, double *coefs, int *coef_ptr)
{
   int    *vlist, the_pt, ncnt, i, j, ind, not_all_zero;
   int    max_vert_per_ele;
   double xy[2];
   double x, y, xdist, ydist, xmax, ymax, xmin, ymin;
   double xwidth, ywidth, coarse_x, coarse_y, local_coef[4];

   /* checking the presence of the set of grid access functions */
 
   if ( gridfcns_basis == NULL )
   {
      printf("Error in compute_basis : no grid functions available. \n");
      exit(0);
   }
 
   /* fetch the vertices (local vertex numbers) for the given coarse  */
   /* element (leng = the number of vertices for the element)         */

   max_vert_per_ele = gridfcns_basis->ML_MaxElmntVert;

   ML_memory_alloc( (void**) &vlist, max_vert_per_ele*sizeof(int), "BAS");

   /* fetch the left-bottom-front and right-top-back vertices, which  */
   /* should give information about its bounds in each dimension      */

   xmax = ymax = -1.0E10;
   xmin = ymin = 1.0E10;
   for (i = 0; i < 4; i++)
   {
      if ( vlist[i] >= 0 )
      {
         the_pt   = vlist[i];
         gridfcns_basis->USR_grid_get_vertex_coordinate(grid, the_pt, xy);
         if (xy[0] > xmax) xmax = xy[0];
         if (xy[0] < xmin) xmin = xy[0];
         if (xy[1] > ymax) ymax = xy[1];
         if (xy[1] < ymin) ymin = xy[1];
      }
   }
   if (xmax == xmin || ymax == ymin )
   {
      printf("Error : get_basis - width = 0. \n");
      exit(-1);
   }
   xwidth = (double) 1.0 / (xmax - xmin);
   ywidth = (double) 1.0 / (ymax - ymin);

   /* Now examine each incoming vertex and determine its relationship */
   /* with the coarse element vertices.                               */

   ncnt   = 0;

   for (i=0; i<ncoord; i++)
   {

      /* fetch the coordinate of the fine vertex */

      x = (double) coord[2*i]; 
      y = (double) coord[2*i+1]; 

      /* for each of the 4 coarse vertices */

      not_all_zero = 0;
      for ( j = 0; j < 4; j++ )
      {

         /* get the coarse vertex coordinate and find its */
         /* distance from the fine vertex in question     */

         ind = vlist[j];
         if ( ind >= 0 )
         {
            gridfcns_basis->USR_grid_get_vertex_coordinate(grid, ind, xy);
            coarse_x = (double) xy[0];
            coarse_y = (double) xy[1];
            xdist    = 1.0 - ML_dabs((x - coarse_x)) * xwidth;
            ydist    = 1.0 - ML_dabs((y - coarse_y)) * ywidth;
            if (xdist > 0.0 && ydist > 0.0 )
            { 
               local_coef[j] = xdist * ydist;
               if (local_coef[j] > 1.0E-6) not_all_zero++;
               else                        local_coef[j] = 0.0;
            } else local_coef[j] = 0.0;
         } else local_coef[j] = 0.0;
      }
      if (not_all_zero > 0)
      {
         for ( j = 0; j < 4; j++ ) coefs[ncnt++] = local_coef[j];
         coef_ptr[i] = 4;
      } 
      else
      {
         coefs[ncnt++] = -1.0;
         coef_ptr[i] = 1;
      }
   }
   ML_memory_free( (void **) &vlist );
   return 0;
}

