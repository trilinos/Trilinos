/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for manipulating ML_ElementAGX objects                     */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : April, 1998                                          */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_elementagx.h"

/* ******************************************************************** */
/* initialization function                                              */
/* -------------------------------------------------------------------- */

int ML_ElementAGX_Create(ML_ElementAGX **elmntp, int nd, int nvert)
{
   ML_ElementAGX *elmnt;

   ML_memory_alloc((void**) elmntp, sizeof(ML_ElementAGX), "EL1");
   elmnt = (*elmntp);
   
   elmnt->ndim = nd;
   elmnt->Nvertices = 0;
   ML_memory_alloc((void**) &(elmnt->vertices), nvert*sizeof(int),"EL2");
   ML_memory_alloc((void**) &(elmnt->x), nvert*sizeof(double),"EL3");
   ML_memory_alloc((void**) &(elmnt->y), nvert*sizeof(double),"EL4");
   if ( nd > 2 ) 
      ML_memory_alloc((void**) &(elmnt->z), nvert*sizeof(double),"EL5");
   else
      elmnt->z = 0;
   return 0;
}

/* ******************************************************************** */
/* destructor function                                                  */
/* -------------------------------------------------------------------- */

int ML_ElementAGX_Destroy(ML_ElementAGX **elmntp)
{
   ML_ElementAGX *elmnt;
  
   elmnt = (*elmntp);
   ML_memory_free( (void **) &(elmnt->vertices) );
   ML_memory_free( (void **) &(elmnt->x) );
   ML_memory_free( (void **) &(elmnt->y) );
   if ( elmnt->z != 0 ) ML_memory_free( (void **) &(elmnt->z) );
   ML_memory_free( (void **) elmntp );
   return 0;
}

/* ******************************************************************** */
/* load element information                                             */
/* -------------------------------------------------------------------- */

int ML_ElementAGX_Load_VertCoordinate(ML_ElementAGX *elmnt, 
                               int nodenum, double x, double y, double z)
{
   int  index;

   index = elmnt->Nvertices++;
   elmnt->vertices[index] = nodenum;
   elmnt->x[index]        = x;
   elmnt->y[index]        = y;
   if ( elmnt->ndim > 2 ) elmnt->z[index] = z;
   return 0;
}

/* ******************************************************************** */
/* reset element parameters for reuse                                   */
/* -------------------------------------------------------------------- */

int ML_ElementAGX_Reuse(ML_ElementAGX *elmnt)
{
   elmnt->Nvertices = 0;
   return 0;
}

/* ******************************************************************** */
/* print element data                                                   */
/* -------------------------------------------------------------------- */

int ML_ElementAGX_Print(ML_ElementAGX *elmnt)
{
   int  i;
  
   printf("ElementAGX : number of vertices = %d \n",elmnt->Nvertices);
   if ( elmnt->ndim == 2 ) 
      for ( i = 0; i < elmnt->Nvertices; i++ ) 
         printf("    node number, x, y = %d %e %e \n",
                elmnt->vertices[i], elmnt->x[i], elmnt->y[i]);
   else
      for ( i = 0; i < elmnt->Nvertices; i++ ) 
         printf("    node number, x, y, z = %d %e %e %e \n",
                elmnt->vertices[i], elmnt->x[i], elmnt->y[i], elmnt->z[i]);
   return 0;
}

/* ******************************************************************** */
/* get information of a vertex from the element                         */
/* -------------------------------------------------------------------- */

int ML_ElementAGX_Get_VertCoordinate(ML_ElementAGX *elmnt, int index, 
                          int *nodenum, double *x, double *y, double *z)
{
   (*nodenum) = elmnt->vertices[index];
   (*x) = elmnt->x[index];
   (*y) = elmnt->y[index];
   if ( elmnt->ndim > 2 ) (*z) = elmnt->z[index];
   else                   (*z) = 0;
   return 0;
} 

/* ******************************************************************** */
/* check whether a list of coordinates are under the umbrella of the    */
/* element                                                              */
/* -------------------------------------------------------------------- */

int ML_ElementAGX_ComposeCandidates(ML_ElementAGX *element, int nvert, 
              double *coord, int *vlist, int *fnode_flag, int *ncand, 
              int *cand_list)
{
   int    i, icnt, ncnt=0, nd;
   double xmax, xmin, xin, ymax, ymin, yin, zmax, zmin, zin = 0; 

   xmax = ymax = zmax = - 1.0E10;
   xmin = ymin = zmin =   1.0E10;
   nd = element->ndim;
   for ( i = 0; i < element->Nvertices; i++ ) {
      xin = (double) element->x[i]; 
      yin = (double) element->y[i]; 
      if ( nd > 2 ) zin = (double) element->z[i]; 
      if (xin < xmin)  xmin = xin;
      if (xin > xmax)  xmax = xin;
      if (yin < ymin)  ymin = yin;
      if (yin > ymax)  ymax = yin;
      if ( nd > 2 ) {
         if ( zin < zmin)  zmin = zin;
         if ( zin > zmax)  zmax = zin;
      }
   }
   xmin = xmin - 0.0000001;
   xmax = xmax + 0.0000001;
   ymin = ymin - 0.0000001;
   ymax = ymax + 0.0000001;
   if ( nd > 2 ) {
      zmin = zmin - 0.0000001;
      zmax = zmax + 0.0000001;
   }
   if ( nd > 2 ) {
      if ((xmin >= xmax) || (ymin >= ymax) || (zmin >= zmax) ) { 
         printf("Error : max,min - %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
                 xmin,xmax,ymin,ymax,zmin,zmax); exit(-1);
      }
   } else {
      if ((xmin >= xmax) || (ymin >= ymax)) { 
         printf("Error : max,min - %14.7e %14.7e %14.7e %14.7e\n",
                 xmin,xmax,ymin,ymax); exit(-1);
      }
   }

   icnt = 0;
   for ( i = 0; i < nvert; i++ ) {
      if ( fnode_flag[vlist[i]] == -1 ) {
         xin = (double) coord[icnt++];
         yin = (double) coord[icnt++];
         if ( nd > 2 ) {
            zin = (double) coord[icnt++];
            if ((zin >= zmin) && (zin <= zmax)) {
               if ((yin >= ymin) && (yin <= ymax)) {
                  if ( (xin >= xmin) && (xin <= xmax)) { 
                     cand_list[ncnt++] = vlist[i];
                  } 
               } 
            } 
         } else {
            if ((yin >= ymin) && (yin <= ymax)) {
               if ( (xin >= xmin) && (xin <= xmax)) { 
                  cand_list[ncnt++] = vlist[i];
               } 
            } 
         }
      } else {
         icnt = icnt + 2;
         if ( nd > 2 ) icnt++;
      }
   }
   (*ncand) = ncnt;
   return 0;
} 

