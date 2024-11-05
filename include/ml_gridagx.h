/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* local (to ML) data structure to hold mesh information given a finite */
/* element mesh (used in automatic grid transfer generation)            */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : April, 1997                                          */
/* ******************************************************************** */

#ifndef __MLGRIDAGX__
#define __MLGRIDAGX__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include <stdio.h>
/* #include <stdlib.h> */
#include "ml_common.h"
#include "ml_defs.h"
#include "ml_memory.h"
#include "ml_intlist.h"
#include "ml_elementagx.h"

/* ******************************************************************** */
/* definition of the grid structure                                     */
/*  Ndim              : no. of spatial dimensions                       */
/*  Nvertices         : no. of vertices residing in local processor     */
/*  Nelements         : no. of elements residing locally                */
/*  ele_nodes         : element to node list (local node no.)           */
/*  x[i],y[i],z[i]    : coordinate of local (and ghost) node i          */
/*  global_element[i] : global element number of the i-th local element */
/*  global_vertex[i]  : global vertex number of the i-th local vertex   */
/*  elmnt_proc_map[i] : processor where element[global_element[i]] is   */
/*  node_proc_map[i]  : processor where global_vertex[i]] is            */
/*                                                                      */
/* Note : x,y,z, and global_vertex arrays can have length longer than   */
/*        Nvertices since they also have to store external vertices     */
/*        that are vertices to local elements.                          */
/* -------------------------------------------------------------------- */

typedef struct ML_GridAGX_Struct
{
   int              ML_id;
   int              Ndim;
   int              Nvertices, Nvertices_expanded;
   int              Nelements;
   ML_IntList       *ele_nodes;
   double           *x, *y, *z;
   ml_big_int       *global_element;
   int              *global_vertex;
   int              *elmnt_proc_map;
   int              *node_proc_map;

} ML_GridAGX;

/* ******************************************************************** */
/* functions to manipulate the grid structure                           */
/* -------------------------------------------------------------------- */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

extern int  ML_GridAGX_Create( ML_GridAGX ** );
extern int  ML_GridAGX_Destroy( ML_GridAGX ** );
extern int  ML_GridAGX_Get_Element(ML_GridAGX *,int,ML_ElementAGX *);
extern int  ML_GridAGX_Print( ML_GridAGX * );

extern int  ML_GridAGX_Get_Dimension( ML_GridAGX * );
extern int  ML_GridAGX_Get_NVert( ML_GridAGX * );
extern int  ML_GridAGX_Get_NElmnts( ML_GridAGX * );
extern ml_big_int  ML_GridAGX_Get_ElmntGlobalNum( ML_GridAGX *, int );
extern int  ML_GridAGX_Get_ElmntNVert( ML_GridAGX *, int );
extern int  ML_GridAGX_Get_ElmntVertList( ML_GridAGX *, int, int * );
extern int  ML_GridAGX_Get_VertGlobalNum( ML_GridAGX *, int );
extern int  ML_GridAGX_Get_VertCoordinate(ML_GridAGX *grid,int,double*);

extern int  ML_GridAGX_Set_Dimension(ML_GridAGX *, int);
extern int  ML_GridAGX_Set_NVert(ML_GridAGX *, int);
extern int  ML_GridAGX_Set_NElmnts(ML_GridAGX *, int, int);
extern int  ML_GridAGX_Load_ElmntGlobalNum(ML_GridAGX *, int, ml_big_int *);
extern int  ML_GridAGX_Load_VertGlobalNum(ML_GridAGX *, int, int *);
extern int  ML_GridAGX_Load_ElmntVertList(ML_GridAGX *, int, int *);
extern int  ML_GridAGX_Load_AllVertCoordinates(ML_GridAGX*,int,double*);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
