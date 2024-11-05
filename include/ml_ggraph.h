/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* data structure to hold grid information in the form of a graph       */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : September, 1998                                      */
/* ******************************************************************** */

#ifndef __MLGRIDG__
#define __MLGRIDG__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include <stdio.h>
/* #include <stdlib.h> */

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_comm.h"
#include "ml_gridfunc.h"
#include "ml_memory.h"
#include "ml_comminfoop.h"
#include "ml_operator.h"
#include "ml_mat_formats.h"

/* ******************************************************************** */
/* definition of the grid graph structure                               */
/* -------------------------------------------------------------------- */

typedef struct ML_GGraph_Struct
{
   int  ML_id;
   int  Npoints, Nselected;
   int  ML_rank;
   int  *row_ptr, *col_ptr;
   int  send_cnt, *send_leng, *send_proc, **send_list;
   int  recv_cnt, *recv_leng, *recv_proc, **recv_list;
   char *bdry_type;
   char *vertex_state;

} ML_GGraph;

/* ******************************************************************** */
/* functions to manipulate the grid graph data structure                */
/* -------------------------------------------------------------------- */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

extern int  ML_GGraph_Create( ML_GGraph ** );
extern int  ML_GGraph_Destroy( ML_GGraph ** );
extern int  ML_GGraph_Print( ML_GGraph * );
extern int  ML_GGraph_Load_BdryTypes( ML_GGraph *, int , char *);
extern int  ML_GGraph_Coarsen(ML_GGraph*, ML_Comm *);
extern int  ML_GGraph_Gen_NodeGraph(ML_GGraph*,void*,void (*func),ML_Comm *);
extern int  ML_GGraph_Get_NodeStates(ML_GGraph*, int *, char **);
extern int  ML_GGraph_Gen_ElementGraph(ML_GGraph*,void*,void (*gf),ML_Comm*);
extern int  ML_GGraph_Gen_Restrictor(ML_GGraph*);
extern int ML_GGraph_CheckMIS( ML_GGraph *ml_gg, ML_Comm *comm );
extern int ML_GGraph_Find_NeighborElements(int leng1, int *list1, int leng2,
					   int *list2, int *vlist3);
extern int ML_GGraph_LabelVertices(int, int *, int, char *, char *, int,
                     int *, int *, int, int **, int, int **, int *, int *,
                     int, int **, int *, int *, int **, int, ML_Comm *);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
