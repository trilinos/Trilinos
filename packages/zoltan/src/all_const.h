/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#ifndef __ALL_CONST_H
#define __ALL_CONST_H

#ifndef lint
static char *cvs_all_const_h = "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "all_allo.h"
#include "par_const.h"

/****************************************************************************
 *  Definitions
 ****************************************************************************/
/*
 *  Graph Iterator Macro -- this macro provides a loop over all vertices in
 *  a graph.  It uses the public graph functions First_Vertex and Next_Vertex.
 *  It is used in place of a "for" statement.
 *  This approach was taken over the more object-oriented approach of 
 *  passing a function (performing the body of the loop) and a variable 
 *  argument list to an iterator function, because in the object-oriented
 *  approach, argument checking could not be done for the body function.
 */

#define FOR_EACH_VERTEX(graph, vertex, i) \
        for (vertex = graph->First_Vertex(graph, &i); vertex != NULL; \
             vertex = graph->Next_Vertex(graph, &i))

/****************************************************************************
 *  Typedefs for important structures.
 ****************************************************************************/

typedef struct Vertex_Struct VERTEX;
typedef struct Graph_Struct GRAPH;
typedef void   GRAPH_DATA;
typedef struct ID_Struct ID;
typedef struct ID_Util_Struct ID_UTIL;
typedef void   ID_NEW_FN(ID *, int);
typedef int    ID_INT_FN(ID *, ID *);
typedef void   ID_VOID_FN(ID *, ID *);
typedef void   ID_PRINT_FN(ID *);
typedef void  *LOOP_CONTROL;

typedef enum Boolean_Type {
  FALSE,
  TRUE 
} BOOLEAN;


/****************************************************************************
 *  Prototypes
 ****************************************************************************/

extern GRAPH *new_graph();
extern GRAPH *build_graph_ala_chaco(int, int *, int *, ID *, 
                                    float *, int, float *, float *, float *);

extern VERTEX *new_vertex(ID *, int, int, ID *, float *, int,
                          double *);

extern void free_vertex(VERTEX **);
#endif
