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
#ifndef lint
static char *cvs_gr_utils_c = "$Id$";
#endif

#include "all_const.h"
#include "vx_const.h"
#include "gr_const.h"
#include "gr_util_const.h"
#include "id_util_const.h"

void LB_print_graph(GRAPH *graph)
{
LOOP_CONTROL i; 
VERTEX *vertex;
int j, num_nbors;

  FOR_EACH_VERTEX(graph, vertex, i) {
    BL_ID_Util.Print_ID(&(vertex->Id));
    printf("  Weight:  %d  ", vertex->Weight);
    if (vertex->Coor != NULL) {
      printf("Coors:  (%f,%f,%f)", 
             vertex->Coor[0], vertex->Coor[1], vertex->Coor[2]);
    }
    printf("\n");
    num_nbors = vertex->Num_Nbors;
    printf("   Edges: ");
    for (j = 0; j < num_nbors; j++) {
      BL_ID_Util.Print_ID(&(vertex->Nbor_List[j]));
      printf("  ");
    }
    printf("\n   Edge Weights: ");
    if (vertex->Edge_Weights != NULL) {
      for (j = 0; j < num_nbors; j++) {
        printf("%f  ", vertex->Edge_Weights[j]);
      }
    }
    printf("\n\n");
  }
}
