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
static char *cvs_vx_c = "$Id$";
#endif

#include "all_const.h"
#include "vx_const.h"
#include "id_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

VERTEX *LB_new_vertex(
  ID *id,                    /* ID to be assigned to the new vertex.         */
  int weight,                /* weight to be assigned to the new vertex.     */
  int num_nbors,             /* degree of the new vertex.                    */
  ID *nbor_list,             /* neighbor nodes of the new vertex.            */
  float *edge_weights,       /* edge weights for each edge, if appropriate.  */
  int num_dim,               /* number of dimensions in the problem; 
                                num_dim == 0 if geometry info is not used.   */
  double *coor               /* coordinates of the vertices.                 */
)
{
VERTEX *vertex;
int i;
static int size_vertex = sizeof(VERTEX);
int size_nbor_list = num_nbors * sizeof(ID);
int size_edge_weights = num_nbors * sizeof(float);
int size_coords = MAX_DIM * sizeof(double);
int size;
char *tmp;

  /*
   *  Compute size of entire vertex structure, including the storage for
   *  the neighbor lists, edge weights (if appropriate) and geometry info
   *  (if appropriate). 
   */

  printf("Adding vertex "); BL_ID_Util.Print_ID(id); printf("\n");

  size = size_vertex + size_nbor_list;
  if (edge_weights != NULL)
    size += size_edge_weights;
  if (num_dim > 0)
    size += size_coords;

  /*
   *  Allocate the storage contiguously and set pointers to appropriate 
   *  fields.
   */

  tmp = (char *) LB_MALLOC(size);
  vertex = (VERTEX *) tmp;
  tmp += size_vertex;
  vertex->Nbor_List = (ID *) (tmp);
  tmp += size_nbor_list;
  if (edge_weights != NULL) {

    /* 
     *  Assign pointer for array of edge weights.
     */

    vertex->Edge_Weights = (float *) (tmp);
    tmp += size_edge_weights;
  }
  else
    vertex->Edge_Weights = NULL;

  if (num_dim > 0) {
   
    /*
     *  Assign pointer for geometry info.
     */
   
    vertex->Coor = (double *) (tmp);
    tmp += size_coords;
  }
  else
    vertex->Coor = NULL;

  /*
   *  Initialize the vertex with the information provided by the user.
   */

  BL_ID_Util.Assign(&(vertex->Id), id);
  vertex->Weight = weight;
  vertex->Num_Nbors = num_nbors;
  for (i = 0; i < num_nbors; i++) {
    BL_ID_Util.Assign(&(vertex->Nbor_List[i]), &(nbor_list[i]));
    if (edge_weights != NULL) 
      vertex->Edge_Weights[i] = edge_weights[i];
  }

  if (num_dim > 0)
    for (i = 0; i < MAX_DIM; i++)
      vertex->Coor[i] = coor[i];

  return(vertex);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_free_vertex(VERTEX **vertex)
{
  LB_FREE(vertex);
}
