#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "ha_const.h"
#include "hilbert_const.h"
#include "sfc_const.h"
#include "sfc.h"


void sfc_get_normed_coords(SFC_VERTEX_PTR sfc_vert_ptr, double min_bounding_box[], 
				  double max_bounding_box[], double normed_coords[], int num_dims)
{
  int i;
  double denominator, numerator;
  for(i=0;i<num_dims;i++) {
    numerator = (sfc_vert_ptr->coord[i]-min_bounding_box[i]);
    denominator = (max_bounding_box[i]-min_bounding_box[i]);
    if(numerator < denominator)
      normed_coords[i] = numerator/denominator;
    else
      normed_coords[i] = 1.0;
  }
  return;
}

/*---------------------------------------------------------------------------*/

void sfc_create_info(
  LB *lb,                       /* The load-balancing structure with info for
                                   the SFC balancer.                         */
  double min_bounding_box[],
  double max_bounding_box[],
  int num_dims,
  int num_local_objects,
  int wgt_dim,
  SFC_VERTEX_PTR sfc_vert_ptr,
  int sfc_keylength

)
{
  int i;
  unsigned unsigned_sfc_keylength = sfc_keylength;
  double normed_coords[3];
  if(num_dims == 2) {
    for(i=0;i<num_local_objects;i++) {
      sfc_get_normed_coords((sfc_vert_ptr+i), min_bounding_box,
			    max_bounding_box, normed_coords, 2);
      LB_fhsfc2d(normed_coords, &unsigned_sfc_keylength, sfc_vert_ptr[i].sfc_key);
    }
  }
  else {  /* if num_dims ==3 */
    for(i=0;i<num_local_objects;i++) {
      sfc_get_normed_coords((sfc_vert_ptr+i), min_bounding_box,
			    max_bounding_box, normed_coords, 3);
      LB_fhsfc3d(normed_coords, &unsigned_sfc_keylength, sfc_vert_ptr[i].sfc_key);      
    }
  }
  
  return;
}
