#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "ha_const.h"
#include "hilbert_const.h"
#include "sfc.h"

/* routine returns the normed coordinates */
void sfc_get_normed_coords(double min_bounding_box[],
			   double max_bounding_box[],
			   double normed_coords[],
			   int num_dims, double my_coords[])
{
  int i;
  double denominator, numerator;
  for(i=0;i<num_dims;i++) {
    numerator = (my_coords[i]-min_bounding_box[i]);
    denominator = (max_bounding_box[i]-min_bounding_box[i]);
    if(numerator < denominator)
      normed_coords[i] = numerator/denominator;
    else
      normed_coords[i] = 1.0;
  }
 
  return;
}

/*---------------------------------------------------------------------------*/
/* routine gets the normed coordinates and then calculates the sfc key of 
   an object */

void sfc_create_info(
  ZZ *zz,                       /* The Zoltan structure with info for
                                   the SFC balancer.                         */
  double min_bounding_box[],
  double max_bounding_box[],
  int num_dims,
  int num_local_objects,
  int wgt_dim,
  SFC_VERTEX_PTR sfc_vert_ptr,
  double coords[]
)
{
  char yo[] = "sfc_create_info";
  int i, j;
  unsigned unsigned_sfc_keylength = SFC_KEYLENGTH;
  double normed_coords[3];

  ZOLTAN_TRACE_ENTER(zz, yo);
  if(num_dims == 2) {
    for(i=0;i<num_local_objects;i++) {
      for(j=0;j<SFC_KEYLENGTH;j++)
	sfc_vert_ptr[i].sfc_key[j] = 0;
      sfc_get_normed_coords(min_bounding_box, max_bounding_box,
			    normed_coords, 2, (coords+i*num_dims));
      Zoltan_SFC_fhsfc2d(normed_coords, &unsigned_sfc_keylength, sfc_vert_ptr[i].sfc_key);
    }
  }
  else {  /* if num_dims ==3 */
    for(i=0;i<num_local_objects;i++) {
      sfc_get_normed_coords(min_bounding_box, max_bounding_box, 
			    normed_coords, 3, (coords+i*num_dims));
      Zoltan_SFC_fhsfc3d(normed_coords, &unsigned_sfc_keylength, sfc_vert_ptr[i].sfc_key);      
    }
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return;
}
