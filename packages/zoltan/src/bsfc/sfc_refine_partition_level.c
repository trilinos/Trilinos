/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"
#include <values.h>
#include <limits.h>
#include "hilbert_const.h"
#include "sfc_const.h"
#include "sfc.h"


int sfc_refine_partition_level(LB* lb, int* local_balanced_flag, int *amount_of_used_bits,
			       int num_vert_in_cut, SFC_VERTEX_PTR vert_in_cut_ptr,
			       int sfc_keylength, int size_of_unsigned, unsigned imax, int wgt_dim,
			       float* wgts_in_cut_ptr, float* work_percent_array,
			       float* total_weight_array, float* global_actual_work_allocated,
			       int number_of_cuts, int* max_cuts_in_bin, int* ll_bins_head)
{
  int i=0, j=0;
  int amount_of_bits;
  float* binned_weight_array, *summed_binned_weight_array;
  float* work_prev_allocated;
  int* bin_proc_array;
  int* ll_prev_bins;
  int ll_counter, ll_location;

  /* amount of sub-bins in a bin, probably want this as a passed in parameter */
  int number_of_bins = 10;

  printf("refining some bins on proc %d!!!\n",lb->Proc);

  /* if there are a lot of cuts in a bin, we want the amount of bins
     to be greater than the amount of cuts */
  if(*max_cuts_in_bin >= number_of_bins)
    number_of_bins = *max_cuts_in_bin + 1;

  /*increase sub-bins so that there is a power of 2 */
  
  while(number_of_bins > pow(2,i))
    i++;
  amount_of_bits = i;
  number_of_bins = pow(2,i);
    
  /* loop over all bins that have a cut in them using linklist to find objects in the cut bins */
  for(ll_counter=0;ll_counter<number_of_cuts;ll_counter++) {

  /* calculate new bin numbers for objects that are in a cut bin */
    ll_location = ll_bins_head[ll_counter];
    while(ll_location != -1) {
      vert_in_cut_ptr[ll_location].my_bin = get_array_location(number_of_bins, amount_of_bits,
							       *amount_of_used_bits, (vert_in_cut_ptr+ll_location), 
							       sfc_keylength, size_of_unsigned, imax);
      ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
    }  
  
    binned_weight_array = (float*) LB_MALLOC(number_of_bins*wgt_dim*sizeof(float));
    for(i=0;i<number_of_bins*wgt_dim;i++)
      binned_weight_array[i] = 0;
    
    /* fill up the weight array */
    ll_location = ll_bins_head[ll_counter];
    while(ll_location != -1) {  
      for(j=0;j<wgt_dim;j++)
	binned_weight_array[vert_in_cut_ptr[ll_location].my_bin+j] += wgts_in_cut_ptr[ll_location*wgt_dim+j];
      ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
    }
    
    summed_binned_weight_array = (float*) LB_MALLOC(sizeof(float) * wgt_dim);
    for(i=0;i<wgt_dim;i++)
      summed_binned_weight_array[i] = 0;
    for(i=0;i<number_of_bins;i++)
      for(j=0;j<wgt_dim;j++)
	summed_binned_weight_array[j] += binned_weight_array[i*wgt_dim+j];
    
    bin_proc_array = (int*) LB_MALLOC(sizeof(int) * number_of_bins);
    bin_proc_array[0] = 0;
    for(i=1;i<number_of_bins;i++)
      bin_proc_array[i] = number_of_bins;
    work_prev_allocated = (float*) LB_MALLOC(sizeof(float) * wgt_dim);
    /* update work previously allocated to include work in all bins with higher keys than this bin */
    for(i=0;i<wgt_dim;i++)
      work_prev_allocated[i] = global_actual_work_allocated[(lb->Proc)*wgt_dim+i] -
	summed_binned_weight_array[i];
    
    LB_FREE(&summed_binned_weight_array);
    /* find new cut(s) in the sub-bins */
    if(wgt_dim==1) 
      single_wgt_calc_partition(wgt_dim, work_prev_allocated[0], total_weight_array, 
				bin_proc_array, lb, binned_weight_array, 
				(work_percent_array+lb->Proc- number_of_cuts),
				global_actual_work_allocated, number_of_bins, max_cuts_in_bin,
				number_of_cuts, 0);
    else {
      /* multi_wgt_dim_calc_partition(); */  
      /*fill in with erik's multi-weight stuff when it is ready, use first weight for now */      
      
    }
    
    ll_prev_bins = (int*) LB_MALLOC(sizeof(int) * number_of_cuts);
    
    /* specify which processor an object belongs to,
       we will know this because we know what bin an object 
       belongs to and we know what processor a bin belongs to */
    
    ll_location = ll_bins_head[ll_counter];
    while(ll_location != -1) {
      /* bad search method but easy, should probably replace */
      j=number_of_cuts;
      while((int) vert_in_cut_ptr[ll_location].my_bin < bin_proc_array[j])
	j--;
      vert_in_cut_ptr[ll_location].destination_proc = lb->Proc+j- number_of_cuts;
      if((int) vert_in_cut_ptr[ll_location].my_bin != bin_proc_array[j])
	vert_in_cut_ptr[ll_location].cut_bin_flag = SFC_NO_CUT;
      else {
	vert_in_cut_ptr[ll_location].cut_bin_flag = SFC_CUT; 
	/*	if(ll_bins_head[j] != -1) 
		vert_in_cut_ptr[ll_prev_bins[j]].next_sfc_vert_index = i;
		else
		ll_bins_head[j] = i;
		
		ll_prev_bins[j] = i;*/
      }
      ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
    }
  
    LB_FREE(&ll_prev_bins);
    LB_FREE(&binned_weight_array);
    LB_FREE(&bin_proc_array);
    LB_FREE(&work_prev_allocated);
  }

  *amount_of_used_bits += amount_of_bits;

  /* check if there are any partitions that are not balanced */
  i = 0;
  *local_balanced_flag = SFC_BALANCED;
  while(*local_balanced_flag == SFC_BALANCED && i < number_of_cuts) {
    *local_balanced_flag=single_wgt_find_imbalance(work_percent_array, 
						   global_actual_work_allocated[(lb->Proc-i)*wgt_dim],
						   total_weight_array[0], lb->Proc-i, lb);
    i++;
  }
  
  return LB_OK;
}


