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
			       int number_of_cuts, int* max_cuts_in_bin, int* ll_bins_head,
			       float* work_prev_allocated) 
{
  int i=0, j=0, k;
  int amount_of_bits;
  float* binned_weight_array, *work_prev_allocated_copy;
  int* bin_proc_array;
  int* ll_prev_bins;
  int ll_counter, ll_location, *ll_bins_head_copy;

  /* amount of sub-bins in a bin, probably want this as a passed in parameter */
  int number_of_bins = 1; /* should equal user specified parameter sub_bins_per_bin */

  printf("refining some bins on proc %d!!!\n",lb->Proc);

  /* if there are a lot of cuts in a bin, we want the amount of bins
     to be greater than the amount of cuts */
  if(*max_cuts_in_bin >= number_of_bins)
    number_of_bins = *max_cuts_in_bin + 1;

  *max_cuts_in_bin = 0;

  /*increase sub-bins so that there is a power of 2 */
  i=0;
  while(number_of_bins > pow(2,i))
    i++;
  amount_of_bits = i;
  number_of_bins = pow(2,i);
    
  ll_prev_bins = (int*) LB_MALLOC(sizeof(int) * (number_of_cuts+1));
  
  ll_bins_head_copy = (int*) LB_MALLOC(sizeof(int) * (number_of_cuts+1));
  for(i=0;i<=number_of_cuts;i++)
    ll_bins_head_copy[i] = -1;
  
  work_prev_allocated_copy = (float*) LB_MALLOC(sizeof(float) *wgt_dim * (number_of_cuts+1));
  
  /* loop over all bins that have a cut in them using linklist to find objects in the cut bins */
  for(ll_counter=0;ll_counter<=number_of_cuts;ll_counter++) 
    if(ll_bins_head[ll_counter] != -1) {
      
      /* calculate new bin numbers for objects that are in a cut bin */
      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {
	vert_in_cut_ptr[ll_location].my_bin = 
	  get_array_location(number_of_bins, amount_of_bits,
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
	  binned_weight_array[vert_in_cut_ptr[ll_location].my_bin+j] += 
	    wgts_in_cut_ptr[ll_location*wgt_dim+j];
	ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
      }

      for(i=0;i<number_of_bins;i++)
	printf("weight[%d] = %e ************* proc = %d\n",i, 
	       binned_weight_array[i], lb->Proc); 
      printf("@@@@@@@@@@ work prev allocated is %e on proc %d for ll_counter %d\n",
	     work_prev_allocated[ll_counter], lb->Proc, ll_counter);  
      
      bin_proc_array = (int*) LB_MALLOC(sizeof(int) * (1+number_of_cuts));
      bin_proc_array[0] = -1;
      for(i=1;i<=number_of_cuts;i++)
	bin_proc_array[i] = number_of_bins;
      
      /* find new cut(s) in the sub-bins */
      if(wgt_dim==1) {
	int temp_max_cuts_in_bin = 0;
	single_wgt_calc_partition(wgt_dim, work_prev_allocated[ll_counter], total_weight_array, 
				  bin_proc_array, lb, binned_weight_array, 
				  (work_percent_array+lb->Proc+ll_counter-2*number_of_cuts),
				  global_actual_work_allocated, number_of_bins, 
				  &temp_max_cuts_in_bin, number_of_cuts, 0);
	if(temp_max_cuts_in_bin > *max_cuts_in_bin)
	  *max_cuts_in_bin = temp_max_cuts_in_bin;
      }
      else {
	/* multi_wgt_dim_calc_partition(); */  
	/*fill in with erik's multi-weight stuff when it is ready, use first weight for now */      
	
      }
      
      
      /* specify which processor an object belongs to,
	 we will know this because we know what bin an object 
	 belongs to and we know what processor a bin belongs to */
      
      /* initialize link list */
      for(i=0;i<=number_of_cuts;i++)
	ll_prev_bins[i] = -1;

      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {
	int non_cut_bin_counter = 0;  /* keeps track of how many bins between number_of_cuts and j
					 do not have a cut in them */
	/* bad search method but easy, should probably replace */
	j=number_of_cuts;
	while((int) vert_in_cut_ptr[ll_location].my_bin < bin_proc_array[j]) {
	  j--;
	  if(bin_proc_array[j] == number_of_bins)
	    non_cut_bin_counter++;
	}
	vert_in_cut_ptr[ll_location].destination_proc += 
	  j-number_of_cuts+non_cut_bin_counter;

	if(vert_in_cut_ptr[ll_location].destination_proc < 0)
	  printf("bad proc right here!!!!!!! i'm proc %d >>>>>>>\n",lb->Proc);

	/* if this object is in a bin with a cut... */
	if((int) vert_in_cut_ptr[ll_location].my_bin == bin_proc_array[j]) {
	  if(ll_bins_head_copy[number_of_cuts-(lb->Proc)+vert_in_cut_ptr[ll_location].destination_proc] != -1) 
	    vert_in_cut_ptr[ll_prev_bins[j]].next_sfc_vert_index = ll_location;
	  else {
	    ll_bins_head_copy[number_of_cuts-(lb->Proc)+vert_in_cut_ptr[ll_location].destination_proc] = ll_location;
	    /* calculate work_prev_allocated for this new partition */
	    for(i=0;i<wgt_dim;i++)
	      work_prev_allocated_copy[number_of_cuts-(lb->Proc)+vert_in_cut_ptr[ll_location].destination_proc*wgt_dim+i] = 
		work_prev_allocated[i+wgt_dim*ll_counter];
	    for(i=bin_proc_array[j]+1;i<number_of_bins;i++)
	      for(k=0;k<wgt_dim;k++)
		work_prev_allocated_copy[number_of_cuts-(lb->Proc)+vert_in_cut_ptr[ll_location].destination_proc*wgt_dim+k] += binned_weight_array[i*wgt_dim+k];
	  }
	  
	  ll_prev_bins[j] = ll_location;
	}
	ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
      }
      /* end the linklists created above */
      for(i=0;i<=number_of_cuts;i++)
	if(ll_prev_bins[i] != -1)
	  vert_in_cut_ptr[ll_prev_bins[i]].next_sfc_vert_index = -1;

      LB_FREE(&binned_weight_array);
      LB_FREE(&bin_proc_array);
    }
  

  
  for(i=0;i<=number_of_cuts;i++) {
    ll_bins_head[i] = ll_bins_head_copy[i];
    for(j=0;j<wgt_dim;j++)
      work_prev_allocated[i*wgt_dim+j] = work_prev_allocated_copy[i*wgt_dim+j];
  }
  
  LB_FREE(&ll_prev_bins);
  LB_FREE(&ll_bins_head_copy);
  LB_FREE(&work_prev_allocated_copy);
  

  *amount_of_used_bits += amount_of_bits;
  
  /* check if there are any partitions that are not balanced */
  i = 0;
  *local_balanced_flag = SFC_NOT_BALANCED;
/*  while(*local_balanced_flag == SFC_BALANCED && i <= number_of_cuts) {
    *local_balanced_flag=single_wgt_find_imbalance(work_percent_array, 
						   global_actual_work_allocated[(lb->Proc-i)*wgt_dim],
						   total_weight_array[0], lb->Proc-i, lb);
    i++;
  }*/
  
  return LB_OK;
}


