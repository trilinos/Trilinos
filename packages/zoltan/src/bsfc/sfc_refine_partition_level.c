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
#include "sfc.h"


int sfc_refine_partition(ZZ *zz, int* local_balanced_flag, 
			 int *amount_of_used_bits, int num_vert_in_cut,
			 SFC_VERTEX_PTR vert_in_cut_ptr, int size_of_unsigned,
			 int wgt_dim, float* wgts_in_cut_ptr,
			 float* work_percent_array, float* total_weight_array,
			 float* global_actual_work_allocated, 
			 int number_of_cuts, int* max_cuts_in_bin,
			 int* ll_bins_head, float* work_prev_allocated,
			 int subbins_per_bin, int* local_balanced_flag_array,
			 int bin_refinement_method) 
{
  char yo[] = "sfc_refine_partition";
  int i=0, j=0, k;
  int amount_of_bits;
  float* binned_weight_array, *work_prev_allocated_copy;
  int* bin_proc_array;
  int* ll_prev_bins;
  int ll_counter, ll_location, *ll_bins_head_copy;

  /* amount of sub-bins in a bin, probably want this as a passed in parameter */
  int number_of_bins = subbins_per_bin;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* check to see that all of the bits of the sfc key 
     have not already been used */
  if(*amount_of_used_bits >= size_of_unsigned * SFC_KEYLENGTH * 8) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "No more refinement is possible.");
    *local_balanced_flag = SFC_BALANCED;
    return(ZOLTAN_OK);
  }
  
  /*  assume initially that all the partitions on this processor are balanced.
      we will check later on whether any are not balanced */
  *local_balanced_flag = SFC_BALANCED;

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
  if(amount_of_bits + *amount_of_used_bits > 8*size_of_unsigned * SFC_KEYLENGTH)
    amount_of_bits = 8*size_of_unsigned * SFC_KEYLENGTH - *amount_of_used_bits;
  number_of_bins = pow(2,i);
    
  ll_prev_bins = (int*) ZOLTAN_MALLOC(sizeof(int) * (number_of_cuts+1));
  if(ll_prev_bins == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }    
  
  ll_bins_head_copy = (int*) ZOLTAN_MALLOC(sizeof(int) * (number_of_cuts+1));
  if(ll_bins_head_copy == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  } 

  for(i=0;i<=number_of_cuts;i++)
    ll_bins_head_copy[i] = -1;
  
  work_prev_allocated_copy = 
    (float*) ZOLTAN_MALLOC(sizeof(float) *wgt_dim * (number_of_cuts+1));
  if(work_prev_allocated_copy == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  } 
  
  /* loop over all bins that have a cut in them using linklist 
     to find objects in the cut bins */
  for(ll_counter=0;ll_counter<=number_of_cuts;ll_counter++) 
    if((bin_refinement_method==1 ||
	local_balanced_flag_array[ll_counter]==SFC_NOT_BALANCED)
       && ll_bins_head[ll_counter] != -1) {

      int temp_max_cuts_in_bin = 0;
      
      /* calculate new bin numbers for objects that are in a cut bin */
      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {
	vert_in_cut_ptr[ll_location].my_bin = 
	  sfc_get_array_location(number_of_bins, amount_of_bits, 
				 *amount_of_used_bits, (vert_in_cut_ptr+ll_location));
	ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
      }  
      
      binned_weight_array = (float*) ZOLTAN_MALLOC(number_of_bins*wgt_dim*sizeof(float));
      if(binned_weight_array == NULL) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
	return(ZOLTAN_MEMERR);
      } 
      for(i=0;i<number_of_bins*wgt_dim;i++)
	binned_weight_array[i] = 0;
      
      /* fill up the weight array */
      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {  
	for(j=0;j<wgt_dim;j++)
	  binned_weight_array[wgt_dim*vert_in_cut_ptr[ll_location].my_bin+j] += 
	    wgts_in_cut_ptr[ll_location*wgt_dim+j];
	ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
      }
      
      bin_proc_array = (int*) ZOLTAN_MALLOC(sizeof(int) * (1+number_of_cuts));
      if(bin_proc_array == NULL) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
	return(ZOLTAN_MEMERR);
      } 
      bin_proc_array[0] = -1;
      for(i=1;i<=number_of_cuts;i++)
	bin_proc_array[i] = number_of_bins;
      
      /* find new cut(s) in the sub-bins */
      if(wgt_dim==1) {
	sfc_single_wgt_calc_partition(wgt_dim, work_prev_allocated[ll_counter],
				      total_weight_array, bin_proc_array, zz,
				      binned_weight_array, 
				      (work_percent_array+zz->Proc+ll_counter-2*number_of_cuts),
				      (global_actual_work_allocated+(zz->Proc-number_of_cuts*wgt_dim)),
				      number_of_bins, &temp_max_cuts_in_bin,
				      number_of_cuts, 0, NULL);

	if(temp_max_cuts_in_bin > *max_cuts_in_bin)
	  *max_cuts_in_bin = temp_max_cuts_in_bin;
      }
      else {
	/* multi_wgt_dim_calc_partition(); */  
	/*fill in with erik's multi-weight stuff when it is ready, use first weight for now */      
	sfc_single_wgt_calc_partition(wgt_dim, work_prev_allocated[ll_counter*wgt_dim],
				      total_weight_array, bin_proc_array, zz, 
				      binned_weight_array, 
				      (work_percent_array+zz->Proc+ll_counter-2*number_of_cuts),
				      (global_actual_work_allocated+(zz->Proc-number_of_cuts*wgt_dim)),
				      number_of_bins, &temp_max_cuts_in_bin, 
				      number_of_cuts, 0, NULL);

	if(temp_max_cuts_in_bin > *max_cuts_in_bin)
	  *max_cuts_in_bin = temp_max_cuts_in_bin;	
      }
      
      
      /* specify which processor an object belongs to,
	 we will know this because we know what bin an object 
	 belongs to and we know what processor a bin belongs to */
      
      /* initialize link list */
      for(i=0;i<=number_of_cuts;i++)
	ll_prev_bins[i] = -1;

      ll_location = ll_bins_head[ll_counter];
      while(ll_location != -1) {
	int non_cut_bin_counter = 0;  /* keeps track of how many bins between
					 number_of_cuts and j do not have 
					 a cut in them */
	/* bad search method but easy, should probably replace */
	j=number_of_cuts;
	while((int) vert_in_cut_ptr[ll_location].my_bin < bin_proc_array[j]) {
	  j--;
	  if(bin_proc_array[j] == number_of_bins)
	    non_cut_bin_counter++;
	}
	vert_in_cut_ptr[ll_location].destination_proc += 
	  j-number_of_cuts; 

	/* if this object is in a bin with a cut... */
	if((int) vert_in_cut_ptr[ll_location].my_bin == bin_proc_array[j]) {
	  if(ll_bins_head_copy[number_of_cuts-(zz->Proc)+vert_in_cut_ptr[ll_location].destination_proc] != -1) 
	    vert_in_cut_ptr[ll_prev_bins[j]].next_sfc_vert_index = ll_location;
	  else {
	    ll_bins_head_copy[number_of_cuts-(zz->Proc)+vert_in_cut_ptr[ll_location].destination_proc] = ll_location;
	    /* calculate work_prev_allocated for this new partition */
	    for(i=0;i<wgt_dim;i++)
	      work_prev_allocated_copy /* array continued on next line */
		[(number_of_cuts-(zz->Proc)+vert_in_cut_ptr[ll_location].destination_proc)*wgt_dim+i] = 
		work_prev_allocated[i+wgt_dim*ll_counter];
	    for(i=bin_proc_array[j]+1;i<number_of_bins;i++)
	      for(k=0;k<wgt_dim;k++) 
		work_prev_allocated_copy /* array continued on next line */
		  [(number_of_cuts-(zz->Proc)+vert_in_cut_ptr[ll_location].destination_proc)*wgt_dim+k] 
		  += binned_weight_array[i*wgt_dim+k];
	  }
	  
	  ll_prev_bins[j] = ll_location;
	}
	ll_location = vert_in_cut_ptr[ll_location].next_sfc_vert_index;
      }
      /* end the linklists created above */
      for(i=0;i<=number_of_cuts;i++)
	if(ll_prev_bins[i] != -1)
	  vert_in_cut_ptr[ll_prev_bins[i]].next_sfc_vert_index = -1;

      ZOLTAN_FREE(&binned_weight_array);
      ZOLTAN_FREE(&bin_proc_array);
    }
  
  for(i=0;i<=number_of_cuts;i++) {
    ll_bins_head[i] = ll_bins_head_copy[i];
    for(j=0;j<wgt_dim;j++)
      work_prev_allocated[i*wgt_dim+j] = work_prev_allocated_copy[i*wgt_dim+j];
  }
  
  ZOLTAN_FREE(&ll_prev_bins);
  ZOLTAN_FREE(&ll_bins_head_copy);
  ZOLTAN_FREE(&work_prev_allocated_copy);
  

  *amount_of_used_bits += amount_of_bits;
  
  /* check which partitions that are not balanced */
  for(i=0;i<=number_of_cuts;i++) {
    if(ll_bins_head[i] != -1 || local_balanced_flag_array[i] != SFC_BALANCED)
      {
	if(wgt_dim == 1) {
	  local_balanced_flag_array[i] =
	    sfc_single_wgt_find_imbalance(work_percent_array, 
					  global_actual_work_allocated[(zz->Proc+i-number_of_cuts)*wgt_dim],
					  total_weight_array[0], zz->Proc+i-number_of_cuts, zz);
	}
	else {
	  /* put in multi-dimensional algorithm to calculate imbalance of the partitions */
	  local_balanced_flag_array[i] =
	    sfc_single_wgt_find_imbalance(work_percent_array, 
					  global_actual_work_allocated[(zz->Proc+i-number_of_cuts)*wgt_dim],
					  total_weight_array[0], zz->Proc+i-number_of_cuts, zz);
	}
      }
  }
  
  /* check if any of the partitions are not balanced */
  *local_balanced_flag = SFC_BALANCED;
  i=0;
  while(*local_balanced_flag == SFC_BALANCED && i<=number_of_cuts) {
    *local_balanced_flag = local_balanced_flag_array[i];
    i++;
  }
  
  /* check the partitions to see if any more improvement can be made on them */
  if(*local_balanced_flag == SFC_NOT_BALANCED) {
    for(i=0;i<=number_of_cuts;i++)
      if(ll_bins_head[i] != -1) {
	/* check if there is only 1 object in this bin. if there is, 
	   no further bin refinement will improve load-balance */
	if(vert_in_cut_ptr[ll_bins_head[i]].next_sfc_vert_index == -1) {
	  ll_bins_head[i] = -1;
	  local_balanced_flag_array[i] = SFC_BALANCED;
	  ZOLTAN_PRINT_WARN(zz->Proc, yo, 
			"Bin refinement cannot improve load balance on this processor.");
	}
	/* check if the objects in the bin have all the same sfc_key */
	else {
	  unsigned sfc_key[SFC_KEYLENGTH];
	  int same_flag = 0;  /* flag to indicate if all the sfc_keys are the same */
	  int amount_of_objects_in_bin = 0;
	  ll_counter = ll_bins_head[i];
	  for(j=0;j<SFC_KEYLENGTH;j++)
	    sfc_key[j] = vert_in_cut_ptr[ll_counter].sfc_key[j];
	  ll_counter = vert_in_cut_ptr[ll_counter].next_sfc_vert_index;
	  while(ll_counter != -1 && same_flag == 0) {
	    for(j=0;j<SFC_KEYLENGTH;j++)
	      if(vert_in_cut_ptr[ll_counter].sfc_key[j] != sfc_key[j])
		same_flag = 1;
	    
	    amount_of_objects_in_bin++;
	    ll_counter = vert_in_cut_ptr[ll_counter].next_sfc_vert_index;
	  }
	  if(ll_counter == -1 && same_flag == 0) {
	    /* 
	       all of the objects in this bin have the same sfc_key, options:
	       1. stop refinement of the bin (no improvement of the partition can be obtained)
	       2. create new sfc_keys for the objects and continue refinement
	       note: acbauer has chosen option 2
	    */
	    unsigned umax;
	    if(size_of_unsigned == sizeof(unsigned))
	      umax = ~(0u);
	    else 
	      umax = pow(2,size_of_unsigned) - 1;
	    ll_counter = ll_bins_head[i];
	    j=1;
	    while(ll_counter != -1) {
	      unsigned new_key = umax*((float) j/(float) 10*amount_of_objects_in_bin);
	      for(k=0;k<SFC_KEYLENGTH;k++)
		vert_in_cut_ptr[ll_counter].sfc_key[k] = new_key;
	      j++;
	      ll_counter = vert_in_cut_ptr[ll_counter].next_sfc_vert_index; 
	    }
	    ZOLTAN_PRINT_WARN(zz->Proc, yo, 
			  "All objects in a bin that needs to be refined have the same sfc key.");
	  }
	}
      }
    /* check again if any of the partitions are not balanced */
    *local_balanced_flag = SFC_BALANCED;
    j=0;
    while(*local_balanced_flag == SFC_BALANCED && j<=number_of_cuts) {
      *local_balanced_flag = local_balanced_flag_array[j];
      j++;
    }
  }
  
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}


