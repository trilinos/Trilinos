#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "lb_util_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "comm_const.h"
#include <values.h>
#include <limits.h>
#include "hilbert_const.h"
#include "sfc_const.h"
#include "sfc.h"

int sfc_refine_overloaded_bins(LB* lb, int max_cuts_in_bin, int max_number_of_cuts,
			       int* number_of_cuts_in_bin, int number_of_bins, int wgt_dim,
			       SFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[],
			       int* amount_of_bits_used, int num_local_objects)
{
  int i, j, k, l, m, ierr, *istore, *istore2, *bins_to_refine;
  int number_of_bits, number_of_bins_to_refine;
  int gl_max_number_of_cuts, *refine_flag;

  int* overloaded_bin_flag = (int*) LB_MALLOC(sizeof(int) * lb->Num_Proc);
  for(i=0;i<lb->Num_Proc;i++)
    overloaded_bin_flag[i] = 0;
  j=0;
  for(i=0;i<number_of_bins;i++)
    if(number_of_cuts_in_bin[i] > max_cuts_in_bin)
      j++;
  if(j>1)
    overloaded_bin_flag[lb->Proc] = -j;
  else if(j==1)
    overloaded_bin_flag[lb->Proc] = max_number_of_cuts;
  
  istore = (int*) LB_MALLOC(sizeof(int) * lb->Num_Proc);
  
  ierr = MPI_Allreduce(overloaded_bin_flag, istore, lb->Num_Proc, MPI_INT, 
		       MPI_SUM, lb->Communicator);
  
  LB_FREE(&overloaded_bin_flag);
  
  number_of_bins_to_refine=0;
  for(i=0;i<lb->Num_Proc;i++) {
    if(istore[i] <= 0)
      number_of_bins_to_refine += -istore[i];
    else
      number_of_bins_to_refine += 1;
  }  
  if(number_of_bins_to_refine==0) {
    /* all of the bins are okay! */
    LB_FREE(&istore);
    return LB_OK;
  }
  /* if we need to refine some coarse bins... */
  istore2 = (int*) LB_MALLOC(sizeof(int) * number_of_bins_to_refine);
  for(i=0;i<number_of_bins_to_refine;i++)
    istore2[i] = 0;
  k=0;
  for(i=0;i<lb->Num_Proc;i++) {
    if(i==lb->Proc && istore[i] != 0) {
      if(istore[i] > 0) {
	l=0;
	while(number_of_cuts_in_bin[l] == 0)
	  l++;
	istore2[k] = lb->Proc*number_of_bins+l;
      }
      else {
	m=0;
	for(l=0;l<number_of_bins;l++)
	  if(number_of_cuts_in_bin[l] != 0) {
	    istore2[k+m] = lb->Proc*number_of_bins+l;
	    m++;
	  }
      }
    }
    if(istore[i] <= 0)
      k += -istore[i];
    else 
      k++;
  }
  /* calculate how many subbins to break a bin into */
  ierr = MPI_Allreduce(&max_number_of_cuts, &gl_max_number_of_cuts, 1,
		       MPI_INT, MPI_SUM, lb->Communicator);

  i=0;
  while(pow(2,i) < (gl_max_number_of_cuts+1))
    i++;
  number_of_bins = pow(2,i);
  number_of_bits = i;

  bins_to_refine = (int*) LB_MALLOC(sizeof(int) * number_of_bins_to_refine);
  ierr = MPI_Allreduce(istore2, bins_to_refine, number_of_bins_to_refine, 
		       MPI_INT, MPI_SUM, lb->Communicator);

  LB_FREE(&istore2);
  refine_flag = (int*) LB_MALLOC(sizeof(int) * num_local_objects);
  for(i=0;i<num_local_objects;i++)
    refine_flag[i] = -1;
  for(i=0;i<num_local_objects;i++)
    for(j=0;j<number_of_bins_to_refine;j++) 
      if((int) sfc_vert_ptr[i].my_bin == bins_to_refine[j])
	refine_flag[i] = bins_to_refine[j];



  for(i=0;i<number_of_bins_to_refine;i++) {

    

  }  
      
  
  LB_FREE(&refine_flag);
  LB_FREE(&istore);
  
  return LB_OK;
}

int refine_coarse_bin(LB* lb, int num_local_objects, SFC_VERTEX_PTR sfc_vert_ptr, 
		      float objs_wgt[], int wgt_dim, int refined_bin, int proc,
		      int sfc_keylength, int number_of_bins, int number_of_bits, 
		      int prev_used_bits, int size_of_unsigned, unsigned imax,
		      int* number_of_cuts, int* refined_flag, int current_loc,
		      float* work_percent_array, float* total_weight_array, 
		      float* actual_work_allocated, float work_prev_allocated,
		      int* refine_flag) 
{
  int i, j, counter = 0, *bin_proc_array, ierr;
  float *binned_wgt_array, *summed_binned_wgt_array = NULL;
  
  binned_wgt_array = (float*) LB_MALLOC(sizeof(float) * number_of_bins * wgt_dim);
  for(i=0;i<number_of_bins*wgt_dim;i++)
    binned_wgt_array[i] = 0;
  for(i=0;i<num_local_objects; i++) {
    if(sfc_vert_ptr[i].my_bin != refined_bin)
      refine_flag[i] = -1;
    else {
      refine_flag[i] = get_array_location(number_of_bins, number_of_bits, prev_used_bits,
					   (sfc_vert_ptr+i), sfc_keylength, size_of_unsigned,
					   imax);
      for(j=0;j<wgt_dim;j++)
	binned_wgt_array[j+wgt_dim*refined_flag[i]] += objs_wgt[i*wgt_dim+j];
    counter++;
    }
  }
  if(proc == lb->Proc) 
    summed_binned_wgt_array = (float*) LB_MALLOC(sizeof(float) * number_of_bins * wgt_dim);
  ierr = MPI_Reduce(binned_wgt_array, summed_binned_wgt_array, number_of_bins*wgt_dim, 
		    MPI_FLOAT, MPI_SUM, proc, lb->Communicator);

  /* only proc figures out the cuts */
  if(proc == lb->Proc) {
    
    if(wgt_dim == 1)
      single_wgt_calc_partition(wgt_dim, work_prev_allocated, total_weight_array,
				bin_proc_array, lb, summed_binned_wgt_array, 
				work_percent_array, actual_work_allocated, number_of_bins, 
				number_of_cuts, current_loc, 0, NULL);
    else {
      /* multi_wgt_dim_calc_partition(); */  
      /*fill in with erik's multi-weight stuff when it is ready, use first weight for now */
      single_wgt_calc_partition(wgt_dim, work_prev_allocated, total_weight_array,
				bin_proc_array, lb, summed_binned_wgt_array, 
				work_percent_array, actual_work_allocated, number_of_bins, 
				number_of_cuts, current_loc, 0, NULL);
    }
      
	


  }


  if(proc == lb->Proc)
    LB_FREE(&summed_binned_wgt_array);

  return LB_OK;
}
