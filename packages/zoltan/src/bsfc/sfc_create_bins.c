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
#include "lb_util_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "comm_const.h"
#include <values.h>
#include <limits.h>
#include "hilbert_const.h"
#include "sfc_const.h"
#include "sfc.h"


 /* global_actual_work_allocated is used to make sure that each processor 
    knows how much extra work every processor has allocated
    (extra work[proc] = global_actual_work_allocated[proc] - work_percent_array[proc]*total_work */
int sfc_create_bins(LB* lb, int num_local_objects, 
		    int wgt_dim, SFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[], int* amount_of_bits_used,
		    int sfc_keylength, int size_of_unsigned, unsigned imax, 
		    float* global_actual_work_allocated, float *work_percent_array, 
		    float* total_weight_array, int* balanced_flag,
		    SFC_VERTEX_PTR *vert_in_cut_ptr, float** wgts_in_cut_ptr, 
		    int* num_vert_in_cut, int* number_of_cuts, int bins_per_proc, 
		    int hashtable_divider, COMM_OBJ **plan, int* num_vert_sent)
{
  char    yo[] = "sfc_create_bins";
  int i, j, number_of_bins, ierr = 0;
  int array_location = 0;
  int comm_tag = 4190; 
  int * proclist;
  int nreturn = 0;
  int off_proc_objects = 0;  /*counter to keep track of how many objects will be off processor*/
  float * binned_weight_array;
  SFC_HASH_OBJ_PTR * sfc_hash_ptr;
  SFC_HASH_OBJ_PTR extra_hash_ptr;
  int hashtable_length;
  int counter = 0;
  SFC_BIN_WEIGHT_PTR send_buffer, rcv_buffer;
  float *extra_float_array;
  float my_work_percent;
  int *bin_proc_array;
  float *scanned_work_prev_allocated; /*scanned_work_prev_allocated is the amount of work allocated to higher ranked procs */
  float *actual_work_allocated;
  int* global_bin_proc_array;
  int amount_of_bits;
  SFC_VERTEX_PTR send_vert_buffer;
  float* send_wgt_buffer;
  int current_proc;
  float* extra_float_array2 = NULL;
  int local_balanced_flag;
  

  binned_weight_array = (float *) LB_MALLOC(sizeof(float) * 2 * bins_per_proc * wgt_dim);
  for(i=0;i<2*bins_per_proc*wgt_dim;i++)
    binned_weight_array[i] = 0;

  /*assume initially that each processor has the same amount of bins*/
  number_of_bins = lb->Num_Proc * bins_per_proc;
  i=0;
  while(number_of_bins > pow(2,i))
    i++;
  amount_of_bits = i;
  number_of_bins = pow(2,i);
  *amount_of_bits_used = amount_of_bits;

  /*hash table */
  if(hashtable_divider<1) /* hashtable_divider must be >= 1 */
    hashtable_divider = 1;
  hashtable_length = number_of_bins/hashtable_divider + 1;  

  sfc_hash_ptr = (SFC_HASH_OBJ_PTR *) LB_MALLOC(sizeof(SFC_HASH_OBJ_PTR) * hashtable_length);
  for(i=0;i<hashtable_length;i++)
    sfc_hash_ptr[i] = NULL;


  for(i=0;i<num_local_objects;i++) {
    sfc_vert_ptr[i].my_bin = 
      get_array_location(number_of_bins, amount_of_bits, 0,
			 (sfc_vert_ptr+i), sfc_keylength, size_of_unsigned, imax);
    sfc_vert_ptr[i].destination_proc = (sfc_vert_ptr[i].my_bin)/(2*bins_per_proc);
    if(sfc_vert_ptr[i].destination_proc != lb->Proc) {
      array_location = LB_Hash(&(sfc_vert_ptr[i].my_bin), 1, hashtable_length);
      ierr = put_in_hashtable(sfc_hash_ptr, array_location, &(sfc_vert_ptr[i]),
			      wgt_dim, &(objs_wgt[i*wgt_dim]));
      if(ierr != LB_OK) {
	LB_PRINT_ERROR(lb->Proc, yo, "Zoltan error in put_in_hashtable function.");
	return(ierr);
      }
    }
    else {
      for(j=0;j<wgt_dim;j++)
	binned_weight_array[(sfc_vert_ptr[i].my_bin % (2*bins_per_proc))*wgt_dim+j] +=
	  objs_wgt[i*wgt_dim+j];
    }
  }
    
  off_proc_objects = 0;
  for(i=0;i<hashtable_length;i++) {
    if(sfc_hash_ptr[i] != NULL) {
      extra_hash_ptr = sfc_hash_ptr[i];
      while(extra_hash_ptr != NULL) {
	off_proc_objects++;
	extra_hash_ptr = extra_hash_ptr->next;
      }
    }
  }
  counter = 0;

  /*create array of processors to send data to and create array of objects to send */
  send_buffer = (SFC_BIN_WEIGHT_PTR) LB_MALLOC(sizeof(SFC_BIN_WEIGHT) * off_proc_objects);
  proclist = (int*) LB_MALLOC(sizeof(int) * off_proc_objects);
  for(i=0;i<hashtable_length;i++) {
    if(sfc_hash_ptr[i] != NULL) {
      extra_hash_ptr = sfc_hash_ptr[i];
      while(extra_hash_ptr != NULL) {
	proclist[counter] = extra_hash_ptr->destination_proc;
	send_buffer[counter].bin = extra_hash_ptr->id;
	send_buffer[counter].weight = extra_hash_ptr->weight_ptr[0];
	counter++;
	extra_hash_ptr = extra_hash_ptr->next;
      }
    }
  }  
  
  ierr = LB_Comm_Create(plan, off_proc_objects, proclist, lb->Communicator, comm_tag, &nreturn);

  rcv_buffer = (SFC_BIN_WEIGHT_PTR) LB_MALLOC(sizeof(SFC_BIN_WEIGHT) * nreturn);


  /*  loop for multiple weight  */
  for(i=0;i<wgt_dim;i++) {
    /* refill the weights in the send buffer for wgt_dim > 1 */
    if(i>1) 
      for(j=0;j<hashtable_length;j++) 
	if(sfc_hash_ptr[j] != NULL) {
	  extra_hash_ptr = sfc_hash_ptr[j];
	  while(extra_hash_ptr != NULL) {
	    send_buffer[counter].bin = extra_hash_ptr->id;
	    send_buffer[counter].weight = extra_hash_ptr->weight_ptr[i];
	    counter++;
	    extra_hash_ptr = extra_hash_ptr->next;
	  }
	}

    ierr = LB_Comm_Do(*plan, comm_tag+1, (char *) send_buffer, 
		      sizeof(SFC_BIN_WEIGHT), (char *) rcv_buffer);
    
    /* put weights from other processors in their bins */
    for(j=0;j<nreturn;j++) 
      binned_weight_array[(rcv_buffer[j].bin % (2*bins_per_proc))*wgt_dim+i] +=
	rcv_buffer[j].weight;
      
  }
  for(i=0;i<2*bins_per_proc;i++)
    printf("weight[%d] on proc %d is %e \n",i+lb->Proc*bins_per_proc*2,
	   lb->Proc, binned_weight_array[i]);

  ierr = LB_Comm_Destroy(plan);
  *plan = NULL; 
  sfc_clear_hashtable(sfc_hash_ptr, hashtable_length);

  LB_FREE(&sfc_hash_ptr);
  LB_FREE(&send_buffer);
  LB_FREE(&rcv_buffer);
  LB_FREE(&proclist);

  /* global distributed array has been created, now perform the scan operation on it */
  
  /* first, must sum up total weight */
  if(wgt_dim < lb->Num_Proc)
    extra_float_array = (float*) LB_MALLOC(sizeof(float) * lb->Num_Proc);
  else
    extra_float_array = (float*) LB_MALLOC(sizeof(float) * wgt_dim);
  for(i=0;i<wgt_dim;i++)
    extra_float_array[i] = 0;
  for(i=0;i<2*bins_per_proc;i++)
    for(j=0;j<wgt_dim;j++)
      extra_float_array[j] += binned_weight_array[i*wgt_dim+j];

  ierr = MPI_Allreduce(extra_float_array, total_weight_array, 
		       wgt_dim, MPI_FLOAT, MPI_SUM, lb->Communicator);
    
  /* put in desired amount of work here, needs to be changed for varying workloads */
  my_work_percent = 1.0/((float) lb->Num_Proc);

  for(i=0;i<lb->Num_Proc;i++)
    extra_float_array[i] = 0.0;
  extra_float_array[lb->Proc] = my_work_percent;

  ierr = MPI_Allreduce(extra_float_array, work_percent_array, 
		       lb->Num_Proc, MPI_FLOAT, MPI_SUM, lb->Communicator);

  for(i=lb->Num_Proc-2;i>=0;i--)
    work_percent_array[i] += work_percent_array[i+1];

  /* each processor needs to know which bins get partitioned into which processor */
  bin_proc_array = (int*) LB_MALLOC(sizeof(int) * lb->Num_Proc); /* lists max bin that a processor should get */
  for(i=0;i<lb->Num_Proc;i++)
    bin_proc_array[i] = 0.;
  
  for(i=0;i<wgt_dim;i++)
    extra_float_array[i] = 0.;
  for(i=0;i<2*bins_per_proc;i++)
    for(j=0;j<wgt_dim;j++) 
      extra_float_array[j] += binned_weight_array[i*wgt_dim+j];

  extra_float_array2 = (float*) LB_MALLOC(sizeof(float) * wgt_dim);
  for(i=0;i<wgt_dim;i++)
    extra_float_array2[i] = extra_float_array[i];

  scanned_work_prev_allocated = (float*) LB_MALLOC(sizeof(float) * wgt_dim);
  actual_work_allocated = (float*) LB_MALLOC(sizeof(float)*(lb->Num_Proc)*wgt_dim);
  for(i=0;i<lb->Num_Proc*wgt_dim;i++)
    actual_work_allocated[i] = 0;
  
  ierr = MPI_Scan(extra_float_array, scanned_work_prev_allocated, wgt_dim,
		  MPI_FLOAT, MPI_SUM, lb->Communicator);
  /* make scan start from proc(num_proc - 1) and finish at proc(0) */
  for(i=0;i<wgt_dim;i++) {
    scanned_work_prev_allocated[i] = total_weight_array[i] - scanned_work_prev_allocated[i];
    printf("proc %d scanned_prev is %e my work is %e total is %e\n", lb->Proc, 
	   scanned_work_prev_allocated[i], extra_float_array2[i], total_weight_array[i]);
  }
  LB_FREE(&extra_float_array);
  LB_FREE(&extra_float_array2);
  


  if(wgt_dim == 1) {
    current_proc = lb->Num_Proc-1;
    while(scanned_work_prev_allocated[0]>work_percent_array[current_proc]*total_weight_array[0] 
	  && current_proc!= 0) 
      current_proc--;
    
    single_wgt_calc_partition(wgt_dim, scanned_work_prev_allocated[0],
			      total_weight_array, bin_proc_array, lb, 
			      binned_weight_array, work_percent_array,
			      actual_work_allocated, 2*bins_per_proc, 
			      number_of_cuts, current_proc, SFC_COARSE_LEVEL_FLAG);
  }
  else {
    /* multi_wgt_dim_calc_partition(); */  
    /*fill in with erik's multi-weight stuff when it is ready, use first weight for now */
    current_proc = lb->Num_Proc-1;
    while(scanned_work_prev_allocated[0]>work_percent_array[current_proc]*total_weight_array[0] 
	  && current_proc!= 0) 
      current_proc--;
    single_wgt_calc_partition(wgt_dim, scanned_work_prev_allocated[0], 
			      total_weight_array, bin_proc_array, lb, 
			      binned_weight_array, work_percent_array,
			      actual_work_allocated, 2*bins_per_proc, 
			      number_of_cuts, current_proc, SFC_COARSE_LEVEL_FLAG);
  }

  ierr = MPI_Allreduce(actual_work_allocated, global_actual_work_allocated, 
		       (lb->Num_Proc)*wgt_dim, MPI_FLOAT, MPI_MAX, lb->Communicator);

  LB_FREE(&actual_work_allocated);
  
  global_bin_proc_array = (int*) LB_MALLOC(sizeof(int)*lb->Num_Proc);
  ierr = MPI_Allreduce(bin_proc_array, global_bin_proc_array, lb->Num_Proc, MPI_INT, 
		       MPI_MAX, lb->Communicator);

  LB_FREE(&bin_proc_array);
  if(lb->Proc == 0)
    for(i=0;i<lb->Num_Proc;i++)
      printf("global_bin_proc_array[%d]= %d\n",i, global_bin_proc_array[i]);


  /* specify which processor an object belongs to,
     we will know this because we know what bin an object 
     belongs to and we know what processor a bin belongs to */
  for(i=0;i<num_local_objects;i++) {
    /* bad search method but easy, should probably replace */
    j=lb->Num_Proc-1;
    while((int) sfc_vert_ptr[i].my_bin < global_bin_proc_array[j])
      j--;
    sfc_vert_ptr[i].destination_proc = j;
    if(sfc_vert_ptr[i].my_bin != global_bin_proc_array[j])
      sfc_vert_ptr[i].cut_bin_flag = SFC_NO_CUT;
    else
      sfc_vert_ptr[i].cut_bin_flag = SFC_CUT;
  }
    
  local_balanced_flag = single_wgt_find_imbalance(work_percent_array,
						  global_actual_work_allocated[lb->Proc*wgt_dim],
						  total_weight_array[0], lb->Proc, lb);

  ierr = MPI_Allreduce(&local_balanced_flag, balanced_flag, 1, MPI_INT, MPI_MAX, lb->Communicator);
  printf("proc %d has balanced_flag = %d\n",lb->Proc, *balanced_flag);
  *balanced_flag = SFC_NOT_BALANCED;   /*change this!!!*/
  LB_FREE(&global_bin_proc_array);
  LB_FREE(&binned_weight_array);
  LB_FREE(&scanned_work_prev_allocated);
  if(*balanced_flag == SFC_BALANCED) {
    return LB_OK;
  }

  /* move the objects that belong to any bin that contains a cut to the proper processor */
  off_proc_objects = 0;

  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].cut_bin_flag == SFC_CUT) 
      off_proc_objects++;  /* actually, this includes objects on this processor as well */

  /*create array of processors to send data to and create array of objects to send */
  send_vert_buffer = (SFC_VERTEX_PTR) LB_MALLOC(sizeof(SFC_VERTEX) * off_proc_objects);
  send_wgt_buffer = (float*) LB_MALLOC(sizeof(float) * wgt_dim * off_proc_objects);
  proclist = (int*) LB_MALLOC(sizeof(int) * off_proc_objects);
  counter = 0;
  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].cut_bin_flag == SFC_CUT)  {
      send_vert_buffer[counter] = sfc_vert_ptr[i];
      for(j=0;j<wgt_dim;j++)
	send_wgt_buffer[counter*wgt_dim+j] = objs_wgt[i*wgt_dim+j];
      proclist[counter] = sfc_vert_ptr[i].destination_proc;	
      counter++;
    }
  comm_tag+=10;  /*  create new comm tag (10 is chosen arbitrarily) */
  *num_vert_sent = off_proc_objects;
  ierr = LB_Comm_Create(plan, off_proc_objects, proclist, lb->Communicator,
			comm_tag, num_vert_in_cut);

  /* send out vertices */
  *vert_in_cut_ptr = (SFC_VERTEX_PTR) LB_MALLOC(sizeof(SFC_VERTEX) * (*num_vert_in_cut));
  ierr = LB_Comm_Do(*plan, comm_tag, (char *) send_vert_buffer, 
		    sizeof(SFC_VERTEX), (char *) *vert_in_cut_ptr);
  LB_FREE(&send_vert_buffer);
    
  /* send out weights of vertices */
  *wgts_in_cut_ptr = (float*) LB_MALLOC(sizeof(float) * (*num_vert_in_cut) *wgt_dim);
  ierr = LB_Comm_Do(*plan, comm_tag+2, (char *) send_wgt_buffer, 
		    sizeof(float)*wgt_dim, (char *) *wgts_in_cut_ptr);
  LB_FREE(&send_wgt_buffer);
  
  LB_FREE(&proclist);

  /* objects that are in a bin that has a cut in it have been sent to 
     their corresponding processors */

  return LB_OK;
}

void single_wgt_calc_partition(int wgt_dim, float work_prev_allocated,
			       float* total_weight_array, int* bin_proc_array, 
			       LB* lb, float* binned_weight_array, 
			       float* work_percent_array, float* actual_work_allocated,
			       int number_of_bins, int* number_of_cuts, int current_loc,
			       int level_flag)
{
  int i;
  int number_of_cuts2 = 0;
  *number_of_cuts = 0;
  
  for(i=number_of_bins-1;i>=0;i--) {
    work_prev_allocated += binned_weight_array[i*wgt_dim];
    if(work_prev_allocated >= total_weight_array[0]*work_percent_array[current_loc]) {
      if(level_flag != SFC_COARSE_LEVEL_FLAG)
	bin_proc_array[current_loc] = i;
      else
	bin_proc_array[current_loc] = number_of_bins*(lb->Proc) + i;
      actual_work_allocated[current_loc*wgt_dim] = work_prev_allocated;
      number_of_cuts2 = 1;
      if(current_loc-number_of_cuts2 < 0)
	printf("problem here !!!\n");
      while(current_loc-number_of_cuts2 >=0 && work_prev_allocated > 
	    total_weight_array[0]*work_percent_array[current_loc-number_of_cuts2] ) {
	actual_work_allocated[(current_loc-number_of_cuts2)*wgt_dim] = work_prev_allocated;
	number_of_cuts2++;
      }
      current_loc--;
    }
    if(*number_of_cuts < number_of_cuts2)
      *number_of_cuts = number_of_cuts2;

  }
  /* make sure that on the coarse level partition, 
     the last proc (proc 0) gets the rest of the work */
  if(level_flag == SFC_COARSE_LEVEL_FLAG)
    bin_proc_array[0] = -1;

  return;

}

void sfc_clear_hashtable(SFC_HASH_OBJ_PTR * sfc_hash_ptr, int hashtable_length)
{
  int i;
  SFC_HASH_OBJ_PTR extra_hash_ptr, another_hash_ptr;

  for(i=0;i<hashtable_length;i++) {
    extra_hash_ptr = sfc_hash_ptr[i];
    while(extra_hash_ptr != NULL) {
      another_hash_ptr = extra_hash_ptr;
      extra_hash_ptr = extra_hash_ptr->next;
      LB_FREE(&(another_hash_ptr->weight_ptr));
      LB_FREE(&another_hash_ptr);
    }
  }

  return;
}

/* will likely speed up this operation if an ordered linklist is used instead of a randomly ordered linklist */
int put_in_hashtable(SFC_HASH_OBJ_PTR * sfc_hash_ptr, int array_location,
		     SFC_VERTEX_PTR sfc_vert_ptr, int wgt_dim, float* obj_wgt)
{
  int i;
  SFC_HASH_OBJ_PTR extra_hash_ptr;

  extra_hash_ptr = sfc_hash_ptr[array_location];

  if(sfc_hash_ptr[array_location] == NULL) {
    sfc_hash_ptr[array_location] = (SFC_HASH_OBJ_PTR) LB_MALLOC(sizeof(SFC_HASH_OBJ));
    (sfc_hash_ptr[array_location])->id = sfc_vert_ptr->my_bin;
    (sfc_hash_ptr[array_location])->destination_proc = sfc_vert_ptr->destination_proc;
    (sfc_hash_ptr[array_location])->next = NULL;
    (sfc_hash_ptr[array_location])->weight_ptr = (float *) LB_MALLOC(sizeof(float) * wgt_dim);
    for(i=0;i<wgt_dim;i++)
      (sfc_hash_ptr[array_location])->weight_ptr[i] = obj_wgt[i];
  }
  else {
    while(extra_hash_ptr->next != NULL && extra_hash_ptr->id != sfc_vert_ptr->my_bin) 
      extra_hash_ptr = extra_hash_ptr->next;

    if(extra_hash_ptr->next == NULL) {
      extra_hash_ptr->next = (SFC_HASH_OBJ_PTR) LB_MALLOC(sizeof(SFC_HASH_OBJ));
      extra_hash_ptr = extra_hash_ptr->next;
      extra_hash_ptr->id = sfc_vert_ptr->my_bin;
      extra_hash_ptr->destination_proc = sfc_vert_ptr->destination_proc;
      extra_hash_ptr->next = NULL;
      extra_hash_ptr->weight_ptr = (float *) LB_MALLOC(sizeof(float) * wgt_dim);
      for(i=0;i<wgt_dim;i++)
	extra_hash_ptr->weight_ptr[i] = obj_wgt[i];
    }
    else {
      extra_hash_ptr = extra_hash_ptr->next;
      for(i=0;i<wgt_dim;i++)
	extra_hash_ptr->weight_ptr[i] += obj_wgt[i];
    }      
  }
  
  return LB_OK;
}
    
/* right_shift is how many bits to ignore on the right
   left_shift is how many to ignore on the left 
   thus, the amount of bits used is size_of_unsigned*8 - right_shift - left_shift
*/

int get_array_location(int number_of_bins, int number_of_bits, int prev_used_bits, 
		       SFC_VERTEX_PTR sfc_vert_ptr, int sfc_keylength, 
		       int size_of_unsigned, unsigned imax)
{
  int counter = 0;
  unsigned ilocation, ilocation2;

  
  
  if(prev_used_bits == 0)
    ilocation = (sfc_vert_ptr->sfc_key[0]) >> (size_of_unsigned*8 - number_of_bits);
  else {
    /* in case prev_used_bits is larger than an unsigned integer */
    while((counter+1)*size_of_unsigned*8 < prev_used_bits)
      counter++;
    prev_used_bits = prev_used_bits - counter*size_of_unsigned*8;

    ilocation2 = (sfc_vert_ptr->sfc_key[counter]) << prev_used_bits;
    ilocation =  ilocation2 >> (size_of_unsigned*8-number_of_bits);
    /* if some of the bits that we want are in the next array value
       this might not be correct!!! */
    if(prev_used_bits+number_of_bits > size_of_unsigned*8) 
      ilocation += ((sfc_vert_ptr->sfc_key[counter+1]) >> 
		    (2*size_of_unsigned*8 - prev_used_bits - number_of_bits));
    
  }
  
  if(ilocation >= number_of_bins) {
    ilocation = number_of_bins - 1;
    printf("possible bad address in get_array_location\n");
  }

  return(ilocation);
}

int single_wgt_find_imbalance(float* work_percent_array, float cumulative_work,
			      float total_work, int which_proc, LB* lb)
/* which_proc is not required to be equal to lb->Proc */
{
  int balanced_flag;
  float my_extra_work;
  float my_ideal_work;

  my_extra_work = cumulative_work - work_percent_array[which_proc]*total_work;

  if(which_proc != lb->Num_Proc - 1) 
    my_ideal_work = (work_percent_array[which_proc] - 
		     work_percent_array[which_proc+1])*total_work;
  else
    my_ideal_work = work_percent_array[which_proc] * total_work;

  if(my_ideal_work == 0) {
    if(my_extra_work != 0) 
      balanced_flag = SFC_NOT_BALANCED;
    else
      balanced_flag = SFC_BALANCED;
  }
  else {
    if(1 + my_extra_work/my_ideal_work > lb->Imbalance_Tol)
      balanced_flag = SFC_NOT_BALANCED;
    else
      balanced_flag = SFC_BALANCED;
  }
  if(balanced_flag==SFC_BALANCED)
    printf("proc %d is balanced!!!!!!\n",which_proc);
  else
    printf("proc %d is not balanced.  tolerance is %e actual is %e\n",which_proc, 
	   lb->Imbalance_Tol,my_extra_work/my_ideal_work
);

  return(balanced_flag);
}

