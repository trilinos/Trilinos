/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <values.h>
#include <limits.h>

#include "zz_const.h"
#include "zz_util_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "hilbert_const.h"
#include "sfc.h"

int Zoltan_BSFC_refine_overloaded_bins(ZZ *zz, int max_cuts_in_bin, 
			       int actual_bins_per_proc,
			       int* number_of_cuts_in_bin, int wgt_dim,
			       BSFC_VERTEX_PTR sfc_vert_ptr,
			       float objs_wgt[], int num_local_objects,
			       int prev_used_bits, int size_of_unsigned,
			       float* work_percent_array,
			       float* total_weight_array,
			       float* actual_work_allocated);

 /* global_actual_work_allocated is used to make sure that each processor 
    knows how much extra work every processor has allocated
    extra work[proc] =
    global_actual_work_allocated[proc] - work_percent_array[proc]*total_work */

int Zoltan_BSFC_create_bins(ZZ *zz, int num_local_objects, int wgt_dim,
		    BSFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[],
		    int* amount_of_bits_used, int size_of_unsigned, 
		    float* global_actual_work_allocated, 
		    float *work_percent_array, float* total_weight_array, 
		    int* balanced_flag, BSFC_VERTEX_PTR *vert_in_cut_ptr,
		    float** wgts_in_cut_ptr, int* num_vert_in_cut,
		    int* number_of_cuts, int bins_per_proc, 
		    int hashtable_divider, ZOLTAN_COMM_OBJ **plan,
		    int* num_vert_sent, int max_cuts_in_bin)
{
  char    yo[] = "Zoltan_BSFC_create_bins";
  int i, j, number_of_bins, ierr = 0;
  int array_location = 0;
  int comm_tag = 4190; 
  int * proclist;
  int nreturn = 0;
  int off_proc_objects = 0;  /*counter to keep track of how 
			       many objects will be off processor*/
  float * binned_weight_array;
  BSFC_HASH_OBJ_PTR * sfc_hash_ptr;
  BSFC_HASH_OBJ_PTR extra_hash_ptr;
  int hashtable_length;
  int counter = 0;
  BSFC_BIN_WEIGHT_PTR send_buffer, rcv_buffer;
  float *extra_float_array;
  float my_work_percent;
  int *bin_proc_array;
  float *scanned_work_prev_allocated; /*scanned_work_prev_allocated is the 
					amount of work allocated to higher
					ranked procs */
  float *actual_work_allocated;
  int* global_bin_proc_array;
  int amount_of_bits;
  BSFC_VERTEX_PTR send_vert_buffer;
  float* send_wgt_buffer;
  int current_proc;
  float* extra_float_array2 = NULL;
  int local_balanced_flag;
  int* number_of_cuts_in_bin;
  
  ZOLTAN_TRACE_ENTER(zz, yo);
  binned_weight_array = 
    (float *) ZOLTAN_MALLOC(sizeof(float) * 2 * bins_per_proc * wgt_dim);

  if(binned_weight_array == NULL && bins_per_proc > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  for(i=0;i<2*bins_per_proc*wgt_dim;i++)
    binned_weight_array[i] = 0;

  /*assume initially that each processor has the same amount of bins*/
  number_of_bins = zz->Num_Proc * bins_per_proc;
  i=0;
  while(number_of_bins > pow(2,i))
    i++;
  amount_of_bits = i;
  /* check to see that we have not used up all of the bits */
  if(amount_of_bits > 8*size_of_unsigned * BSFC_KEYLENGTH)
    amount_of_bits = 8*size_of_unsigned * BSFC_KEYLENGTH;
  number_of_bins = pow(2,i);
  *amount_of_bits_used = amount_of_bits;

  /*hash table */
  if(hashtable_divider<1) /* hashtable_divider must be >= 1 */
    hashtable_divider = 1;
  hashtable_length = number_of_bins/hashtable_divider + 1;  /* hashtable length must be
							       greater than 0 */

  sfc_hash_ptr = (BSFC_HASH_OBJ_PTR *)
    ZOLTAN_MALLOC(sizeof(BSFC_HASH_OBJ_PTR) * hashtable_length);

  if(sfc_hash_ptr == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  for(i=0;i<hashtable_length;i++)
    sfc_hash_ptr[i] = NULL;

  for(i=0;i<num_local_objects;i++) {
    sfc_vert_ptr[i].my_bin = 
      Zoltan_BSFC_get_array_location(number_of_bins, amount_of_bits, 0, (sfc_vert_ptr+i));
    sfc_vert_ptr[i].destination_proc = 
      (sfc_vert_ptr[i].my_bin)/(2*bins_per_proc);
    if(sfc_vert_ptr[i].destination_proc != zz->Proc) {
      array_location = Zoltan_Hash(&(sfc_vert_ptr[i].my_bin), 1, hashtable_length);
      ierr = Zoltan_BSFC_put_in_hashtable(zz, sfc_hash_ptr, array_location, 
				  &(sfc_vert_ptr[i]), wgt_dim, (objs_wgt+i*wgt_dim));
      if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
		       "Zoltan error in put_in_hashtable function.");
	return(ierr);
      }
    }
    else {
      for(j=0;j<wgt_dim;j++)
	binned_weight_array[((sfc_vert_ptr[i].my_bin % 
			      (2*bins_per_proc))*wgt_dim)+j] +=
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

  /* create array of processors to send data to
     and create array of objects to send */
  send_buffer = (BSFC_BIN_WEIGHT_PTR) 
    ZOLTAN_MALLOC(sizeof(BSFC_BIN_WEIGHT) * off_proc_objects);

  if(send_buffer == NULL && off_proc_objects > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }    
  proclist = (int*) ZOLTAN_MALLOC(sizeof(int) * off_proc_objects);
  if(proclist == NULL && off_proc_objects > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  counter = 0;
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
  
  ierr = Zoltan_Comm_Create(plan, off_proc_objects, proclist, 
			zz->Communicator, comm_tag, &nreturn);

  if(ierr == ZOLTAN_WARN) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Create.");
  }
  else if(ierr == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Create.");
    return(ierr);
  }      
  else if(ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Create.");
    return(ierr);
  }

  rcv_buffer = (BSFC_BIN_WEIGHT_PTR) ZOLTAN_MALLOC(sizeof(BSFC_BIN_WEIGHT) * nreturn);
  if(rcv_buffer == NULL && nreturn > 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
  }

  /*  loop for multiple weight  */
  for(i=0;i<wgt_dim;i++) {
    /* refill the weights in the send buffer for wgt_dim > 1 */
    if(i>=1) {
      counter = 0;
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
    }

    ierr = Zoltan_Comm_Do(*plan, comm_tag+1+i, (char *) send_buffer, 
		      sizeof(BSFC_BIN_WEIGHT), (char *) rcv_buffer);
    if(ierr == ZOLTAN_WARN) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Do.");
    }
    else if(ierr == ZOLTAN_FATAL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Do.");
      return(ierr);
    }      
    else if(ierr == ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Do.");
      return(ierr);
    }    
    /* put weights from other processors in their bins */
    for(j=0;j<nreturn;j++) 
      binned_weight_array[(rcv_buffer[j].bin %
			   (2*bins_per_proc))*wgt_dim+i] +=
	rcv_buffer[j].weight;   
  }

  ierr = Zoltan_Comm_Destroy(plan);
  if(ierr == ZOLTAN_WARN) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Destory.");
  }
  else if(ierr == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Destroy.");
    return(ierr);
  }      
  else if(ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Destroy.");
    return(ierr);
  }
  *plan = NULL; 
  Zoltan_BSFC_clear_hashtable(sfc_hash_ptr, hashtable_length);

  ZOLTAN_FREE(&sfc_hash_ptr);
  ZOLTAN_FREE(&send_buffer);
  ZOLTAN_FREE(&rcv_buffer);
  ZOLTAN_FREE(&proclist);

  /* global distributed array has been created,
     now perform the scan operation on it */
  
  /* first, must sum up total weight */
  if(wgt_dim < zz->Num_Proc)
    extra_float_array = (float*) ZOLTAN_MALLOC(sizeof(float) * zz->Num_Proc);
  else
    extra_float_array = (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
  if(extra_float_array == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }   
  for(i=0;i<wgt_dim;i++)
    extra_float_array[i] = 0;
  for(i=0;i<2*bins_per_proc;i++)
    for(j=0;j<wgt_dim;j++)
      extra_float_array[j] += binned_weight_array[i*wgt_dim+j];
  
  ierr = MPI_Allreduce(extra_float_array, total_weight_array, 
		       wgt_dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
  
  /* put in desired amount of work here, needs to
     be changed for varying workloads */
  my_work_percent = 1.0/((float) zz->Num_Proc);
  
  for(i=0;i<zz->Num_Proc;i++)
    extra_float_array[i] = 0.0;
  extra_float_array[zz->Proc] = my_work_percent;
  /* make sure that proc 0 gets all of the rest of the work percent */
  if(zz->Proc == 0)
    extra_float_array[0] = 1.1;    
  
  ierr = MPI_Allreduce(extra_float_array, work_percent_array, 
		       zz->Num_Proc, MPI_FLOAT, MPI_SUM, zz->Communicator);
  
  for(i=zz->Num_Proc-2;i>=0;i--)
    work_percent_array[i] += work_percent_array[i+1];
  
  for(i=0;i<wgt_dim;i++)
    extra_float_array[i] = 0.;
  for(i=0;i<2*bins_per_proc;i++)
    for(j=0;j<wgt_dim;j++) 
      extra_float_array[j] += binned_weight_array[i*wgt_dim+j];
  
  extra_float_array2 = (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
  if(extra_float_array2 == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  for(i=0;i<wgt_dim;i++)
    extra_float_array2[i] = extra_float_array[i];
  
  scanned_work_prev_allocated = (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
  if(scanned_work_prev_allocated == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  actual_work_allocated = 
    (float*) ZOLTAN_MALLOC(sizeof(float)*(zz->Num_Proc)*wgt_dim);

  if(actual_work_allocated == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  for(i=0;i<zz->Num_Proc*wgt_dim;i++)
    actual_work_allocated[i] = 0;
  
  ierr = MPI_Scan(extra_float_array, scanned_work_prev_allocated, 
		  wgt_dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
  /* make scan start from proc(num_proc - 1) and finish at proc(0) */
  for(i=0;i<wgt_dim;i++) 
    scanned_work_prev_allocated[i] = 
      total_weight_array[i] - scanned_work_prev_allocated[i];
  
  ZOLTAN_FREE(&extra_float_array);
  ZOLTAN_FREE(&extra_float_array2);
  
  /* each processor needs to know which bins get partitioned
     into which processor, bin_proc_array lists max bin that
     a processor should get */
  bin_proc_array = (int*) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
  if(bin_proc_array == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }    
  number_of_cuts_in_bin = (int*) ZOLTAN_MALLOC(sizeof(int) * 2*bins_per_proc);
  if(number_of_cuts_in_bin == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  for(i=0;i<zz->Num_Proc;i++)
    bin_proc_array[i] = 2*number_of_bins;

  if(wgt_dim == 1) {
    current_proc = zz->Num_Proc-1;
    while(current_proc!= 0 && scanned_work_prev_allocated[0] > 
	  work_percent_array[current_proc]*total_weight_array[0]) 
      current_proc--;
    
    Zoltan_BSFC_single_wgt_calc_partition(wgt_dim, scanned_work_prev_allocated[0],
				  total_weight_array, bin_proc_array, zz, 
				  binned_weight_array, work_percent_array,
				  actual_work_allocated, 2*bins_per_proc, 
				  number_of_cuts, current_proc, 
				  BSFC_COARSE_LEVEL_FLAG, number_of_cuts_in_bin);
  }
  else {
    /* multi_wgt_dim_calc_partition(); */  
    /*fill in with erik's multi-weight stuff when it is ready,
      use first weight for now */
    current_proc = zz->Num_Proc-1;
    while(current_proc != 0 && scanned_work_prev_allocated[0] > 
	  work_percent_array[current_proc]*total_weight_array[0] ) 
      current_proc--;
    Zoltan_BSFC_single_wgt_calc_partition(wgt_dim, scanned_work_prev_allocated[0], 
			      total_weight_array, bin_proc_array, zz, 
			      binned_weight_array, work_percent_array,
			      actual_work_allocated, 2*bins_per_proc, 
			      number_of_cuts, current_proc, 
			      BSFC_COARSE_LEVEL_FLAG, number_of_cuts_in_bin);
  }

  /* -1 is used to make sure that the last bin does not have a cut in it */
  bin_proc_array[0] = -1;

  ierr = MPI_Allreduce(actual_work_allocated, global_actual_work_allocated, 
		       (zz->Num_Proc)*wgt_dim, MPI_FLOAT, 
		       MPI_MAX, zz->Communicator);

  ZOLTAN_FREE(&actual_work_allocated);
  
  global_bin_proc_array = (int*) ZOLTAN_MALLOC(sizeof(int)*zz->Num_Proc);
  if(global_bin_proc_array == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  ierr = MPI_Allreduce(bin_proc_array, global_bin_proc_array, zz->Num_Proc,
		       MPI_INT, MPI_MIN, zz->Communicator);

  ZOLTAN_FREE(&bin_proc_array);

  /* specify which processor an object belongs to,
     we will know this because we know what bin an object 
     belongs to and we know what processor a bin belongs to */
  for(i=0;i<num_local_objects;i++) {
    /* bad search method but easy, should probably replace */
    j=zz->Num_Proc-1;
    while((int) sfc_vert_ptr[i].my_bin < global_bin_proc_array[j])
      j--;
    sfc_vert_ptr[i].destination_proc = j;
    if((int) sfc_vert_ptr[i].my_bin != global_bin_proc_array[j])
      sfc_vert_ptr[i].cut_bin_flag = BSFC_NO_CUT;
    else
      sfc_vert_ptr[i].cut_bin_flag = BSFC_CUT;
  }

  /* check to see if any cut-bin has too many objects in it and refine it. 
     the problem with too many objects in a cut-bin is that the processor 
     that gets assigned to that bin will get swamped with communication
     and might not have enough memory to hold all of the information.  we
     detect the cut-bins with too many objects by how many cuts are in that 
     bin.  zz->Num_Proc - 1 is the amount of cuts and if this is less than
     or equal to max_cuts_in_bin then there is no possibility that a
     coarse bin is overloaded */
  if(zz->Num_Proc - 1 > max_cuts_in_bin)
    ierr = Zoltan_BSFC_refine_overloaded_bins(zz, max_cuts_in_bin, 2*bins_per_proc, 
				      number_of_cuts_in_bin, wgt_dim,
				      sfc_vert_ptr, objs_wgt,num_local_objects,
				      amount_of_bits, size_of_unsigned,
				      work_percent_array, total_weight_array, 
				      global_actual_work_allocated);
  
  ZOLTAN_FREE(&number_of_cuts_in_bin);    

  if(wgt_dim == 1) {
    local_balanced_flag = 
      Zoltan_BSFC_single_wgt_find_imbalance(work_percent_array,
				    global_actual_work_allocated[zz->Proc*wgt_dim],
				    total_weight_array[0], zz->Proc, zz);
  }
  else {
    /* put in routine here to calculate imbalance for multiple weights */
    local_balanced_flag = 
      Zoltan_BSFC_single_wgt_find_imbalance(work_percent_array,
				    global_actual_work_allocated[zz->Proc*wgt_dim],
				    total_weight_array[0], zz->Proc, zz);
  }

  ierr = MPI_Allreduce(&local_balanced_flag, balanced_flag, 1,
		       MPI_INT, MPI_MAX, zz->Communicator);

  ZOLTAN_FREE(&global_bin_proc_array);
  ZOLTAN_FREE(&binned_weight_array);
  ZOLTAN_FREE(&scanned_work_prev_allocated);
  /* if the current partitioning is acceptable, the algorithm is finished */
  if(*balanced_flag == BSFC_BALANCED) {
    return ZOLTAN_OK;
  }

  /* if the size of an unsigned integer is different on different processors,
     need to 'shrink' down the key size of unsigned integers on all processors
     to have the same amount of bits but only need to do this for sfc objects
     that are in a cut bin */
  if(size_of_unsigned != sizeof(unsigned)) {
    int k;
    int my_size_of_unsigned = sizeof(unsigned);
    int difference = sizeof(unsigned) - size_of_unsigned;
    unsigned copy_key[BSFC_KEYLENGTH];
    for(i=0;i<num_local_objects;i++) 
      if(sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT) {
	copy_key[0] = sfc_vert_ptr[i].sfc_key[0] >> (8*difference); 
	for(j=1;j<BSFC_KEYLENGTH;j++) {
	  k=0;
	  while(size_of_unsigned*(j+1) > my_size_of_unsigned*(k+1))
	    k++;
	  
	  /* all needed bits in one sfc_key */
	  if(size_of_unsigned*(j+1) > my_size_of_unsigned*k) {
	    copy_key[j] = sfc_vert_ptr[i].sfc_key[k] << 
	      (8*(size_of_unsigned*j-my_size_of_unsigned*k));
	    copy_key[j] = copy_key[j] >> (8*difference);
	  }
	  /* needed bits in 2 different sfc_keys */
	  else { 
	    /* first key */
	    copy_key[j] = sfc_vert_ptr[i].sfc_key[k-1] <<  
	      (8*(size_of_unsigned*j-my_size_of_unsigned*(k-1)));
	    copy_key[j] = copy_key[j] >> (8*difference);
	    /* second key */
	    copy_key[j] += sfc_vert_ptr[i].sfc_key[k] >>
	      (8*2*(my_size_of_unsigned*k-size_of_unsigned*j));
	  }
	}
	for(j=0;j<BSFC_KEYLENGTH;j++)
	  sfc_vert_ptr[i].sfc_key[j] = copy_key[j];
      }
  }

  /* move the sfc objects that belong to any bin that contains a cut
     to the proper processor */
  off_proc_objects = 0;

  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT) 
      off_proc_objects++;  /* actually, this includes objects on this
			      processor as well, this is done to simplify 
			      the code */

  /*create array of processors to send data to and
    create array of objects to send */
  send_vert_buffer = 
    (BSFC_VERTEX_PTR) ZOLTAN_MALLOC(sizeof(BSFC_VERTEX) * off_proc_objects);
  if(send_vert_buffer == NULL && off_proc_objects > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  send_wgt_buffer = 
    (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim * off_proc_objects);
  if(send_wgt_buffer == NULL && off_proc_objects > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  proclist = (int*) ZOLTAN_MALLOC(sizeof(int) * off_proc_objects);
  if(proclist == NULL && off_proc_objects > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  counter = 0;
  for(i=0;i<num_local_objects;i++) 
    if(sfc_vert_ptr[i].cut_bin_flag == BSFC_CUT)  {
      send_vert_buffer[counter] = sfc_vert_ptr[i];
      for(j=0;j<wgt_dim;j++)
	send_wgt_buffer[counter*wgt_dim+j] = objs_wgt[i*wgt_dim+j];
      proclist[counter] = sfc_vert_ptr[i].destination_proc;	
      counter++;
    }
  comm_tag+=10;  /* create new comm tag (10 is chosen arbitrarily) */
  *num_vert_sent = off_proc_objects;
  ierr = Zoltan_Comm_Create(plan, off_proc_objects, proclist, 
			zz->Communicator, comm_tag, num_vert_in_cut);
  if(ierr == ZOLTAN_WARN) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Create.");
  }
  else if(ierr == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Create.");
    return(ierr);
  }      
  else if(ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Create.");
    return(ierr);
  }
  /* send out vertices */
  *vert_in_cut_ptr = 
    (BSFC_VERTEX_PTR) ZOLTAN_MALLOC(sizeof(BSFC_VERTEX) * (*num_vert_in_cut));
  if(*vert_in_cut_ptr == NULL && *num_vert_in_cut > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  ierr = Zoltan_Comm_Do(*plan, comm_tag, (char *) send_vert_buffer, 
		    sizeof(BSFC_VERTEX), (char *) *vert_in_cut_ptr);
  if(ierr == ZOLTAN_WARN) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Do.");
  }
  else if(ierr == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Do.");
    return(ierr);
  }      
  else if(ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Do.");
    return(ierr);
  }
  ZOLTAN_FREE(&send_vert_buffer);
  
  /* send out weights of vertices */
  *wgts_in_cut_ptr = 
    (float*) ZOLTAN_MALLOC(sizeof(float) * (*num_vert_in_cut) * wgt_dim);
  if(*wgts_in_cut_ptr == NULL && *num_vert_in_cut > 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  ierr = Zoltan_Comm_Do(*plan, comm_tag+2, (char *) send_wgt_buffer, 
		    sizeof(float)*wgt_dim, (char *) *wgts_in_cut_ptr);
  if(ierr == ZOLTAN_WARN) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning from Zoltan_Comm_Do.");
  }
  else if(ierr == ZOLTAN_FATAL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Fatal error in Zoltan_Comm_Do.");
    return(ierr);
  }      
  else if(ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error in Zoltan_Comm_Do.");
    return(ierr);
  }
  ZOLTAN_FREE(&send_wgt_buffer);
  
  ZOLTAN_FREE(&proclist);

  /* objects that are in a bin that has a cut in it have been sent to 
     their corresponding processors */
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}
/*  done sfc_create_bins routine */

/* routine calculates what bins the cuts are in, how many cuts are in 
   each bin and the actual work allocated up to that cut.  routine
   will work for wgt_dim > 1 but will only use the first weight 
   of an object */
void Zoltan_BSFC_single_wgt_calc_partition(int wgt_dim, float work_prev_allocated,
				   float* total_weight_array, int* bin_proc_array,
				   ZZ *zz, float* binned_weight_array, 
				   float* work_percent_array, 
				   float* actual_work_allocated, int number_of_bins,
				   int* number_of_cuts, int current_loc,
				   int level_flag, int* number_of_cuts_in_bin)
{
  int i;
  int number_of_cuts2 = 0;
  *number_of_cuts = 0;

  if(level_flag == BSFC_COARSE_LEVEL_FLAG)
    for(i=0;i<number_of_bins;i++)
      number_of_cuts_in_bin[i] = 0;

  for(i=number_of_bins-1;i>=0;i--) {
    work_prev_allocated += binned_weight_array[i*wgt_dim];
    if(work_prev_allocated >= 
       total_weight_array[0] * work_percent_array[current_loc]) {
      actual_work_allocated[current_loc*wgt_dim] = work_prev_allocated;
      if(level_flag != BSFC_COARSE_LEVEL_FLAG)
	bin_proc_array[current_loc] = i;
      else
	bin_proc_array[current_loc] = number_of_bins*(zz->Proc) + i;
      number_of_cuts2 = 1;
      while(current_loc-number_of_cuts2 >=0 && work_prev_allocated > 
	    total_weight_array[0]*work_percent_array[current_loc-number_of_cuts2]) {
	actual_work_allocated[(current_loc-number_of_cuts2)*wgt_dim] =
	  work_prev_allocated;
	number_of_cuts2++;
      }
      current_loc=current_loc-number_of_cuts2;
    
      if(number_of_cuts_in_bin != NULL)
	number_of_cuts_in_bin[i] = number_of_cuts2;

      if(*number_of_cuts < number_of_cuts2)
	*number_of_cuts = number_of_cuts2;
    }
  }

  /* make sure that the first bin gets the rest of the work */
  bin_proc_array[0] = -1;
  if(level_flag != BSFC_COARSE_LEVEL_FLAG && current_loc >= 0)
    bin_proc_array[current_loc] = -1;

  return;
}

/* free all of the memory allocated in the hashtable (including
   the linklist used for collisions in the hashtable) */
void Zoltan_BSFC_clear_hashtable(BSFC_HASH_OBJ_PTR * sfc_hash_ptr, 
			 int hashtable_length)
{
  int i;
  BSFC_HASH_OBJ_PTR extra_hash_ptr, another_hash_ptr;

  for(i=0;i<hashtable_length;i++) {
    extra_hash_ptr = sfc_hash_ptr[i];
    while(extra_hash_ptr != NULL) {
      another_hash_ptr = extra_hash_ptr;
      extra_hash_ptr = extra_hash_ptr->next;
      ZOLTAN_FREE(&(another_hash_ptr->weight_ptr));
      ZOLTAN_FREE(&another_hash_ptr);
    }
  }

  return;
}

/* this routine puts adds up the weights of the objects in the coarse
   bins and stores the information in a hashtable.  the hashtable is 
   used because otherwise each processor would need to alloocate an
   array of size (2 * bins_per_proc * zz->Num_Proc), with the
   hashtable, only an array of size 
   (2 * bins_per_proc * Num_Proc / HASHTABLE_DIVIDER) needs to be
   allocated.  if the destination bin has not been put in the 
   hashtable yet, memory is allocated and the object weight(s) is put
   in.  if the destination bin has already been allocated, the object
   weight(s) is added in. */

int Zoltan_BSFC_put_in_hashtable(ZZ *zz, BSFC_HASH_OBJ_PTR * sfc_hash_ptr, 
			 int array_location, BSFC_VERTEX_PTR sfc_vert_ptr,
			 int wgt_dim, float* obj_wgt)
{
  int i;
  BSFC_HASH_OBJ_PTR extra_hash_ptr;
  char    yo[] = "Zoltan_BSFC_put_in_hashtable";

  extra_hash_ptr = sfc_hash_ptr[array_location];
  /* if this location has not been filled yet */
  if(sfc_hash_ptr[array_location] == NULL) { 
    sfc_hash_ptr[array_location] = 
      (BSFC_HASH_OBJ_PTR) ZOLTAN_MALLOC(sizeof(BSFC_HASH_OBJ));
    if(sfc_hash_ptr[array_location] == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    (sfc_hash_ptr[array_location])->id = sfc_vert_ptr->my_bin;
    (sfc_hash_ptr[array_location])->destination_proc = 
      sfc_vert_ptr->destination_proc;
    (sfc_hash_ptr[array_location])->next = NULL;
    (sfc_hash_ptr[array_location])->prev = NULL;
    (sfc_hash_ptr[array_location])->weight_ptr = 
      (float *) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
    if((sfc_hash_ptr[array_location])->weight_ptr == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    for(i=0;i<wgt_dim;i++)
      (sfc_hash_ptr[array_location])->weight_ptr[i] = obj_wgt[i];
  }
  /* location is at the beginning of the linklist */
  else if(extra_hash_ptr->id > sfc_vert_ptr->my_bin) { 
    sfc_hash_ptr[array_location] = 
      (BSFC_HASH_OBJ_PTR) ZOLTAN_MALLOC(sizeof(BSFC_HASH_OBJ));
    if(sfc_hash_ptr[array_location] == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    (sfc_hash_ptr[array_location])->id = sfc_vert_ptr->my_bin;
    (sfc_hash_ptr[array_location])->destination_proc = 
      sfc_vert_ptr->destination_proc;
    (sfc_hash_ptr[array_location])->next = extra_hash_ptr;
    extra_hash_ptr->prev = sfc_hash_ptr[array_location];
    (sfc_hash_ptr[array_location])->prev = NULL;
    (sfc_hash_ptr[array_location])->weight_ptr = 
      (float *) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
    if((sfc_hash_ptr[array_location])->weight_ptr == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    for(i=0;i<wgt_dim;i++)
      (sfc_hash_ptr[array_location])->weight_ptr[i] = obj_wgt[i];
  }
  else {
    while(extra_hash_ptr->next != NULL &&
	  extra_hash_ptr->id < sfc_vert_ptr->my_bin) 
      extra_hash_ptr = extra_hash_ptr->next;

    if(extra_hash_ptr->id == sfc_vert_ptr->my_bin) {
      for(i=0;i<wgt_dim;i++)
	extra_hash_ptr->weight_ptr[i] += obj_wgt[i];
    }      
    else if(extra_hash_ptr->next == NULL) {
      BSFC_HASH_OBJ_PTR another_hash_ptr = extra_hash_ptr;
      extra_hash_ptr->next =
	(BSFC_HASH_OBJ_PTR) ZOLTAN_MALLOC(sizeof(BSFC_HASH_OBJ));
      if(extra_hash_ptr->next == NULL) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
	return(ZOLTAN_MEMERR);
      }
      extra_hash_ptr = extra_hash_ptr->next;
      extra_hash_ptr->id = sfc_vert_ptr->my_bin;
      extra_hash_ptr->destination_proc = sfc_vert_ptr->destination_proc;
      extra_hash_ptr->next = NULL;
      extra_hash_ptr->prev = another_hash_ptr;
      extra_hash_ptr->weight_ptr = 
	(float *) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
      if(extra_hash_ptr->weight_ptr == NULL) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
	return(ZOLTAN_MEMERR);
      }
      for(i=0;i<wgt_dim;i++)
	extra_hash_ptr->weight_ptr[i] = obj_wgt[i];
    }
    else {
      BSFC_HASH_OBJ_PTR another_hash_ptr = extra_hash_ptr;
      extra_hash_ptr = extra_hash_ptr->prev;
      extra_hash_ptr->next =
	(BSFC_HASH_OBJ_PTR) ZOLTAN_MALLOC(sizeof(BSFC_HASH_OBJ));
      if(extra_hash_ptr->next == NULL) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
	return(ZOLTAN_MEMERR);
      }
      extra_hash_ptr = extra_hash_ptr->next;
      extra_hash_ptr->prev = another_hash_ptr->prev;
      another_hash_ptr->prev = extra_hash_ptr;
      extra_hash_ptr->id = sfc_vert_ptr->my_bin;
      extra_hash_ptr->destination_proc = sfc_vert_ptr->destination_proc;
      extra_hash_ptr->next = another_hash_ptr;
      extra_hash_ptr->weight_ptr = 
	(float *) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
      if(extra_hash_ptr->weight_ptr == NULL) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
	return(ZOLTAN_MEMERR);
      }
      for(i=0;i<wgt_dim;i++)
	extra_hash_ptr->weight_ptr[i] = obj_wgt[i];
    }
  }
  
  return ZOLTAN_OK;
}
    
/* 
   routine calculates the new bin number of an object using its
   sfc_key.  prev_used_bits is how many bits have already been
   used and number_of_bits is how many bits to use to calculate
   the key.  the output is the bin number and will be a value
   between 0 and (2^number_of_bits - 1)
*/

int Zoltan_BSFC_get_array_location(int number_of_bins, int number_of_bits, 
			   int prev_used_bits, BSFC_VERTEX_PTR sfc_vert_ptr)
{
  int counter = 0;
  unsigned ilocation = 0;
  unsigned ilocation2 = 0;
  int pub = prev_used_bits;
  int size_of_unsigned = sizeof(unsigned);

  if(prev_used_bits == 0)
    ilocation = 
      (sfc_vert_ptr->sfc_key[0]) >> (size_of_unsigned*8 - number_of_bits);
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
		    (2*size_of_unsigned*8-prev_used_bits-number_of_bits));  
  }
  if(ilocation >= number_of_bins)  {
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    printf("OH DAMN OH on proc %d, ilocation is %d and # bins is %d pub is %d \n", myid, ilocation, number_of_bins, pub);
    ilocation = number_of_bins - 1;
  }

  return(ilocation);
}


/* routine finds the imbalance for processor which_proc given the
   cumulative work assigned to processors which_proc through the
   last processor (zz->Num_Proc-1).  routine returns value which 
   indicates whether this imbalance is beyond a specified 
   tolerance.  routine works for wgt_dim > 1 but only calculates
   imbalance from first weight */
int Zoltan_BSFC_single_wgt_find_imbalance(float* work_percent_array, 
				  float cumulative_work, 
				  float total_work,
				  int which_proc, ZZ *zz)
/* NOTE:  which_proc is not required to be equal to zz->Proc */
{
  int balanced_flag;
  float my_extra_work;
  float my_ideal_work;

  my_extra_work = 
    cumulative_work - work_percent_array[which_proc]*total_work;

  if(which_proc != zz->Num_Proc - 1) 
    my_ideal_work = (work_percent_array[which_proc] - 
		     work_percent_array[which_proc+1])*total_work;
  else
    my_ideal_work = work_percent_array[which_proc] * total_work;

  /* if processor which_proc is not supposed to have any work,
     it is imbalanced if it has any amount of work greater
     than 0 */
  if(my_ideal_work == 0) {
    if(my_extra_work != 0) 
      balanced_flag = BSFC_NOT_BALANCED;
    else
      balanced_flag = BSFC_BALANCED;
  }
  else {
    if(1 + my_extra_work/my_ideal_work > zz->LB.Imbalance_Tol)
      balanced_flag = BSFC_NOT_BALANCED;
    else
      balanced_flag = BSFC_BALANCED;
  }

  return(balanced_flag);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
