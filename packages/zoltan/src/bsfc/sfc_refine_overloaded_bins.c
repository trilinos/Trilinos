/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * This software is distributed under the GNU Lesser General Public License. *
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
#include "zz_const.h"
#include "params_const.h"
#include "timer_const.h"
#include <values.h>
#include <limits.h>
#include "hilbert_const.h"
#include "sfc.h"
#include "all_allo_const.h"


int Zoltan_BSFC_refine_coarse_bin(ZZ *zz, int num_local_objects, 
			  BSFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[],
			  int wgt_dim, unsigned* refine_key,
			  unsigned*, int proc, int my_bins_per_proc,
			  int number_of_bits, int prev_used_bits, 
			  int number_of_cuts, float* work_percent_array, 
			  float* total_weight_array, 
			  float* actual_work_allocated,
			  int max_cuts_in_bin, int level_flag);

/* routine figures out if any coarse bins have too many cuts in 
   them and if the do, it will refine them until they do not
   have too many cuts in them */
int Zoltan_BSFC_refine_overloaded_bins(ZZ *zz, int max_cuts_in_bin,
			       int actual_bins_per_proc,
			       int* number_of_cuts_in_bin, int wgt_dim,
			       BSFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[],
			       int num_local_objects, int prev_used_bits, 
			       int size_of_unsigned,
			       float* work_percent_array,
			       float* total_weight_array,
			       float* actual_work_allocated)
{
  char yo[] = "Zoltan_BSFC_refine_overloaded_bins";
  int i, j, k, *istore, *istore2, *bins_to_refine, number_of_bins, ierr;
  int number_of_bits, number_of_bins_to_refine;
  int gl_max_cuts;
  float* fstore;
  int* overloaded_bin_flag;

  if(size_of_unsigned != sizeof(unsigned)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "The size of an unsigned integer is smaller on another processor. Cannot use this routine. Try increasing max_cuts_in_bin.");
    return(ZOLTAN_FATAL);
  }      

  overloaded_bin_flag = (int*) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
  if(overloaded_bin_flag == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }  
  for(i=0;i<zz->Num_Proc;i++)
    overloaded_bin_flag[i] = 0;
  j=0;
  for(i=0;i<actual_bins_per_proc;i++)
    if(number_of_cuts_in_bin[i] > max_cuts_in_bin) 
      j++;
  
  overloaded_bin_flag[zz->Proc] = j;
    
  istore = (int*) ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
  if(istore == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }  
  
  i = MPI_Allreduce(overloaded_bin_flag, istore, zz->Num_Proc, 
		    MPI_INT, MPI_SUM, zz->Communicator);
  
  ZOLTAN_FREE(&overloaded_bin_flag);
  
  number_of_bins_to_refine=0;
  for(i=0;i<zz->Num_Proc;i++) 
    number_of_bins_to_refine += istore[i];
    
  if(number_of_bins_to_refine==0) {
    /* all of the bins are okay! */
    ZOLTAN_FREE(&istore);
    return(ZOLTAN_OK);
  }

  if(istore[zz->Proc])
    ZOLTAN_PRINT_WARN(0, yo,
      "A coarse bin has too many cuts in it and it is being refined globally.");
  
  /* since we need to refine some coarse bins... */

  /* calculate how many subbins to break a bin into */
  j = 0;
  for(i=0;i<actual_bins_per_proc;i++)
    if(number_of_cuts_in_bin[i] > j)
      j = number_of_cuts_in_bin[i];

  ierr = MPI_Allreduce(&j, &gl_max_cuts, 1, MPI_INT,
		       MPI_MAX, zz->Communicator);

  i=0;
  while(pow(2,i) < (gl_max_cuts+1))
    i++;
  number_of_bins = pow(2,i);
  number_of_bits = i;
  
  /* figure out which coarse bins to refine */
  istore2 = (int*) ZOLTAN_MALLOC(sizeof(int) * number_of_bins_to_refine);
  if(istore2 == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }  
  for(i=0;i<number_of_bins_to_refine;i++)
    istore2[i] = 0;
  k=0;
  for(i=0;i<zz->Num_Proc;i++) {
    if(i==zz->Proc) {
      for(j=0;j<actual_bins_per_proc;j++)
	if(number_of_cuts_in_bin[j] > max_cuts_in_bin) {
	  istore2[k] = zz->Proc*actual_bins_per_proc+j;
	  k++;
	}
    }
    else 
      k += istore[i]; 
  }
  
  bins_to_refine = (int*) ZOLTAN_MALLOC(sizeof(int) * number_of_bins_to_refine);
  if(bins_to_refine == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }  
  i = MPI_Allreduce(istore2, bins_to_refine, number_of_bins_to_refine, 
		    MPI_INT, MPI_SUM, zz->Communicator);
  
  ZOLTAN_FREE(&istore2);
  
  /* now go through all of the bins separately and refine them until 
     all bins are less than max_cuts_in_bin */
  k=0;
  for(i=0;i<number_of_bins_to_refine;i++) {
    /* first BSFC_KEYLENGTH parts of the array is for the compare key 
       and the second is for the AND_operator_array.  if prev_used_bits
       of the sfc key array are the same as the compare key, objects
       belong to the same bin.  the AND_operator array is used to look 
       at only the first prev_used_bits of an sfc key */
    unsigned refine_key[2*BSFC_KEYLENGTH], refine_key_gl[2*BSFC_KEYLENGTH];

    int refine_proc = bins_to_refine[i]/actual_bins_per_proc;
    int refine_bin = bins_to_refine[i];

    /* find an object which belongs to refine_bin and then calculate
       compare_key and AND_operator array */
    k=0;
    while(k<num_local_objects && refine_bin >= 0) {
      if((int) sfc_vert_ptr[k].my_bin == refine_bin) {
	ierr = Zoltan_BSFC_create_compare_key(zz, sfc_vert_ptr[k].sfc_key, refine_key, 
				      (refine_key+BSFC_KEYLENGTH), prev_used_bits);
	if(ierr == ZOLTAN_FATAL)
	  return ZOLTAN_FATAL;
	
	refine_bin = -1; /* we have found an object that is in refine_bin */
      }
      k++;
    }
    ierr = MPI_Allreduce(refine_key, refine_key_gl, 2*BSFC_KEYLENGTH,
			 MPI_UNSIGNED, MPI_MAX, zz->Communicator);
    ierr = Zoltan_BSFC_refine_coarse_bin(zz, num_local_objects, sfc_vert_ptr,
				 objs_wgt, wgt_dim, refine_key_gl, 
				 (refine_key_gl+BSFC_KEYLENGTH), 
				 refine_proc, number_of_bins,  
				 number_of_bits, prev_used_bits, 
				 gl_max_cuts, work_percent_array,
				 total_weight_array, actual_work_allocated, 
				 max_cuts_in_bin, 0);
    if(ierr!= ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_BSFC_refine_coarse_bin function.");
      return(ierr);
    }
  }

  ZOLTAN_FREE(&istore);
  ZOLTAN_FREE(&bins_to_refine);
  
  /* make sure that every processor has a correct copy of the work allocated
     to all processors */
  
  fstore = (float*) ZOLTAN_MALLOC(sizeof(float) * zz->Num_Proc);
  for(i=0;i<zz->Num_Proc;i++)
    fstore[i] = actual_work_allocated[i];
  
  i = MPI_Allreduce(fstore, actual_work_allocated, zz->Num_Proc, MPI_FLOAT,
		    MPI_MIN, zz->Communicator);
  
  ZOLTAN_FREE(&fstore);
  
  return(ZOLTAN_OK);
}

/* routine does the actual refinement of a coarse bin into subbins.  if
   a refined bin still has too many cuts in it, it will call this 
   routine again until all subbins have at most max_cuts_in_bin cuts
   in them.  we know an object belongs to a bin with too many cuts
   in it if its first prev_used_bits are the same as the compare key's
   first prev_used_bits */
int Zoltan_BSFC_refine_coarse_bin(ZZ *zz, int num_local_objects, 
			  BSFC_VERTEX_PTR sfc_vert_ptr, float objs_wgt[], 
			  int wgt_dim, unsigned* refine_key, 
			  unsigned* AND_operator_array, int proc,
			  int my_bins_per_proc, int number_of_bits, 
			  int prev_used_bits, int number_of_cuts,
			  float* work_percent_array,
			  float* total_weight_array, 
			  float* actual_work_allocated, 
			  int max_cuts_in_bin, int level_flag) 
{
  char yo[] = "Zoltan_BSFC_refine_coarse_bin";
  int i, j, k, l, *bin_proc_array;
  int *number_of_cuts_in_bin, destination_proc = -1;
  int new_number_of_cuts, ierr;
  float *binned_wgt_array, *summed_binned_wgt_array = NULL, work_prev_allocated;

  unsigned check_bin[BSFC_KEYLENGTH+1], gl_check_bin[BSFC_KEYLENGTH+1]; 
  unsigned umax = ~(0u);

  /* check that we haven't used all of the bits already */
  if(prev_used_bits >= 8*sizeof(unsigned) * BSFC_KEYLENGTH) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "No more refinement is possible.");
    return(ZOLTAN_OK);
  }
  if(prev_used_bits + number_of_bits > 8*sizeof(unsigned) * BSFC_KEYLENGTH) 
    number_of_bits = 8*sizeof(unsigned) * BSFC_KEYLENGTH - prev_used_bits;

  /* check_bin is used to check if there is only 1 object
     in bin or if all objects in the bin have the same sfc key */
  for(i=0;i<=BSFC_KEYLENGTH;i++)
    check_bin[i] = 0;
  
  binned_wgt_array = (float*) ZOLTAN_MALLOC(sizeof(float)*my_bins_per_proc*wgt_dim); 
  if(binned_wgt_array == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }
  for(i=0;i<my_bins_per_proc*wgt_dim;i++)
    binned_wgt_array[i] = 0;
  for(i=0;i<num_local_objects;i++) {
    if(Zoltan_BSFC_check_refine(sfc_vert_ptr[i].sfc_key, refine_key, AND_operator_array)) {
      sfc_vert_ptr[i].my_bin = Zoltan_BSFC_get_array_location(my_bins_per_proc, number_of_bits,
						      prev_used_bits, (sfc_vert_ptr+i));
      destination_proc = sfc_vert_ptr[i].destination_proc;
      for(j=0;j<wgt_dim;j++)
	binned_wgt_array[j+wgt_dim*sfc_vert_ptr[i].my_bin] += objs_wgt[i*wgt_dim+j];
      
      for(j=0;j<BSFC_KEYLENGTH;j++)
	check_bin[j] = sfc_vert_ptr[i].sfc_key[j];
      check_bin[BSFC_KEYLENGTH] += 1;	
    } 
  } 
  /* check if refining the bin will help at all */
  j = (int) check_bin[BSFC_KEYLENGTH];
  i = MPI_Allreduce(check_bin, gl_check_bin, BSFC_KEYLENGTH+1, MPI_UNSIGNED, 
		    MPI_MAX, zz->Communicator);
  if((int) gl_check_bin[BSFC_KEYLENGTH] == 1) {
    ierr = MPI_Allreduce(&j, &i, 1, MPI_INT, MPI_SUM, zz->Communicator);
    if(i == 1) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, 
         "A coarse bin has too many cuts in it and the cuts are all from one object.");
      ZOLTAN_FREE(&binned_wgt_array);
      return(ZOLTAN_OK);
    }
  }
  check_bin[BSFC_KEYLENGTH] = 0;
  i=0;
  while((int) check_bin[BSFC_KEYLENGTH] == 0 && i < num_local_objects) {
    if(Zoltan_BSFC_check_refine(sfc_vert_ptr[i].sfc_key, refine_key, AND_operator_array)) {
      for(j=0;j<BSFC_KEYLENGTH;j++)
	if(check_bin[j] != sfc_vert_ptr[i].sfc_key[j])
	  check_bin[BSFC_KEYLENGTH] = 1;
    }
    i++;
  }
  i = MPI_Allreduce((check_bin+BSFC_KEYLENGTH), (gl_check_bin+BSFC_KEYLENGTH), 1,
		    MPI_UNSIGNED, MPI_SUM, zz->Communicator);
  if((int) gl_check_bin[BSFC_KEYLENGTH] == 0) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, 
	"All objects have the same sfc key in a coarse bin with too many cuts in it.");
    /* go through all of the objects with the same key and give them 
       random, non-repeating keys for the bits that haven't been used yet */
    for(i=0;i<num_local_objects;i++) 
      if(Zoltan_BSFC_check_refine(sfc_vert_ptr[i].sfc_key, refine_key, AND_operator_array)) { 
	for(j=0;j<BSFC_KEYLENGTH;j++)
	  sfc_vert_ptr[i].sfc_key[j] = (sfc_vert_ptr[i].sfc_key[i] & AND_operator_array[i]) + 
	    (((umax/(num_local_objects*zz->Num_Proc))*i*(zz->Proc+1)) & ~(AND_operator_array[i]));   
      }
  }
    
  if(proc == zz->Proc) {
    summed_binned_wgt_array = 
      (float*) ZOLTAN_MALLOC(sizeof(float)*my_bins_per_proc*wgt_dim);
    if(summed_binned_wgt_array == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }  
  }
  
  i = MPI_Reduce(binned_wgt_array, summed_binned_wgt_array, 
		 my_bins_per_proc*wgt_dim, MPI_FLOAT, MPI_SUM, 
		 proc, zz->Communicator);
  
  bin_proc_array = (int*) ZOLTAN_MALLOC(sizeof(int)*(number_of_cuts+1));
  if(bin_proc_array == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }  
  for(i=0;i<=number_of_cuts;i++)
    bin_proc_array[i] = my_bins_per_proc;
  number_of_cuts_in_bin = (int*) ZOLTAN_MALLOC(sizeof(int) * my_bins_per_proc);
  if(number_of_cuts_in_bin == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }  
  for(i=0;i<my_bins_per_proc;i++)
    number_of_cuts_in_bin[i] = 0;

  i = destination_proc;
  j = MPI_Allreduce(&i, &destination_proc, 1, MPI_INT, MPI_MAX, zz->Communicator);

  /* only proc finds the new cuts */
  if(proc == zz->Proc) {
    float* summed_wgts = (float*) ZOLTAN_MALLOC(sizeof(float) * wgt_dim);
    for(i=0;i<wgt_dim;i++)
      summed_wgts[i] = 0;
    for(i=0;i<my_bins_per_proc;i++)
      for(j=0;j<wgt_dim;j++)
	summed_wgts[j] += summed_binned_wgt_array[i*wgt_dim+j];

    work_prev_allocated = actual_work_allocated[destination_proc]-summed_wgts[0]; 
    if(wgt_dim == 1)
      Zoltan_BSFC_single_wgt_calc_partition(wgt_dim, work_prev_allocated, total_weight_array,
				    bin_proc_array, zz, summed_binned_wgt_array, 
				    work_percent_array+destination_proc-number_of_cuts,
				    actual_work_allocated+wgt_dim*(destination_proc-number_of_cuts),
				    my_bins_per_proc, &new_number_of_cuts, number_of_cuts,
				    0, number_of_cuts_in_bin);
    else {
      /* multi_wgt_dim_calc_partition(); */  
      /*fill in with erik's multi-weight stuff when it is ready, use first weight for now */
      Zoltan_BSFC_single_wgt_calc_partition(wgt_dim, work_prev_allocated, total_weight_array,
				    bin_proc_array, zz, summed_binned_wgt_array, 
				    work_percent_array+destination_proc-number_of_cuts,
				    actual_work_allocated+wgt_dim*(destination_proc-number_of_cuts),
				    my_bins_per_proc, &new_number_of_cuts, number_of_cuts,
				    0, number_of_cuts_in_bin);
    }
    ZOLTAN_FREE(&summed_wgts);
  }

  /* send out bin_proc_array from proc to the rest of the processors */
  ierr = MPI_Bcast(bin_proc_array, number_of_cuts+1, MPI_INT, proc, zz->Communicator); 
  /* send out how many cuts in each bin from proc to the rest of the processors */
  ierr = MPI_Bcast(number_of_cuts_in_bin, my_bins_per_proc, MPI_INT, 
		   proc, zz->Communicator);

  /* specify which processor an object belongs to,
     we will know this because we know what bin an object 
     belongs to and we know what processor a bin belongs to */
  for(i=0;i<num_local_objects;i++) 
    if(Zoltan_BSFC_check_refine(sfc_vert_ptr[i].sfc_key, refine_key, AND_operator_array)) {
      j=number_of_cuts;
      while((int) sfc_vert_ptr[i].my_bin < bin_proc_array[j])
	j--;
      sfc_vert_ptr[i].destination_proc = destination_proc+j-number_of_cuts;
      if((int) sfc_vert_ptr[i].my_bin != bin_proc_array[j]) 
	sfc_vert_ptr[i].cut_bin_flag = BSFC_NO_CUT;
      else 
	sfc_vert_ptr[i].cut_bin_flag = BSFC_CUT;	
    }
  /* check to see if any of the refined bins still have too many
     cuts in them and then refine them if they do */
  j = 0;
  for(i=my_bins_per_proc-1;i>=0;i--) {
    if(number_of_cuts_in_bin[i] > max_cuts_in_bin) {
      unsigned new_refine_key[BSFC_KEYLENGTH];
      unsigned new_AND_operator_array[BSFC_KEYLENGTH];
      /* calculate new_refine_key and new_AND_operator_array */
      k=0;
      l=-1;
      while(l < 0 && k < num_local_objects ) {
	if(Zoltan_BSFC_check_refine(sfc_vert_ptr[k].sfc_key, refine_key, AND_operator_array)) {
	  if((int) sfc_vert_ptr[k].my_bin == i)
	    l=k;
	}
	k++;
      }
      ierr = Zoltan_BSFC_create_compare_key(zz, sfc_vert_ptr[l].sfc_key, new_refine_key, 
				    new_AND_operator_array, prev_used_bits+number_of_bits);
      if(ierr == ZOLTAN_FATAL) 
	return ZOLTAN_FATAL;
      /* call this routine again */
      ierr=Zoltan_BSFC_refine_coarse_bin(zz, num_local_objects, sfc_vert_ptr, objs_wgt,
				 wgt_dim, new_refine_key, new_AND_operator_array,
				 proc, my_bins_per_proc, number_of_bits,
				 prev_used_bits+number_of_bits,
				 number_of_cuts_in_bin[i], 
				 work_percent_array, total_weight_array, 
				 actual_work_allocated, max_cuts_in_bin, level_flag+1);
      if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
	ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_BSFC_refine_coarse_bin.");
	return(ZOLTAN_FATAL);
      } 
    } 
    j += number_of_cuts_in_bin[i];
  }
      
  ZOLTAN_FREE(&number_of_cuts_in_bin);
  ZOLTAN_FREE(&binned_wgt_array);
  ZOLTAN_FREE(&bin_proc_array);
  ZOLTAN_FREE(&summed_binned_wgt_array);

  return ZOLTAN_OK;
}

/*
  takes an sfc_key and makes a new key which contains the same first prev_used_bits
  and 0 for the rest of the bits.  
  output can be used to determine if objects belong to the same bin because all 
  previously used bits are the same. 
*/
int Zoltan_BSFC_create_compare_key(ZZ *zz, unsigned sfc_key[], unsigned compare_key[], 
			   unsigned AND_operator_array[], int prev_used_bits)
{
  int i;
  unsigned umax = ~(0u);
  int size_of_unsigned = sizeof(unsigned);
  char yo[] = "Zoltan_BSFC_create_compare_key";

  if(prev_used_bits/(size_of_unsigned*8) >= BSFC_KEYLENGTH) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Too many previously used bits.");
    return(ZOLTAN_FATAL);
  }
    
  for(i=0;i<prev_used_bits/(size_of_unsigned*8);i++) 
    AND_operator_array[i] = umax;

  AND_operator_array[prev_used_bits/(size_of_unsigned*8)] =
    (umax >> (8*size_of_unsigned-prev_used_bits)) << (8*size_of_unsigned-prev_used_bits);

  for(i=1+prev_used_bits/(size_of_unsigned*8);i<BSFC_KEYLENGTH;i++) 
    AND_operator_array[i] = 0;

  for(i=0;i<BSFC_KEYLENGTH;i++)
    compare_key[i] = sfc_key[i] & AND_operator_array[i];

  return ZOLTAN_OK;
}

/* 
   checks a key with compare_key to see if an object is in this bin
   returns 1 if it is in the bin and 0 if it is in another bin
*/
int Zoltan_BSFC_check_refine(unsigned* sfc_key, unsigned* compare_key,
		     unsigned* AND_operator_array)
{
  int i=0, iflag = 1;
  
  while(i<BSFC_KEYLENGTH && iflag == 1) {
    if(compare_key[i] != (unsigned) (sfc_key[i] & AND_operator_array[i])) 
      iflag = 0;

    i++;
  }

  return iflag;
}
