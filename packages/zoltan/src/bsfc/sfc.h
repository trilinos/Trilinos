#ifndef _LB_SFC_H
#define _LB_SFC_H
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"
#include <values.h>

struct sfc_vertex {         /* SFC vertex  */
  double coord[3];      /* Global coordinates */
  unsigned sfc_key[SFC_KEYLENGTH];   /* space-filling curve key */
  int destination_proc;
  unsigned my_bin;
  int cut_bin_flag;  /* =SFC_NO_CUT (0) if this object does not belong to a bin with a cut in it
			=SFC_CUT (1) if this object belongs to a bin with a cut in it */
  int next_sfc_vert_index; /* used for creating a linked list of vertices 
			      (linklist is fortran style) */
  
};
typedef struct sfc_vertex SFC_VERTEX;
typedef struct sfc_vertex * SFC_VERTEX_PTR; 

struct sfc_hash_obj {
  unsigned id;
  int destination_proc;
  struct sfc_hash_obj * next;
  float* weight_ptr;
};

typedef struct sfc_hash_obj SFC_HASH_OBJ;
typedef struct sfc_hash_obj * SFC_HASH_OBJ_PTR;

struct sfc_bin_weight {
  int bin;
  float weight;
};

typedef struct sfc_bin_weight SFC_BIN_WEIGHT;
typedef struct sfc_bin_weight * SFC_BIN_WEIGHT_PTR;

/* declare functions in the sfc routines */
int single_wgt_find_imbalance(float* work_percent_array, float cumulative_work,
			      float total_work, int which_proc, LB* lb);

void sfc_clear_hashtable(SFC_HASH_OBJ_PTR * sfc_hash_ptr, int hashtable_length);

void single_wgt_calc_partition(int wgt_dim, float work_allocated,
			       float* total_weight_array, int* bin_proc_array, 
			       LB* lb, float* binned_weight_array, 
			       float* work_percent_array, float* actual_work_allocated,
			       int number_of_bins, int* number_of_cuts, int current_proc,
			       int level_flag);

int put_in_hashtable(SFC_HASH_OBJ_PTR * sfc_hash_ptr, int array_location,
		     SFC_VERTEX_PTR sfc_vert_ptr, int wgt_dim, float* obj_wgt);

int get_array_location(int number_of_bins, int number_of_bits, int prev_used_bits, 
		       SFC_VERTEX_PTR sfc_vert_ptr, int sfc_keylength, 
		       int size_of_unsigned, unsigned imax);

void sfc_get_normed_coords(SFC_VERTEX_PTR sfc_vert_ptr, double min_bounding_box[], 
				  double max_bounding_box[], double normed_coords[], int num_dims);

void sfc_create_info(LB *lb, double min_bounding_box[], double max_bounding_box[], 
		     int num_dims, int num_local_objects, int wgt_dim, 
		     SFC_VERTEX_PTR sfc_vert_ptr, int sfc_keylength);

int sfc_refine_partition_level(LB* lb, int* local_balanced_flag, int *amount_of_used_bits,
			       int num_vert_in_cut, SFC_VERTEX_PTR vert_in_cut_ptr,
			       int sfc_keylength, int size_of_unsigned, unsigned imax, int wgt_dim,
			       float* wgts_in_cut_ptr, float* work_percent_array,
			       float* total_weight_array, float* global_actual_work_allocated,
			       int number_of_cuts, int* max_cuts_in_bin, int* ll_bins_head,
			       float* work_prev_allocated, int subbins_per_bin);


#endif /* _LB_SFC_H */
