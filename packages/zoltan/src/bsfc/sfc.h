#ifndef _ZOLTAN_SFC_H
#define _ZOLTAN_SFC_H
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "sfc_const.h"
#include <values.h>

/* define some constants that are used in multiple files */
#define SFC_KEYLENGTH 3
#define SFC_NO_CUT 0
#define SFC_CUT 1
#define SFC_NOT_BALANCED 1
#define SFC_BALANCED 0
#define SFC_COARSE_LEVEL_FLAG 2

struct sfc_vertex {         /* SFC vertex  */
  unsigned sfc_key[SFC_KEYLENGTH];   /* space-filling curve key */
  int destination_proc;
  unsigned my_bin;
  int cut_bin_flag;  /* =SFC_NO_CUT (0) if this object does not belong to
			a bin with a cut in it
			=SFC_CUT (1) if this object belongs to a bin with
			a cut in it */
  int next_sfc_vert_index; /* used for creating a linked list of vertices 
			      (linklist is fortran style) */
  
};
typedef struct sfc_vertex SFC_VERTEX;
typedef struct sfc_vertex * SFC_VERTEX_PTR; 

struct sfc_hash_obj {
  unsigned id;
  int destination_proc;
  struct sfc_hash_obj * next;
  struct sfc_hash_obj * prev;
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
int sfc_single_wgt_find_imbalance(float* work_percent_array, 
				  float cumulative_work, float total_work,
				  int which_proc, ZZ *zz);

void sfc_clear_hashtable(SFC_HASH_OBJ_PTR * sfc_hash_ptr,
			 int hashtable_length);

void sfc_single_wgt_calc_partition(int wgt_dim, float work_allocated,
				   float* total_weight_array, 
				   int* bin_proc_array, ZZ *zz, 
				   float* binned_weight_array, 
				   float* work_percent_array, 
				   float* actual_work_allocated,
				   int number_of_bins, int* number_of_cuts,
				   int current_proc, int level_flag, int*);

int sfc_put_in_hashtable(ZZ *zz, SFC_HASH_OBJ_PTR * sfc_hash_ptr, 
			 int array_location, SFC_VERTEX_PTR sfc_vert_ptr, 
			 int wgt_dim, float* obj_wgt);

int sfc_get_array_location(int number_of_bins, int number_of_bits, 
			   int prev_used_bits, SFC_VERTEX_PTR sfc_vert_ptr);

void sfc_get_normed_coords(double min_bounding_box[], 
			   double max_bounding_box[], 
			   double normed_coords[],
			   int num_dims, double my_coords[]);

void sfc_create_info(ZZ *zz, double min_bounding_box[], 
		     double max_bounding_box[], int num_dims,
		     int num_local_objects, int wgt_dim, 
		     SFC_VERTEX_PTR sfc_vert_ptr, 
		     double* coords);

int sfc_refine_partition(ZZ *zz, int* local_balanced_flag,
			 int *amount_of_used_bits, int num_vert_in_cut,
			 SFC_VERTEX_PTR vert_in_cut_ptr,
			 int size_of_unsigned,  
			 int wgt_dim, float* wgts_in_cut_ptr, 
			 float* work_percent_array, 
			 float* total_weight_array,
			 float* global_actual_work_allocated,
			 int number_of_cuts, int* max_cuts_in_bin, 
			 int* ll_bins_head, float* work_prev_allocated,
			 int subbins_per_bin, int* local_balanced_flag_array,
			 int bin_refinement_method);

int sfc_create_compare_key(ZZ *zz, unsigned sfc_key[], unsigned compare_key[], 
			   unsigned AND_operator_array[], int prev_used_bits);

int sfc_check_refine(unsigned* sfc_key, unsigned* compare_key,
		     unsigned* AND_operator_array);

#endif /* _ZOLTAN_SFC_H */
