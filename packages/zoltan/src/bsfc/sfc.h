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
/*  int current_proc; /*necessary ????*/
  int destination_proc;
  unsigned my_bin;
/*  int id_location; /* the array location for the GID and LID on the original proc */  /*necessary ????*/
  int cut_bin_flag;  /* =SFC_NO_CUT (0) if this object does not belong to a bin with a cut in it
		    =SFC_CUT (1) if this object belongs to a bin with a cut in it */
  int next_sfc_vert_index; /* used for creating a linked list of vertices (linklist is fortran style) */
  
};
typedef struct sfc_vertex SFC_VERTEX;
typedef struct sfc_vertex * SFC_VERTEX_PTR;

struct sfc_hash_obj {
  int id;
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

#endif /* _LB_SFC_H */
