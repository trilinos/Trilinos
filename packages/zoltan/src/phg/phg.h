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

#ifndef __ZOLTAN_PHG_H
#define __ZOLTAN_PHG_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "hg.h"
    
#include "phg_comm.h"
#include "phg_const.h"
#include "phg_util.h"
#include "params_const.h"
#include "zoltan_comm.h"

/********************************************************************
 * Data structure for Zoltan's base hypergraph. Includes Zoltan IDs 
 * corresponding to local objects (vertices) and a HGraph as used 
 * by the partitioning algorithms. 
 ********************************************************************/

struct Zoltan_PHGraph {
  int nObj;                 /* Number of on-processor objects. */
  ZOLTAN_ID_PTR Global_IDs; /* Global IDs for on-processor objects.  */
  ZOLTAN_ID_PTR Local_IDs;  /* Local IDs for on-processor objects.   */
  int *Parts;               /* Initial partition #s for on-processor objects */
                            /* KDD In parallel version Part may be in HG */
  ZOLTAN_COMM_OBJ *VtxPlan; /* Communication plan mapping GIDs to GNOs 
                               within row communicators. */
  HGraph PHG;              /* Hypergraph for initial objects.       */
};
typedef struct Zoltan_PHGraph ZPHG;

/************************************************/
/* Mappings supporting the 2D data distribution */
/************************************************/

/* Mapping of global number to local number           */
/* Code should call VTX_GNO_TO_LNO or EDGE_GNO_TO_LNO */

#define GNO_TO_LNO(gno, dist, myblock) \
    (gno) - (dist)[(myblock)]
#define VTX_GNO_TO_LNO(phg, gno) \
    GNO_TO_LNO(gno, (phg)->dist_x, (phg)->comm->myProc_x)
#define EDGE_GNO_TO_LNO(phg, gno) \
    GNO_TO_LNO(gno, (phg)->dist_y, (phg)->comm->myProc_y)


/* Mapping of local number to global number           */
/* Code should call VTX_LNO_TO_GNO or EDGE_LNO_TO_GNO */

#define LNO_TO_GNO(lno, dist, myblock) \
    (lno) + (dist)[(myblock)]
#define VTX_LNO_TO_GNO(phg, lno) \
    LNO_TO_GNO(lno, (phg)->dist_x, (phg)->comm->myProc_x)
#define EDGE_LNO_TO_GNO(phg, lno) \
    LNO_TO_GNO(lno, (phg)->dist_y, (phg)->comm->myProc_y)


/* Mapping of global number to processor block.     */
/* Code should call EDGE_TO_PROC_Y or VTX_TO_PROC_X */

#define EDGE_TO_PROC_Y(phg, gno) \
    Zoltan_PHG_Gno_To_Proc_Block((gno), (phg)->dist_y, (phg)->comm->nProc_y)

#define VTX_TO_PROC_X(phg, gno) \
    Zoltan_PHG_Gno_To_Proc_Block(gno, (phg)->dist_x, (phg)->comm->nProc_x)


/* Function types for options to hypergraph partitioning */
struct PHGPartParamsStruct;  /* Forward declaration */

typedef int ZOLTAN_PHG_MATCHING_FN(ZZ*, HGraph*, Matching);
typedef int ZOLTAN_PHG_COARSEPARTITION_FN(ZZ*, HGraph*, int, Partition,
                                          struct PHGPartParamsStruct*);
typedef int ZOLTAN_PHG_REFINEMENT_FN(ZZ*, HGraph*, int, Partition,
                                     struct PHGPartParamsStruct*, float);


    
/******************************************/
/* Parameters to the hypergraph functions */
/******************************************/
struct PHGPartParamsStruct {
  float bal_tol;                 /* Balance tolerance in % of average */
  int kway;                      /* 1 -> direct kway, 0->recursive bisection */
  int redl;                      /* Reduction limit (constant). */
  char redm_str[MAX_PARAM_STRING_LEN];   /* Reduction method string. */
  ZOLTAN_PHG_MATCHING_FN *matching;      /* Pointers to Matching function */
  ZOLTAN_PHG_MATCHING_FN *matching_opt;  /* Pointers to Matching opt function */
  int ews;                               /* type of hyperedge weight scaling */
  char coarsepartition_str[MAX_PARAM_STRING_LEN]; 
                                         /* Coarse partitioning string */
  ZOLTAN_PHG_COARSEPARTITION_FN *CoarsePartition;
                                         /* pointer to coarse partitioning fn */
  char refinement_str[MAX_PARAM_STRING_LEN]; /* Refinement string and */
  ZOLTAN_PHG_REFINEMENT_FN *Refinement;      /* pointer to refinement fn */

  int fm_loop_limit;    /* Number of FM loops if the refinement is FM */
  int fm_max_neg_move;  /* Maximum number of vertex moves with negative gain;
                           an early stop condition. <0 means try all moves */
    
  int check_graph;      /* Flag indicating whether the input hypergraph should 
                         * be checked for errors. */
  int output_level;     /* Flag indicating amount of output from HG algorithms.
                         * See levels PHG_DEBUG_* below.  */
  PHGComm comm;   /* UVCUVC: although this is not a paramater; we'll keep it here
                     for now; later we can move it out */
};

typedef struct PHGPartParamsStruct PHGPartParams;

    
/*****************************/
/* Hypergraph output levels: */
/*****************************/
#define PHG_DEBUG_NONE 0
#define PHG_DEBUG_LIST 1
#define PHG_DEBUG_ALL  2
#define PHG_DEBUG_PRINT 3
#define PHG_DEBUG_PLOT 4


/*********************/
/* Scale Edge Weight */
/*********************/
int Zoltan_PHG_Scale_HGraph_Weight (ZZ*, HGraph*, float*, int);

/**********************/
/* Matching functions */
/**********************/
int Zoltan_PHG_Matching (ZZ*, HGraph*, Matching, PHGPartParams*);
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams*);

/**************/
/* Coarsening */
/**************/
int Zoltan_PHG_Coarsening(ZZ*, HGraph*, Matching, HGraph*, int*);

/*********************************/
/* Coarse Partitioning functions */
/*********************************/
extern int Zoltan_PHG_Gather_To_All_Procs(ZZ*, HGraph*, PHGComm*, HGraph**);
extern int Zoltan_PHG_CoarsePartition(ZZ*, HGraph*, int, Partition, 
                                      PHGPartParams*);
ZOLTAN_PHG_COARSEPARTITION_FN *Zoltan_PHG_Set_CoarsePartition_Fn(char*);

/************************/
/* Refinement functions */ 
/************************/
int Zoltan_PHG_Refinement (ZZ*, HGraph*, int, Partition, PHGPartParams*);
ZOLTAN_PHG_REFINEMENT_FN *Zoltan_PHG_Set_Refinement_Fn(char*);

/*****************/
/* Communication */
/*****************/
extern int Zoltan_PHG_gather_slice_root(
        int, int, int, int, int, char*, int*, char**, MPI_Comm*, int);
extern int Zoltan_PHG_gather_slice(
        int, int, int, int, int, char*, int*, char**, MPI_Comm*, int);

/*****************************/
/* Other Function Prototypes */
/*****************************/

extern int Zoltan_PHG_rdivide (int,  int, Partition, ZZ *, HGraph *,
                               PHGPartParams *, int);
    
extern int Zoltan_PHG_Set_Part_Options(ZZ*, PHGPartParams*);
extern int Zoltan_PHG_HPart_Lib(ZZ*, HGraph*, int, Partition, PHGPartParams*, 
                                int);
extern int Zoltan_PHG_HPart_Info(ZZ*, HGraph*, int, Partition, PHGPartParams*);
extern double Zoltan_PHG_hcut_size_total(PHGComm*, HGraph*, Partition, int);
extern double Zoltan_PHG_hcut_size_links(PHGComm*, HGraph*, Partition, int);    
extern double Zoltan_PHG_HPart_balance(ZZ*, HGraph*, int, Partition);

extern int Zoltan_PHG_Build_Hypergraph(ZZ*, ZPHG**, PHGPartParams*);
extern void Zoltan_PHG_HGraph_Print(ZZ*, ZPHG*,  HGraph*, FILE*);
extern void Zoltan_PHG_Plot(int, int, int, int*, int*, int*, char*);
extern void Zoltan_PHG_Plot_2D_Distrib(ZZ*, HGraph*);

extern int Zoltan_PHG_Gno_To_Proc_Block(int gno, int*, int);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_H */
