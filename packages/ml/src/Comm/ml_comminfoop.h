/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the New stuff                                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#ifndef __MLCOMMINFOOP__
#define __MLCOMMINFOOP__

#include "ml_common.h"
#include "ml_comm.h"

typedef struct ML_CommInfoOP_Struct ML_CommInfoOP;
typedef struct ML_NeighborList_Struct ML_NeighborList;
typedef struct ML_Comm_Envelope_Struct ML_Comm_Envelope;

/* ******************************************************************** */
/* ******************************************************************** */
/* Internal data struction definition                                   */
/* ******************************************************************** */
/* ******************************************************************** */

/* ******************************************************************** */
/* data definition for the ML_Operator Class                            */
/* ******************************************************************** */
/* -------------------------------------------------------------------- */
/* 'NeighborList' stores the send and receive information of a given    */
/* neighbor.  This data structure is used primarily in the construction */
/* of Galerkin coarse grid operator.                                    */
/* -------------------------------------------------------------------- */

struct ML_NeighborList_Struct 
{
   int ML_id;
   int N_rcv, N_send;
   int *rcv_list, *send_list;
};

/* -------------------------------------------------------------------- */
/* 'CommInfoOP' records communication information: no. of neighboring   */
/* processors, an array of pointers to the neighbors, and 'add_rcfd'    */
/* which indicates whether or not received information is added to      */
/* existing values or overwrites existing values.  That is, if we       */
/* receive a value that should be stored in the kth location of a       */
/* vector, do we add it with the value already in the kth location or   */
/* do we over-write. Additionally, if remap != NULL, the locations of   */
/* the values are effectively permuted after communication. In          */
/* particular, remap[i] = -1 indicates that the ith value of the vector */
/* is no longer used while remap[i] = k indicates that the ith value of */
/* the vector is actually stored in the kth index of the array.         */
/* -------------------------------------------------------------------- */

struct ML_CommInfoOP_Struct {
   int             N_neighbors;
   ML_NeighborList *neighbors;
   int             add_rcvd;
   int             *remap;
   int             total_rcv_length;
   int             minimum_vec_size, remap_length, remap_max;
};

/* -------------------------------------------------------------------- 
   'Comm_Envelope' holds message envelope information (e.g., message    
   tag). The motivation is to avoid asynchronous communication behavior.        

   tag              the message tag
   -------------------------------------------------------------------- */

struct ML_Comm_Envelope_Struct {
   int             tag;
};

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern int  ML_CommInfoOP_Generate(ML_CommInfoOP **comm_info,
                   int (*user_comm)(double *, void *), void *user_data, 
                   ML_Comm *ml_comm, int N_cols, int Nghost );
extern int  ML_CommInfoOP_Clone(ML_CommInfoOP **newone, ML_CommInfoOP *oldone);

extern ML_CommInfoOP *ML_CommInfoOP_Create(void);
extern void ML_CommInfoOP_Destroy(ML_CommInfoOP **comm_info);

extern int  ML_CommInfoOP_Get_Nneighbors(ML_CommInfoOP *c_info);
extern int *ML_CommInfoOP_Get_neighbors( ML_CommInfoOP *c_info);
extern int *ML_CommInfoOP_Get_sendlist(  ML_CommInfoOP *c_info, int neighbor);
extern int  ML_CommInfoOP_Get_Nsendlist( ML_CommInfoOP *c_info, int neighbor);
extern int  ML_CommInfoOP_Get_Nrcvlist(  ML_CommInfoOP *c_info, int neighbor);
extern int *ML_CommInfoOP_Get_rcvlist(   ML_CommInfoOP *c_info, int neighbor);
extern int  ML_CommInfoOP_Set_neighbors(ML_CommInfoOP  **c_info,int N_neighbors,
                   int *neighbors, int add_or_not, int *remap, int remap_leng);
extern int  ML_CommInfoOP_Set_exch_info(ML_CommInfoOP *comm_info, int k,
                   int N_rcv, int *rcv_list, int N_send, int *send_list);

extern int ML_CommInfoOP_Compute_TotalRcvLength(ML_CommInfoOP *comm_info);

extern int  ML_CommInfoOP_Print(ML_CommInfoOP *c_info, char *label);

extern int ML_CommInfoOP_TransComm(ML_CommInfoOP *pre_comm, 
				   ML_CommInfoOP **post_comm,
				   int invec_leng);

extern void ML_create_unique_col_id(int Ncols, int **map, ML_CommInfoOP *,
                                    int *max_per_proc, ML_Comm *comm);
extern void ML_create_unique_id(int N_local, int **map, ML_CommInfoOP *, ML_Comm *comm);
extern void ML_cheap_exchange_bdry(double dtemp[], ML_CommInfoOP *comm_info,
                                   int start_location, int total_send,
                                   ML_Comm *comm);
extern void ML_exchange_bdry(double dtemp[], ML_CommInfoOP *comm_info,
                             int start_location, ML_Comm *comm,
                             int overwrite_or_add, ML_Comm_Envelope *envelope);
extern int ML_reverse_exchange(double *, ML_CommInfoOP *, int, ML_Comm *comm);

extern void ML_transposed_exchange_bdry(double x[], ML_CommInfoOP *comm_info,
                     int start_location, ML_Comm *comm, int overwrite_or_add); 

extern int ML_Comm_Envelope_Create(ML_Comm_Envelope**);
extern int ML_Comm_Envelope_Init(ML_Comm_Envelope* envelope);
extern int ML_Comm_Envelope_Destroy(ML_Comm_Envelope*);
extern int ML_Comm_Envelope_Clean(ML_Comm_Envelope*);
extern int ML_Comm_Envelope_Get_Tag(ML_Comm_Envelope*, int* tag);
extern int ML_Comm_Envelope_Set_Tag(ML_Comm_Envelope*, int, int);
extern int ML_Comm_Envelope_Increment_Tag(ML_Comm_Envelope*);
#define ML_Comm_Envelope_Print

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif

