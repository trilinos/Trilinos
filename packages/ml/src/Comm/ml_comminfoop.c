/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* functions for creating communicators for the ML_Operator class       */
/* -------------------------------------------------------------------- */

#include <stdlib.h>
#include "ml_comminfoop.h"
#include "ml_utils.h"
#include "ml_memory.h"
#include "ml_xyt.h"

/* ******************************************************************** */
/* Fill in the pre-communication struction of an ML_Operator's getrow   */
/* by using a communication routine supplied by the user.               */
/* -------------------------------------------------------------------- */

int ML_CommInfoOP_Generate(ML_CommInfoOP **comm_info, 
	int (*user_comm)(double *,void *), void *user_data, ML_Comm *ml_comm, 
	int N_cols, int Nghost) 
{
   double *data_vec;
   int    *procs, *tempo;
   int    i, j, index, N_rcv_procs, *rcv_neighbors, *send_neighbors, **rcv_list;
   int    *rcv_length, *send_length, N_send_procs;
   int    proc_id;
   int *send_list, *collect, *rindex, N_neighbors, partner;
   USR_REQ         *request;
   int type, largest_itemp, *itemp;

   proc_id  = ml_comm->ML_mypid;

   data_vec = (double *) ML_allocate((N_cols + Nghost + 1)*sizeof(double));
   procs    = (int    *) ML_allocate(ml_comm->ML_nprocs * sizeof(int));
   tempo    = (int    *) ML_allocate(ml_comm->ML_nprocs * sizeof(int));
   
   for (i = 0; i < N_cols+Nghost; i++) data_vec[i] = (double ) proc_id;
   user_comm(data_vec, user_data);

   /* Compute the number of elements recvd from the ith */
   /* processor in procs[i] and store the total number  */
   /* of processors from which we rcv in 'N_rcv_procs'. */

   for (i = 0; i < ml_comm->ML_nprocs; i++) procs[i] = 0;
   for (i = 0; i < N_cols+Nghost ; i++) {
      procs[ (int) data_vec[i] ]++;
   }
   procs[ proc_id ] = 0;

   N_rcv_procs = 0;
   for (i = 0; i < ml_comm->ML_nprocs; i++) {
      if ( procs[i] > 0) N_rcv_procs++; 
   }

   /* Store the neighbors in rcv_neighbors[k] and the number of */
   /* elements rcvd in rcv_length[k]. Allocate an array to hold */
   /* the rcv list and finally store the index k in procs[i].   */

   rcv_neighbors  = (int * ) ML_allocate( (N_rcv_procs+1)*sizeof(int));
   rcv_length     = (int * ) ML_allocate( (N_rcv_procs+1)*sizeof(int));
   rcv_list       = (int **) ML_allocate( (N_rcv_procs+1)*sizeof(int *));
   N_rcv_procs    = 0;
   for (i = 0; i < ml_comm->ML_nprocs; i++) 
   {
      if ( procs[i] > 0) 
      { 
         rcv_neighbors[N_rcv_procs] = i;
         rcv_list[N_rcv_procs] = (int *) ML_allocate(procs[i] * sizeof(int) );
         procs[i] = N_rcv_procs++;
      }
   }

   /* store the rcv list */

   for (i = 0; i < N_rcv_procs; i++) rcv_length[i] = 0;

   for (i = 0; i < N_cols+Nghost; i++) 
   {
      j = (int) data_vec[i];
      if ( j != proc_id ) 
      {
          index = procs[j];
          rcv_list[   index    ][ rcv_length[index]++ ] = i;
      }
   }

   /* figure out the number of neighbors that we send to */

   for (i = 0; i < N_rcv_procs; i++) procs[rcv_neighbors[i]] = 1;
   ML_gsum_vec_int(&procs, &tempo, ml_comm->ML_nprocs, ml_comm);
   N_send_procs = procs[proc_id];
   ML_free(tempo);
   ML_free(procs);

   /* figure out to whom we send and how many we send to them */

   i = N_send_procs + N_rcv_procs + 1;
   send_neighbors  = (int *    ) ML_allocate( (i)*sizeof(int));
   send_length     = (int *    ) ML_allocate( (i)*sizeof(int));
   request         = (USR_REQ *) ML_allocate((N_send_procs+N_rcv_procs+1)*
                                              sizeof(USR_REQ    ));

   type = 4901;
   for (i = 0; i < N_send_procs ; i++) 
   {
     send_neighbors[i] = -1; /* receive from anyone */
     ml_comm->USR_irecvbytes((void *) &(send_length[i]), sizeof(int) ,
                &(send_neighbors[i]), &type, ml_comm->USR_comm, request+i);
   }
   for (i = 0; i < N_rcv_procs; i++) 
   {
      ml_comm->USR_sendbytes((void *) &(rcv_length[i]), sizeof(int), 
                          rcv_neighbors[i], type, ml_comm->USR_comm);
   }
   for (i = 0; i < N_send_procs ; i++) 
   {
     ml_comm->USR_cheapwaitbytes((void *) &(send_length[i]), sizeof(int) ,
                &(send_neighbors[i]), &type, ml_comm->USR_comm, request+i);
   }
   ML_az_sort( send_neighbors, N_send_procs , send_length, NULL);

   /* Fill in the send list */

   j = 1;
   for (i = 0; i < N_send_procs; i++) j += send_length[i];
   send_list = (int *) ML_allocate(sizeof(int)*j);
   j = 1;
   for (i = 0; i < N_rcv_procs; i++) 
      if (j < rcv_length[i]) j = rcv_length[i];
   collect = (int *) ML_allocate(sizeof(int)*j);


   for (i = 0; i < N_cols+Nghost ; i++) data_vec[i] = (double) i;
   user_comm(data_vec, user_data);

   type++;
   j = 0;
   for (i = 0; i < N_send_procs ; i++) 
   {
      ml_comm->USR_irecvbytes((void *) &(send_list[j]), sizeof(int)*
			       send_length[i], &(send_neighbors[i]), &type, 
			       ml_comm->USR_comm, request+i);
      j += send_length[i];
   }
   for (i = 0; i < N_rcv_procs; i++) 
   {
      for (j = 0; j < rcv_length[i]; j++) 
         collect[j] = (int) data_vec[ rcv_list[i][j] ];
      ml_comm->USR_sendbytes((void *) collect, rcv_length[i]*sizeof(int), 
                          rcv_neighbors[i], type, ml_comm->USR_comm);
   }
   j = 0;
   for (i = 0; i < N_send_procs ; i++) 
   {
      ml_comm->USR_cheapwaitbytes((void *) &(send_list[j]), sizeof(int)*
			     send_length[i], &(send_neighbors[i]), &type, 
			     ml_comm->USR_comm, request+i);
      j += send_length[i];
   }
   ML_free(collect);
   ML_free(request);
   ML_free(data_vec);
  
   /* Finally, put all of this stuff into ML */
     
   N_neighbors = N_send_procs;
   j = 0;
   rindex = (int *) ML_allocate(sizeof(int)*(N_send_procs+N_rcv_procs+1));
   for (i = 0; i < N_send_procs+N_rcv_procs; i++) rindex[i] = -1;

   for (i = 0; i < N_rcv_procs; i++) 
   {
      while (  (j < N_neighbors)&&(send_neighbors[j] < rcv_neighbors[i]))
         j++;

      if ( (j == N_neighbors) || (send_neighbors[j] != rcv_neighbors[i])) 
      {
         rindex[N_neighbors] = i;
         send_neighbors[N_neighbors++] = rcv_neighbors[i];
      }
      else rindex[j++] = i;
   }
   ML_free(rcv_neighbors);

   /* kludge in something so that the neighbors are in Aztec order. In */
   /* particular, make sure that those neighbors with lower rcv_lists  */
   /* appear first. Some of the aggregation code seems to be assuming  */
   /* this.                                                            */

   itemp = (int *) ML_allocate(sizeof(int)*(N_neighbors+1));
   rcv_neighbors = (int *) ML_allocate(sizeof(int)*(N_neighbors+1));
   if (rcv_neighbors == NULL) 
   {
      printf("ML_CommInfoOP_Generate: Not enough space\n");
      exit(1);
   }
   for (i = 0; i < N_neighbors; i++) itemp[i] = -1;
   for (i = 0; i < N_neighbors; i++) rcv_neighbors[i] = send_neighbors[i];
   largest_itemp = 100;
   for (i = 0; i < N_neighbors; i++) 
   {
      partner = rindex[i];
      if (partner != -1) 
      {
         if (rcv_length[partner] != -1) 
         {
            itemp[i] = rcv_list[partner][0];
            if (largest_itemp < itemp[i]) largest_itemp = itemp[i];
         }
      }
   }
   for (i = 0; i < N_neighbors; i++) 
      if (itemp[i] == -1) itemp[i] = largest_itemp+10;
   ML_az_sort(itemp, N_neighbors, rcv_neighbors,NULL);

   ML_CommInfoOP_Set_neighbors(comm_info, N_neighbors, 
			       rcv_neighbors, ML_OVERWRITE, NULL, 0);
/*
			       send_neighbors, ML_OVERWRITE, NULL, 0);
*/
   ML_free(rcv_neighbors); ML_free(itemp);

   j = 0;
   for (i = 0; i < N_send_procs; i++) 
   {
      partner = rindex[i];
      if (partner != -1)
         ML_CommInfoOP_Set_exch_info(*comm_info, send_neighbors[i],
		      rcv_length[partner], rcv_list[partner],send_length[i],
		      &(send_list[j]));
      else 
         ML_CommInfoOP_Set_exch_info(*comm_info, send_neighbors[i],
				 0, NULL,send_length[i],&(send_list[j]));
      j += send_length[i];
   }
   for (i = N_send_procs; i < N_neighbors; i++) 
   {
      partner = rindex[i];
      ML_CommInfoOP_Set_exch_info(*comm_info, send_neighbors[i], 
				  rcv_length[partner],rcv_list[partner],0,NULL);
   }
   ML_free(rindex);
   ML_free(send_list);
   ML_free(send_length);
   ML_free(send_neighbors);
   for (i = 0; i < N_rcv_procs; i++) ML_free(rcv_list[i]);
   ML_free(rcv_list);
   ML_free(rcv_length);
   return 1;
}

/* ******************************************************************** */
/* constructor                                                          */
/* -------------------------------------------------------------------- */

ML_CommInfoOP *ML_CommInfoOP_Create()
{
   ML_CommInfoOP *comm_info;

   comm_info = (ML_CommInfoOP *) ML_allocate(sizeof(ML_CommInfoOP));
   comm_info->remap = NULL;
   comm_info->neighbors = NULL;
   comm_info->N_neighbors = 0;
   comm_info->add_rcvd = 0;
   comm_info->total_rcv_length = 0;
   comm_info->minimum_vec_size = 0;
   comm_info->remap_length = 0;
   comm_info->remap_max = 0;
   return (comm_info);
}

/* ******************************************************************** */
/* destructor                                                           */
/* -------------------------------------------------------------------- */

void ML_CommInfoOP_Destroy(ML_CommInfoOP **comm_info)
{
   int i;
   ML_CommInfoOP *c_info;

   c_info = *comm_info;
   if (c_info != NULL) 
   {
      if (c_info->remap != NULL) ML_free(c_info->remap);
      for (i = 0; i < c_info->N_neighbors; i++)
	  {
         if (c_info->neighbors[i].rcv_list != NULL)
            ML_free(c_info->neighbors[i].rcv_list);
         if (c_info->neighbors[i].send_list != NULL)
            ML_free(c_info->neighbors[i].send_list);
      }
      if (c_info->neighbors != NULL)
         ML_free(c_info->neighbors);
      ML_free(c_info);
      c_info = NULL;
      *comm_info = NULL;
   }
}

/* ******************************************************************** */

int ML_CommInfoOP_Get_Nneighbors(ML_CommInfoOP *c_info)
{
   if (c_info == NULL) return(0);
   else return(c_info->N_neighbors);
}

/* ******************************************************************** */

int *ML_CommInfoOP_Get_neighbors(ML_CommInfoOP *c_info)
{
   int *itemp, i;
   
   if (c_info == NULL) return(NULL);
   else 
   {
      itemp = (int *) ML_allocate(sizeof(int)*c_info->N_neighbors);
      if ((itemp == NULL) && (c_info->N_neighbors != 0))
         pr_error("ML_CommInfoOP_Get_neighbors: no space\n");
      
      for (i = 0; i < c_info->N_neighbors; i++) 
         itemp[i] = c_info->neighbors[i].ML_id;
      return(itemp);
   }
}

/* ******************************************************************** */

int *ML_CommInfoOP_Get_sendlist(ML_CommInfoOP *c_info, int neighbor)
{
   int *itemp, i, j;
   
   if (c_info == NULL) return(NULL);
   else 
   {
      for (i = 0; i < c_info->N_neighbors; i++) 
         if ( c_info->neighbors[i].ML_id == neighbor) break;
      if (i == c_info->N_neighbors) return(NULL);

      itemp = (int *) ML_allocate(sizeof(int)*c_info->neighbors[i].N_send);
      for (j = 0; j < c_info->neighbors[i].N_send; j++) 
         itemp[j] = c_info->neighbors[i].send_list[j];
      return(itemp);
   }
}

/* ******************************************************************** */

int ML_CommInfoOP_Get_Nsendlist(ML_CommInfoOP *c_info, int neighbor)
{
   int i;

   if (c_info == NULL) return(0);
   else 
   {
      for (i = 0; i < c_info->N_neighbors; i++) 
         if ( c_info->neighbors[i].ML_id == neighbor) break;
      if (i == c_info->N_neighbors) return(0);
      else return(c_info->neighbors[i].N_send);
   }
}

/* ******************************************************************** */

int ML_CommInfoOP_Clone(ML_CommInfoOP **newone, ML_CommInfoOP *oldone)
{
   int i, *neighbors;

   if (oldone == NULL) 
   {
      *newone = NULL;
      return(0);
   }
   neighbors =(int *) ML_allocate((oldone->N_neighbors+1)*sizeof(int));
   if ( neighbors == NULL) 
   {
      printf("Not enough space in ML_CommInfoOP_Clone\n"); exit(1);
   }
   for (i = 0; i < oldone->N_neighbors; i++)
      neighbors[i] = oldone->neighbors[i].ML_id;

   ML_CommInfoOP_Set_neighbors(newone, oldone->N_neighbors,
      neighbors, oldone->add_rcvd, oldone->remap, oldone->remap_length);
   ML_free(neighbors);

   for (i = 0; i < oldone->N_neighbors; i++) 
   {
      ML_CommInfoOP_Set_exch_info(*newone, oldone->neighbors[i].ML_id,
                                  oldone->neighbors[i].N_rcv,
                                  oldone->neighbors[i].rcv_list,
                                  oldone->neighbors[i].N_send,
                                  oldone->neighbors[i].send_list);
   }
   return(1);
}

/* ******************************************************************** */

int ML_CommInfoOP_Print(ML_CommInfoOP *c_info, char *label)
{
   int i,j;

   if (c_info == NULL) return(0);
   printf("%s :: Number of neighbors = %d\n",label,c_info->N_neighbors);
   for (i = 0; i < c_info->N_neighbors; i++) 
   {
       printf("%s :: %dth neighbor = %4d (N_send = %4d, N_rcv = %4d)\n",
              label, i+1, c_info->neighbors[i].ML_id, 
	      c_info->neighbors[i].N_send, c_info->neighbors[i].N_rcv);
       for (j = 0; j < c_info->neighbors[i].N_send; j++) 
       {
          printf("%s ::      send(%d) = %d\n",label,j,
			c_info->neighbors[i].send_list[j]);
       }
       if (c_info->neighbors[i].rcv_list != NULL) 
       {
          for (j = 0; j < c_info->neighbors[i].N_rcv; j++) 
          {
             printf("%s ::      rcv(%d) = %d\n",label,j,
			c_info->neighbors[i].rcv_list[j]);
          }
       }
   }
   return(1);
}

/* ******************************************************************** */

int ML_CommInfoOP_Get_Nrcvlist(ML_CommInfoOP *c_info, int neighbor)
{
   int i;

   if (c_info == NULL) return(0);
   else 
   {
      for (i = 0; i < c_info->N_neighbors; i++) 
         if ( c_info->neighbors[i].ML_id == neighbor) break;
      if (i == c_info->N_neighbors) return(0);
      else return(c_info->neighbors[i].N_rcv);
   }
}

/* ******************************************************************** */

int *ML_CommInfoOP_Get_rcvlist(ML_CommInfoOP *c_info, int neighbor)
{
   int *itemp, i, j;
   
   if (c_info == NULL) return(NULL);
   else 
   {
      for (i = 0; i < c_info->N_neighbors; i++) 
         if ( c_info->neighbors[i].ML_id == neighbor) break;
      if (i == c_info->N_neighbors) return(NULL);

      if (NULL == c_info->neighbors[i].rcv_list) return(NULL);
      itemp = (int *) ML_allocate(sizeof(int)*c_info->neighbors[i].N_rcv);
      for (j = 0; j < c_info->neighbors[i].N_rcv; j++) {
         itemp[j] = c_info->neighbors[i].rcv_list[j];
      }
      return(itemp);
   }
}

/* ******************************************************************** */
/*
 * Set the communication data structure to reflect the number of
 * neighbors with which to communicate and whether or not received
 * information is to be added into existing values or whether it should
 * overwrite existing values.
 *
 * Parameters
 * ==========
 *
 * matrix         On input, matrix to which communication information
 *                will be associated. On output, basic neighbor
 *                information will be stored with matrix object.
 *
 * N_neighbor     On input, number of processors this node must
 *                communicate with to exchange vector components.
 *
 * neighbors      On input, processor ids of neighbors.
 *
 * add_or_not     On input, indicates whether received values add or
 *                overwrite existing vector components. Normally, for
 *                matrices distributed by rows overwriting is
 *                appropriate (ML_OVERWRITE). For column distributed
 *                matrices adding is often necessary (ML_ADD).
 *
 * remap          On input, defines a mapping for vector components
 *                after communication. Normally, remap is set to NULL
 *                (indicating no re-mapping) for matrices distributed
 *                by rows. However for column distributed matrices,
 *                re-mapping is often done to indicate values that
 *                are no longer needed.
 *
 */

int ML_CommInfoOP_Set_neighbors(ML_CommInfoOP **c_info, int N_neighbors,
	int *neighbors, int add_or_not, int *remap, int remap_length)
{
  int i;
  ML_CommInfoOP *comm_info;

  if (*c_info != NULL) 
  {
     printf("ML_CommInfoOP_Set_neighbors: c_info not NULL! Does communication structure already exist?\n");
     exit(1);
  }

  comm_info = (ML_CommInfoOP *) ML_allocate(sizeof(ML_CommInfoOP));
  *c_info   = comm_info;
  comm_info->total_rcv_length   = 0;
  comm_info->minimum_vec_size   = 0;
  comm_info->remap_length       = 0;
  comm_info->remap_max          = -1;

  comm_info->N_neighbors = N_neighbors;
  if (N_neighbors ==0) comm_info->neighbors = NULL;
  else
    comm_info->neighbors = (ML_NeighborList *) ML_allocate(N_neighbors*
                                              sizeof(ML_NeighborList));
  if ((N_neighbors != 0) && (comm_info->neighbors == NULL)) 
  {
      printf("Out ot memory in ML_CommInfoOP_Set_neighbors\n");
      exit(1);
  }

  for (i = 0; i < N_neighbors; i++ ) 
  {
      comm_info->neighbors[i].ML_id  = neighbors[i];
      comm_info->neighbors[i].N_send = 0;
      comm_info->neighbors[i].N_rcv  = 0;
      comm_info->neighbors[i].send_list = NULL;
      comm_info->neighbors[i].rcv_list  = NULL;
  }

  if ((add_or_not != ML_OVERWRITE) && (add_or_not != ML_ADD)) 
  {
      printf("ML_CommInfoOP_Set_neighbors: Invalid value for 'add_or_not'\n");
      exit(1);
  }

  comm_info->add_rcvd = add_or_not;
  if (remap != NULL) 
  {
      comm_info->remap_length = remap_length;
      comm_info->remap = (int *) ML_allocate(sizeof(int)*(remap_length+1));
      if (comm_info->remap == NULL) 
      {
          printf("Error: Not enough space for remap vector of length = %d\n",
                  remap_length);
          exit(1);
      }
      for (i = 0; i < remap_length; i++) 
      {
         comm_info->remap[i] = remap[i];
         if (remap[i] > comm_info->remap_max)
            comm_info->remap_max = remap[i];
      }
  }
  else comm_info->remap = NULL;
  return 1;
}

/* ******************************************************************** */
/*
 * Fill in the communication information (Number of doubles to send and
 * receive, the list of indices to send, and the list of indices where
 * where received information is stored) for the neighbor with id = k.
 * NOTE: the rcv_list can be NULL in which case received information
 * is appended to the end of the vector.
 *
 * Parameters
 * ==========
 *   comm_info     On input, comm_info->N_neighbors is set to the number of
 *                 neighbors and comm_info->neighbors[0..N_neighbors-1].id' is
 *                 set to the neighbor ids.
 *                 On output, comm_info->neighbors[] is filled for
 *                 the neighbor id corresponding to k.
 *
 *   k             On input, id of neighbor to be filled in.
 *
 *   N_rcv         On input, number of doubles to receive.
 *
 *   rcv_list      On input, list of indices where received information from
 *                 processor k is stored. NOTE: if NULL, received information
 *                 is appended to the end of the vector.
 *
 *   N_send        On input, number of doubles to send.
 *
 *   send_list     On input, list of indices to send to processor k.
 *
 *   pre_or_post   On input, == ML_PRE  ==> matrix->pre_comm is filled
 *                           == ML_POST ==> matrix->post_comm is filled
 */

int ML_CommInfoOP_Set_exch_info(ML_CommInfoOP *comm_info, int k,
        int N_rcv, int *rcv_list, int N_send, int *send_list)
{
   int oldone, i, j, *list;

   if (comm_info == NULL) 
   {
      printf("ML_CommInfoOP_Set_exch_info: communication structure \
              does not exist.\n");
      exit(1);
   }
   for (i = 0; i < comm_info->N_neighbors; i++ ) 
   {
      if ( comm_info->neighbors[i].ML_id == k) break;
   }

   if (i >= comm_info->N_neighbors) 
   {
      printf("Error: neighbor (%d) not found\n",k);
      exit(1);
   }
   oldone = comm_info->neighbors[i].N_rcv;
   comm_info->total_rcv_length   += (N_rcv - comm_info->neighbors[i].N_rcv);
   comm_info->neighbors[i].N_rcv  = N_rcv;
   comm_info->neighbors[i].N_send = N_send;

   if (N_send > 0) 
   {
      list = (int *) ML_allocate(N_send*sizeof(int));
      for (j = 0; j < N_send; j++) 
      {
         if (send_list[j] >= comm_info->minimum_vec_size)
            comm_info->minimum_vec_size = send_list[j] + 1;
         list[j] = send_list[j];
      }
      if (comm_info->neighbors[i].send_list != NULL) 
         ML_free(comm_info->neighbors[i].send_list);
      comm_info->neighbors[i].send_list = list;
   }
   else comm_info->neighbors[i].send_list = NULL;

   if ((N_rcv > 0) && (rcv_list != NULL)) 
   {
      list = (int *) ML_allocate(N_rcv*sizeof(int));
      for (j = 0; j < N_rcv; j++) 
      {
         if (rcv_list[j] >= comm_info->minimum_vec_size)
            comm_info->minimum_vec_size = rcv_list[j] + 1;
         list[j] = rcv_list[j];
      }
      if (comm_info->neighbors[i].rcv_list != NULL) 
         ML_free(comm_info->neighbors[i].rcv_list);
      comm_info->neighbors[i].rcv_list = list;
   }
   else 
   {
      comm_info->neighbors[i].rcv_list = NULL;
      comm_info->minimum_vec_size   += (N_rcv - oldone);
   }
   return 1;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ML_create_unique_id(int N_local, int **map, 
                ML_CommInfoOP *comm_info, ML_Comm *comm)
{
/* Create a map between local variables on this processor and a unique
 * global number where local variables on different processors which
 * correspond to the same global variable have the same unique global number.
 *
 * Parameters
 * ==========
 *   N_local       On input, number of local variables assigned to this node.
 * 
 *   map           On output, map[k] is the unique global id of the kth local
 *                 variable. Note: if the kth local variable on processor P0
 *                 corresponds to the jth local variable on processor P1, then
 *                 map[k] on P0 is equal to map[j] on P1.
 *   
 *   comm_info     On input, communcation information (see ml_rap.h) which
 *                 indicates which local variables are sent to other processors
 *                 and where received information is stored locally.
 *
 *   max_per_proc  On output, the maximum number of local variables on any one
 *                 processor.
 */

   int i, j, count, N_rcvd, N_send, offset, flag = 0;
   double *dtemp;

   /* compute the number of variables to receive and send */

   N_rcvd = 0;
   N_send = 0;
   if (comm_info != NULL) 
   {
      for (i = 0; i < comm_info->N_neighbors; i++)  
      {
         N_rcvd += (comm_info->neighbors)[i].N_rcv;
         N_send += (comm_info->neighbors)[i].N_send;
         if (  ((comm_info->neighbors)[i].N_rcv != 0) &&
            ((comm_info->neighbors)[i].rcv_list != NULL) )  flag = 1;
      }
   }

   dtemp  = (double *) ML_allocate((N_local + N_rcvd + 1)*sizeof(double));
   if (dtemp == NULL) 
   {
     printf("out of space in ML_create_unique_col_ids\n");
     exit(1);
   }

   /* Set the N_local components of 'map' and 'dtemp' */
   /* to unique numbers on each processor.            */


   offset = ML_gpartialsum_int(N_local, comm);

   *map = (int    *) ML_allocate((N_local + N_rcvd + 1) * sizeof(int));
   if (map == NULL) 
   {
      printf("out of space in ML_create_unique_col_ids\n");
      exit(1);
   }
   for (i = 0 ; i < N_local; i++ ) 
   {
      (*map)[i]    = offset + i;
      dtemp[i] = (double) (*map)[i];
   }

   /* exchange these global ids with the neighbors, appending */
   /* received information (starting at dtemp[N_local])       */

   if (comm_info != NULL)
   {
      ML_cheap_exchange_bdry(dtemp, comm_info, N_local, N_send, comm);
   }

   if (flag == 1) 
   {
      count = N_local;
      for (i = 0; i < comm_info->N_neighbors; i++) 
      {
         for (j = 0; j < comm_info->neighbors[i].N_rcv; j++) 
         {
            (*map)[comm_info->neighbors[i].rcv_list[j]] = (int) dtemp[count++];
         }
      }
   }
   else 
      for (i = N_local; i < N_local + N_rcvd; i++ ) (*map)[i] = (int) dtemp[i];

   ML_free(dtemp);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ML_create_unique_col_id(int N_local, int **map, 
                ML_CommInfoOP *comm_info, int *max_per_proc, ML_Comm *comm)
{
/* Create a map between local variables on this processor and a unique
 * global number where local variables on different processors which
 * correspond to the same global variable have the same unique global number.
 *
 * Parameters
 * ==========
 *   N_local       On input, number of local variables assigned to this node.
 * 
 *   map           On output, map[k] is the unique global id of the kth local
 *                 variable. Note: if the kth local variable on processor P0
 *                 corresponds to the jth local variable on processor P1, then
 *                 map[k] on P0 is equal to map[j] on P1.
 *   
 *   comm_info     On input, communcation information (see ml_rap.h) which
 *                 indicates which local variables are sent to other processors
 *                 and where received information is stored locally.
 *
 *   max_per_proc  On output, the maximum number of local variables on any one
 *                 processor.
 */

   int i, j, count, N_rcvd, N_send, offset, flag = 0;
   double *dtemp;

   /* compute the number of variables to receive and send */

   N_rcvd = 0;
   N_send = 0;
   if (comm_info != NULL) 
   {
      for (i = 0; i < comm_info->N_neighbors; i++)  
      {
         N_rcvd += (comm_info->neighbors)[i].N_rcv;
         N_send += (comm_info->neighbors)[i].N_send;
         if (  ((comm_info->neighbors)[i].N_rcv != 0) &&
            ((comm_info->neighbors)[i].rcv_list != NULL) )  flag = 1;
      }
   }

   dtemp  = (double *) ML_allocate((N_local + N_rcvd + 1)*sizeof(double));
   if (dtemp == NULL) 
   {
     printf("out of space in ML_create_unique_col_ids\n");
     exit(1);
   }

   /* Set the N_local components of 'map' and 'dtemp' */
   /* to unique numbers on each processor.            */

   *max_per_proc = ML_gmax_int(N_local, comm);
   offset       = *max_per_proc*(comm->ML_mypid);

   *map = (int    *) ML_allocate((N_local + N_rcvd + 1) * sizeof(int));
   if (map == NULL) 
   {
      printf("out of space in ML_create_unique_col_ids\n");
      exit(1);
   }
   for (i = 0 ; i < N_local; i++ ) 
   {
      (*map)[i]    = offset + i;
      dtemp[i] = (double) (*map)[i];
   }

   /* exchange these global ids with the neighbors, appending */
   /* received information (starting at dtemp[N_local])       */

   if (comm_info != NULL)
   {
      ML_cheap_exchange_bdry(dtemp, comm_info, N_local, N_send, comm);
   }

   if (flag == 1) 
   {
      count = N_local;
      for (i = 0; i < comm_info->N_neighbors; i++) 
      {
         for (j = 0; j < comm_info->neighbors[i].N_rcv; j++) 
         {
            (*map)[comm_info->neighbors[i].rcv_list[j]] = (int) dtemp[count++];
         }
      }
   }
   else 
      for (i = N_local; i < N_local + N_rcvd; i++ ) (*map)[i] = (int) dtemp[i];

   ML_free(dtemp);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ML_cheap_exchange_bdry(double x[], ML_CommInfoOP *comm_info, 
                      int start_location, int total_send, ML_Comm *comm)

/*******************************************************************************

  Routine to locally exchange components of the vector "x". This routine
  gathers necessary components of the vector and sends the required
  messages. Messages which are received are placed contiguously starting
  from location x[start_location]. 

  Author:          Ray Tuminaro, SNL, 9222
  =======

  Parameter list:
  ===============

  x                On input, vector of unknowns defined on current processor.
                   On output, x is appended with information received from
                   other processors (dictated by 'comm_info').

  comm_info        On input, comm_info contains information on what needs to be
                   exchanged with other processors. See ml_rap.h.

  start_location   On input, starting location in 'x' where received information
                   will be placed. 

  total_send       On input, total number of components in 'x' to be sent.
   
*******************************************************************************/

{
  register double *ptrd;
  double          *ptr_send_list, *ptr_recv_list, *orig_ptr;
  int              type, N_neighbors, *temp, i, j, rtype;
  USR_REQ         *request;
  ML_NeighborList *neighbor;

  /**************************** execution begins ******************************/

  N_neighbors = comm_info->N_neighbors;
  if (N_neighbors == 0) return;

  /* Set up send messages: Gather send unknowns from "x" vector */

  request = (USR_REQ *)  ML_allocate(N_neighbors*sizeof(USR_REQ));
  ptrd = (double *) ML_allocate( (total_send+1)*sizeof(double));
  if (ptrd == NULL) 
  {
     printf("Out of space in ML_cheap_exchange_bdry\n");
     exit(1);
  }
  ptr_send_list = ptrd;
  orig_ptr      = ptrd;

  for (i = 0; i < N_neighbors; i++) 
  {
     temp = comm_info->neighbors[i].send_list;
     for (j = 0; j < comm_info->neighbors[i].N_send; j++) 
     {
        *ptrd++ = x[ temp[j] ];
     }
  }

  type = 1991;

  /* post receives for all messages */

  ptr_recv_list = &x[start_location];
  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_rcv;
    comm->USR_irecvbytes((void *) ptr_recv_list, (unsigned int)j, 
		&(neighbor->ML_id), &rtype, comm->USR_comm, request+i);
    ptr_recv_list += neighbor->N_rcv;
  }

  /* write out all messages */

  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    j = sizeof(double)* neighbor->N_send;
    comm->USR_sendbytes((void *) ptr_send_list, (unsigned) j, neighbor->ML_id, 
                          rtype, comm->USR_comm);
    ptr_send_list += neighbor->N_send;
  }

  /* wait for all messages */

  ptr_recv_list = &x[start_location];
  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_rcv;
    comm->USR_cheapwaitbytes((void *) ptr_recv_list, j, &(neighbor->ML_id),
                        &rtype, comm->USR_comm, request+i);
    ptr_recv_list += neighbor->N_rcv;
  }

  ML_free(orig_ptr);
  ML_free(request);

} /* ML_cheap_exchange_bdry */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*******************************************************************************

  Routine to locally exchange components of the vector "x". This routine
  gathers necessary components of the vector and sends the required
  messages. Messages which are received are placed contiguously starting
  from location x[start_location]. 

  Author:          Ray Tuminaro, SNL, 9222
  =======

  Parameter list:
  ===============

  x                On input, vector of unknowns defined on current processor.
                   On output, x is appended with information received from
                   other processors (dictated by 'comm_info').

  comm_info        On input, comm_info contains information on what needs to be
                   exchanged with other processors. See ml_rap.h.

  start_location   On input, starting location in 'x' where received information
                   will be placed. 

  envelope         On input, contains MPI message tag.

*******************************************************************************/

void ML_exchange_bdry(double x[], ML_CommInfoOP *comm_info, int start_location, 
	ML_Comm *comm, int overwrite_or_add, ML_Comm_Envelope *envelope)
{
  double          *send_buf, **rcv_buf, *tempv;
  int              type, N_neighbors, *temp, i, j, k, rtype;
  USR_REQ         *request;
  ML_NeighborList *neighbor;

  /**************************** execution begins ******************************/
  if (comm_info == NULL) return;
  N_neighbors              = comm_info->N_neighbors;
  if (N_neighbors == 0) return;

  /* Set up send messages: Gather send unknowns from "x" vector */

  if ( N_neighbors > 0 )
  {
     request = (USR_REQ  *)  ML_allocate(N_neighbors*sizeof(USR_REQ ));
     rcv_buf = (double  **)  ML_allocate(N_neighbors*sizeof(double *));
  } else { request = NULL; rcv_buf = NULL;}

  /*
  type = 1991;
  */
  if ( ML_Comm_Envelope_Get_Tag(envelope,&type) )
     type = 1991;

  /* post receives for all messages */

  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_rcv;
    rcv_buf[i] = (double *)  ML_allocate(j);
    comm->USR_irecvbytes((void *) rcv_buf[i], (unsigned int)j, 
		&(neighbor->ML_id), &rtype, comm->USR_comm, request+i);
  }

  /* write out all messages */

  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    j = sizeof(double)* neighbor->N_send;
    send_buf = (double *)  ML_allocate(j);
    temp = comm_info->neighbors[i].send_list;
    for (k = 0; k < neighbor->N_send; k++) 
    {
        send_buf[k] = x[ temp[k] ];
    }
    comm->USR_sendbytes((void *) send_buf, (unsigned) j, neighbor->ML_id, 
                          rtype, comm->USR_comm);
    if (send_buf != NULL) ML_free(send_buf);
  }

  /* wait for all messages */

  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_rcv;
    comm->USR_cheapwaitbytes((void *) rcv_buf[i], j, &(neighbor->ML_id),
                        &rtype, comm->USR_comm, request+i);
    temp = comm_info->neighbors[i].rcv_list;
    if (temp == NULL) 
    {
       if (overwrite_or_add == ML_ADD) 
       {
          for (k = 0; k < neighbor->N_rcv; k++) 
          {
             x[ start_location ] = rcv_buf[i][k] + x[start_location];
             start_location++;
          }
       }
       else 
       {
          for (k = 0; k < neighbor->N_rcv; k++) 
          {
             x[ start_location ] = rcv_buf[i][k];
             start_location++;
          }
       }
    }
    else 
    {
       if (overwrite_or_add == ML_ADD) 
       {
          for (k = 0; k < neighbor->N_rcv; k++) 
          {
             x[ temp[k] ] = rcv_buf[i][k] + x[temp[k]];
          }
       }
       else 
       {
          for (k = 0; k < neighbor->N_rcv; k++) 
          {
             x[ temp[k] ] = rcv_buf[i][k];
          }
       }
    }
    if (rcv_buf[i] != NULL) ML_free(rcv_buf[i]);
  }
  if ( N_neighbors > 0 ) ML_free(rcv_buf);
  if ( N_neighbors > 0 ) ML_free(request);

  if (comm_info->remap != NULL) 
  {
     tempv = (double *) ML_allocate((comm_info->remap_max+1)*sizeof(double));
     for (k = 0; k < comm_info->remap_length; k++) 
     {
        j = comm_info->remap[k];
        if (j >= 0) tempv[ j ] = x[k];
     }
     for (i = 0; i < comm_info->remap_max; i++) x[i] = tempv[i];
     ML_free(tempv);
  }
} /* ML_exchange_bdry */

int ML_CommInfoOP_Compute_TotalRcvLength(ML_CommInfoOP *comm_info)
{
  int i;

  if (comm_info == NULL) return 0;

  comm_info->total_rcv_length = 0;
  for (i = 0; i < comm_info->N_neighbors; i++ ) {
   comm_info->total_rcv_length   += comm_info->neighbors[i].N_rcv;
  }
  return comm_info->total_rcv_length;
} 

/******************************************************************************/

void ML_transposed_exchange_bdry(double x[], ML_CommInfoOP *comm_info,
    int start_location, ML_Comm *comm, int overwrite_or_add)
{
  double          *send_buf, **rcv_buf;
  int              type, N_neighbors, *temp, i, j, k, rtype;
  USR_REQ         *request;
  ML_NeighborList *neighbor;

  /**************************** execution begins ******************************/

  N_neighbors              = comm_info->N_neighbors;
  if (N_neighbors == 0) { 
    /* bogus code to avoid warning */
    ML_avoid_unused_param((void *) &start_location);
    return;
  }

  /* Set up send messages: Gather send unknowns from "x" vector */

  if ( N_neighbors > 0 )
  {
     request = (USR_REQ  *)  ML_allocate(N_neighbors*sizeof(USR_REQ ));
     rcv_buf = (double  **)  ML_allocate(N_neighbors*sizeof(double *));
  } else { request = NULL; rcv_buf = NULL;}

  type = 2001;

  /* post receives for all messages */

  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_send;
    rcv_buf[i] = (double *)  ML_allocate(j);
    comm->USR_irecvbytes((void *) rcv_buf[i], (unsigned int)j, 
		&(neighbor->ML_id), &rtype, comm->USR_comm, request+i);
  }

  /* write out all messages */

  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    j = sizeof(double)* neighbor->N_rcv;
    send_buf = (double *)  ML_allocate(j);
    temp = comm_info->neighbors[i].rcv_list;
	if (temp == NULL && j != 0)
	{
	   printf("In ML_transposed_exchange_bdry: comm_info->neighbors[i].rcv_list cannot be NULL\n");
	   exit(1);
    }
    for (k = 0; k < neighbor->N_rcv; k++) 
    {
        send_buf[k] = x[ temp[k] ];
    }
    comm->USR_sendbytes((void *) send_buf, (unsigned) j, neighbor->ML_id, 
                          rtype, comm->USR_comm);
    if (send_buf != NULL) ML_free(send_buf);
  }

  /* wait for all messages */

  for (i = 0; i < N_neighbors; i++) 
  {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_send;
    comm->USR_cheapwaitbytes((void *) rcv_buf[i], j, &(neighbor->ML_id),
                        &rtype, comm->USR_comm, request+i);
    temp = comm_info->neighbors[i].send_list;
    if (overwrite_or_add == ML_ADD) 
    {
      for (k = 0; k < neighbor->N_send; k++) 
        x[ temp[k] ] += rcv_buf[i][k];
    }
    else {
      for (k = 0; k < neighbor->N_send; k++) 
        x[ temp[k] ] = rcv_buf[i][k];
    }
    if (rcv_buf[i] != NULL) ML_free(rcv_buf[i]);
  }
  if ( N_neighbors > 0 ) ML_free(rcv_buf);
  if ( N_neighbors > 0 ) ML_free(request);

  if (comm_info->remap != NULL) 
  {
     printf("comm_info->remap != NULL\n");
     exit(1);
  }
} /* ML_transposed_exchange_bdry */

/*******************************************************************************
 Create data structure for message tags.
*******************************************************************************/
int ML_Comm_Envelope_Create(ML_Comm_Envelope** envelope)
{
   *envelope = (ML_Comm_Envelope*) ML_allocate( sizeof(ML_Comm_Envelope) );
   ML_Comm_Envelope_Init(*envelope);
   return 0;
}

/*******************************************************************************
 Initialize data structure for message tags.
*******************************************************************************/
int ML_Comm_Envelope_Init(ML_Comm_Envelope* envelope)
{
   envelope->tag = 0;
   return 0;
}

/*******************************************************************************
 Destroy data structure for message tags.
*******************************************************************************/
int ML_Comm_Envelope_Destroy(ML_Comm_Envelope* envelope)
{
   if (envelope != NULL)
      ML_free(envelope);
   return 0;
}

/*******************************************************************************
 Reset data structure for message tag but don't free it.
*******************************************************************************/
int ML_Comm_Envelope_Clean(ML_Comm_Envelope* envelope)
{
   if (envelope != NULL)
      envelope->tag = 0;
   return 0;
}

/*******************************************************************************
 Fetch message tag.
*******************************************************************************/
int ML_Comm_Envelope_Get_Tag(ML_Comm_Envelope* envelope, int* tag)
{
   if (envelope != NULL)
   {
      *tag = envelope->tag; 
      return 0;
   }
   else
      return 1;
}

/*******************************************************************************
 Set message tag.  Start with a base value of ML_TAG_BASE.  Each level will
 have a range of values, starting at ML_TAG_BASE +  sm->mylevel*100.
 This range is further broken into subranges for pre- and post-smoothing.
*******************************************************************************/
int ML_Comm_Envelope_Set_Tag(ML_Comm_Envelope* envelope,
                             int level, int pre_or_post) 
{
   envelope->tag = ML_TAG_BASE +  level*200 + pre_or_post;
   return 0;
}

/*******************************************************************************
 Increment message tag.
*******************************************************************************/
int ML_Comm_Envelope_Increment_Tag(ML_Comm_Envelope* envelope)
{
   envelope->tag++;
   return 0;
}

/********************************************************************
 * Take a standard 'pre_comm' communication object and transpose it
 * so that we can do post communication. This routine is used when
 * transposing a matrix. 
 *    On input:
 *        pre_comm     An existing communication structure that
 *                     will be transposed.
 *        invec_leng   The length of the local vector for which the
 *                     pre  communication is defined.
 *
 *    On output
 *        post_comm    A new communication structure that corresponds
 *                     to taking a matrix and transposing it (but 
 *                     keeping the transpose of the data local).
 *
 *    Returns the number of ghost nodes in the new communicator. 
 *
 ********************************************************************/
int ML_CommInfoOP_TransComm(ML_CommInfoOP *pre_comm, ML_CommInfoOP **post_comm,
			    int invec_leng)
{
   int osize, Nneighbors, remap_leng, Nrcv, Nsend, i, j;
   int *neigh_list, *send_list, *rcv_list, *remap;
   int Nghost = 0, Nghost2 = 0;

   osize = invec_leng;

   /* transpose communication list. This means that PRE communication */
   /* is replaced by POST, ML_OVERWRITE is replaced by ML_ADD, and the send  */
   /* send and receive lists are swapped.                                    */

   Nneighbors = ML_CommInfoOP_Get_Nneighbors(pre_comm);
   neigh_list = ML_CommInfoOP_Get_neighbors(pre_comm);
   remap_leng = osize;
   Nrcv = 0;
   Nsend = 0;
   for (i = 0; i < Nneighbors; i++) {
      Nrcv  += ML_CommInfoOP_Get_Nrcvlist (pre_comm, neigh_list[i]);
      Nsend += ML_CommInfoOP_Get_Nsendlist(pre_comm, neigh_list[i]);
   }
   remap_leng = osize + Nrcv + Nsend;
   remap = (int *) ML_allocate( remap_leng*sizeof(int));
   for (i = 0; i < osize; i++) remap[i] = i;
   for (i = osize; i < osize+Nrcv+Nsend; i++) 
      remap[i] = -1;
 
   ML_CommInfoOP_Set_neighbors(post_comm, Nneighbors,
 			      neigh_list,ML_ADD,remap,remap_leng);
   ML_free(remap);
   for (i = 0; i < Nneighbors; i++) {
      Nsend      = ML_CommInfoOP_Get_Nsendlist(pre_comm, neigh_list[i]);
      send_list  = ML_CommInfoOP_Get_sendlist (pre_comm, neigh_list[i]);
      Nrcv       = ML_CommInfoOP_Get_Nrcvlist (pre_comm, neigh_list[i]);
      Nghost    += Nrcv;
      rcv_list   = ML_CommInfoOP_Get_rcvlist(pre_comm, neigh_list[i]);
      /* handle empty rows ... i.e. ghost variables not used */
      if (rcv_list != NULL) {
         for (j = 0; j < Nrcv; j++) {
            if (rcv_list[j] > Nghost2 + osize - 1)
               Nghost2 = rcv_list[j] - osize + 1;
         }
      }
 
      ML_CommInfoOP_Set_exch_info(*post_comm, neigh_list[i], Nsend, send_list,
 				 Nrcv,rcv_list);
      if (send_list != NULL) ML_free(send_list);
      if ( rcv_list != NULL) ML_free( rcv_list);
   }
   if (neigh_list != NULL) ML_free(neigh_list);
   if (Nghost2 > Nghost) Nghost = Nghost2;

   return Nghost;
}

int ML_reverse_exchange(double *x_over, ML_CommInfoOP 
			*nonOverlapped_2_Overlapped,
			int nn_over , ML_Comm *comm)
  
{
  int i, j, *itemp;

  /* flip send and rcv lists */

  for (i = 0; i < nonOverlapped_2_Overlapped->N_neighbors; i++ ) {
    j = nonOverlapped_2_Overlapped->neighbors[i].N_rcv;
    itemp = nonOverlapped_2_Overlapped->neighbors[i].rcv_list;
    nonOverlapped_2_Overlapped->neighbors[i].N_rcv = 
      nonOverlapped_2_Overlapped->neighbors[i].N_send;
    nonOverlapped_2_Overlapped->neighbors[i].rcv_list = 
      nonOverlapped_2_Overlapped->neighbors[i].send_list;
    nonOverlapped_2_Overlapped->neighbors[i].N_send = j;
    nonOverlapped_2_Overlapped->neighbors[i].send_list = itemp;
  }
  ML_exchange_bdry(x_over,nonOverlapped_2_Overlapped,nn_over,comm,ML_ADD,NULL);

  /* flip back */

  for (i = 0; i < nonOverlapped_2_Overlapped->N_neighbors; i++ ) {
    j = nonOverlapped_2_Overlapped->neighbors[i].N_rcv;
    itemp = nonOverlapped_2_Overlapped->neighbors[i].rcv_list;
    nonOverlapped_2_Overlapped->neighbors[i].N_rcv = 
      nonOverlapped_2_Overlapped->neighbors[i].N_send;
    nonOverlapped_2_Overlapped->neighbors[i].rcv_list = 
      nonOverlapped_2_Overlapped->neighbors[i].send_list;
    nonOverlapped_2_Overlapped->neighbors[i].N_send = j;
    nonOverlapped_2_Overlapped->neighbors[i].send_list = itemp;
  }

  return 0;

}
