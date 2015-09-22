/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "az_aztec.h"

#define ALLOCATE      1
#define FREE_IT       2
#define UPDATE_IT     3
#define V_SET         4
#define B_SET         5
#define LASTUSED_SET  6
#define ESTIMATED_SET 7
 
 
#define OUT_OF_BOUNDS 0
#define IN_V          1
#define IN_B          2
 
#define SET_VAL(a,b)   AZ_allocate_or_free((void *) (a),(unsigned int)(b),V_SET)
#define SET_BINDX(a,b) AZ_allocate_or_free((void *)(a),(unsigned int)(b),B_SET)
#define SET_USAGE(a,b) AZ_allocate_or_free(NULL,(unsigned int)(a),LASTUSED_SET); \
                       AZ_allocate_or_free(NULL,(unsigned int) (b),ESTIMATED_SET)
#define BV_ALLOC(a)    AZ_allocate_or_free(NULL,(a),ALLOCATE)
#define BV_FREE(a)     AZ_allocate_or_free((void *) (a),(long int) NULL,FREE_IT)
#define RESET_USAGE(a) AZ_allocate_or_free(NULL,(unsigned int)(a),LASTUSED_SET);

extern void AZ_setup_sendlist(int , int *, int *, int *, int **,
                              int *, int , int , int *);

/******************************************************************************/
/*******************************************************************************

  Setup subroutine for overlapping domain decomposition preconditioner
  (currently for matrices with MSR format only).

  Author:          Charles H. Tong (8117)
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  Input :

    N_rows:        number of matrix rows residing locally
    Rownum:        an array containing the row numbers of the local rows
    bindx :        an index array for the DMSR sparse matrix storage 
    val:           an array containing the nonzero entries of the matrix 
    olap_size:     the degree of overlap requested
    proc_config:   Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

Output :
    New_N_rows:    the number of matrix rows for the enlarged matrix
    bindx:         an index array for the DMSR storage of the enlarged matrix
    val:           an array containing the nonzeros of the enlarged matrix


*******************************************************************************/

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

void AZ_setup_dd_olap_msr(int N_rows, int *New_N_rows, int *bindx, double *val,
   int olap_size, int *proc_config, int *data_org[], int **map3,
   int bindx_length,int name, int *prev_data_org, int estimated_requirements,
   struct context *context)
{
  /* ----------------- variable and function declaration ---------------- */
  /* local variables */

  int    i,j,k,jj,m,ndist,ind,ncnt,nzeros,ret_index,off_set,type;
  int    nprocs,status,num_neigh,*proc_neigh;
  int    enlarged_N, *enlarged_Rownum, *sorted_New_Rownum;

  int    *send_leng, *actual_send_leng, **send_rownum, *send_proc;
  int    *recv_leng, *actual_recv_leng, *recv_proc, *send_rowcnt;
  int    *start_send;

  int    N_ext_node, *ext_nodelist, *sorted_exts, *externals;
  int    *Rownum, *tRownum, net_gain_N;
  int    offset, max_per_proc;

  double *comm_buffer = NULL, *dptr = NULL, *dbl_grows;
  int    comm_size = 0, *iptr = NULL;

  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */

char str[80];

  /* ------------------------- subroutine begins ------------------------ */

  /* Initialize overlap memory scheme */

  enlarged_N = bindx[0] - 1;
  nzeros     = bindx[enlarged_N];
  if (bindx_length < nzeros) {
      AZ_printf_err("Not enough memory allocated for overlapping\n");
      exit(1);
  }
  SET_VAL(val    , bindx_length);
  SET_BINDX(bindx, bindx_length);
  if (estimated_requirements > nzeros) {
     SET_USAGE(nzeros,estimated_requirements);
  }
  else { SET_USAGE(nzeros,bindx_length);}

  /* tranform the local matrices to a global matrix */
  /* by putting in fake numbers                     */

  max_per_proc = AZ_gmax_int(N_rows,proc_config);
  max_per_proc++;
  offset       = max_per_proc*proc_config[AZ_node];
  
  /* determine the global row number corresponding to */ 
  /* external variables using AZ_exchange_bdry().     */

  Rownum    = (int *)    BV_ALLOC((N_rows+1) * sizeof(int));
  externals = (int *)    BV_ALLOC((prev_data_org[AZ_N_external]+1)*sizeof(int));
  dbl_grows = (double *) BV_ALLOC((N_rows + prev_data_org[AZ_N_external] + 1)* 
                                   sizeof(double) );
     
  for (i = 0 ; i < N_rows; i++ ) {
     Rownum[i]    = offset + i;
     dbl_grows[i] = (double) Rownum[i];
  }
  AZ_exchange_bdry(dbl_grows, prev_data_org, proc_config);

  for (i = 0 ; i < prev_data_org[AZ_N_external] ; i++ )
     externals[i] = (int) dbl_grows[i + N_rows];

  /* change matrix columns to reflect global numbers */

  for ( i = N_rows+1; i < bindx[N_rows] ; i++ ) {
     if ( bindx[i] < N_rows) bindx[i] += offset;
     else bindx[i] = externals[bindx[i] - N_rows];
  }

  BV_FREE((char *) externals);
  BV_FREE((char *) dbl_grows);

  /*  fetch processor information */

  nprocs = proc_config[AZ_N_procs];


  enlarged_Rownum = (int*)    BV_ALLOC((enlarged_N+1) * sizeof(int));
  if (enlarged_Rownum == NULL) {
    AZ_printf_err("Error in allocating memory for enlarged_Rownum.\n");
    exit(-1);
  }
  for (i=0; i<enlarged_N; i++) enlarged_Rownum[i] = Rownum[i];

  /* if no overlap is requested, just return the original matrix */

  if (olap_size <= 0) {
    (*New_N_rows) = enlarged_N;
    return;
  }

  /* Allocate temporary storage space for further processing 
   *   - send_proc : a processor flag array to indicate send processors
   *   - send_leng : data length (in doubles) for each send-to processor
   *   - send_rownum : row numbers of matrix to the send-to processors
   *   - send_rowcnt : number of matrix rows to the send-to processors
   *   - recv_proc : a processor flag array to indicate recv-from processors
   *   - recv_leng : data length (in doubles) for each recv-from processor
   *   - proc_neigh : a compressed list of neighbors for communication
   *   - actual_send_leng : compressed send_leng array
   *   - actual_recv_leng : compressed recv_leng array
   *   - start_send : pointing to subarrays for each send processor
   */

  send_rownum      = (int**) BV_ALLOC(    nprocs *sizeof(int*) +
                                         (9* nprocs)*sizeof(int));
  send_proc        = (int *) &(send_rownum[nprocs]);
  send_leng        = &(send_proc[nprocs]);
  send_rowcnt      = &(send_leng[nprocs]);
  recv_proc        = &(send_rowcnt[nprocs]);
  recv_leng        = &(recv_proc[nprocs]);
  proc_neigh       = &(recv_leng[nprocs]);
  actual_send_leng = &(proc_neigh[nprocs]);
  actual_recv_leng = &(actual_send_leng[nprocs]);
  start_send       = &(actual_recv_leng[nprocs]);

  for (i=0; i<nprocs; i++) send_rownum[i] = NULL; 


  /* duplicate the sorted Rownum list for use later in the loop */

  sorted_New_Rownum = (int*) BV_ALLOC((enlarged_N+1) * sizeof(int));
  if (sorted_New_Rownum == NULL) {
    AZ_printf_err("Error in allocating memory for sorted_New_Rownum.\n");
    exit(-1);
  }
  for (i=0; i<enlarged_N; i++) sorted_New_Rownum[i] = Rownum[i];

  /* enlarging the local matrix with a distance of one at a time using
     breadth first search */

  for (ndist=0; ndist<olap_size; ndist++) {

    /* Starting with the enlarged system, compose the external node list. 
     *  Input  : enlarged_N,sorted_New_Rownum,bindx 
     *  Output : N_ext_node, ext_nodelist
     */

    PAZ_compose_external(enlarged_N,sorted_New_Rownum, bindx, 
                         &N_ext_node, &ext_nodelist);

    BV_FREE(sorted_New_Rownum);

    /* use the external node list to construct the send information 
     *  Input  : N_ext_node, ext_nodelist, Rownum (must be sorted)
     *  Output : send_proc, send_rowcnt, send_rownum
     */

     AZ_setup_sendlist(N_ext_node, ext_nodelist, send_proc, send_rowcnt,
                       send_rownum, proc_config, max_per_proc, N_rows,Rownum);

    /* Now generate receive information data  */

    for (i=0; i<nprocs; i++)     recv_proc[i] = 0;
    for (i=0; i<N_ext_node; i++) recv_proc[ext_nodelist[i]/max_per_proc]++;
    BV_FREE(ext_nodelist);

    /* form a compressed processor list for communication
     */

    num_neigh = 0;
    for (i=0; i<nprocs; i++)
      if (send_proc[i] > 0 || recv_proc[i] > 0) proc_neigh[num_neigh++] = i;

    /* compute the size that the linear system will  */
    /* be when the new rows are incorporated         */

    net_gain_N = 0;
    for (i = 0 ; i < num_neigh; i++ )
       net_gain_N  += recv_proc[proc_neigh[i]];

    /* Pass new matrix size to the memory allocator */

    nzeros = bindx[enlarged_N] + net_gain_N;
    if (bindx_length < nzeros) {
       AZ_printf_err("Not enough memory allocated for overlapping\n");
       exit(1);
    }
    RESET_USAGE(nzeros);


    /* Make room for the new matrix entries by shifting    */
    /* the off-diagonals by 'net_gain_N' storage locations */

    for (i = bindx[enlarged_N]-1 ; i >= bindx[0] ; i--) {
        k        = i + net_gain_N;
        bindx[k] = bindx[i];
        val[k]   = val[i];
    }
    for (i=0; i<=enlarged_N; i++) bindx[i] += net_gain_N;


    /* Figure out the total number of off-diagonals */
    /* that we need to send out to update the rows. */

    k = 0;
    for (i = 0 ; i < num_neigh; i++) {
       jj = proc_neigh[i];
       for (j = 0 ; j < send_rowcnt[jj] ; j++ ) {
          k += (bindx[send_rownum[jj][j]+1] - 
                bindx[send_rownum[jj][j]]);
       }
    }
    if (k+1 > comm_size) {
       if (comm_buffer != NULL) { BV_FREE(comm_buffer); comm_buffer = NULL; }
       comm_size = k+1;
       comm_buffer = (double *) BV_ALLOC( comm_size * sizeof(double));
       if (comm_buffer == NULL) {
	 AZ_printf_err("AZ_setup_dd_olap_msr: Not enough memory for comm_buff\n");
	 exit(1);
       }
       iptr = (int *) comm_buffer;
       dptr = (double *) comm_buffer;
    }







    /* expand the Rownum list: */
    /*   1) allocate space     */
    /*   2) post receives      */
    /*   3) send messages      */
    /*   4) wait to receive    */

    tRownum  = (int*)    BV_ALLOC((net_gain_N + enlarged_N+1) * sizeof(int));
    if (tRownum == NULL) {
      AZ_printf_err("Error in allocating memory for tRownum (%d).\n",enlarged_N);
      exit(-1);
    }
    for (i=0; i<enlarged_N; i++) tRownum[i] = enlarged_Rownum[i];
    BV_FREE(enlarged_Rownum); 
    enlarged_Rownum = tRownum;


    type    =proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + 
                                      AZ_MSG_TYPE;
    off_set = enlarged_N;
    for (i = 0 ; i < num_neigh; i++ ) {
       jj = proc_neigh[i];
       (void) mdwrap_iread((void *) &(enlarged_Rownum[off_set]),
                         recv_proc[jj]*sizeof(int), &jj, &type, &(request[i]));
       off_set += recv_proc[jj];
    }
    for (i = 0 ; i < num_neigh; i++) {
       jj = proc_neigh[i];
       for (k = 0 ; k < send_rowcnt[jj] ; k++ )
             iptr[k] = enlarged_Rownum[send_rownum[jj][k] ] ;
       mdwrap_write((void*)iptr,send_rowcnt[jj]*sizeof(int),jj,type,&status);
    }

    off_set = enlarged_N;
    for (i = 0 ; i < num_neigh; i++ ) {
       jj = proc_neigh[i];
       mdwrap_wait((char*)&(enlarged_Rownum[off_set]),
                     recv_proc[jj]*sizeof(int),&jj, &type, &status, request+i);
       off_set += recv_proc[jj];
    }


    /* Obtain the number of columns in each new row */

    type            =proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + 
                                    AZ_MSG_TYPE;
    off_set = enlarged_N+1;
    for (i = 0 ; i < num_neigh; i++ ) {
       jj = proc_neigh[i];
       (void) mdwrap_iread((void *) &(bindx[off_set]),
                         recv_proc[jj]*sizeof(int), &jj, &type, &(request[i]));
       off_set += recv_proc[jj];
    }

    for (i = 0 ; i < num_neigh; i++) {
       jj = proc_neigh[i];
       for (k = 0 ; k < send_rowcnt[jj] ; k++ )
          iptr[k] = bindx[send_rownum[jj][k] + 1] - bindx[send_rownum[jj][k]];
       mdwrap_write((void*)iptr,send_rowcnt[jj]*sizeof(int),jj,type,&status);
    }

    off_set = enlarged_N + 1;
    for (i = 0 ; i < num_neigh; i++ ) {
       jj = proc_neigh[i];
       mdwrap_wait((char*) &(bindx[off_set]), recv_proc[jj]*sizeof(int),&jj, 
                    &type, &status, request+i);
       off_set += recv_proc[jj];
    }

    /* Now convert the 'row size' into the proper bindx entry */

    for (i = enlarged_N+1 ; i <= enlarged_N+net_gain_N ; i++ )
       bindx[i] += bindx[i-1];


    /* Pass new matrix size (including off-diags) to the memory allocator */

    nzeros = bindx[enlarged_N + net_gain_N];
    if (bindx_length < nzeros) {
       AZ_printf_err("Not enough memory allocated for overlapping\n");
       exit(1);
    }
    RESET_USAGE(nzeros);


    /* Obtain the diagonal entry of the new rows    */

    type =proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + 
                              AZ_MSG_TYPE;

    off_set = enlarged_N;
    for (i = 0 ; i < num_neigh; i++ ) {
       jj = proc_neigh[i];
       (void) mdwrap_iread((void*)&(val[off_set]),recv_proc[jj]*sizeof(double),
                            &jj, &type, &(request[i]));
       off_set += recv_proc[jj];
    }
   
    for (i = 0 ; i < num_neigh; i++) {
       jj = proc_neigh[i];
       for (k = 0 ; k < send_rowcnt[jj] ; k++ )
          dptr[k] = val[send_rownum[jj][k]];
       mdwrap_write((void*)dptr,send_rowcnt[jj]*sizeof(double),jj,
                     type,&status);
    }

    off_set = enlarged_N;
    for (i = 0 ; i < num_neigh; i++ ) {
       jj = proc_neigh[i];
       mdwrap_wait((char*) &(val[off_set]), recv_proc[jj]*sizeof(double),
                    &jj, &type, &status, request+i);
       off_set += recv_proc[jj];
    }

    /* Obtain the off-diagonal column indices */

    type            =proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +
                                   AZ_MSG_TYPE;
    off_set = enlarged_N;
    for (i = 0 ; i < num_neigh; i++ ) {
       jj = proc_neigh[i];
       actual_recv_leng[i] = bindx[off_set + recv_proc[jj]] - bindx[off_set];
       off_set += recv_proc[jj];
    }
   
    off_set = 0;
    for (i = 0 ; i < num_neigh; i++) {
       jj = proc_neigh[i];
       actual_send_leng[i] = 0;
       start_send[i]       = off_set;
       for (k = 0 ; k < send_rowcnt[jj] ; k++ ) {
          actual_send_leng[i] += (bindx[send_rownum[jj][k]+1] - 
                                     bindx[send_rownum[jj][k]]);
          for (j= bindx[send_rownum[jj][k]];j < bindx[send_rownum[jj][k]+1];j++)
             iptr[off_set++] = bindx[j];
       }
    }
    AZ_splitup_big_msg(num_neigh,(char *)iptr,(char *)&bindx[bindx[enlarged_N]],
		       sizeof(int), start_send, actual_send_leng, 
		       actual_recv_leng, proc_neigh,type,&ncnt, proc_config);

    /* Obtain the off-diagonal values */

    type            =proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +AZ_MSG_TYPE;
   
    off_set = 0;
    for (i = 0 ; i < num_neigh; i++) {
       jj = proc_neigh[i];
       for (k = 0 ; k < send_rowcnt[jj] ; k++ ) {
          for (j =bindx[send_rownum[jj][k]];j < bindx[send_rownum[jj][k]+1];j++)
             dptr[off_set++] = val[j];
       }
    }
    AZ_splitup_big_msg(num_neigh,(char *) dptr,(char *) &val[bindx[enlarged_N]],
		       sizeof(double), start_send, actual_send_leng, 
		       actual_recv_leng, proc_neigh,type,&ncnt, proc_config);



    enlarged_N += net_gain_N;

    /* Sort the new row number index array for further processing
     *   - allocate index array for new sorted row indices (enlarged matrix)
     *   - call AZ_sort to order the update index array 
     *   (enlarged_Rownum) ==> sorted_New_Rownum
     */ 

    sorted_New_Rownum = (int*) BV_ALLOC((enlarged_N+1) * sizeof(int));
    if (sorted_New_Rownum == NULL) {
      AZ_printf_err("Error in allocating memory for sorted_New_Rownum.\n");
      exit(-1);
    }
    for (i=0; i<enlarged_N; i++) sorted_New_Rownum[i] = enlarged_Rownum[i];
    AZ_sort(sorted_New_Rownum, enlarged_N, NULL, NULL);

    /* Prune the links that are not in the current enlarged node list
     *  - do this only at the last iteration
     *  - use PAZ_sorted_search to identify the outcasts
     */

    if (ndist == olap_size-1) {
      k   = bindx[0];
      ind = bindx[0];
      for (i=0; i<enlarged_N; i++) {
        m = bindx[i+1];
        for (j=k; j<m; j++) {
          ret_index = PAZ_sorted_search(bindx[j],enlarged_N,
                                        sorted_New_Rownum);
          if (ret_index >= 0) {
            bindx[ind] = bindx[j];
            val[ind++] = val[j];
          }   
        }
        k = m;
        bindx[i+1] = ind;
      }
    }

    /* clean up */

    for (i=0; i<nprocs; i++) { 
      if (send_rownum[i] != NULL) {
         BV_FREE(send_rownum[i]);  send_rownum[i] = NULL;
      }
    }
  }
  BV_FREE(comm_buffer);
  BV_FREE(sorted_New_Rownum);
  BV_FREE(send_rownum);

  nzeros = bindx[enlarged_N];

  /* extract a list of external nodes -> ext_nodelist, and
   * re-order the input matrix according to the local ordering 
   *  - at exit, bindx has been reordered with the local
   *    indices 0->N_rows-1 and the external indices 
   *    N_rows->N_ext_node+N_rows-1.  The global indices for the
   *    external list are recorded in ext_nodelist  
   */

  N_ext_node = enlarged_N - N_rows;
  sorted_exts = (int *) BV_ALLOC((N_ext_node+1)*sizeof(int));
  for (i = 0 ; i < N_ext_node ; i++ ) 
     sorted_exts[i] = enlarged_Rownum[i + N_rows];

  sprintf(str,"over_map %s",context->tag); 
  *map3 = (int *) AZ_manage_memory((N_ext_node+1)*sizeof(int),AZ_ALLOC,name,
				   str,&i);

  for (i = 0 ; i < N_ext_node; i++) (*map3)[i] = i + N_rows;
  AZ_sort(sorted_exts, N_ext_node, *map3, NULL);

  PAZ_find_local_indices(N_rows,bindx,Rownum,sorted_exts,
                        N_ext_node,*map3);

  /* use all current information to construct the data_org structure */

  PAZ_set_message_info(N_ext_node,N_rows,sorted_exts, Rownum,proc_config,
                      NULL, data_org, AZ_MSR_MATRIX,name, max_per_proc,context);

  /* put the enlarged matrix in the output parameters */

  (*New_N_rows) = enlarged_N;

  /* clean up */

  BV_FREE(sorted_exts);
  BV_FREE(enlarged_Rownum);
  BV_FREE(Rownum);

  (*data_org)[AZ_matrix_type] = AZ_MSR_MATRIX;
  (*data_org)[AZ_N_internal]  = N_rows;
  (*data_org)[AZ_N_border]    = 0;
  (*data_org)[AZ_N_external]  = N_ext_node;
  (*data_org)[AZ_N_int_blk]   = N_rows;
  (*data_org)[AZ_N_bord_blk]  = 0;
  (*data_org)[AZ_N_ext_blk]   = N_ext_node;
  
} /* AZ_dd_overlap */
 
/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Given a local matrix with N_local number of rows (N_local, bindx) and
 * a list of the local nodes, construct the external node list 
 *
 * Input :
 *
 *   N_local       : the number of local nodes
 *   sorted_Rownum : a list of row indices (sorted) for the local rows
 *   bindx         : the connectivity of the matrix in MSR
 *
 * Output :
 *
 *   N_ext_node : the number of external nodes found
 *   ext_nodelist : a list of node numbers for the external nodes
 * ----------------------------------------------------------------------
 */

void PAZ_compose_external(int N_local, int *sorted_Rownum, int *bindx, 
                          int *N_ext_node, int **ext_nodelist)
{  
  int  i, ret_index, enode_leng, *enode_list;
  int k,kk;

  enode_leng = 0;
  for (i=N_local+1; i<bindx[N_local]; i++) {
    ret_index = PAZ_sorted_search(bindx[i],N_local,sorted_Rownum);
    if (ret_index < 0) enode_leng++;
  }

  /* allocate initially certain amount of memory for holding the list */

  enode_list = (int*) BV_ALLOC((1+enode_leng) * sizeof(int));
  if (enode_list == NULL) 
    AZ_perror("Error in allocating memory for enode_list.\n");

  /* fetch the column indices in bindx to find out if they are in the
   * local list (sorted_Rownum).  If not, record them as external nodes.
   * In order not to duplicate the external nodes on the external node
   * list, we call PAZ_insert_list to do the job.  We also keep track
   * that the list does not overflow, and if an overflow is detected,
   * re-allocate more memory to hold the additional nodes.
   */    

  enode_leng = 0;
  for (i=N_local+1; i<bindx[N_local]; i++) {
    ret_index = PAZ_sorted_search(bindx[i],N_local,sorted_Rownum);
    if (ret_index < 0) 
       enode_list[enode_leng++] = bindx[i];
  }
  AZ_sort(enode_list, enode_leng, NULL, NULL);

  /* remove duplicates */
 
  kk = 0;
  for (k = 1; k < enode_leng ; k++) {
    if (enode_list[kk] != enode_list[k]) {
      kk++;
      enode_list[kk] = enode_list[k];
    }
  }
 
  if (enode_leng != 0) kk++;
  enode_leng = kk;

  

  /* clean up */

  (*N_ext_node)   = enode_leng;
  (*ext_nodelist) = enode_list; 
}
/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Given a sorted list of indices and the key, find the position of the
 * key in the list.  If not found, return where it would have been stored.
 *
 * author : Charles Tong (8117)
 *------------------------------------------------------------------------
 */

int PAZ_sorted_search(int key, int nlist, int *list)
{
  int  nfirst, nlast, nmid, found, index = 0;

  if (nlist <= 0) return -1;
  nfirst = 0;  
  nlast  = nlist-1; 
  if (key > list[nlast])  return -(nlast+1);
  if (key < list[nfirst]) return -(nfirst+1);
  found = 0;
  while ((found == 0) && ((nlast-nfirst)>1)) {
    nmid = (nfirst + nlast) / 2; 
    if (key == list[nmid])     {index  = nmid; found = 1;}
    else if (key > list[nmid])  nfirst = nmid;
    else                        nlast  = nmid;
  }
  if (found == 1)               return index;
  else if (key == list[nfirst]) return nfirst;
  else if (key == list[nlast])  return nlast;
  else                          return -(nfirst+1);
}

/******************************************************************************/
/******************************************************************************/

void PAZ_order_ele(int extern_index[], int N_update, int externs[], 
                   int N_external, int *m1, int *m3, int max_per_proc)

/*******************************************************************************

  Find an ordering for the elements in 'update' and 'external' (if option ==
  AZ_EXTERNS, we only order 'external'). We first order 'update' by numbering
  the internal elements first followed by the border elements. On output
  update_index[] contains the numbers 0->N_update-1 ordered such that
  update_index[k] < update_index[j] if update[k] is internal and update[j] is
  border. Next, we order the external elements such that elements that are
  updated by the same processor are contiguous.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  update_index:    On output, update_index[i] gives the local numbering of
                   global point 'update[i]'. See User's Guide.

  extern_index:    On output, extern_index[i] gives the local numbering of
                   global point 'external[i]'. See User's Guide.

  internal:        On output, equals the number of internal elements on this
                   processor.

  border:          On output, equals the number of border elements on this
                   processor.

  N_update:        Number of elements updated on this processor.

  bpntr, bindx:    Nonzero matrix information (see User's Guide).
                   Note: for MSR arrays the user invokes this routine by calling
                   AZ_order_ele(...,bindx,bindx, ....). That is both 'bpntr' and
                   'bindx' point to 'bindx'.

  externs:     

  N_external:      Number of external elements on this processor.

  option:          option = AZ_ALL ==> order internal, border and extern ele.
                   option = AZ_EXTERNS ==> order only external elements.
  m_type:          On input, m_type = AZ_MSR_MATRIX
                     or      m_type = AZ_VBR_MATRIX

*******************************************************************************/

{

  /* local variables */

  /**************************** execution begins ******************************/

  int  i, j, count;

  /*
   * Sift through the external elements. For each newly encountered external
   * point assign it the next index in the sequence. Then look for other
   * external elements who are update by the same node and assign them the next
   * set of index numbers in the sequence (ie. elements updated by the same node
   * have consecutive indices).
   */

  count = N_update;
  for (i = 0; i < N_external; i++) extern_index[i] = -1;

  for (i = 0; i < N_external; i++) {
    if (extern_index[i] == -1) {
      extern_index[i] = count++;
      m3[ extern_index[i] - N_update ] = m1[i];

      for (j = i + 1; j < N_external; j++) {
        if ((externs[j]/max_per_proc) == (externs[i]/max_per_proc)) {
           extern_index[j] = count++;
           m3[ extern_index[j] - N_update ] = m1[j];
        }
      }
    }
  }

} /* PAZ_order_ele */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void PAZ_find_local_indices(int N_update, int bindx[], int update[],
                           int *sorted_ext, int N_external, int map[])

/*******************************************************************************

  Given the global column indices for the matrix and a list of elements updated
  on this processor, compute the external elements and store them in the list
  'external' and change the global column indices to local column indices. In
  particular,

  On input, the column index bindx[k] is converted to j on output where

          update[j] = bindx[k]
   or
          external[j - N_update] = bindx[k]

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  bindx:           MSR or VBR column indices. On input, they refer to global
                   column indices. On output, they refer to local column indices
                   as described above. See User's Guide for more information.

  update:          List (global indices) of elements updated on this node.

  external:        List (global indices) of external elements on this node.

  N_external:      Number of external elements on this processor.

  mat_type:        Indicates whether this is an MSR (= AZ_MSR_MATRIX) or a
                   VBR (= AZ_VBR_MATRIX).

  bnptr:           'bpntr[N_update]' indicates the location of the
                   last VBR nonzero.

*******************************************************************************/

{

  /* local variables */

  int  j, k;
  int *bins,shift;
  int  start,end;

  /**************************** execution begins ******************************/

  /* set up some bins so that we will be able to use AZ_quick_find() */

  bins = (int *) BV_ALLOC((N_update / 4 + 10)*sizeof(int));
  if  (bins == NULL) {
    (void) AZ_printf_err( "ERROR: Not enough temp space\n");
    exit(-1);
  }

  AZ_init_quick_find(update, N_update, &shift, bins);

  /*
   * Compute the location of the first and last column index that is stored in
   * the bindx[].
   */

  start = bindx[0]; end = bindx[bindx[0]-1]; 
  

  /*
   * Estimate the amount of space we will need by counting the number of
   * references to columns not found among 'update[]'.  At the same time replace
   * column indices found in update[] by the appropriate index into update[].
   * Add N_update to columns not found in 'update[]' (this effectively marks
   * them as external elements).
   *
   * Note the space estimate will be an over-estimate as we do not take into
   * account that there will be duplicates in the external list.
   */

  for (j = start; j < end; j++) {
    k = AZ_quick_find(bindx[j], update, N_update,shift,bins);

    if (k != -1) bindx[j] = k;
    else {
       k = AZ_find_index(bindx[j], sorted_ext,N_external);
       if (k != -1) bindx[j] = map[k];
       else {
        (void) AZ_printf_err( "ERROR: column number not found %d\n",
                       bindx[j]);
        exit(-1);
      }
    }
  }

  BV_FREE((char *) bins);

} /* PAZ_find_local_indices */


/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* Given a local matrix with N_local number of rows (N_local, bindx) and
 * a list of the local nodes, construct the external node list 
 *
 * Input :
 *
 *   N_local       : the number of local nodes
 *   sorted_Rownum : a list of row indices (sorted) for the local rows
 *   bindx         : the connectivity of the matrix in MSR
 *
 * Output :
 *
 *   N_ext_node : the number of external nodes found
 *   ext_nodelist : a list of node numbers for the external nodes
 * ----------------------------------------------------------------------
 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void PAZ_set_message_info(int N_external, int N_update,
                         int external[], int update[], int proc_config[], 
			 int cnptr[], int *data_org[], int mat_type, int name,
                         int max_per_proc, struct context *context)

/*******************************************************************************

  Perform initializations so that local communication can occur to update the
  external elements. This initialization includes:

   1) determine the number of neighbors to which we send information as well as
      the number of neighbors from which we receive information.  These two
      should be the same. If they are not, we will set up things so that we send
      0 length messages to processors from which we receive (but have no
      information to send them).
   2) determine the total number of unknowns that we must send out.  Using this
      information we allocate space for data_org[].
   3) Initialize data_org[] (see User's Guide) so that it contains the number of
      messages to be sent/received, the node numbers of the processors to which
      we must send and receive, the length of the messages that we must send and
      that we expect to receive from each of the neighbors, and finally, a list
      of the indices of the elements that will be send to other processors (in
      the order that they will be sent).

      NOTE: Implicitly the neighbors are numbered using the ordering of the
      external elements. In particular, the external elements are ordered such
      that those elements updated by the same processor are contiguous. In this
      way, the ordering of the external elements defines an ordering for the
      neighbors.

 The information stored in 'data_org[]' is as follows:

     data_org[AZ_neighbors]         = node number of first neighbor
     data_org[AZ_neighbors+1]       = node number of second neighbor
     data_org[AZ_neighbors+2]       = node number of third  neighbor
         . . .
     data_org[AZ_rec_length]        = # of values to receive from 1st neighbor
     data_org[AZ_rec_length+1]      = # of values to receive from 2nd neighbor
     data_org[AZ_rec_length+2]      = # of values to receive from 3rd neighbor
     . . .
     data_org[AZ_send_length]       = # of values to send to 1st neighbor
     data_org[AZ_send_length+1]     = # of values to send to 2nd neighbor
     data_org[AZ_send_length+2]     = # of values to send to 3rd neighbor
     . . .
     data_org[AZ_send_list]         = elements to be sent to neighbors (starting
                                      with elements to the first neighbor [in
                                      the order that they appear on that
                                      neighbor], followed by elements to the
                                      second neighbor, etc).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_external:      Number of external elements on this processor.

  extern_index:    On output, extern_index[i] gives the local numbering of
                   global point 'external[i]'. See User's Guide.

  N_update:        Number of elements updated on this processor.

  external:        List (global indices) of external elements on this node.

  extern_proc:     extern_proc[i] is updating processor of external[i].

  update:          List (global indices) of elements updated on this node.

  update_index:    On output, update_index[i] gives the local numbering of
                   global point 'update[i]'. See User's Guide.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

  cnptr:           VBR column pointer array (see User's Guide).

  data_org:        Array use to specifiy communication information. See User's
                   Guide.
                   On output, this array contains neighbor and message.

  mat_type:        Indicates whether this is an MSR (= AZ_MSR_MATRIX) or a
                   VBR (= AZ_VBR_MATRIX).

*******************************************************************************/

{

  /* local variables */

  int   i, j, ii, type, start, cflag, partner, newlength;
  int   total_to_be_sent, found;
  int  *new_external, *new_extern_proc;
  int  *neighbors, *tempneigh;
  int  *recv_list, *send_list;
  int   tj;
  int   num_recv_neighbors;
  int   num_send_neighbors;
  int  *bins, shift;
  int  *send_ptr, *lens;
  int   proc,nprocs;
  int   firstone,current;
  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */
  char  str[80];
  unsigned int length;

  char *yo = "PAZ_set_message_info: ";

  /*---------------------- execution begins -----------------------------*/

  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  neighbors = (int *) BV_ALLOC(nprocs*sizeof(int));

  tempneigh = (int *) BV_ALLOC(nprocs*sizeof(int));
  for (i = 0 ; i < nprocs ; i++ ) neighbors[i] = 0;

  /* Produce a list of the external updating processors corresponding */
  /* to each external point in the order given by 'extern_index[]'    */

  new_extern_proc = (int *) BV_ALLOC((N_external+1)*sizeof(int));
  if (new_extern_proc == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  for (i = 0; i < N_external; i++)
    new_extern_proc[i] = external[i]/max_per_proc;

  /*
   * Count the number of neighbors from which we receive information to update
   * our external elements. Additionally, fill the array neighbors in the
   * following way:
   *      neighbors[i] = 0   ==>  No external elements are updated by
   *                              processor i.
   *      neighbors[i] = x   ==>  (x-1)/nprocs elements are updated from
   *                              processor i.
   */

  num_recv_neighbors = 0;
  length             = 1;

  for (i = 0; i < N_external; i++) {
    if (neighbors[new_extern_proc[i]] == 0) {
      num_recv_neighbors++;
      neighbors[new_extern_proc[i]] = 1;
    }

    if (mat_type != AZ_VBR_MATRIX)
      neighbors[new_extern_proc[i]] += nprocs;
    else
      neighbors[new_extern_proc[i]] += (nprocs * (cnptr[i + 1 + N_update]
                                                  - cnptr[i + N_update]));
  }

  /*
   * Make a list of the neighbors that will send information to update our
   * external elements (in the order that we will receive this information).
   */

  recv_list = (int *) BV_ALLOC(AZ_MAX_NEIGHBORS*sizeof(int));
  if (recv_list == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  j = 0;
  recv_list[j++] = new_extern_proc[0];
  for (i = 1; i < N_external; i++) {
    if (new_extern_proc[i - 1] != new_extern_proc[i]) {
      recv_list[j++] = new_extern_proc[i];
    }
  }

  /* sum over all processors all the neighbors arrays */

  AZ_gsum_vec_int(neighbors, tempneigh, proc_config[AZ_N_procs], proc_config);

  /* decode the combined 'neighbors' array from all the processors */

  num_send_neighbors = neighbors[proc] % nprocs;

  /* decode 'neighbors[proc] to deduce total number of elements we must send */

  total_to_be_sent = (neighbors[proc] - num_send_neighbors) / nprocs;

  BV_FREE((char *) neighbors);
  BV_FREE((char *) tempneigh);

  /* send a 0 length message to each of our recv neighbors */

  send_list = (int *) BV_ALLOC((num_send_neighbors+1)*sizeof(int));
  if (send_list == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }
  for (i = 0 ; i < num_send_neighbors; i++ ) send_list[i] = 0;


  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +AZ_MSG_TYPE;


  /* first post receives */

  length = 0;
  for (i = 0; i < num_send_neighbors; i++) {
    partner = -1;
    mdwrap_iread((void *) external, length, &partner, &type, request+i);
  }

  /* send messages */

  for (i = 0; i < num_recv_neighbors; i++) {
    mdwrap_write((void *) external, 0, recv_list[i], type, &cflag);
  }

  /*
   * Receive message from each send neighbor to construct 'send_list'.
   */

  length = 0;
  for (i = 0; i < num_send_neighbors; i++) {
    partner = -1;
    mdwrap_wait((void *) external, length, &partner, &type, &cflag, request+i);
    if (partner != -1) send_list[i] = partner;
  }



  /*
   * Compare the two lists. In most cases they should be the same.  However, if
   * they are not add new entries to the recv list that are in the send list
   * (but not already in the recv list).
   */

  for (j = 0; j < num_send_neighbors; j++) {
    found = 0;
    for (i = 0; i < num_recv_neighbors; i++) {
      if (recv_list[i] == send_list[j])
        found = 1;
    }

    if (found == 0) {
      recv_list[num_recv_neighbors] = send_list[j];
      (num_recv_neighbors)++;
    }
  }

  BV_FREE((char *) send_list);
  num_send_neighbors = num_recv_neighbors;

  /* create data_org array */

  sprintf(str,"padded d_org %s",context->tag); 
  *data_org = (int *) AZ_manage_memory((total_to_be_sent + AZ_send_list)*
                                   sizeof(int), AZ_ALLOC,name,str, &i);

  (*data_org)[AZ_total_send] = total_to_be_sent;
  send_ptr = &((*data_org)[AZ_send_list]);

  /*
   * Create 'new_external' which explicitly put the external elements in the
   * order given by 'extern_index'
   */

  new_external = (int *) BV_ALLOC((N_external+1)*sizeof(int));
  if (new_external == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  for (i = 0; i < N_external; i++) {
    new_external[i] = external[i];
  }

  /*
   * Send each processor the global index list of the external elements in the
   * order that I will want to receive them when updating my external elements
   */



  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  lens = (int *) BV_ALLOC((num_recv_neighbors+1)*sizeof(int));
  if (lens == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  for (i = 0; i < num_recv_neighbors; i++) {
    length     = sizeof(int);
    partner    = recv_list[i];
    mdwrap_iread((void *) &(lens[i]), length, &partner, &type, request+i);
  }


  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    start  = j;
    newlength = 0;

    /* go through list of external elements until updating processor changes */

    while ((j < N_external) && (new_extern_proc[j] == recv_list[i])) {

      if (mat_type != AZ_VBR_MATRIX) newlength++;
      else newlength += (cnptr[j + 1 + N_update] - cnptr[j + N_update]);

      j++;
      if (j == N_external) break;
    }

    (*data_org)[AZ_rec_length+i] = newlength;
    (*data_org)[AZ_neighbors+i]  = recv_list[i];

    length = j - start;
    mdwrap_write((void *) &length, sizeof(int), recv_list[i], type, &cflag);
  }

  /* receive from each neighbor the global index list of external ele */

  for (i = 0; i < num_recv_neighbors; i++) {
    length     = sizeof(int);
    partner    = (*data_org)[AZ_neighbors+i];
    mdwrap_wait((void *) &(lens[i]), length, &partner, &type, &cflag,request+i);
    (*data_org)[AZ_send_length + i] = lens[i];
  }
  BV_FREE(lens);








  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;


  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length     = (  (*data_org)[AZ_send_length + i] ) * sizeof(int);
    partner    = (*data_org)[AZ_neighbors+i];
    mdwrap_iread((void *) &(send_ptr[j]), length, &partner, &type, request+i);
    j += (*data_org)[AZ_send_length + i];
  }

  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    start  = j;
    newlength = 0;

    /* go through list of external elements until updating processor changes */

    while ((j < N_external) && (new_extern_proc[j] == recv_list[i])) {

      if (mat_type != AZ_VBR_MATRIX) newlength++;
      else newlength += (cnptr[j + 1 + N_update] - cnptr[j + N_update]);

      j++;
      if (j == N_external) break;
    }
    mdwrap_write((void *) &(new_external[start]), (j-start)* sizeof(int),
             recv_list[i], type, &cflag);
  }

  /* receive from each neighbor the global index list of external ele */

  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length     = (  (*data_org)[AZ_send_length + i] ) * sizeof(int);
    partner    = (*data_org)[AZ_neighbors+i];
    mdwrap_wait((void *) &(send_ptr[j]), length, &partner, &type, &cflag,
		 request+i);
    j += (*data_org)[AZ_send_length + i];
  }



  tj = 0;
  if (mat_type == AZ_VBR_MATRIX) {
    (*data_org)[AZ_matrix_type] = AZ_VBR_MATRIX;

    /* move the information to the back of the send_ptr */

    for (i = j - 1; i >= 0; i--) send_ptr[total_to_be_sent+i-j] = send_ptr[i];
    tj = total_to_be_sent - j;
  }
  else (*data_org)[AZ_matrix_type] = AZ_MSR_MATRIX;

  total_to_be_sent = j;

  /* replace global indices by local indices */

  ii = 0;
  firstone  = AZ_send_length;
  current   = firstone;

  while ((current-firstone < num_recv_neighbors) && ((*data_org)[current] == 0))
    current++;

  bins = (int *) BV_ALLOC((N_update / 4 + 10)*sizeof(int));
  if (bins == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  AZ_init_quick_find(update, N_update, &shift, bins);
  for (i = tj; i < total_to_be_sent + tj; i++) {

    /*
     * j = AZ_find_index(send_ptr[i], update, N_update);
     */

    j = AZ_quick_find(send_ptr[i], update, N_update, shift, bins);
    if (j == -1) {
      (void) AZ_printf_err( "%sERROR: A point (%d) requested from processor "
                     "%d\nwas not found among the update ele of proc %d\n", yo,
                     send_ptr[i], (*data_org)[AZ_neighbors+current-firstone],
                     proc);
      exit(-1);
    }

    send_ptr[ii] = j;
    ii +=  1;
  }

  BV_FREE((char *) bins);
  BV_FREE((char *) recv_list);
  BV_FREE((char *) new_external);
  BV_FREE((char *) new_extern_proc);

  (*data_org)[AZ_N_neigh] = num_send_neighbors;

} /* PAZ_set_message_info */


/**************************************************************/
/**************************************************************/
/**************************************************************/
char *AZ_allocate_or_free(void *ptr, unsigned int input_size, int action)
{
/*
 * Simple memory allocated that tries to allocate space off the
 * back end of bindx and val. The general scheme is as follows:
 *   1) we give this routine the length of bindx and val.
 *   2) we give this routine an estimated value of how much
 *      space we will use in bindx and val
 *   3) we first try to allocate space from bindx or val if there
 *      is additional memory beyond the estimated requirements.
 *   4) if 3) is unsuccessful, we try and malloc space.
 *   5) if 4) is unsuccessful, we try to allocate space
 *      from bindx and val that might violate the estimates.
*/

    long int           size;
    static long int    b_estimated_need, v_estimated_need;
    static long int    v_freelist, b_freelist;
    static long int    b_lastused, v_lastused;
    static double      *v,*b;
    static long int    v_end, b_end;
    static long int    v_smallest_free, b_smallest_free;

    double *dptr;
    char *t_ptr;
    long int msize, where; 
    long int current, prev;

    size = (long int) input_size;

    if (action == ALLOCATE) {
       msize = size;
       size += sizeof(double) - size%sizeof(double);
       size = size/sizeof(double);
       size = size + 1;

       /* check the v free list */

       current = v_freelist;
       prev    = -1;
       while (current != -1) {
          if ( (long int) v[current-1] >= size) break;
          prev = current;
          current = (long int) v[current];
       }
       if ( current != -1) {
          if (prev == -1) v_freelist = (long int) v[current];
          else v[prev] = v[current];
          return( (char *) &(v[current])   );
       }

       /* check the b free list */

       current = b_freelist;
       prev    = -1;
       while (current != -1) {
          if ( (long int) b[current-1] >= size) break;
          prev = current;
          current = (long int) b[current];
       }
       if ( current != -1) {
          if (prev == -1) b_freelist = (long int) b[current];
          else b[prev] = b[current];
          return( (char *) &(b[current])   );
       }

       /* check (conservative) if there is  */
       /* additional space in either v or b */
     
       if (v_smallest_free - v_estimated_need >= b_smallest_free -
                                                   b_estimated_need) {
          if (v_smallest_free - size > v_estimated_need) {
             v_smallest_free    -= size;
             v[v_smallest_free]  = (double) size;
             return( (char *) &(v[v_smallest_free+1]) );
          }
       }
       else {
          if (b_smallest_free - size > b_estimated_need ) {
             b_smallest_free    -= size;
             b[b_smallest_free]  = (double) size;
             return( (char *) &(b[b_smallest_free+1]) );
          }
       }

       /* check if there is malloc space */
    
       t_ptr = (char *) AZ_allocate((unsigned int) msize);
       if (t_ptr != NULL) return(t_ptr);


       /* check if there is any space  */
       /* in either v or b             */
     
       if (v_smallest_free - v_lastused >= b_smallest_free - b_lastused){
          if (v_smallest_free - size > v_lastused ) {
             v_smallest_free -= size;
             v[v_smallest_free] = size;
             return( (char *) &(v[v_smallest_free+1]) );
          }
       }
       else {
          if (b_smallest_free - size > b_lastused) {
             b_smallest_free -= size;
             b[b_smallest_free] = size;
             return( (char *) &(b[b_smallest_free+1]) );
          }
       }
       return(NULL);

    }
    else if (action == FREE_IT) {

       dptr = (double *) ptr;

       where = OUT_OF_BOUNDS;
       if ( (dptr >= v) && (dptr < &(v[v_end]) ) ) where = IN_V;

       if ( (dptr >= b) && (dptr < &(b[b_end]) ) ) where = IN_B;

       if (where == OUT_OF_BOUNDS) {
          AZ_free(ptr);
          return(NULL);
       }

       dptr--;
       size = (int) dptr[0];

       if (where == IN_V) {
          dptr[1] = (int) v_freelist;
          dptr++;
          v_freelist = ((long int) dptr - (long int) v)/sizeof(double);
       }
       else {
          dptr[1] = (long int) b_freelist;
          dptr++;
          b_freelist = ((long int) dptr - (long int) b)/sizeof(double);
       }
    }
    else if (action == V_SET) {
        v          = (double *) ptr;
        v_freelist = -1;
        v_end      = size;
        v_smallest_free = v_end;
    }
    else if (action == B_SET) {
        b          = (double *) ptr;
        b_freelist = -1;
        size *= sizeof(int);
        size -= size%sizeof(double);
        size /= sizeof(double);
        b_end      = size;
        b_smallest_free = b_end;
    }
    else if (action == LASTUSED_SET) {
        size *= sizeof(int);
        if (size%sizeof(double) != 0) 
           size += (sizeof(double) - size%sizeof(double));
        
        b_lastused = size/sizeof(double);
        v_lastused = size/sizeof(int);
     
        if (b_lastused > b_smallest_free) {
          AZ_printf_err("Error: Out of space due to poor estimate of memory needed\n");
          AZ_printf_err("       for overlapping.\n");
          exit(1);
        }
        if (v_lastused > v_smallest_free) {
          AZ_printf_err("Error: Out of space due to poor estimate of memory needed\n");
          AZ_printf_err("       for overlapping.\n");
        }
    }
    else if (action == ESTIMATED_SET) {
        size *= sizeof(int);
        if (size%sizeof(double) != 0) 
           size += (sizeof(double) - size%sizeof(double));

        b_estimated_need = size/sizeof(double);
        v_estimated_need = size/sizeof(int);
    }
    else if (action == -43) {
       AZ_printf_out("v_list: ");
       current = v_freelist;
       prev    = -1;
       while (current != -1) {
          AZ_printf_out("(%d, %d) ",(int) current,(int) v[current-1]);
          prev = current;
          current = (long int) v[current];
       }
       AZ_printf_out("\n");
       AZ_printf_out("b_list: ");
       current = b_freelist;
       prev    = -1;
       while (current != -1) {
          AZ_printf_out("(%d, %d) ",(int) current,(int) b[current-1]);
          prev = current;
          current = (long int) b[current];
       }
       AZ_printf_out("\n\n");

    }
    else {
        AZ_printf_err("Error: unknown option (%d) for allocate_or_free\n",action);
        exit(1);
    }
    return((char *) NULL);

}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */

/* This subroutine uses the external node list input (ext_nodelist) to
 * construct the send-to-processor information (send_proc, send_rowcnt,
 * and send_rownum.  A sorted local node list and the auxilliary list
 * (the one obtained in AZ_sort) are also required
 *
 * Input :
 *
 *    N_ext             : the length of the external node list
 *    externs           : the actual external node list. This list must
 *                        be sorted.
 *    N_rows            : the number of local rows
 *    Rownum            : the indices of the local rows (MUST BE SORTED)
 *    proc_config       : parallel machine information
 *
 * Output :
 *
 *    send_proc         : if send_proc[i] = 1, then there is something to
 *                        send to processor i. (Note : this array should be
 *                        allocated as nprocs integers before entering)
 *    send_rowcnt       : send_rowcnt[i] holds the number of rows to be
 *                        sent to processor i. (Again, this array should be
 *                        allocated as nprocs integers before entering)
 *    send_rownum       : send_rownum[i] holds the actual row numbers (local)
 *                        of the rows to be sent to processor i. (Again, this
 *                        array should be allocated as nprocs int pointers
 *                        before entering)
 */



#define proc(i)  (externs[i]/max_per_proc)
void AZ_setup_sendlist(int N_ext, int externs[], int send_proc[],
		       int send_rowcnt[], int *send_rownum[], 
                       int proc_config[], int max_per_proc, int N_rows,
                       int Rownum[])
{
   int node, nprocs;
   int N_send, N_recv;                     /* Number of messages processor */
                                           /* must send and recv           */
   MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle               */
   int type,st;
   int length;                             /* message send length          */
   int *temp1,*temp2;
   int i, j, start, processor, *ttt;


   nprocs    = proc_config[AZ_N_procs];
   node      = proc_config[AZ_node];
   temp1     = send_proc;
   temp2     = send_rowcnt;

   /* determine how many messages I need to send  */
   /* to update external rows of other processors */

   for (i = 0 ; i < nprocs; i++) temp1[i] = 0;
   for (i = 0 ; i < N_ext ; i++) temp1[proc(i)] = 1;
   AZ_gsum_vec_int(temp1,temp2, nprocs,proc_config);
   N_send = temp1[node];
   
   /* exchange information so that each processor knows: */
   /*    1) processor to which it must send rows         */
   /*    2) number of rows to be sent to each processor  */

   type            =proc_config[AZ_MPI_Tag];
   proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

   for (i = 0 ; i < N_send; i++ ) {
      send_proc[i] = -1;
      (void) mdwrap_iread((void *) &(send_rowcnt[i]), sizeof(int), 
                           &(send_proc[i]), &type, &(request[i]));
   }

   N_recv = 0;
   length = 1;
   for (i = 1 ; i < N_ext ; i++ ) {
      if (proc(i) != proc(i-1)) {
         (void) mdwrap_write((void *) &length, sizeof(int),proc(i-1),type,&st);
         length = 0;
         N_recv++;
      }
      length++;
   }
   if (N_ext != 0) {
      (void) mdwrap_write((void *)&length, sizeof(int),proc(N_ext-1),type,&st);
      N_recv++;
   }

   for (i = 0 ; i < N_send; i++ ) {
       (void) mdwrap_wait((void *) &(send_rowcnt[i]), sizeof(int), 
                           &(send_proc[i]), &type, &st, &(request[i]));
   }

   AZ_sort(send_proc,N_send,send_rowcnt,NULL);   /* sort the processor info */

   /* now receive the actual row numbers */

   type            =proc_config[AZ_MPI_Tag];
   proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;


   for (i = 0 ; i < N_send; i++ ) {
      send_rownum[i] = (int *) BV_ALLOC((send_rowcnt[i]+1) * sizeof(int));
      (void) mdwrap_iread((void *) send_rownum[i], send_rowcnt[i]*sizeof(int), 
                           &(send_proc[i]), &type, &(request[i]));
   }
   start = 0;
   length = 1;
   for (i = 1 ; i < N_ext ; i++ ) {
      if (proc(i) != proc(i-1)) {
         (void) mdwrap_write((void *) &(externs[start]), 
                               length*sizeof(int),proc(i-1),type,&st);
         start += length;
         length = 0;
      }
      length++;
   }
   if (N_ext != 0) {
      (void) mdwrap_write((void *) &(externs[start]), length*sizeof(int),
                           proc(N_ext-1),type,&st);
   }

   for (i = 0 ; i < N_send; i++ ) {
       (void) mdwrap_wait((void *) send_rownum[i], send_rowcnt[i]*sizeof(int), 
                           &(send_proc[i]), &type, &st, &(request[i]));
   }


   for (i = N_send ; i < nprocs ; i++ ) {
      send_proc[i] = 0;
      send_rowcnt[i]  = 0;
   }

   for (i = N_send-1 ; i >= 0 ; i-- ) {
      processor               = send_proc[i];
      j                       = send_rowcnt[i];
      ttt                     = send_rownum[i];
      if (i < processor) {send_proc[i] = 0;       /* The 'if' statement is  */
                          send_rowcnt[i] = 0;     /* needed.                */
                          send_rownum[i] = NULL;}
      send_rowcnt[processor]  = send_rowcnt[i];
      send_rowcnt[processor]  = j;
      send_rownum[processor]  = ttt;

      for (j = 0 ; j < send_rowcnt[processor] ; j++ ) 
         send_rownum[processor][j] =PAZ_sorted_search(send_rownum[processor][j],
                                                      N_rows,Rownum);

      send_proc[processor]    = 1;
   }
}
