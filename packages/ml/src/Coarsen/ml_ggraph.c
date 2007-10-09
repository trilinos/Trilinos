/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions to communicate between processors                          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : Nov, 1999                                            */
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include "ml_ggraph.h"

/* ******************************************************************** */
/* Constructor                                                          */
/* -------------------------------------------------------------------- */

int ML_GGraph_Create( ML_GGraph **gg )
{
   ML_memory_alloc( (void **) gg, sizeof(ML_GGraph), "gg0" );
  (*gg)->ML_id     = ML_ID_GGRAPH;
  (*gg)->ML_rank   = -1;
  (*gg)->Npoints   = 0;
  (*gg)->send_cnt  = 0;
  (*gg)->recv_cnt  = 0;
  (*gg)->row_ptr   = NULL;
  (*gg)->col_ptr   = NULL;
  (*gg)->recv_list = NULL;
  (*gg)->send_list = NULL;
  (*gg)->recv_leng = NULL;
  (*gg)->send_leng = NULL;
  (*gg)->recv_proc = NULL;
  (*gg)->send_proc = NULL;
  (*gg)->bdry_type = NULL;
  (*gg)->vertex_state = NULL;
  return 0;
}

/* ******************************************************************** */
/* destructor                                                           */
/* -------------------------------------------------------------------- */

int ML_GGraph_Destroy( ML_GGraph **gg )
{
   int i;

   if ( (*gg)->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_Destroy : wrong object. \n");
      exit(1);
   }
   if ( (*gg)->row_ptr   != NULL ) ML_free( (*gg)->row_ptr );
   if ( (*gg)->col_ptr   != NULL ) ML_free( (*gg)->col_ptr );
   if ( (*gg)->send_proc != NULL ) ML_free( (*gg)->send_proc );
   if ( (*gg)->recv_proc != NULL ) ML_free( (*gg)->recv_proc );
   if ( (*gg)->send_leng != NULL ) ML_free( (*gg)->send_leng );
   if ( (*gg)->recv_leng != NULL ) ML_free( (*gg)->recv_leng );
   if ( (*gg)->send_list != NULL )
   {
      for ( i = 0; i < (*gg)->send_cnt; i++ )
         if ( (*gg)->send_list[i] != NULL ) ML_free((*gg)->send_list[i]);
      ML_free((*gg)->send_list);
      (*gg)->send_list = NULL;
   }
   if ( (*gg)->recv_list != NULL )
   {
      for ( i = 0; i < (*gg)->recv_cnt; i++ )
         if ( (*gg)->recv_list[i] != NULL ) ML_free((*gg)->recv_list[i]);
      ML_free((*gg)->recv_list);
      (*gg)->recv_list = NULL;
   }
   if ( (*gg)->bdry_type    != NULL ) ML_free( (*gg)->bdry_type );
   if ( (*gg)->vertex_state != NULL ) ML_free( (*gg)->vertex_state );
   ML_memory_free( (void **) gg );
   return 0;
}

/* ******************************************************************** */
/* print the grid graph data structure                                  */
/* -------------------------------------------------------------------- */

int ML_GGraph_Print( ML_GGraph *ml_gg )
{
   int i;

   if ( ml_gg->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_Print : wrong object. \n");
      exit(1);
   }

   printf(" ************* ML_GGraph Data Structure ************* \n");
   printf(" Number of vertices = %d \n", ml_gg->Npoints);
   if ( ml_gg->bdry_type != NULL )
   {
      for ( i = 0; i < ml_gg->Npoints; i++ )
         printf("    Boundary type %d = %c \n", i, ml_gg->bdry_type[i] );
   }
   printf(" Number of edges    = %d \n", ml_gg->row_ptr[ml_gg->Npoints]);
  
   printf(" Number of points selected = %d \n", ml_gg->Nselected);
   for ( i = 0; i < ml_gg->Npoints; i++ )
      printf(" vertex state %d = %c \n", i, ml_gg->vertex_state[i]);
   
   /*
   printf(" Number of send destinations = %d \n", ml_gg->send_cnt);
   for ( i = 0; i < ml_gg->send_cnt; i++ )
   {
      printf("    destination %2d (length = %d) \n", ml_gg->send_proc[i],
                                                   ml_gg->send_leng[i]);
      for ( j = 0; j < ml_gg->send_leng[i]; j++ )
         printf("      index %2d = %d \n", j, ml_gg->send_list[i][j]);
   }
   printf(" Number of recv sources = %d \n", ml_gg->recv_cnt);
   for ( i = 0; i < ml_gg->recv_cnt; i++ )
   {
      printf("    source %2d (length = %d) \n", ml_gg->recv_proc[i],
                                            ml_gg->recv_leng[i]);
      for ( j = 0; j < ml_gg->recv_leng[i]; j++ )
         printf("      index %2d = %d \n", j, ml_gg->recv_list[i][j]);
   }
   */
   return 0;
}

/* ******************************************************************** */
/* print the grid graph data structure                                  */
/* -------------------------------------------------------------------- */

int ML_GGraph_Load_BdryTypes( ML_GGraph *ml_gg, int n, char *btype )
{
   int  i;

   /* ------------------------------------------------------------- */
   /* initial error checking (for proper ML_GGraph structure)       */
   /* ------------------------------------------------------------- */
   if ( ml_gg->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_Load_BdryTypes : wrong object. \n");
      exit(1);
   }
   if ( (ml_gg->Npoints != 0 && ml_gg->Npoints != n) || n <= 0 )
   {
      printf("ML_GGraph_LoadBdryType : wrong length. \n");
      exit(1);
   }

   /* ------------------------------------------------------------- */
   /* allocate storage and load the boundary types                  */
   /* ------------------------------------------------------------- */

   ml_gg->Npoints = n;
   ml_gg->bdry_type = (char *) ML_allocate( n * sizeof( char ) );
   for ( i = 0; i < n; i++ )
   {
      if ( btype[i] != ML_BDRY_INSIDE && btype[i] != ML_BDRY_RIDGE &&
           btype[i] != ML_BDRY_FACE   && btype[i] != ML_BDRY_CORNER )
      {
         printf("ML_GGraph_LoadBdryType : wrong boundary type. \n");
         exit( 0 );
      }
      ml_gg->bdry_type[i] = btype[i]; 
   }
   return 0;
}
   
/* ******************************************************************** */
/* generate node graph => row_ptr, col_ptr                              */
/* -------------------------------------------------------------------- */

int ML_GGraph_Gen_NodeGraph(ML_GGraph *ml_gg,void *grid,void (*gf),
                            ML_Comm *comm) 
{
   int         i, j, k, m, count, ncount, index, mypid, nprocs, status;
   int         total_count, total_nodes, nelements, nvertices;
   int         vlength, *vlist, *proc_array, *inttmp, *adjacency_cnt;
   int         *templist, *remote_list, *global_list, *node_ia, *node_ja;
   int         send_cnt, *send_leng, *send_proc, **send_list, msgtype;
   int         recv_cnt, *recv_leng, *recv_proc, tot_recv_leng;
   int         fromproc, *intarray, vlengmax, **recv_list, **send_list2;
   int         max_remote_index;
   ML_GridFunc *grid_fcns;
   USR_REQ     *Request;
   unsigned int nbytes;

   /* ------------------------------------------------------------- */
   /* initial error checking (for proper ML_GGraph structure)       */
   /* ------------------------------------------------------------- */

   if ( ml_gg->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_Gen_NodeGraph : wrong object. \n");
      exit(1);
   }

   /* ------------------------------------------------------------- */
   /* probe the incoming grid structure to figure out about the     */
   /* storage requirement to store the grid graph                   */
   /* ------------------------------------------------------------- */

   grid_fcns = (ML_GridFunc *) gf;
   mypid     = comm->ML_mypid;
   nprocs    = comm->ML_nprocs;
   nvertices = grid_fcns->USR_grid_get_nvertices( grid );
   nelements = grid_fcns->USR_grid_get_nelements( grid );
   if ( nvertices > 0 )
   {
      adjacency_cnt  = (int *) ML_allocate( nvertices * sizeof(int) );
   }
   else
   {
      printf("%d : ML_GGraph_NodeGraph : nvertices <= 0\n", mypid);
      adjacency_cnt = NULL;
   }
   for ( i = 0; i < nvertices; i++ ) adjacency_cnt[i] = 0;
   vlengmax = 20;
   vlist  = (int *) ML_allocate( vlengmax * sizeof(int) );
   for ( i = 0; i < nelements; i++ )
   {
      vlength = grid_fcns->USR_grid_get_element_nvertices( grid, i );
      if ( vlength > vlengmax ) 
      {
         vlengmax = vlength;
         ML_free( vlist );
         vlist = (int *) ML_allocate( vlengmax * sizeof(int) );
      }
      grid_fcns->USR_grid_get_element_vlist( grid, i, vlist );
      for ( j = 0; j < vlength; j++ )
      {
         index = vlist[j];
         if ( index < nvertices )
         {
            for ( k = 0; k < vlength; k++ )
               if ( vlist[k] != index ) adjacency_cnt[vlist[j]]++;
         }
      }
   }
   ML_free( vlist );
   total_nodes = 0;
   for ( i = 0; i < nvertices; i++ ) total_nodes += adjacency_cnt[i];

   /* ------------------------------------------------------------- */
   /* allocate the node_ia and node_ja array for the graph          */
   /* ------------------------------------------------------------- */

   node_ia = (int *) ML_allocate( (nvertices + 1) * sizeof(int) );
   if ( total_nodes > 0 )
   {
      node_ja = (int *) ML_allocate( total_nodes * sizeof(int) );
   }
   else
   {
      node_ja = NULL;
   }
   node_ia[0] = 0;
   for ( i = 1; i <= nvertices; i++ ) 
      node_ia[i] = node_ia[i-1] + adjacency_cnt[i-1];

   /* ------------------------------------------------------------- */
   /* put the node connectivity into the node_ia and node_ja array  */
   /* ------------------------------------------------------------- */

   if ( nvertices > 0 ) adjacency_cnt[0] = 0;
   for ( i = 1; i < nvertices; i++ ) adjacency_cnt[i] = node_ia[i]; 
   vlist   = (int *) ML_allocate( vlengmax * sizeof(int) );
   for ( i = 0; i < nelements; i++ )
   {
      vlength = grid_fcns->USR_grid_get_element_nvertices( grid, i );
      grid_fcns->USR_grid_get_element_vlist( grid, i, vlist );
      for ( j = 0; j < vlength; j++ )
      {
         index = vlist[j];
         if ( index < nvertices )
         {
            for ( k = 0; k < vlength; k++ )
               if ( vlist[k] != index ) 
                  node_ja[adjacency_cnt[index]++] = vlist[k];
         }
      }
   }
   ML_free( vlist );
   
   /* ------------------------------------------------------------- */
   /* sort the array and take out repeated indices                  */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < nvertices; i++ ) {
      count = adjacency_cnt[i] - node_ia[i];
      ML_sort( count, &(node_ja[node_ia[i]]) );
      for ( j = 0; j < count-1; j++ ) 
      {
         index = node_ja[node_ia[i]+j];
         for ( k = j+1; k < count; k++ ) 
         {
            if ( node_ja[node_ia[i]+k] != index ) break;
         } 
         ncount = k - j - 1;
         if ( ncount > 0 )
         {
            index = j + 1;
            for ( m = k; m < count; m++ ) 
            {
               node_ja[node_ia[i]+index] = node_ja[node_ia[i]+m]; 
               index++;
            }
            count = count - ncount;
         }
      }
      adjacency_cnt[i] = count;
   }

   /* ------------------------------------------------------------- */
   /* finally compress the array into CSR format                    */ 
   /* ------------------------------------------------------------- */

   index = 0;
   k = node_ia[0];
   for ( i = 0; i < nvertices; i++ ) {
      count = adjacency_cnt[i];
      for ( j = 0; j < count; j++ ) node_ja[index+j] = node_ja[k+j];
      index += count;
      k = node_ia[i+1];
      node_ia[i+1] = index; 
   }

   ML_free( adjacency_cnt );

   /* ------------------------------------------------------------- */
   /* initialize the ML_GGraph data structure                       */
   /* ------------------------------------------------------------- */

   ml_gg->Npoints   = nvertices;
   ml_gg->ML_rank   = mypid;
   ml_gg->row_ptr   = node_ia;
   ml_gg->col_ptr   = node_ja;

   /* ------------------------------------------------------------- */
   /* form the remote portion of the parallel graph                 */
   /* ------------------------------------------------------------- */

   if ( nprocs > 1 ) 
   {

      /* ---------------------------------------------------------- */
      /* count the number of references to remote nodes             */
      /* ---------------------------------------------------------- */

      count = 0;
      nbytes = node_ia[nvertices] * sizeof(int);
      ML_memory_alloc( (void **) &templist, nbytes, "gg1" );
      max_remote_index = nvertices - 1;
      for ( i = 0; i < node_ia[nvertices]; i++ )
      {
         index = node_ja[i];
         if ( index >= nvertices )
         {
            if ( index > max_remote_index ) max_remote_index = index;
            m = grid_fcns->USR_grid_get_vertex_global_num( grid, index );
            ML_search_insert_sort( m, templist, &count, 0 );
         }
      }

      /* ---------------------------------------------------------- */
      /* tabulate total number of remote vertices on all processors */
      /* and put it in total_count, while proc_array is an 'ia'     */
      /* array which indicates which portion of vertex numbers      */
      /* (global no) in remote_list will belong (are remote         */
      /* vertices) to which processor                               */
      /* ---------------------------------------------------------- */

      nbytes = ( nprocs + 1 ) * sizeof( int );
      ML_memory_alloc( (void **) &proc_array, nbytes, "gg2" );
      nbytes = nprocs * sizeof( int );
      ML_memory_alloc( (void **) &inttmp, nbytes, "gg3" );
      for ( i = 0; i < nprocs; i++ ) proc_array[i] = 0;
      proc_array[mypid] = count;
      ML_gsum_vec_int( &proc_array, &inttmp, nprocs, comm );
      ML_memory_free( (void **) &inttmp );
      total_count = 0;
      for ( i = 0; i < nprocs; i++ )
      {
         index  = proc_array[i];
         proc_array[i] = total_count;
         total_count += index;
      }
      proc_array[nprocs] = total_count;

      /* ---------------------------------------------------------- */
      /* next form a global array containing all remote references  */
      /* on all processors                                          */
      /* ---------------------------------------------------------- */

      nbytes = total_count * sizeof( int );
      ML_memory_alloc( (void **) &remote_list, nbytes, "gg4" );
      for ( i = 0; i < count; i++ ) remote_list[i] = templist[i];
      ML_memory_free( (void **) &templist );
      ML_Comm_GappendInt( comm, remote_list, &count, total_count );

      /* ---------------------------------------------------------- */
      /* load the global_list with global numbers of local vertices */
      /* ---------------------------------------------------------- */

      nbytes = nvertices * sizeof( int );
      ML_memory_alloc( (void **) &global_list, nbytes, "gg5" );
      for ( i = 0; i < nvertices; i++ )
      {
         global_list[i] = 
            grid_fcns->USR_grid_get_vertex_global_num( grid, i );
      }

      /* ---------------------------------------------------------- */
      /* construct the send_list based on data in remote_list and   */
      /* proc_array                                                 */
      /* ---------------------------------------------------------- */
   
      nbytes = nprocs * sizeof( int );
      ML_memory_alloc( (void **) &templist, nbytes, "gg6" );
      for ( i = 0; i < nprocs; i++ ) templist[i] = 0;
      for ( i = 0; i < nprocs; i++ )
      {
         if ( i != mypid ) 
         {
            for ( j = proc_array[i]; j < proc_array[i+1]; j++ )
            {
               index = remote_list[j];
               status = ML_sorted_search( index, nvertices, global_list );
               if ( status >= 0 ) templist[i]++;
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* templist[i] != 0 means that my processor has some nodes    */
      /* referenced by processor i.  Next construct the send list.  */
      /* ---------------------------------------------------------- */

      send_cnt = 0;
      for ( i = 0; i < nprocs; i++ )
      {
         if ( templist[i] > 0 ) send_cnt++;
      }
      if ( send_cnt > 0 ) 
      {
         nbytes = send_cnt * sizeof( int );
         ML_memory_alloc( (void **) &send_proc, nbytes, "gg7" );
         ML_memory_alloc( (void **) &send_leng, nbytes, "gg8" );
         nbytes = send_cnt * sizeof( int* );
         ML_memory_alloc( (void **) &send_list, nbytes, "gg9" );
         ML_memory_alloc( (void **) &send_list2, nbytes, "gga" );
      }
      else
      {
         send_proc  = NULL;
         send_leng  = NULL;
         send_list  = NULL;
         send_list2 = NULL;
      }
      index = 0;
      for ( i = 0; i < nprocs; i++ )
      {
         if ( templist[i] > 0 ) 
         {
            send_leng[index]  = templist[i];
            nbytes = templist[i] * sizeof( int );
            ML_memory_alloc( (void **) &send_list[index],nbytes,"ggb" );
            ML_memory_alloc( (void **) &send_list2[index],nbytes,"ggc");
            send_proc[index]  = i;
            index++;
         }
      }

      for ( i = 0; i < send_cnt; i++ ) send_leng[i] = 0;
      for ( i = 0; i < nprocs; i++ )
      {
         if ( i != mypid ) 
         {
            for ( k = 0; k < send_cnt; k++ )
               if ( send_proc[k] == i ) break;
            for ( j = proc_array[i]; j < proc_array[i+1]; j++ )
            {
               index = remote_list[j];
               status = ML_sorted_search( index, nvertices, global_list );
               if ( status >= 0 )
               {
                   send_list[k][send_leng[k]]    = status;
                   send_list2[k][send_leng[k]++] = global_list[status];
               }
            }
         }
      }
      ML_memory_free( (void **) &global_list );
      ML_memory_free( (void **) &remote_list );
      
      /* ---------------------------------------------------------- */
      /* each processor notifies the processors it will send data   */
      /* to about its intention to send and the length of data sent */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < nprocs; i++ ) proc_array[i] = 0;
      for ( i = 0; i < send_cnt; i++ )
      { 
          if ( send_leng[i] != 0 ) proc_array[send_proc[i]] = 1;
      }
      nbytes = nprocs * sizeof( int );
      ML_memory_alloc( (void **) &inttmp, nbytes, "ggd" );
      ML_gsum_vec_int( &proc_array, &inttmp, nprocs, comm );
      recv_cnt = proc_array[mypid];
      msgtype = 539;
      for ( i = 0; i < nprocs; i++ ) inttmp[i] = 0;
      if ( recv_cnt > 0 )
      {
         nbytes = recv_cnt * sizeof( USR_REQ );
         ML_memory_alloc( (void **) &Request, nbytes, "gge" );
         nbytes = recv_cnt * sizeof( int );
         ML_memory_alloc( (void **) &intarray, nbytes, "ggf" );
      }   
      for ( i = 0; i < recv_cnt; i++ )
      {
         fromproc = -1;
         comm->USR_irecvbytes((char*) &intarray[i], sizeof(int), &fromproc,
#ifdef ML_CPP
                     &msgtype, comm->USR_comm, &Request[i] );
#else
                     &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      }
      for ( i = 0; i < send_cnt; i++ )
      {
         comm->USR_sendbytes((void*) &send_leng[i], sizeof(int), 
                             send_proc[i], msgtype, comm->USR_comm );
      }
      for ( i = 0; i < recv_cnt; i++ )
      {
         fromproc = -1;
         comm->USR_cheapwaitbytes((char*) &intarray[i], sizeof(int), &fromproc,
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i] );
#else
                           &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
         inttmp[fromproc] = intarray[i];
      }
      if ( recv_cnt > 0 )
      {
         nbytes = recv_cnt * sizeof( int );   
         ML_memory_alloc( (void **) &recv_proc, nbytes, "ggg" );
         ML_memory_alloc( (void **) &recv_leng, nbytes, "ggh" );
         nbytes = recv_cnt * sizeof( int * );   
         ML_memory_alloc( (void **) &recv_list, nbytes, "ggi" );
         recv_proc = (int *) ML_allocate( recv_cnt * sizeof( int ) );   
         recv_leng = (int *) ML_allocate( recv_cnt * sizeof( int ) );   
      }   
      recv_cnt = tot_recv_leng = 0;
      for ( i = 0; i < nprocs; i++ )
      {
         if ( inttmp[i] > 0 ) {
            recv_proc[recv_cnt]   = i;
            recv_leng[recv_cnt]   = inttmp[i];
            nbytes = recv_leng[recv_cnt] * sizeof( int );
            recv_list[recv_cnt++] = (int*) ML_allocate( nbytes );
            tot_recv_leng        += inttmp[i];
         }
      }
      ML_memory_free( (void **) &inttmp );

      /* ---------------------------------------------------------- */
      /* each processor sends its own send_list to the receiving    */
      /* processors so that they can reconstruct its communication  */ 
      /* pattern                                                    */
      /* ---------------------------------------------------------- */

      msgtype = 541;
      for ( i = 0; i < recv_cnt; i++ )
      {
         fromproc = recv_proc[i];
         nbytes = recv_leng[i] * sizeof( int );
         comm->USR_irecvbytes((char*) recv_list[i], nbytes, &fromproc,
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      }
      for ( i = 0; i < send_cnt; i++ )
      {
         nbytes = send_leng[i] * sizeof( int );
         comm->USR_sendbytes((void*) send_list2[i], nbytes,
                             send_proc[i], msgtype, comm->USR_comm );
      }
      for ( i = 0; i < recv_cnt; i++ )
      {
         fromproc = recv_proc[i];
         nbytes = recv_leng[i] * sizeof( int );
         comm->USR_cheapwaitbytes((char*) recv_list[i], nbytes, &fromproc,
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      }

      /* ---------------------------------------------------------- */
      /* convert the indices in recv_list to local nodes            */
      /* ---------------------------------------------------------- */

      ML_memory_free( (void**) &templist );
      m = max_remote_index - nvertices + 1;
      if ( m > 0 )
      {
         nbytes = m * sizeof( int );
         ML_memory_alloc( (void **) &templist, nbytes, "ggj" );
      }
      for ( i = nvertices; i <= max_remote_index; i++ )
      {
         k = grid_fcns->USR_grid_get_vertex_global_num( grid, i );
         templist[i-nvertices] = k;   
      }
      for ( i = 0; i < recv_cnt; i++ ) 
      {
         for ( j = 0; j < recv_leng[i]; j++ ) 
         {
            index = recv_list[i][j];
            for ( k = 0; k < m; k++ ) 
            {
               if ( templist[k] == index ) break;
            }
            recv_list[i][j] = k + nvertices;
         }
      }

      /* ---------------------------------------------------------- */
      /* fill in the ML_GGraph structure                            */
      /* ---------------------------------------------------------- */

      ml_gg->send_cnt  = send_cnt;
      ml_gg->send_list = send_list;
      ml_gg->send_leng = send_leng;
      ml_gg->send_proc = send_proc;
      ml_gg->recv_cnt  = recv_cnt;
      ml_gg->recv_list = recv_list;
      ml_gg->recv_leng = recv_leng;
      ml_gg->recv_proc = recv_proc;

      /* ---------------------------------------------------------- */
      /* clean up                                                   */
      /* ---------------------------------------------------------- */

      if ( recv_cnt > 0 )
      {
         ML_memory_free( (void **) &Request );
         ML_memory_free( (void **) &intarray );
      }   
      if ( send_list2 != NULL )
      {
         for ( i = 0; i < send_cnt; i++ ) 
            ML_memory_free( (void**) &send_list2[i] );
         ML_memory_free( (void**) &send_list2 );
      }
      ML_memory_free( (void **) &templist );
      ML_memory_free( (void **) &proc_array );
   }
   return 0;
}

/* ******************************************************************** */
/* graph coarsening                                                     */
/* -------------------------------------------------------------------- */

int ML_GGraph_Coarsen(ML_GGraph *ml_gg, ML_Comm *comm) 
{
   int   i, j, m, myrank, nvertices, vertex_cnt, index, *short_list;
   int   short_cnt, *templist, **proclist, *rptr, *cptr, ext_nvertices;
   int   send_cnt, *send_leng, **send_list, *send_proc;
   int   recv_cnt, *recv_leng, **recv_list, *recv_proc, msgtype;
   int   char_leng, fproc, offset, nselected;
   int   **send_buf, **recv_buf;
   char  *vertex_state, *btypes, *vertex_type, *send_carray, *recv_carray;
   USR_REQ *Request;
   unsigned int nbytes;

   /* ------------------------------------------------------------- */
   /* error checking (for proper ML_GGraph structure)               */
   /* ------------------------------------------------------------- */
   if ( ml_gg->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_Coarsen : wrong object. \n");
      exit(1);
   }
   if ( ml_gg->Npoints <= 0 )
   {
      printf("ML_GGraph_Coarsen : Npoints <= 0. \n");
      exit(1);
   }
   if ( ml_gg->bdry_type == NULL )
   {
      printf("ML_GGraph_Coarsen : no boundary types. \n");
      exit(1);
   }
   if ( ml_gg->row_ptr == NULL || ml_gg->col_ptr == NULL )
   {
      printf("ML_GGraph_Coarsen : no graph. \n");
      exit(1);
   }

   /* ------------------------------------------------------------- */
   /* fetch the data in the ML_GGraph data structure                */
   /* ------------------------------------------------------------- */
  
   myrank            = ml_gg->ML_rank;
   nvertices         = ml_gg->Npoints;
   btypes            = ml_gg->bdry_type;
   rptr              = ml_gg->row_ptr;
   cptr              = ml_gg->col_ptr;
   send_cnt          = ml_gg->send_cnt;
   send_leng         = ml_gg->send_leng;
   send_list         = ml_gg->send_list;
   send_proc         = ml_gg->send_proc;
   recv_cnt          = ml_gg->recv_cnt;
   recv_leng         = ml_gg->recv_leng;
   recv_list         = ml_gg->recv_list;
   recv_proc         = ml_gg->recv_proc;

   /* ------------------------------------------------------------- */
   /* construct the proclist array, which contains information      */
   /* about where the vertex information has to be sent to for      */
   /* updates.  If the vertex is a remote vertex, it contains info  */
   /* about where the receive processor and the receive index       */
   /* ------------------------------------------------------------- */

   /* search for the maximum index of the external vertices */

   ext_nvertices = nvertices - 1;
   for ( i = 0; i < rptr[nvertices]; i++ )
   {
      if ( cptr[i] > ext_nvertices ) ext_nvertices = cptr[i];
   }
   ext_nvertices++;

   /* find out how many other processors each of my local vertices */
   /* needs to update in the case of changing states               */

   nbytes = nvertices * sizeof( int );
   ML_memory_alloc( (void**) &templist, nbytes, "ggk" );
   for ( i = 0; i < nvertices; i++ ) templist[i] = 0;
   for ( i = 0; i < send_cnt; i++ )
   {
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[i][j];     
         if ( index >= nvertices || index < 0 )
         {
            printf("%d : Error : in Coarsening.\n", myrank);
            exit(0);
         }
         templist[index]++;
      }
   }

   /* Allocate proclist to record the processors and indices each of */
   /* my local vertices are to send.  The first element of the array */
   /* is a counter of how many processors, followed by a number of   */
   /* processor and index pairs.                                     */

   nbytes = ext_nvertices * sizeof( int *);
   ML_memory_alloc( (void**) &proclist, nbytes, "gf2" );
   for ( i = 0; i < nvertices; i++ ) 
   {
      proclist[i] = (int *) ML_allocate( (2*templist[i]+1) * sizeof( int ) ); 
      proclist[i][0] = 0;
      templist[i] = 0;
   }
   for ( i = 0; i < send_cnt; i++ )
   {
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[i][j];     
         proclist[index][templist[index]+1] = i;
         proclist[index][templist[index]+2] = j;
         templist[index] += 2;
         proclist[index][0]++;
      }
   }
   for ( i = nvertices; i < ext_nvertices; i++ ) 
   {
      proclist[i] = (int *) ML_allocate( sizeof( int ) ); 
   }
   for ( i = 0; i < recv_cnt; i++ )
   {
      for ( j = 0; j < recv_leng[i]; j++ )
      {
         index = recv_list[i][j];     
         proclist[index][0] = recv_proc[i];
      }
   }
   ML_memory_free( (void **) &templist );

   /* ------------------------------------------------------------- */
   /* fill in the initial state of all local and remote vertices    */
   /* (initial vertex state = 'F' - free)                           */
   /* ------------------------------------------------------------- */
   
   vertex_cnt   = ext_nvertices;
   nbytes = vertex_cnt * sizeof( char );
   ML_memory_alloc( (void **) &vertex_state, nbytes, "ggl" );
   ml_gg->vertex_state = vertex_state;
   for ( i = 0; i < vertex_cnt; i++ ) vertex_state[i] = 'F';

   /* ------------------------------------------------------------- */
   /* fill in the boundary type of local vertices                   */
   /* ------------------------------------------------------------- */

   nbytes = vertex_cnt * sizeof( char );
   ML_memory_alloc( (void **) &vertex_type, nbytes, "ggm" );
   for ( i = 0; i < nvertices; i++ ) vertex_type[i] = btypes[i];

   /* ------------------------------------------------------------- */
   /* fill in the boundary type of remote vertices                  */
   /* ------------------------------------------------------------- */

   char_leng = 0;
   for ( i = 0; i < recv_cnt; i++ ) char_leng += recv_leng[i];
   if ( char_leng > 0 )
   {
      nbytes = char_leng * sizeof( char );
      ML_memory_alloc( (void **) &recv_carray, nbytes, "gf5" );
   }
   msgtype = 679;
   if ( recv_cnt > 0 )
   {
      nbytes = recv_cnt * sizeof( USR_REQ );
      ML_memory_alloc( (void **) &Request, nbytes, "gf6" );
   }  
   offset = 0;
   for ( i = 0; i < recv_cnt; i++ )
   {
      fproc = recv_proc[i];
      nbytes = recv_leng[i] * sizeof(char);   
      comm->USR_irecvbytes((char*) &recv_carray[offset], nbytes, &fproc, 
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      offset += recv_leng[i];
   }
   for ( i = 0; i < send_cnt; i++ )
   {
      nbytes = send_leng[i] * sizeof( char );
      send_carray = (char *) ML_allocate( nbytes );
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[i][j];
         send_carray[j] = btypes[index];
      }
      comm->USR_sendbytes((void*) send_carray, nbytes, send_proc[i], 
                          msgtype, comm->USR_comm );
      ML_free( send_carray );
   }
   offset = 0;
   for ( i = 0; i < recv_cnt; i++ )
   {
      fproc = recv_proc[i];
      nbytes = recv_leng[i] * sizeof(char);   
      comm->USR_cheapwaitbytes((char*) &recv_carray[offset], nbytes, &fproc, 
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      for ( j = 0; j < recv_leng[i]; j++ )
      {
         index = recv_list[i][j];
         vertex_type[index] = recv_carray[offset+j];
      }
      offset += recv_leng[i];
   }
   if ( recv_cnt > 0 ) ML_memory_free( (void**) &Request );
   if ( char_leng > 0 ) ML_memory_free( (void**) &recv_carray );
   
   /* ------------------------------------------------------------- */
   /* now everything is almost ready                                */
   /*   proclist contains send information for vertex updates       */
   /*   vertex_type contains boundary type info for all vertices    */
   /*   vertex_state contains current state info for all vertices   */
   /* before moving on, make sure the communication buffers have    */
   /* been sufficiently allocated                                   */ 
   /* ------------------------------------------------------------- */

   nbytes = nvertices * sizeof( int );
   ML_memory_alloc( (void**) &short_list, nbytes, "gf7" );
   short_cnt = 0;
   if ( send_cnt > 0 )
   {
      nbytes = send_cnt * sizeof( int * );
      ML_memory_alloc( (void**) &send_buf, nbytes, "gf8" );
      for ( i = 0; i < send_cnt; i++ )
      {
         nbytes = (send_leng[i] + 1 ) * sizeof( int );
         ML_memory_alloc( (void**) &send_buf[i], nbytes, "gf9" );
         for ( j = 0; j <= send_leng[i]; j++ ) send_buf[i][j] = 0;
      }
   }
   if ( recv_cnt > 0 )
   {
      nbytes = recv_cnt * sizeof( int * );
      ML_memory_alloc( (void**) &recv_buf, nbytes, "gfa" );
      for ( i = 0; i < recv_cnt; i++ )
      {
         nbytes = (recv_leng[i] + 1) * sizeof( int );
         ML_memory_alloc( (void**) &recv_buf[i], nbytes, "gfb" );
      }
      nbytes = recv_cnt * sizeof( USR_REQ );
      ML_memory_alloc( (void**) &Request, nbytes, "gfc" );
   }

   /* ------------------------------------------------------------- */
   /* now do coarsening of the corner vertex                        */
   /* ------------------------------------------------------------- */

   nselected = 0;
   short_cnt = 0;
   for ( i = 0; i < nvertices; i++ )
   {
      if ( btypes[i] == 'C' ) short_list[short_cnt++] = i;
   }

   m = ML_GGraph_LabelVertices(short_cnt, short_list, 'C', vertex_state,
             vertex_type, nvertices, rptr, cptr, myrank, proclist, 
             send_cnt, send_buf, send_proc, send_leng, recv_cnt, 
             recv_buf, recv_proc, recv_leng, recv_list, 9413, comm);
   nselected += m;

   /* ------------------------------------------------------------- */
   /* now do coarsening of the ridge vertices                       */
   /* ------------------------------------------------------------- */

   short_cnt = 0;
   for ( i = 0; i < nvertices; i++ )
   {
      if ( btypes[i] == 'R' && vertex_state[i] == 'F' ) 
         short_list[short_cnt++] = i;
   }

   m = ML_GGraph_LabelVertices(short_cnt, short_list, 'R', vertex_state,
             vertex_type, nvertices, rptr, cptr, myrank, proclist, 
             send_cnt, send_buf, send_proc, send_leng, recv_cnt, 
             recv_buf, recv_proc, recv_leng, recv_list, 9414, comm);
   nselected += m;

   /* ------------------------------------------------------------- */
   /* now do coarsening of the face vertices                        */
   /* ------------------------------------------------------------- */

   short_cnt = 0;
   for ( i = 0; i < nvertices; i++ )
   {
      if ( btypes[i] == 'F' && vertex_state[i] == 'F' ) 
         short_list[short_cnt++] = i;
   }

   m = ML_GGraph_LabelVertices(short_cnt, short_list, 'F', vertex_state,
             vertex_type, nvertices, rptr, cptr, myrank, proclist, 
             send_cnt, send_buf, send_proc, send_leng, recv_cnt, 
             recv_buf, recv_proc, recv_leng, recv_list, 9415, comm);
   nselected += m;

   /* ------------------------------------------------------------- */
   /* finally do coarsening of the interior vertices                */
   /* ------------------------------------------------------------- */

   short_cnt = 0;
   for ( i = 0; i < nvertices; i++ )
   {
      if ( btypes[i] == 'I' && vertex_state[i] == 'F' ) 
         short_list[short_cnt++] = i;
   }

   m = ML_GGraph_LabelVertices(short_cnt, short_list, 'I', vertex_state,
             vertex_type, nvertices, rptr, cptr, myrank, proclist, 
             send_cnt, send_buf, send_proc, send_leng, recv_cnt, 
             recv_buf, recv_proc, recv_leng, recv_list, 9416, comm);
   nselected += m;

   /*
   for ( j = 0; j < nvertices; j++ )
   {
      printf("%d : vertex %d = %c \n", myrank, j, vertex_state[j]);
   }
   */
   ml_gg->Nselected = nselected;
   nselected = ML_Comm_GsumInt( comm, nselected );
   printf("%d : total selected = %d \n", myrank, nselected );

   /* ------------------------------------------------------------- */
   /* Finally, clean up                                             */
   /* ------------------------------------------------------------- */
   if ( recv_cnt > 0 ) ML_memory_free( (void**) &Request );
   ML_memory_free( (void**) &short_list );
   if ( send_cnt > 0 )
   {
      for ( i = 0; i < send_cnt; i++ ) ML_memory_free((void**)&send_buf[i]);
      ML_memory_free( (void**) &send_buf );
   }
   if ( recv_cnt > 0 )
   {
      for ( i = 0; i < recv_cnt; i++ ) ML_memory_free((void**)&recv_buf[i]);
      ML_memory_free( (void**) &recv_buf );
   }
   if ( nvertices > 0 )
   {
      for (i = 0; i < nvertices; i++) ML_free( proclist[i] );
      ML_memory_free( (void **) &proclist );
   }
   ML_memory_free( (void **) &vertex_type );
   return 0;
}

/* ******************************************************************** */
/* A subroutine to get the states of the nodes                          */
/* -------------------------------------------------------------------- */
int ML_GGraph_Get_NodeStates(ML_GGraph *ml_gg, int *size, char **states)
{
   (*size) = ml_gg->Npoints;
   (*states) = ml_gg->vertex_state;
   return 0;
}

/* ******************************************************************** */
/* A subroutine to label vertices of a particular type                  */
/* -------------------------------------------------------------------- */
int ML_GGraph_LabelVertices(int vlist_cnt, int *vlist, int Vtype,
                           char *vertex_state, char *vertex_type,
                           int nvertices, int *rptr, int *cptr, 
                           int myrank, int **proclist, int send_cnt, 
                           int **send_buf, int *send_proc, int *send_leng,
                           int recv_cnt, int **recv_buf, int *recv_proc, 
                           int *recv_leng, int **recv_list, int msgtype, 
                           ML_Comm *comm)
{
  int     i, j, k, m, temp_cnt, index, select_flag, fproc, col;
   int     loop_flag, change_flag, *proc_flag, nselected, *tlist, pref_cnt;
   int     *pref_list, col2, loop_cnt;
   int     pref_flag, pref_index, delete_flag;
   USR_REQ *Request;
   unsigned int nbytes;

   /* ---------------------------------------------------------- */
   /* give the vertices adjacent to deleted vertices preferences */
   /* ---------------------------------------------------------- */

   if ( vlist_cnt > 0 )
   {
      nbytes = vlist_cnt * sizeof( int );
      ML_memory_alloc((void**) &tlist, nbytes, "ggn" );
      for ( i = 0; i < vlist_cnt; i++ ) tlist[i] = vlist[i];
      for ( i = 0; i < vlist_cnt; i++ )
      {
         index = tlist[i];
         for ( j = rptr[index]; j < rptr[index+1]; j++ )
         {
            col = cptr[j];
            if ( vertex_state[col] == 'D' )
            {
               tlist[i] = - index;
               break;
            }
         }
      }
      m = 0;
      for ( i = 0; i < vlist_cnt; i++ )
      {
         if ( tlist[i] < 0 ) vlist[m++] = - tlist[i];
      }
      for ( i = 0; i < vlist_cnt; i++ )
      {
         if ( tlist[i] >= 0 ) vlist[m++] = tlist[i];
      }
      ML_memory_free( (void**) &tlist );
   }
   if ( nvertices > 0 )
   {
      nbytes = nvertices * sizeof( int );
      ML_memory_alloc((void**) &pref_list, nbytes, "ggo" );
   }   
   pref_cnt = 0;
   
   /* -------------------------------------------------------- */
   /* get ready for the coarsening                             */
   /* -------------------------------------------------------- */

   temp_cnt = vlist_cnt;
   if ( recv_cnt > 0 )
   {
      nbytes = recv_cnt * sizeof( USR_REQ );
      ML_memory_alloc((void**) &Request, nbytes, "ggp" );
      nbytes = recv_cnt * sizeof( int );
      ML_memory_alloc((void**) &proc_flag, nbytes, "ggq" );
      for ( i = 0; i < recv_cnt; i++ ) proc_flag[i] = 0;
   }
   nselected = 0;

   /* -------------------------------------------------------- */
   /* let's actually do coarsening                             */
   /* -------------------------------------------------------- */

   loop_flag = recv_cnt;
   change_flag = 1;
   loop_cnt = 0;
   pref_index = 0;     /* pointer to a stack of vertex numbers */

   do {
      /* loop_cnt is to monitor the performance of coarsening */

      loop_cnt++;

      /* reset all buffers to zero only if it has been changed */

      if ( change_flag == 1 )
      {
         for ( j = 0; j < send_cnt; j++ )
            for ( k = 0; k <= send_leng[j]; k++ ) send_buf[j][k] = 0;
         change_flag = 0;
      }

      /* examine the vertices in vlist */

      for ( i = 0; i < vlist_cnt; i++ )
      {

         /* handle the preference list first, if there is any */
         /* Note : we want to fetch the pref_list from the    */
         /*        front                                      */

         index = vlist[i];
         pref_flag = 0;
         if ( pref_cnt > 0 )
         {
            index = pref_list[pref_index];    
            for (j = pref_index+1; j < pref_cnt; j++) 
               pref_list[j-1] = pref_list[j];
            pref_cnt--;
            pref_flag = 1;
            i--;
         }

         /* if the vertex in question has not been considered F(ree) */

         delete_flag = 0;
         if ( vertex_state[index] == 'F' )
         {
            select_flag = 1;
            for ( j = rptr[index]; j < rptr[index+1]; j++ )
            {
               /* if its neighbor is selected, delete this vertex */

               col = cptr[j];
               if ( vertex_state[col] == 'S' )
               {
                  vertex_state[index] = 'D';
                  delete_flag = 1;
                  temp_cnt--;  
                  select_flag = 0;
                  break;
               }
               
               /* If its neighbor is of the same type and not been   */
               /* considered. Furthermore, if it is a remote vertex  */
               /* and its owner processor has rank smaller than mine,*/
               /* my processor should wait(thus turn off select_flag)*/

               else if ( vertex_type[col] == Vtype && 
                         vertex_state[col] == 'F')
               {
                  if ( col >= nvertices )
                  {
                     if ( proclist[col][0] < myrank )
                     {
                        /*
                        printf("%d : %d not selected due to N\n",myrank,index);
                        */
                        select_flag = 0;
                        break;
                     }
                  }
               }
            }

            /* if the vertex in question is not any of those considered */
            /* above, select this vertex.                               */

            if ( select_flag == 1 )
            {
               /* printf("%d : %d selected\n", myrank, index); */
               vertex_state[index] = 'S';
               nselected++;
               temp_cnt--;

               /* set the flag that this vertex has been selected in */
               /* the buffer which is to be sent to other processors */

               for ( k = 0; k < proclist[index][0]; k++ )
               {
                  fproc = proclist[index][2*k+1];
                  m     = proclist[index][2*k+2];
                  send_buf[fproc][m] = 1;
                  change_flag = 1;
               }

               /* delete all vertices adjacent to this vertex and */
               /* indicate that also in the communication buffer  */

               for ( j = rptr[index]; j < rptr[index+1]; j++ )
               {
                  col = cptr[j];
                  vertex_state[col] = 'D';
                  if ( col < nvertices )
                  {
                     for ( k = 0; k < proclist[col][0]; k++ )
                     {
                        fproc = proclist[col][2*k+1];
                        m     = proclist[col][2*k+2];
                        send_buf[fproc][m] = 2;
                        change_flag = 1;
                     }

                     /* also, put the next set of vertices into the  */
                     /* preference list (try to mimic the sequential */
                     /* maximally independent set algorithm          */

                     for ( k = rptr[col]; k < rptr[col+1]; k++ )
                     {
                        col2 = cptr[k];
                        if (col2 < nvertices && 
                            vertex_state[col2] == 'F' &&
                            vertex_type[col2] == Vtype )
                        {
                           ML_search_insert_sort(col2,pref_list,&pref_cnt,0);
                        }
                        if (pref_cnt >= nvertices)
                           printf("%d : warning : overflow. \n", myrank);
                     }
                  }
               }
            } 

            /* if not selected due to it has just been deleted, */
            /* still have to tell the other processors          */

            else if ( select_flag == 0 && delete_flag == 1 )
            {
               for ( k = 0; k < proclist[index][0]; k++ )
               {
                  fproc = proclist[index][2*k+1];
                  m     = proclist[index][2*k+2];
                  send_buf[fproc][m] = 2;
                  change_flag = 1;
               }
            }
         }

         /* if after the steps above, the vertex is still not */
         /* selected.  Well, do something about it.           */

         if ( vertex_state[index] == 'F' )
         {
            /* if a vertex in the pref_list has been considered */
            /* but not selected, need to put the vertex back to */
            /* the list, and move on to consider the next one   */
            /* (i.e. advance pref_index)                        */

            if ( pref_flag == 1 )
            {
               for (j = pref_index; j < pref_cnt; j++) 
                  pref_list[j+1] = pref_list[j];
               pref_list[pref_index] = index;
               pref_index++;
               if ( pref_index >= pref_cnt ) pref_index = 0;
            }
         } else if ( pref_flag != 1 ) temp_cnt--;
      }

      /* update the states to/from other processors */

      for ( j = 0; j < recv_cnt; j++ )
      {
         if ( proc_flag[j] == 0 )
         {
            fproc = recv_proc[j];
            nbytes = (recv_leng[j] + 1) * sizeof( int );
            comm->USR_irecvbytes((char*) recv_buf[j], nbytes, &fproc,
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
         }
      }
      for ( j = 0; j < send_cnt; j++ )
      {
         nbytes = (send_leng[j] + 1) * sizeof( int );
         if ( temp_cnt <= 0 ) send_buf[j][send_leng[j]] = 1;
         comm->USR_sendbytes((void*) send_buf[j], nbytes,
                             send_proc[j], msgtype, comm->USR_comm );
      }
      for ( j = 0; j < recv_cnt; j++ )
      {
         if ( proc_flag[j] == 0 )
         {
            fproc = recv_proc[j];
            nbytes = (recv_leng[j] + 1) * sizeof( int );
            comm->USR_cheapwaitbytes((char*) recv_buf[j], nbytes, &fproc,
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
            for ( k = 0; k < recv_leng[j]; k++ )
            {
               index = recv_list[j][k];
               if      (recv_buf[j][k] == 1)
               {
                  vertex_state[index] = 'S';
                  /* printf("%d : incoming state %d = S \n",myrank,index); */
               }
               else if (recv_buf[j][k] == 2)
               {
                  vertex_state[index] = 'D';
                  /* printf("%d : incoming state %d = D \n",myrank,index); */
               }
            }
            if ( recv_buf[j][recv_leng[j]] == 1 )
            {
               proc_flag[j] = 1;
               loop_flag--;
            }
         }
      }
   } while ( loop_flag > 0 || temp_cnt > 0 );

   /* printf("%d : loop_count = %d \n", myrank, loop_cnt ); */
   if ( recv_cnt > 0 )
   {
      ML_memory_free( (void **) &proc_flag );
      ML_memory_free( (void **) &Request );
   }
   if ( nvertices > 0 ) ML_memory_free( (void **) &pref_list );
   return nselected;
}

/* ******************************************************************** */
/* Check to see if the coarse graph is a maximally independent set      */
/* -------------------------------------------------------------------- */

int ML_GGraph_CheckMIS( ML_GGraph *ml_gg, ML_Comm *comm )
{
   int  myrank, *cptr, *rptr, send_cnt, recv_cnt, ext_nvertices;
   int  *send_leng, *send_proc, **send_list, nvertices, index; 
   int  *recv_leng, *recv_proc, **recv_list, msgtype, num_faults;
   int  i, j, char_leng, offset, fproc, nselected, fault_flag;
   char *vertex_state_here, *recv_carray, *send_carray;
   USR_REQ *Request;
   unsigned int nbytes;

   /* ------------------------------------------------------------- */
   /* error checking (for proper ML_GGraph structure)               */
   /* ------------------------------------------------------------- */
   if ( ml_gg->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_CheckMIS : wrong object. \n");
      exit(1);
   }
   if ( ml_gg->vertex_state == NULL )
   {
      printf("Warning : Graph not coarsened yet. \n");
      return -1;
   }
   
   /* ------------------------------------------------------------- */
   /* fetch the data in the ML_GGraph data structure                */
   /* ------------------------------------------------------------- */
  
   myrank            = ml_gg->ML_rank;
   nvertices         = ml_gg->Npoints;
   rptr              = ml_gg->row_ptr;
   cptr              = ml_gg->col_ptr;
   send_cnt          = ml_gg->send_cnt;
   send_leng         = ml_gg->send_leng;
   send_list         = ml_gg->send_list;
   send_proc         = ml_gg->send_proc;
   recv_cnt          = ml_gg->recv_cnt;
   recv_leng         = ml_gg->recv_leng;
   recv_list         = ml_gg->recv_list;
   recv_proc         = ml_gg->recv_proc;

   /* ------------------------------------------------------------- */
   /* search for the maximum index of the external vertices         */
   /* ------------------------------------------------------------- */

   ext_nvertices = nvertices - 1;
   for ( i = 0; i < rptr[nvertices]; i++ )
   {
      if ( cptr[i] > ext_nvertices ) ext_nvertices = cptr[i];
   }
   ext_nvertices++;

   /* ------------------------------------------------------------- */
   /* Create a local vertex state array to accommodate also vertex  */
   /* states from other processors. Then put the local vertex states*/
   /* into this local array                                         */
   /* ------------------------------------------------------------- */

   nbytes = ext_nvertices * sizeof( char );
   ML_memory_alloc( (void **) &vertex_state_here, nbytes, "gh1" );
   for ( i = 0; i < nvertices; i++ )
      vertex_state_here[i] = ml_gg->vertex_state[i];

   /* ------------------------------------------------------------- */
   /* fill in the vertex states of the remote vertices              */
   /* ------------------------------------------------------------- */

   char_leng = 0;
   for ( i = 0; i < recv_cnt; i++ ) char_leng += recv_leng[i];
   if ( char_leng > 0 )
   {
      nbytes = char_leng * sizeof( char );
      ML_memory_alloc( (void **) &recv_carray, nbytes, "gh2" );
   }
   msgtype = 23945;
   if ( recv_cnt > 0 )
   {
      nbytes = recv_cnt * sizeof( USR_REQ );
      ML_memory_alloc( (void **) &Request, nbytes, "gh3" );
   }  
   offset = 0;
   for ( i = 0; i < recv_cnt; i++ )
   {
      fproc = recv_proc[i];
      nbytes = recv_leng[i] * sizeof(char);   
      comm->USR_irecvbytes((char*) &recv_carray[offset], nbytes, &fproc, 
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      offset += recv_leng[i];
   }
   for ( i = 0; i < send_cnt; i++ )
   {
      nbytes = send_leng[i] * sizeof( char );
      send_carray = (char *) ML_allocate( nbytes );
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[i][j];
         send_carray[j] = vertex_state_here[index];
      }
      comm->USR_sendbytes((void*) send_carray, nbytes, send_proc[i], 
                          msgtype, comm->USR_comm );
      ML_free( send_carray );
   }
   offset = 0;
   for ( i = 0; i < recv_cnt; i++ )
   {
      fproc = recv_proc[i];
      nbytes = recv_leng[i] * sizeof(char);   
      comm->USR_cheapwaitbytes((char*) &recv_carray[offset], nbytes, &fproc, 
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      for ( j = 0; j < recv_leng[i]; j++ )
      {
         index = recv_list[i][j];
         vertex_state_here[index] = recv_carray[offset+j];
      }
      offset += recv_leng[i];
   }
   if ( recv_cnt > 0 ) ML_memory_free( (void**) &Request );
   if ( char_leng > 0 ) ML_memory_free( (void**) &recv_carray );
   
   /* ------------------------------------------------------------- */
   /* first check that all vertex states are either 'S' or 'D'      */
   /* ------------------------------------------------------------- */

   fault_flag = 0;
   for ( i = 0; i < nvertices; i++ )
   {
      if ( vertex_state_here[i] != 'S' && vertex_state_here[i] != 'D' )
         fault_flag++;
   }
   printf("%d : ML_GGraph_CheckMIS : %d vertices are mislabeled.\n",
                                     myrank, fault_flag);

   /* ------------------------------------------------------------- */
   /* then check that all vertices labeled 'S' do not have any      */
   /* neighbors with label 'S'                                      */
   /* ------------------------------------------------------------- */

   nselected = 0;
   fault_flag = 0;
   for ( i = 0; i < nvertices; i++ )
   {
      if ( vertex_state_here[i] == 'S' )
      {
         nselected++;
         for ( j = rptr[i]; j < rptr[i+1]; j++ )
            if ( vertex_state_here[cptr[j]] == 'S' ) fault_flag = 1; 
      }
   }
   printf("%d : ML_GGraph_CheckMIS : nselected = %d\n", myrank, nselected);
   nselected = ML_Comm_GsumInt( comm, nselected );
   if ( myrank == 0 )
      printf("%d : ML_GGraph_CheckMIS : TOTAL SELECTED = %d\n",
                                        myrank, nselected);
   if ( fault_flag == 1 )
      printf("%d : ML_GGraph_CheckMIS : FAILED independent subset test.\n",
                                        myrank);
   else
      printf("%d : ML_GGraph_CheckMIS : PASSED independent subset test.\n",
                                        myrank);

   /* ------------------------------------------------------------- */
   /* finally check that all vertices labeled 'D' do have some      */
   /* neighbors with label 'S'                                      */
   /* ------------------------------------------------------------- */

   num_faults = 0;
   for ( i = 0; i < nvertices; i++ )
   {
      fault_flag = 0;
      if ( vertex_state_here[i] == 'D' )
      {
         fault_flag = 1;
         for ( j = rptr[i]; j < rptr[i+1]; j++ )
            if ( vertex_state_here[cptr[j]] == 'S' ) fault_flag = 0; 
      }
      if ( fault_flag == 1 ) num_faults++;
   }
   num_faults = ML_Comm_GsumInt( comm, num_faults );
   if ( num_faults == 0 )
   {
      printf("%d : ML_GGraph_CheckMIS : PASSED maximality test. \n", myrank);
   }
   else
   {
      printf("%d : ML_GGraph_CheckMIS : FAILED maximality test. \n", myrank);
      if ( myrank == 0 )
         printf("%d : ML_GGraph_CheckMIS : total no. of faults = %d \n",
                                           myrank, num_faults);
   }

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free( (void**) &vertex_state_here );
   return 0;
}

/* ******************************************************************** */
/* compare two element lists to see if they have a common edge          */
/* -------------------------------------------------------------------- */

int ML_GGraph_Find_NeighborElements(int leng1, int *list1, int leng2,
                                   int *list2, int *vlist3)
{
   int match_count, cnt1, cnt2;
 
   cnt1 = cnt2 = match_count = 0;
   while ( cnt1 < leng1 && cnt2 < leng2 ) {
      if ( list1[cnt1] == list2[cnt2] ) {
         vlist3[match_count] = list1[cnt1];
         match_count++;
         cnt1++; cnt2++;
      } else if ( list1[cnt1] > list2[cnt2] ) cnt2++;
      else   cnt1++;
   }
   return match_count;
}

/* ******************************************************************** */
/* construct element graph                                              */
/* -------------------------------------------------------------------- */

int ML_GGraph_Gen_ElementGraph(ML_GGraph *ml_gg,void *grid,void (*gf),
                               ML_Comm *comm) 
{
   int         i, j, nelements, nvertices, *vlist, *vlist2, **v2elem_map;
   int         *egraph_ia, *egraph_ja, egraph_cnt, vlength, index;
   int         count, elem_index, elem_num2, vindex, vertnum, vcount;
   int         eindex2, vlength2, num_common_edges, egraph_n, mypid;
   int         process_flag, index2, vlist3[3];
   char        *vertex_state;
   ML_GridFunc *grid_fcns;

   /* ------------------------------------------------------------- */
   /* initial error checking (for proper ML_GGraph structure)       */
   /* ------------------------------------------------------------- */

   printf("ML_GGraph_Gen_ElementGraph : this is sequential for now. \n");
   if ( ml_gg->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_Gen_ElementGraph : wrong object. \n");
      exit(1);
   }

   /* ------------------------------------------------------------- */
   /* find, for each of the vertices, the elements where the vertex */
   /* belongs to (v2elem_map[0] contains the count, and the rest    */
   /* contains the actual element number)                           */
   /* ------------------------------------------------------------- */

   grid_fcns = (ML_GridFunc *) gf;
   mypid     = comm->ML_mypid;
   nvertices = grid_fcns->USR_grid_get_nvertices( grid );
   nelements = grid_fcns->USR_grid_get_nelements( grid );
   if ( nvertices > 0 ) {
      v2elem_map = (int **) ML_allocate( nvertices * sizeof(int) );
      for ( i = 0; i < nvertices; i++ ) {
         v2elem_map[i] = (int *) ML_allocate( 7 * sizeof(int) );
         v2elem_map[i][0] = 0;
      }
   } else {
      printf("%d : ML_GGraph_Gen_ElementGraph : nvertices <= 0\n", mypid);
      return -1;
   }
   vlist = (int *) ML_allocate( 100 * sizeof(int) );
   for ( i = 0; i < nelements; i++ ) {
      vlength = grid_fcns->USR_grid_get_element_vlist( grid, i, vlist );
      for ( j = 0; j < vlength; j++ ) {
         index = vlist[j];
         count = ++v2elem_map[index][0];
         if ( count >= 7 ) {
            printf("ML_GGraph_Gen_ElementGraph : error - \n"); 
            printf("    not enough local space, tune the code to fix it.\n"); 
            exit(1);
         }
         v2elem_map[index][count] = i;
      }
   }

   /* ------------------------------------------------------------- */
   /* now start with element 0, look for its neighboring elements   */
   /* ------------------------------------------------------------- */

   vlist2 = (int *) ML_allocate( 100 * sizeof(int) );
   egraph_ia = (int *) ML_allocate( ( nelements + 1 ) * sizeof(int));
   egraph_ja = (int *) ML_allocate( 6 * nelements * sizeof(int));
   egraph_ia[0] = 0;
   egraph_cnt = 0;
   egraph_n = nelements;
   vertex_state = ml_gg->vertex_state;
   for ( elem_index = 0; elem_index < nelements; elem_index++ ) {
      vlength = grid_fcns->USR_grid_get_element_vlist(grid,elem_index,vlist);
      ML_sort( vlength, vlist );
      for ( vindex = 0; vindex < vlength; vindex++ ) {
         vertnum = vlist[vindex];
         vcount = v2elem_map[vertnum][0];
         for ( eindex2 = 0; eindex2 < vcount; eindex2++ ) {
            elem_num2 = v2elem_map[vertnum][eindex2+1];
            process_flag = 1;
            if ( elem_index == elem_num2) process_flag = 0; 
            if ( process_flag == 1 ) {
               for ( j = egraph_ia[elem_index]; j < egraph_cnt; j++ )
                  if ( egraph_ja[j] == elem_num2 )
                     { process_flag = 0; break;}
            }
            if ( process_flag == 1 ) {
               vlength2 = grid_fcns->USR_grid_get_element_vlist(grid, 
                                         elem_num2, vlist2 );
               ML_sort( vlength2, vlist2 );
               num_common_edges = ML_GGraph_Find_NeighborElements(vlength, 
                                           vlist, vlength2, vlist2, vlist3);
               if ( num_common_edges == 2 ) {
                  index = vlist3[0];
                  index2 = vlist3[1];
                  if (vertex_state[index] != 'S' && vertex_state[index2] != 'S')
                     egraph_ja[egraph_cnt++] = elem_num2; 
               } 
            } 
         } 
      } 
      egraph_ia[elem_index+1] = egraph_cnt;
   }
 
   for ( i = 0; i < egraph_n; i++ ) {
      for ( j = egraph_ia[i]; j < egraph_ia[i+1]; j++ ) {
         printf("row %5d : column = %5d \n", i, egraph_ja[j]);
      }
   }
   ML_free( vlist );
   ML_free( vlist2 );
   for ( i = 0; i < nvertices; i++ ) ML_free( v2elem_map[i] );
   ML_free(v2elem_map);
   return 0;
}

/* ******************************************************************** */
/* Generate a restriction operator, after the ML_GGraph_Coarsen has     */
/* been called.                                                         */
/* -------------------------------------------------------------------- */

/*
int ML_GGraph_Gen_Restrictor(ML_GGraph *ml_gg, struct ML_CSR_MSRdata *mat)
                             ML_Comm *comm) 
{
*/
   /* ------------------------------------------------------------- */
   /* initial error checking (for proper ML_GGraph structure)       */
   /* ------------------------------------------------------------- */

/*
   if ( ml_gg->ML_id != ML_ID_GGRAPH )
   {
      printf("ML_GGraph_Gen_Restrictor : wrong object. \n");
      exit(1);
   }
   if (ml_gg->Npoints == 0 || ml_gg->row_ptr == NULL || ml_gg->col_ptr == NULL)
   {
      printf("ML_GGraph_Gen_Restrictor : graph not given. \n");
      exit(1);
   }
   if (ml_gg->vertex_state == NULL)
   {
      printf("ML_GGraph_Gen_Restrictor : coarsening not done. \n");
      exit(1);
   }
*/

   /* ------------------------------------------------------------- */
   /* ------------------------------------------------------------- */

/*
   Nfine = ml_gg->Npoints;
   Ncoarse = ml_gg->Nselected;
   vertex_state = ml_gg->vertex_state;
   row_ptr = ml_gg->row_ptr;
   col_ptr = ml_gg->col_ptr;
   new_row_ptr = (int *) ML_allocate( Ncoarse * sizeof(int) );
   total_length = 0;
 
   for ( fine_index = 0; fine_index < Nfine; fine_index++ ) {
      if ( vertex_state[fine_index] == 'S' ) {
         total_length += ( row_ptr[fine_index+1] - row_ptr[fine_index] );
*/

