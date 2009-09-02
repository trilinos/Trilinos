/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators using METIS to create the      */
/* aggregate.                                                                */
/* ************************************************************************* */
/* Author        : Michael Gee (SNL)                                         */
/* Date          : November 2004                                             */
/* ************************************************************************* */
/* Local Function :                                                          */
/*    ML_DecomposeGraph_with_VBMETIS                                         */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_agg_METIS.h"
#include "ml_agg_VBMETIS.h"
#include "ml_viz_stats.h"
#include "ml_agg_info.h"

#if defined(OUTPUT_AGGREGATES) || defined(INPUT_AGGREGATES) || (ML_AGGR_INAGGR) || (ML_AGGR_OUTAGGR) || (ML_AGGR_MARKINAGGR)
#ifndef MAXWELL
#ifndef ALEGRA
        extern int *update_index, *update, *extern_index, *external;
#endif
#else
        extern int *reordered_glob_nodes, *global_node_inds,
                   *reordered_node_externs, *global_node_externs;
#endif /*ifdef MAXWELL */
extern int ML_gpartialsum_int(int val, ML_Comm *comm);
#endif

/* ********************************************************************** */
/*! \file ml_agg_VBMETIS.c: functions to decompose a local graph using METIS

   \note metis.h is required to properly define idxtype, and to declare the     
   function `METIS_PartGraphRecursive' or `METIS_PartGraphKway'.          
   By default, idxtype is defined as int, so in principle you can         
   compile also without including metis.h. However, to be sure            
   that use has specified include and library, I require ML_METIS         
   to be defined at compilation time.                                     

   \warning without metis this function compiles, but does not do any  job.  
 ********************************************************************** */

#ifdef HAVE_ML_METIS
#ifdef __cplusplus
extern "C" {
#endif
#include "metis.h"
#ifdef __cplusplus
}
#endif
#else
#define idxtype int
#endif

#define ML_AGGR_VBMETIS  10  /* complete list of aggregation defines is in ml_aggregate.c */

static int ML_DecomposeGraph_with_VBMETIS( ML_Operator *Amatrix,
					   int N_parts,
    					   int graph_decomposition[],
					   char bdry_nodes[],
					   int local_or_global,
					   int offsets[],
					   int reorder_flag,
					   int current_level, int *total_nz);
#define OPTIMAL_VALUE (27*27)
#define EPS9 (1.0e-09)

#ifdef EXTREME_DEBUGGING
static int MyPID_ = 0;
void set_print(int MyPID ) 
{
  MyPID_ = MyPID;
}
#include <stdarg.h>
void print(char * str, ...) 
{
  
  va_list ArgList;
  va_start(ArgList, str);
  printf("===%d=== ", MyPID_);
  vprintf(str,ArgList);
  return;
}
#endif

#ifdef LATER
/* ======================================================================== */
/*!
 \brief find position of the maximum element in an integer vector

*/
/* ------------------------------------------------------------------------ */
static int find_max(int length, int vector[] ) 
{
  int max = -1;
  int pos = -1, i;
  
  for( i=0 ; i<length ; i++ )
    if( vector[i]>max ) {
      max = vector[i];
      pos = i;
    }
  return pos;
} /* find_max */

/* ======================================================================== */
/*!
 \brief find \c key in an integer vector

 \c I cannot use AZ_find_index because my vectors are not necessarely
 ordered in ascending order.

*/
/* ------------------------------------------------------------------------ */
static int find_index( int key, int list[], int N )
{
  int i;
  for( i=0 ; i<N ; i++ )
    if( list[i] == key )
      return i;
  return -1;
} /* find_index */
#endif


/* ------------------------------------------------------------------------ */
/*!
 \brief sets the coarsening scheme to variable block METIS 
 <pre>
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param ag                 ML_Aggregate* (input/output)  The aggregate object
 
 \sa ML_Aggregate_CoarsenVBMETIS ML_Aggregate_Set_CoarsenScheme_VBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS
*/
/* ------------------------------------------------------------------------ */
int ML_Aggregate_Set_CoarsenScheme_VBMETIS( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Set_CoarsenScheme_METIS : wrong object. \n");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      exit(EXIT_FAILURE);
   }
   ag->coarsen_scheme = ML_AGGR_VBMETIS;
   return 0;
}
/* ------------------------------------------------------------------------ */
/*!
 \brief function to pass variable block data to ML
 <pre>
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param ag                 ML_Aggregate* (input/output)  The aggregate object
 \param level              const int     (input)         the level to set the block data for
 \param N_levels           const int     (input)         maximum number of levels
 \param nblocks            const int     (input)         number of variable blocks
 \param blocks             const int *   (input)         row i belongs to block blocks[i]
 \param block_pde          const int *   (input)         row i is the block_pde[i] pde equation
 \param block_dim          const int     (input)         dimension of blocks and block_pde

 \warning block numbers in blocks have to be ascending and c-style numbered
 
 \sa ML_Aggregate_CoarsenVBMETIS ML_Aggregate_Set_CoarsenScheme_VBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS
*/
/* ------------------------------------------------------------------------ */
int ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS( ML_Aggregate *ag, 
                                                    const int level,
                                                    const int N_levels, 
                                                    const int nblocks, 
                                                    const int *blocks, 
                                                    const int *block_pde,
                                                    const int block_dim)
{
   struct aggr_vblock* temp = NULL;
   int i,nbyte;

   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      fprintf(stderr,"ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS : wrong object. \n");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   if (nblocks<=0)
   {
      fprintf(stderr,"ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS: number of blocks <= 0");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   if (blocks==NULL || block_pde==NULL)
   {
      fprintf(stderr,"ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS: no blocks supplied");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }

   if (ag->vblock_data==NULL)
   {
      nbyte = N_levels;
      if (nbyte==0)
      {
         fprintf(stderr,"*ML*ERR** Apply ML_Create() prior to \n"
                        "*ML*ERR** ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS\n"
                        "%s:%d\n",__FILE__,__LINE__);
         fflush(stderr); exit(EXIT_FAILURE);
      }
      temp = (struct aggr_vblock*)ML_allocate(nbyte*sizeof(struct aggr_vblock));
      if (temp==NULL)
      {
         fprintf(stderr,"ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS: not enough space\n");
         fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
         fflush(stderr); exit(EXIT_FAILURE);
      }   

      for (i=0; i<N_levels; i++)
      {
      temp[i].nblocks         = 0;
      temp[i].block_dim       = 0;
      temp[i].blocks          = NULL;
      temp[i].block_pde       = NULL;
      temp[i].old_invec_leng  = 0;
      temp[i].old_outvec_leng = 0;
      }

   }
   else
      temp = (struct aggr_vblock*)(ag->vblock_data);

   if (level < 0 || level >= N_levels)
   {
     fprintf(stderr,"*ML*ERR** given level %d is out of range in \n*ML*ERR** ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS\n%s:%d\n",
	     level, __FILE__, __LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   
   temp[level].nblocks         = nblocks;
   temp[level].block_dim       = block_dim;
   temp[level].old_invec_leng  = 0;
   temp[level].old_outvec_leng = 0;

   if (temp[level].blocks==NULL)
   {
      temp[level].blocks = (int*)ML_allocate(block_dim*sizeof(int));
   }
   else
   {
      ML_free(temp[level].blocks);
      temp[level].blocks = (int*)ML_allocate(block_dim*sizeof(int));
   }
   if (temp[level].blocks==NULL)
   {
      fprintf(stderr,"*ML*ERR** not enough memory to allocate blocks in \n*ML*ERR** ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS\n%s:%d\n",
	      __FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   
   if (temp[level].block_pde==NULL)
   {
      temp[level].block_pde = (int*)ML_allocate(block_dim*sizeof(int));
   }
   else
   {
      ML_free(temp[level].block_pde);
      temp[level].block_pde = (int*)ML_allocate(block_dim*sizeof(int));
   }
   if (temp[level].block_pde==NULL)
   {
      fprintf(stderr,"*ML*ERR** not enough memory to allocate block_pde in \n*ML*ERR** ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS\n%s:%d\n",
	      __FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   
   for (i=0; i<block_dim; i++)
   {
      temp[level].blocks[i]          = blocks[i];
      temp[level].block_pde[i]       = block_pde[i];
   }
   
   
   ag->vblock_data = (void*)(temp);

   return 0;
}


/* ------------------------------------------------------------------------ */
/*!
 \brief Modify matrix so that it uses a getrow wrapper
 <pre>
    
    Modify matrix so that it uses a getrow wrapper that will effectively 
    drop small values and will collapse several rows into a block row.   
    this routine is derived from ML_Operator_AmalgamateAndDropWeak_VBlocks
    but  can handle variable block sizes                                 
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param Amat               ML_Operator * (input/output)  The matrix
 \param block_size         int           (input)         the MAXIMUM blocks size
 \param drop_tolerance     double        (input)         drop tolerance
 \param nblocks            const int     (input)         number of variable blocks
 \param blocks             const int *   (input)         row i belongs to block blocks[i]

 \warning block numbers in blocks have to be ascending and c-style numbered
 
 \sa ML_Aggregate_CoarsenVBMETIS ML_Aggregate_Set_CoarsenScheme_VBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS ML_Operator_UnAmalgamateAndDropWeak_Vblocks
*/
/* ------------------------------------------------------------------------ */
int ML_Operator_AmalgamateAndDropWeak_VBlocks(ML_Operator *Amat, int block_size, 
                                double drop_tolerance, int nblocks, int* blocks)
{
   struct amalg_drop  *new_data;
   struct aggr_vblock *temp_data=NULL;
   int Nneigh, *neighbors, sendleng, rcvleng, *newsend, *newrcv, i, j, k;
   int sendcount, rcvcount, temp, row_length;
   int allocated, *bindx, Nghost, Nrows, block_count, t2, current;
   double *val, *scaled_diag, *dtemp;
   ML_Comm *comm;

   /* create a new widget to hold the amalgamation and drop information */

  comm = Amat->comm;

  if (nblocks<=0 || blocks==NULL)
  {
        fprintf(stderr,"ML_Operator_AmalgamateAndDropWeak_VBlocks: no block information\n");
        fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
        fflush(stderr); exit(EXIT_FAILURE);
  }
  if ( (block_size > 1) || (drop_tolerance >= 0.0)) {
     new_data = (struct amalg_drop *) ML_allocate( sizeof(struct amalg_drop) );
     if (new_data == NULL) {
        fprintf(stderr,"ML_Operator_AmalgamateAndDropWeak_VBlocks: out of space\n");
        fflush(stderr); exit(EXIT_FAILURE);
     }
     Nrows                      = Amat->getrow->Nrows;
     new_data->original_data    = Amat->data;
     new_data->original_getrow  = Amat->getrow;
     new_data->scaled_diag      = NULL;
     new_data->block_size       = block_size;
     new_data->drop_tolerance   = drop_tolerance;
     new_data->Amat             = Amat;
     temp_data = (struct aggr_vblock*)ML_allocate( sizeof(struct aggr_vblock) );
     /* to do ML_Operator_UnAmalgamateAndDropWeak_VBlocks later, we have to*/
     /* store original Amat->invec_leng and Amat->outvec_leng              */
     if (temp_data == NULL) {
        fprintf(stderr,"ML_Operator_AmalgamateAndDropWeak_VBlocks: out of space\n");
        fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
        fflush(stderr); exit(EXIT_FAILURE);}
     temp_data->nblocks         = nblocks;
     temp_data->blocks          = blocks;
     temp_data->old_invec_leng  = Amat->invec_leng;
     temp_data->old_outvec_leng = Amat->outvec_leng;
     new_data->vblock_data      = (void*)temp_data;
     
     
     /* figure out the block indices (need communication for ghost points) */
     /* and store these in new_data->blk_inds[]                            */


     i = Amat->invec_leng + 1;
     if (Amat->getrow->pre_comm != NULL) {
        i += Amat->getrow->pre_comm->total_rcv_length;
     }
     new_data->blk_inds   = (int    *) ML_allocate(sizeof(int)* i );
     dtemp                = (double *) ML_allocate(sizeof(double)* i );
     if (dtemp == NULL) 
        pr_error("ML_Operator_AmalgamateAndDropWeak_VBlocks: out of space\n");
                                        
     for (i = 0; i < Amat->invec_leng; i++)
        dtemp[i] = blocks[i];
        
     if (Amat->getrow->pre_comm != NULL) {
       ML_exchange_bdry(dtemp,Amat->getrow->pre_comm, Amat->invec_leng,
                        comm, ML_OVERWRITE,NULL);
     }

/*
renumber the blocks consecutively starting from 0 on each proc 
including the ghost blocks
*/
     block_count = 0;
     for (i=0; i<Amat->invec_leng; i++)
     {
        current = (int)dtemp[i];
        if (current>=0)
        {
           new_data->blk_inds[i] = block_count;
           for (k=i; k<Amat->invec_leng; k++)
           {
              t2 = (int)dtemp[k];
              if (current==t2)
              {
                 dtemp[k]=-1.;
                 new_data->blk_inds[k] = block_count;
              }
           }
           block_count++;
        }
     }
    
     Nneigh    = ML_CommInfoOP_Get_Nneighbors(Amat->getrow->pre_comm);
     neighbors = ML_CommInfoOP_Get_neighbors(Amat->getrow->pre_comm);

     for (i = 0; i < Nneigh; i++) {
       rcvleng = ML_CommInfoOP_Get_Nrcvlist(Amat->getrow->pre_comm,
                                                neighbors[i]);
       newrcv = ML_CommInfoOP_Get_rcvlist(Amat->getrow->pre_comm, 
                                              neighbors[i]);
       
       for (j = 0; j < rcvleng; j++) {
          current = (int) dtemp[ newrcv[j] ];
          if (current >= 0) {
             new_data->blk_inds[newrcv[j]] = block_count;
             for (k = j; k < rcvleng; k++) {
                t2 = (int) dtemp[ newrcv[k] ];
                if (current == t2) {
                   dtemp[ newrcv[k] ] = -1.;
                   new_data->blk_inds[newrcv[k]] = block_count;
                }
             }
             block_count++;
          }
       }
       ML_free(newrcv);
    }
    ML_free(dtemp);


     /* we need to get the matrix diagonal, scale it by drop_tolerance, */
     /* and store it */

     if ( drop_tolerance >= 0.0) {

        Nghost = 0;
        for (i = 0; i < Nneigh; i++) {
           rcvleng = ML_CommInfoOP_Get_Nrcvlist(Amat->getrow->pre_comm,
                                                neighbors[i]);
           newrcv = ML_CommInfoOP_Get_rcvlist(Amat->getrow->pre_comm, 
                                              neighbors[i]);
           for (j = 0; j < rcvleng; j++) {
              if (newrcv[j] > Nghost + Nrows - 1)
                 Nghost = newrcv[j] - Nrows + 1;
           }
           ML_free(newrcv);
        }
        ML_free(neighbors);

        allocated = 100;
        scaled_diag = (double *) ML_allocate((Nrows+Nghost)*sizeof(double));
        bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
        val   = (double *)  ML_allocate( allocated*sizeof(double));
        if (val == NULL) {
           fprintf(stderr,"ML_Operator_AmalgamateAndDropWeak_VBlocks: out of space\n");
           fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
           fflush(stderr); exit(EXIT_FAILURE);
        }

        for (i = 0 ; i < Nrows; i++) {
           ML_get_matrix_row(Amat,1,&i,&allocated,&bindx,&val,&row_length,0);
           for (j = 0; j < row_length; j++) 
              if (bindx[j] == i) break;

           scaled_diag[i] = 0.0;
           if (j != row_length)  scaled_diag[i] = val[j];
           scaled_diag[i] *= drop_tolerance;
           if (scaled_diag[i] < 0.) scaled_diag[i] = -scaled_diag[i];
        }
        ML_free(val);
        ML_free(bindx);
      
        if ( Amat->getrow->pre_comm != NULL )
           ML_exchange_bdry(scaled_diag,Amat->getrow->pre_comm,Nrows, comm, 
                            ML_OVERWRITE,NULL);

        new_data->scaled_diag = scaled_diag;
     }

     /* We need to create a new getrow structure */
     /* containing a getrow wrapper              */


     Amat->num_PDEs     = 1;

     Amat->invec_leng   = nblocks;
     Amat->outvec_leng  = nblocks;

     Amat->data         = new_data;
     ML_memory_alloc((void**)&(Amat->getrow),sizeof(ML_GetrowFunc),"OF2");
     Amat->getrow->ML_id            = ML_EMPTY;
     Amat->getrow->Nrows            = 0;
     Amat->getrow->pre_comm         = NULL;
     Amat->getrow->post_comm        = NULL;
     Amat->getrow->func_ptr         = NULL;
     Amat->getrow->data             = NULL;
     Amat->getrow->use_loc_glob_map = ML_NO;
     Amat->getrow->loc_glob_map     = NULL;
     Amat->getrow->row_map          = NULL;

     ML_Operator_Set_Getrow(Amat,nblocks,ML_amalg_drop_getrow_VBlocks);

     /* amalgamation needs a new communication structure. Let's create a new */
     /* communication object and modify it if we are doing amalgmation.     */                   
     ML_CommInfoOP_Clone( &(Amat->getrow->pre_comm),  
                          new_data->original_getrow->pre_comm);

     if (block_size > 1) {
        Nneigh    = ML_CommInfoOP_Get_Nneighbors(Amat->getrow->pre_comm);
        neighbors = ML_CommInfoOP_Get_neighbors(Amat->getrow->pre_comm);


        for (i = 0; i < Nneigh; i++) {
           sendleng = ML_CommInfoOP_Get_Nsendlist(Amat->getrow->pre_comm, 
                                                  neighbors[i]);
           newsend = ML_CommInfoOP_Get_sendlist(Amat->getrow->pre_comm, 
                                                neighbors[i]);
           sendcount = 0;
           for (j = 0 ; j < sendleng; j++) {
              temp = new_data->blk_inds[newsend[j]];

              /* search to see if it is already in the list */
              for (k = 0; k < sendcount; k++) 
                 if ( newsend[k] == temp) break;

              if (k == sendcount) newsend[sendcount++] = temp;
           }
           
           
           rcvleng = ML_CommInfoOP_Get_Nrcvlist(Amat->getrow->pre_comm, 
                                                neighbors[i]);
           newrcv = ML_CommInfoOP_Get_rcvlist(Amat->getrow->pre_comm, neighbors[i]);
           rcvcount = 0;
           for (j = 0 ; j < rcvleng; j++) {
              temp = new_data->blk_inds[newrcv[j]];

              /* search to see if it is already in the list */
              for (k = 0; k < rcvcount; k++) 
                 if ( newrcv[k] == temp) break;

              if (k == rcvcount) newrcv[rcvcount++] = temp;
           }
           ML_CommInfoOP_Set_exch_info(Amat->getrow->pre_comm, neighbors[i],
                      rcvcount, newrcv,sendcount, newsend);
           ML_free(newrcv); 
           ML_free(newsend); 
        }
        if (neighbors != NULL) ML_free(neighbors);
     }
  }
  return 0;
}

/* ------------------------------------------------------------------------ */
/*!
 \brief Restores a matrix that has been modified via ML_Operator_AmalgamateAndDropWeak_Vblocks
 <pre>
    
    Restores a matrix that has been modified via                         
    ML_Operator_AmalgamateAndDropWeak_Vblocks back to its original form
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param Amat               ML_Operator * (input/output)  The matrix
 \param block_size         int           (input)         the MAXIMUM blocks size
 \param drop_tolerance     double        (input)         drop tolerance
 
 \sa ML_Aggregate_CoarsenVBMETIS ML_Aggregate_Set_CoarsenScheme_VBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS ML_Operator_AmalgamateAndDropWeak_VBlocks
*/
/* ------------------------------------------------------------------------ */

int ML_Operator_UnAmalgamateAndDropWeak_Vblocks(ML_Operator *Amat, int block_size,
	double drop_tolerance)
{
   struct amalg_drop  *temp;
   struct aggr_vblock *temp_data;
 
   if ( (block_size > 1) || (drop_tolerance >= 0.0)) {
      temp      = (struct amalg_drop *) Amat->data;
      temp_data = (struct aggr_vblock*)(temp->vblock_data);
      
      ML_CommInfoOP_Destroy(&(Amat->getrow->pre_comm));
      ML_memory_free((void**)&(Amat->getrow));
      Amat->data         = temp->original_data;
      Amat->getrow       = temp->original_getrow;
      Amat->invec_leng   = temp_data->old_invec_leng;
      Amat->outvec_leng  = temp_data->old_outvec_leng;
      Amat->num_PDEs     = temp->block_size;
      if (temp->blk_inds != NULL) ML_free(temp->blk_inds);
      if (temp->scaled_diag != NULL) ML_free(temp->scaled_diag);
      ML_free(temp_data);
      ML_free(temp);
   }
   return 0;
}


/* ------------------------------------------------------------------------ */
/*!
 \brief Getrow function for amalgamteded variable block row matrix
 <pre>
    
    Getrow function that is used to drop matrix elements and to collapse 
    several rows into a block. It is assumed that                        
    ML_Operator_AmalgamateAndDropWeak_VBlocks() was previously called to 
    properly set up the data structure (data).                           
    This function is derived from ML_amalg_drop_getrow                   
    of which it is the variable block version                            
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param data               ML_Operator * (input)         The matrix
 \param N_requested_rows   int           (input)         number of rows to receive
 \param requested_rows     int[]         (input)         row number requested
 \param allocated_space    int           (input)         length of columns & values
 \param columns            int[]         (output)        column indizes of requested row
 \param values             double[]      (output)        values of requested row
 \param row_lengths        int[]         (output)        row_lengths[i] is length of requested row requested_rows[i]
 
 \warning currently implemented only N_requested_rows=1
 
 \sa ML_Operator_AmalgamateAndDropWeak_VBlocks ML_Operator_UnAmalgamateAndDropWeak_Vblocks
     
*/
/* ------------------------------------------------------------------------ */

int ML_amalg_drop_getrow_VBlocks(ML_Operator *data, int N_requested_rows, 
   int requested_rows[],int allocated_space, int columns[], double values[], 
   int row_lengths[])
{
   struct amalg_drop  *temp;
   struct aggr_vblock *block_data=NULL;
   int    block_size, row, size, i, j, k, tcol, count;
   int    *tcolumns, tallocated_space;
   double *tvalues, *scaled_diag;
   int offset, status = 1;
   struct ML_GetrowFunc_Struct *amalg_getrow;
   ML_Operator *Amat;
   int nblocks,*blocks,index1,index2,nrows,bsize;
 
   if (N_requested_rows > 1) {
      fprintf(stderr,"ML_amalg_drop_getrow_VBlocks: Not implemented for > 1 row at a time\n");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   Amat = (ML_Operator *) data;
   temp = (struct amalg_drop *) ML_Get_MyGetrowData(Amat);
   Amat = temp->Amat;
   /* get the blocked data */
   block_data   = (struct aggr_vblock*)(temp->vblock_data);
   block_size   = temp->block_size;
   amalg_getrow = Amat->getrow;
   scaled_diag  = temp->scaled_diag;
   nblocks      = block_data->nblocks;
   blocks       = temp->blk_inds;
   /* set the matrix back to original data and getrow */
   Amat->data         = temp->original_data;
   Amat->getrow       = temp->original_getrow;
   Amat->invec_leng   = block_data->old_invec_leng;
   Amat->outvec_leng  = block_data->old_outvec_leng;
   nrows              = Amat->getrow->Nrows;

   tallocated_space = allocated_space*block_size*block_size + 1;
   tcolumns     = (int    *) ML_allocate(sizeof(int)*tallocated_space);
   tvalues      = (double *) ML_allocate(sizeof(double)*tallocated_space);
   if (tvalues == NULL) {
      if (tcolumns != NULL) ML_free(tcolumns);
      Amat->data         = temp;
      Amat->getrow       = amalg_getrow;
      Amat->invec_leng   = nblocks;
      Amat->outvec_leng  = nblocks;
      return(0);
   }
   
   /* determine start rows in block and blocksize */
   index1 = ML_find_index(requested_rows[0],blocks,nrows);
   if (index1==-1)
   {
      ML_free(tvalues); ML_free(tcolumns);
      Amat->data         = temp;
      Amat->getrow       = amalg_getrow;
      Amat->invec_leng   = nblocks;
      Amat->outvec_leng  = nblocks;
      return(0);
   }
   if (requested_rows[0] != nblocks-1)
   {
      index2 = ML_find_index(requested_rows[0]+1,blocks,nrows);
      if (index2==-1)
      {
         ML_free(tvalues); ML_free(tcolumns);
         Amat->data         = temp;
         Amat->getrow       = amalg_getrow;
         Amat->invec_leng   = nblocks;
         Amat->outvec_leng  = nblocks;
         return(0);
      }
      bsize = index2 - index1;
   }
   else 
      bsize = nrows - index1;
      
   if (bsize==0)
   {
      fprintf(stderr,"requested a variable block row of size 0\n");
      fprintf(stderr,"check the vblock input to ML and/or the code\n");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   offset = 0;
   for (i = 0; i < bsize; i++) {
      row = index1+i;
      status = ML_Operator_Getrow(Amat, N_requested_rows, &row, 
                                  tallocated_space, &(tcolumns[offset]), 
				  &(tvalues[offset]), &size );
      if (status == 0) {
         ML_free(tvalues); ML_free(tcolumns);
         Amat->data         = temp;
         Amat->getrow       = amalg_getrow;
         Amat->invec_leng   = nblocks;
         Amat->outvec_leng  = nblocks;
         return(status);
      }
      
      
      
      if (scaled_diag != NULL) {
         count = 0;
         for (j = offset; j < offset + size; j++) {
            tcol = tcolumns[j];
            if (tvalues[j] != 0.0) {
              if (tvalues[j]*tvalues[j] >= scaled_diag[row]*scaled_diag[tcol]) {
                 tcolumns[offset+count]  = tcolumns[j];
                 tvalues[offset+count++] = tvalues[j];
              }
            }
         }
         size = count;
      }
      tallocated_space -= size;
      offset += size;
   }

   row_lengths[0] = 0;

   for (j = 0; j < offset; j++) {
      tcol = temp->blk_inds[tcolumns[j]];
      for (k = 0; k < row_lengths[0]; k++) 
         if (tcol == columns[k]) break;

      if (k == row_lengths[0]) {
         if ( allocated_space == row_lengths[0]) {
            ML_free(tvalues); ML_free(tcolumns);
            Amat->data         = temp;
            Amat->getrow       = amalg_getrow;
            Amat->invec_leng   = nblocks;
            Amat->outvec_leng  = nblocks;
            return(0);
         }
         values[row_lengths[0]] = 1.;
         columns[row_lengths[0]++] = tcol;
      }
   }
   Amat->data         = temp;
   Amat->getrow       = amalg_getrow;
   Amat->invec_leng   = nblocks;
   Amat->outvec_leng  = nblocks;
   ML_free(tvalues); ML_free(tcolumns);
   return(status);
}


/* ------------------------------------------------------------------------ */
/*!
 \brief function to obtain variable block data of a certain level
 <pre>
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param ag                 const ML_Aggregate* (input)   The aggregate object
 \param level              const int           (input)   the level to get the block data from
 \param N_levels           const int           (input)   maximum number of levels
 \param nblocks            int*                (output)  number of variable blocks
 \param blocks             int**               (output)  on output, row i is in blocks[i]
 \param block_pde          int**               (output)  on output, row i is the block_pde[i] pde equation
 
 \warning block-numbers have to be ascending and C-style numbering
 
 \sa ML_Aggregate_CoarsenVBMETIS ML_Aggregate_Set_CoarsenScheme_VBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS
*/
/* ------------------------------------------------------------------------ */
int ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS( const ML_Aggregate *ag, 
                                                    const int level,
                                                    const int N_levels, 
                                                          int *nblocks, 
                                                          int **blocks, 
                                                          int **block_pde)
{
   struct aggr_vblock* temp = NULL;
   int i,start;

   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      fprintf(stderr,"ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS : wrong object. \n");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   /* block data has not been set before, print warning */
   if (ag->vblock_data==NULL)
   {
      fprintf(stdout,"*ML*WRN** ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS:\n");
      fprintf(stdout,"*ML*WRN** no block data set in ML_Aggregate *ag, use\n");
      fprintf(stdout,"*ML*WRN** ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS to\n");
      fprintf(stdout,"*ML*WRN** set block data.\n");
      fflush(stdout);
      *nblocks   = 0;
      *blocks    = NULL;
      *block_pde = NULL;
      return 0;
   }
   /* check range of input data */
   if (level<0 || level >= N_levels)
   {
      fprintf(stderr,"*ML*ERR** ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS:\n*ML*ERR** level %d out of range ( 0 - %d )\n%s:%d\n",
	      level,N_levels,__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   if (level >= ag->max_levels)
   {
      fprintf(stdout,"*ML*WRN** ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS:\n");
      fprintf(stdout,"*ML*WRN** requested level=%d >= ag->max_levels=%d\n",level,ag->max_levels);
      fflush(stdout);
      *nblocks   = 0;
      *blocks    = NULL;
      *block_pde = NULL;
      return 0;
   }
   temp = (struct aggr_vblock*)(ag->vblock_data);
   temp = &(temp[level]);
   if (temp->blocks==NULL || temp->block_pde==NULL)
   {
      fprintf(stdout,"*ML*WRN** ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS:\n");
      fprintf(stdout,"*ML*WRN** no blocks on level %d\n",level);
      fflush(stdout);
      *nblocks   = 0;
      *blocks    = NULL;
      *block_pde = NULL;
      return 0;
   }
   else
   {
      *nblocks     = temp->nblocks;
      (*blocks)    = (int*)ML_allocate(temp->block_dim*sizeof(int));
      (*block_pde) = (int*)ML_allocate(temp->block_dim*sizeof(int));
      if ((*block_pde)==NULL)
      {
         fprintf(stderr,"*ML*ERR** ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS:\n"
                        "*ML*ERR** not enough space\n"
                        "%s:%d\n",__FILE__,__LINE__);
         fflush(stderr); exit(EXIT_FAILURE);
      }
      start = temp->blocks[0];
      for (i=0; i<temp->block_dim; i++)
      {
         (*blocks)[i] = temp->blocks[i]-start;
         (*block_pde)[i] = temp->block_pde[i];
      }
   }

   return 0;
}
/* ------------------------------------------------------------------------ */
/*!
 \brief destroy variable block data of a certain level
 <pre>
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param ag                 const ML_Aggregate* (input)   The aggregate object
 \param level              const int           (input)   the level to destroy the block data
 
 \sa ML_Aggregate_CoarsenVBMETIS ML_Aggregate_Set_CoarsenScheme_VBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS
*/
/* ------------------------------------------------------------------------ */
int ML_Aggregate_Destroy_Vblocks_CoarsenScheme_VBMETIS( const ML_Aggregate *ag, 
                                                        const int level)
{
   struct aggr_vblock* temp = NULL;

   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      fprintf(stderr,"ML_Aggregate_Destroy_Vblocks_CoarsenScheme_VBMETIS : wrong object. \n");
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stderr); exit(EXIT_FAILURE);
   }
   if (ag->vblock_data==NULL)
      return 0;
   /* check range of input data */
   if (level < 0 || level >= ag->max_levels)
   {
      fprintf(stdout,"*ML*WRN** ML_Aggregate_Destroy_Vblocks_CoarsenScheme_VBMETIS:\n");
      fprintf(stdout,"*ML*WRN** requested level=%d >= ag->max_levels=%d\n",level,ag->max_levels);
      fprintf(stderr,"file %s, line %d\n",__FILE__,__LINE__);
      fflush(stdout);
      return 0;
   }
   temp = (struct aggr_vblock*)(ag->vblock_data);
   temp = &(temp[level]);
   if (temp->blocks==NULL || temp->block_pde==NULL)
      return 0;
   else
   {
      ML_free(temp->blocks);
      ML_free(temp->block_pde);
      temp->nblocks = 0;
      temp->block_dim = 0;
      temp->blocks = NULL;
      temp->block_pde = NULL;
      temp->old_invec_leng = 0;
      temp->old_outvec_leng = 0;
   }

   return 0;
}
/* ------------------------------------------------------------------------ */
/*!
 \brief create non-smoothed aggregates using METIS
 <pre>
    In order to use this function, the user has to define the number of aggregate 
    (or the # of nodes in each aggregate) using the functions 
    \c ML_Aggregate_Set_LocalNumber or \c ML_Aggregate_Set . 
    
    Difference between ML_Aggregate_CoarsenMETIS and this routine is that 
    this one respects variable sized blocks in the aggregation process. 
    It guarantees that all rows of a variable block belong to the same 
    aggregate.
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
 </pre>
 \param ml_ag                ML_Aggregate* (input/output)   The aggregate object
 \param Amatrix              ML_Operator*  (input)          the operator of the current level
 \param Pmatrix              ML_Operator** (output)         the prolongation operator
 \param comm                 ML_Comm *     (input)          parallel layout
 
 \sa  ML_Aggregate_Set_CoarsenScheme_VBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS ML_Aggregate_Destroy_Vblocks_CoarsenScheme_VBMETIS

 \note this function is derived from ML_Aggregate_CoarsenMETIS

*/
/* ------------------------------------------------------------------------ */
int ML_Aggregate_CoarsenVBMETIS( ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
			         ML_Operator **Pmatrix, ML_Comm *comm)
{
  unsigned int nbytes, length;
   int     i, j,  k, Nrows, exp_Nrows;
   int     diff_level;
   int     aggr_count, index = 0, mypid, num_PDE_eqns;
   int     *aggr_index = NULL, nullspace_dim;
   int     Ncoarse, count;
   int     *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     exp_Ncoarse;
   int     *aggr_cnt_array = NULL;
   int     level, index3, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info;
   double  *new_val = NULL, epsilon;
   double  *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null = NULL;
   ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
   struct  ML_CSR_MSRdata *csr_data;
   int     total_nz = 0;
   char str[80];

   struct   aggr_vblock *vb_data;
   int      nblocks,*blocks=NULL,oldNrows,*block_pde=NULL;
   int      new_nrows, new_nblocks, *new_blocks=NULL, *new_block_pde=NULL;
   int      index1,index2,printflag=0;
      
   int reorder_flag;
   /*   int kk, old_upper, nnzs, count2, newptr; */
#ifdef CLEAN_DEBUG
   char tlabel[80];
#endif

#ifdef DDEBUG
   int curagg,myagg,*good,*bad, kk;
#endif

#if defined(OUTPUT_AGGREGATES) || defined(DDEBUG) || defined(INPUT_AGGREGATES) || (ML_AGGR_INAGGR) || (ML_AGGR_OUTAGGR) || (ML_AGGR_MARKINAGGR)
   FILE *fp;
   char fname[80];
   static int level_count = 0;
   double *d2temp,*dtemp;
   int agg_offset, vertex_offset;
#endif
   int * graph_decomposition = NULL;
   ML_Aggregate_Viz_Stats * aggr_viz_and_stats;
   ML_Aggregate_Options * aggr_options;
   int mod, Nprocs;
   int optimal_value;
   char * unamalg_bdry = NULL;
 
#ifdef EXTREME_DEBUGGING
 set_print(comm->ML_mypid);
#endif
 
 /* ------------------- execution begins --------------------------------- */

 sprintf( str, "VBMETIS (level %d) :", ml_ag->cur_level );

 /* ============================================================= */
 /* get the machine information and matrix references             */
 /* ============================================================= */
 
 mypid                   = comm->ML_mypid;
 Nprocs                  = comm->ML_nprocs;
 epsilon                 = ml_ag->threshold;
 num_PDE_eqns            = ml_ag->num_PDE_eqns;
 nullspace_dim           = ml_ag->nullspace_dim;
 nullspace_vect          = ml_ag->nullspace_vect;
 vb_data                 = (struct aggr_vblock*)(ml_ag->vblock_data);
 if (vb_data==NULL)
 {
    fprintf( stderr,
             "*ML*ERR* no blocks for variable block coarsening\n"
             "*ML*ERR* use ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS to\n"
             "*ML*ERR* set user supplied variable block information\n"
             "*ML*ERR* (file %s, line %d)\n",__FILE__,__LINE__ );
    fflush(stderr); exit( EXIT_FAILURE );
 }
 
 nblocks                 = vb_data[ml_ag->cur_level].nblocks;
 blocks                  = vb_data[ml_ag->cur_level].blocks;
 block_pde               = vb_data[ml_ag->cur_level].block_pde;
 Nrows                   = Amatrix->outvec_leng;
 
 if (nblocks==0 || blocks==NULL || block_pde==NULL)
 {
    fprintf( stderr,
             "*ML*ERR* no blocks for variable block coarsening\n"
             "*ML*ERR* use ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS to\n"
             "*ML*ERR* set user supplied block information\n"
             "*ML*ERR* (file %s, line %d)\n",__FILE__,__LINE__ );
    fflush(stderr); exit( EXIT_FAILURE );
 }

 if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s variable blocks & num PDE eqns = %d\n",
	    str,
	    num_PDE_eqns);
 }
 
   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */

#ifdef EXTREME_DEBUGGING
   print("# rows orig = %d, # PDE eqns = %d\n", Nrows, num_PDE_eqns);
#endif
   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && 7 < ML_Get_PrintLevel())
   {
      printf("%s current eps = %e\n",
	     str,
	     epsilon);
      if( epsilon != 0.0 ) {
	fprintf( stderr,
		 "%s (note that VBMETIS may not work with dropping)\n",
		 str );
      }
      
   }
   /*
     epsilon = epsilon * epsilon;
   */
   oldNrows = Amatrix->getrow->Nrows;
   
   ML_Operator_AmalgamateAndDropWeak_VBlocks(Amatrix, num_PDE_eqns, epsilon,
                                             nblocks,blocks);
   Nrows     = nblocks;
   exp_Nrows = nblocks;
   
   /* ********************************************************************** */
   /* allocate memory for aggr_index, which will contain the METIS decomp.   */
   /* ********************************************************************** */

   nbytes = (oldNrows) * sizeof(int);

   if ( nbytes > 0 ) {
     ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
     if( aggr_index == NULL ) {
       fprintf( stderr,
		"*ML*ERR* not enough memory for %d bytes\n"
		"*ML*ERR* (file %s, line %d)\n",
		nbytes,
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
     }
   }
   else              aggr_index = NULL;

   for( i=0 ; i<oldNrows ; i++ ) aggr_index[i] = -1;
   
#ifdef EXTREME_DEBUGGING
   print("aggr_index has size %d bytes\n", nbytes);
   print("aggr_options pointer is %x\n", aggr_options);
   print("Amatrix->to->Grid->Grid pointer is  %x\n",
	 Amatrix->to->Grid->Grid );
#endif

   /* ********************************************************************** */
   /* retrive the pointer to the ML_Aggregate_Options, which contains the    */
   /* number of aggregates (or their size), as well as few options for the   */
   /* constructions.                                                         */
   /* ********************************************************************** */

   aggr_options = (ML_Aggregate_Options *)ml_ag->aggr_options;

   if( aggr_options == NULL ) {

     if( mypid == 0 && 8 < ML_Get_PrintLevel() ) {
       printf("%s Using default values\n",
	      str);
     }

     /* this value is hardwired in ml_agg_METIS.c, and can be set by      */
     /* the user with `ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate' */
     
     optimal_value = ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate();
     
     aggr_count = Nrows/optimal_value;
     if( aggr_count < 1 ) aggr_count = 1;

     reorder_flag = ML_NO;
     
   } else {

     if( aggr_options[ml_ag->cur_level].id != ML_AGGREGATE_OPTIONS_ID ) {
       fprintf( stderr,
		"*ML*ERR* `ML_Aggregate_CoarsenVBMETIS' : wrong object\n"
		"*ML*ERR* (file %s, line %d)\n",
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
     }

     /* ******************************************************************** */
     /* Retrive the user's defined choice to define the number of aggregates */
     /* For local number, it is ok.                                          */
     /* If global number of aggregates or nodes per aggregate have been      */
     /* specified, compute the local one (evenly dividing this global number)*/
     /* For those two latter cases, I suppose that the input value is the    */
     /* same on all processes (but I don't check ... )                       */
     /* ******************************************************************** */

     switch( aggr_options[ml_ag->cur_level].choice ) {

     case ML_NUM_LOCAL_AGGREGATES:

       aggr_count = aggr_options[ml_ag->cur_level].Naggregates_local;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Objective : %d local (block) aggregates (on proc 0)\n",
		 str,
		 aggr_count );
       }
       break;

     case ML_NUM_GLOBAL_AGGREGATES:
       
       aggr_count = aggr_options[ml_ag->cur_level].Naggregates_global;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Objective : %d global aggregates\n",
		 str,
		 aggr_count );
       }
       
       if( aggr_count < Nprocs ) {
	 if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	   fprintf( stderr,
		    "*ML*WRN* In CoarsenVBMETIS, %d global (block) aggregates are required,\n"
		    "*ML*WRN* but you have only %d processes. METIS requires at\n"
		    "*ML*WRN* one aggregate per process. Otherwise, you can use ParMETIS\n"
		    "*ML*WRN* as coarsen scheme. Now proceeding with 1 local (block) aggregate\n"
		    "*ML*WRN* (file %s, line %d)\n",
		    aggr_count,
		    Nprocs,
		    __FILE__,
		    __LINE__ );
	 }
	 aggr_count = 1;

       } else { 
       
	 mod = aggr_count % Nprocs;
     
	 aggr_count /= Nprocs;
	 if( mypid == 0 ) {
	   aggr_count += mod;
	 }
	 
	 if( aggr_count < 1 ) {
	   fprintf( stderr,
		    "*ML*WRN* something weird happened... Check the code !!\n"
		    "*ML*WRN* (file %s, line %d)\n",
		    __FILE__,
		    __LINE__ );
	   aggr_count = 1;
	 }
       }
       
       break;
       
     case ML_NUM_NODES_PER_AGGREGATE:
       
       aggr_count = aggr_options[ml_ag->cur_level].Nnodes_per_aggregate;

       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Objective : %d nodes per aggregate\n",
		 str,
		 aggr_count );
       }
       
       if( aggr_count >= Nrows) {

	 i = aggr_count;
	 
	 aggr_count = Nrows/OPTIMAL_VALUE;
	 if( aggr_count == 0 ) aggr_count = 1;
	 
	 if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	   fprintf( stderr,
		    "*ML*WRN* # (block) nodes per (block) aggregate (%d) > # (block) nodes (%d)\n"
		    "*ML*WRN* (on proc 0). Now proceeding with %d aggregates\n",
		    i,
		    Nrows,
		    aggr_count);
	 }
	 
       } else {

	 aggr_count = (Nrows/aggr_count);

	 if( aggr_count == 0 ) aggr_count = 1;
	 
       }

#ifdef ML_MPI
       MPI_Reduce( &aggr_count, &i, 1, MPI_INT, MPI_SUM, 0,
		   comm->USR_comm);
#else
       i = aggr_count;
#endif
       
       if ( mypid == 0 && 7 < ML_Get_PrintLevel() )  {
	 printf("%s avg %f (block) aggr/process\n",
		str,
		1.0*i/Nprocs );
       }
       
       break;
       
     } /* switch */
       
     reorder_flag = aggr_options[ml_ag->cur_level].reordering_flag;
     
   } /* if( aggr_options == NULL )*/
   
   if( aggr_count<=0 ) {
     fprintf( stderr,
	      "*ML*ERR* on proc %d, value of aggr_count not correct (%d)\n"
	      "*ML*ERR* Set this value using ML_Aggregate_Set_LocalNumber\n"
	      "*ML*ERR* or ML_Aggregate_Set_NodesPerAggr or ML_Aggregate_Set_GlobalNumber\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      mypid,
	      aggr_count,
	      __FILE__,
	      __LINE__ );
     fflush(stderr); exit( EXIT_FAILURE );
   }

   /* ********************************************************************** */
   /* to call METIS, we have to create a graph using CSR data format         */
   /* Essentially, this requires a matrix-vector product (on the Amalgamated */
   /* matrix, so with dropped elements)                                      */
   /* ********************************************************************** */

   unamalg_bdry = (char *) ML_allocate( sizeof(char) * (Nrows+1) );

   if( unamalg_bdry == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* on proc %d, not enough space for %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      mypid,
	      (int)sizeof(char) * Nrows,
	      __FILE__,
	      __LINE__ );
     fflush(stderr); exit( EXIT_FAILURE );
   }

#ifdef EXTREME_DEBUGGING
   print("# requested local aggregates = %d, # rows_METIS = %d\n",
	 aggr_count,
	 Nrows );
#endif


   aggr_count = ML_DecomposeGraph_with_VBMETIS( Amatrix,aggr_count,
					        aggr_index, unamalg_bdry,
					        ML_LOCAL_INDICES, NULL,
					        reorder_flag, ml_ag->cur_level,
					        &total_nz);


#ifdef ML_MPI
   MPI_Allreduce( &Nrows, &i, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
   MPI_Allreduce( &aggr_count, &j, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
#else
   i = Nrows;
   j = aggr_count;
#endif

   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
     printf("%s Using %d (variable block) aggregates (globally)\n",
	    str,
	    j );
     printf("%s # (vblock) aggre/ # (vblock) rows = %8.5f %% ( = %d / %d)\n",
	    str,
	    100.0*j/i,
	    j, i);
   }

   j = ML_gsum_int( aggr_count, comm );
   if ( mypid == 0 && 7 < ML_Get_PrintLevel() )  {
     printf("%s %d (variable block) aggregates (globally)\n",
	    str,
	    j );
   }   
   
   /* ********************************************************************** */
   /* I allocate room to copy aggr_index and pass this value to the user,    */
   /* who will be able to analyze and visualize this after the construction  */
   /* of the levels. This way, the only price we have to pay for stats and   */
   /* viz is essentially a little bit of memory.                             */
   /* this memory will be cleaned with the object ML_Aggregate ml_ag.        */
   /* I set the pointers using the ML_Aggregate_Info structure. This is      */
   /* allocated using ML_Aggregate_Info_Setup(ml,MaxNumLevels)               */
   /* ********************************************************************** */

   if( Amatrix->to->Grid->Grid != NULL ) {
     fprintf( stderr,
              "%s *ML*WRN* visualization of variable-block-aggregates not tested!\n%s *ML*WRN* (file %s, line %d)\n",str,str,__FILE__,__LINE__ );
     fflush(stderr);
     graph_decomposition = (int *)ML_allocate(sizeof(int)*(Nrows+1));
     if( graph_decomposition == NULL ) {
       fprintf( stderr,
		"*ML*ERR* Not enough memory for %d bytes\n"
		"*ML*ERR* (file %s, line %d)\n",
		(int)sizeof(int)*Nrows,
		__FILE__,
	        __LINE__ );
       fflush(stderr); exit( EXIT_FAILURE );
     }

     for( i=0 ; i<Nrows ; i++ ) graph_decomposition[i] = aggr_index[i];

     aggr_viz_and_stats = (ML_Aggregate_Viz_Stats *) (Amatrix->to->Grid->Grid);
     aggr_viz_and_stats->graph_decomposition = graph_decomposition;
     aggr_viz_and_stats->Nlocal = Nrows;
     aggr_viz_and_stats->Naggregates = aggr_count;
     aggr_viz_and_stats->local_or_global = ML_LOCAL_INDICES;
     aggr_viz_and_stats->is_filled = ML_YES;
     aggr_viz_and_stats->Amatrix = Amatrix;
   }

   /* ********************************************************************** */
   /* take the decomposition as created by METIS and form the aggregates     */
   /* ********************************************************************** */
   
   total_nz = ML_Comm_GsumInt( comm, total_nz);
   i = ML_Comm_GsumInt( comm, Nrows);

   if ( mypid == 0 && 7 < ML_Get_PrintLevel())
     printf("%s Total (vblock) nnz = %d ( = %5.2f/(variable block)row)\n",
	    str,
	    total_nz,1.0*total_nz/i);
   
   if ( ml_ag->operator_complexity == 0.0 ) {
      ml_ag->fine_complexity = total_nz;
      ml_ag->operator_complexity = total_nz;
   }
   else ml_ag->operator_complexity += total_nz;

   /* fix aggr_index for num_PDE_eqns > 1 and with variable blocks    */
   /* Nrows is number of variale block rows                           */
   /* aggr_index[i] is number of aggregate containing i-th vblock row */
   /* oldNrows is real number of scalar rows                          */
   /* nblocks  is number of variable block rows (== Nrows)            */
   /* blocks (length oldNrows) is blockindizes of all oldNrows rows   */
   index1 = index2 = oldNrows-1;
   for (i = Nrows - 1; i >= 0; i-- )
   {
      while (blocks[index2]==blocks[index1] && index2 >= -1)
         index2--;
      for (j=index1; j>index2; j--)
      {
         if (j<i && i != 0)
         {
            fprintf( stderr,
		     "*ML*ERR* check the code for blockrow %d\n*ML*ERR* (file %s, line %d)\n",
		     i,__FILE__,__LINE__ ); 
                      fflush(stderr); exit( EXIT_FAILURE );
         }
         aggr_index[j] = aggr_index[i];
      }
      index1 = index2;
   }

   if ( mypid == 0 && 8 < ML_Get_PrintLevel())
   {
      printf("%s Calling ML_Operator_UnAmalgamateAndDropWeak_Vblocks\n",str);
      fflush(stdout);
   }

   ML_Operator_UnAmalgamateAndDropWeak_Vblocks(Amatrix, num_PDE_eqns, epsilon);

#ifdef EXTREME_DEBUGGING
   print("After `ML_Operator_UnAmalgamateAndDropWeak_Vblocks'\n");
#endif   
   
   Nrows      = oldNrows;
   exp_Nrows  = oldNrows;

   
   /* count the size of each aggregate */

   aggr_cnt_array = (int *) ML_allocate(sizeof(int)*(aggr_count+1));
   if (aggr_cnt_array==NULL)
   {
      fprintf( stderr,
               "*ML*ERR* not enough space\n*ML*ERR* (file %s, line %d)\n",
	       __FILE__,__LINE__ ); 
                fflush(stderr); exit( EXIT_FAILURE );
   }
   for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
   for (i = 0; i < exp_Nrows; i++) {
     if (aggr_index[i] >= 0) {
       if( aggr_index[i] >= aggr_count ) {
	 fprintf( stderr,
		  "*ML*WRN* on process %d, something weird happened...\n*ML*WRN* node %d belong to aggregate %d (#aggr = %d)\n*ML*WRN* (file %s, line %d)\n",
		  comm->ML_mypid,i,aggr_index[i],aggr_count,__FILE__,__LINE__ );
       } else {
	 aggr_cnt_array[aggr_index[i]]++;
       }
     }
   }

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

#ifdef EXTREME_DEBUGGING
   print("Form tentative prolongator\n");
#endif

   Ncoarse = aggr_count; /* number of aggregates, not number of rows! */
   
   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = (Nrows+1) * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGl");
   count = aggr_count;
   for ( i = 0; i < Nrows; i++ ) 
   {
      if ( aggr_index[i] >= 0 )
      {
            ml_ag->aggr_info[level][i] = aggr_index[i];
         if (aggr_index[i] >= count) 
         {
            fprintf(stderr,
		    "**WRNG** %d : aggr_index[%d]=%d >= count=%d\n**WRNG** (file %s, line %d)\n",
		    mypid,i,aggr_index[i],count,__FILE__,__LINE__ );
            count = aggr_index[i] + 1;
         }
      }
      /* this seems to be obsolete, but I'm not sure mwgee 10/05 */
      else
      {
          fprintf(stderr,"%d : CoarsenVBMETIS error : aggr_index[%d] < 0\nsomething wrong with exlusion of Dirichlet rows\n",
                          mypid,i);
          fflush(stderr); exit(EXIT_FAILURE);
      }
   }
   ml_ag->aggr_count[level] = count; /* for relaxing boundary points */ 
   
   /* ============================================================= */
   /* set up the new operator                                       */
   /* ------------------------------------------------------------- */

#ifdef EXTREME_DEBUGGING
   print("nullspace_dim = %d\n", nullspace_dim );
   print("Nrows = %d\n", Nrows );
#endif
   
   new_Nrows   = Nrows;
   exp_Ncoarse = Nrows;
   
   for ( i = 0; i < new_Nrows; i++ ) 
   {
      if ( aggr_index[i] >= exp_Ncoarse ) 
      {
         fprintf(stderr,"*ML*WRN* index out of bound %d = %d(%d)\n",
		i, aggr_index[i], 
                exp_Ncoarse); fflush(stderr);
      }
   }
   nbytes = ( new_Nrows+1 ) * sizeof(int); 
   ML_memory_alloc((void**)&(new_ia), nbytes, "AIA");
   nbytes = ( new_Nrows+1)  * nullspace_dim * sizeof(int); 
   ML_memory_alloc((void**)&(new_ja), nbytes, "AJA");
   nbytes = ( new_Nrows+1)  * nullspace_dim * sizeof(double); 
   ML_memory_alloc((void**)&(new_val), nbytes, "AVA");
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;
   
   /* ------------------------------------------------------------- */
   /* set up the space for storing the new null space               */
   /* ------------------------------------------------------------- */

   nbytes = (Ncoarse+1) * nullspace_dim * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
   
   for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++) 
      new_null[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */
/*
   for (i = 0; i <= Nrows; i++) new_ia[i] = i * nullspace_dim;
*/
   /* trying this when a Dirichlet row is taken out */
   /* at the moment, assuming that the dimension of the */
   /* nullspace is a fixed constant                     */
   j = 0;
   new_ia[0] = 0;
   for (i = 0; i < Nrows; i++) {
      if (aggr_index[i] != -1) j += nullspace_dim;
      new_ia[i+1] = j;
   }

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"MLs");
   for (i = 0; i < aggr_count; i++) {
     nbytes = aggr_cnt_array[i]+1;
     rows_in_aggs[i] = (int *) ML_allocate(nbytes*sizeof(int));
     if (rows_in_aggs[i] == NULL)  {
       fprintf(stderr,
	       "*ML*ERR* couldn't allocate memory in CoarsenMETIS\n*ML*ERR* (file %s, line %d)\n",
	       __FILE__,__LINE__);
       fflush(stderr); exit(EXIT_FAILURE);
     }
     aggr_cnt_array[i] = 0;
   }

   for (i = 0; i < exp_Nrows; i++) 
   {
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count)
      {
          index = aggr_cnt_array[aggr_index[i]];
          aggr_cnt_array[aggr_index[i]]++;
          rows_in_aggs[aggr_index[i]][index] = i;
      }
   }
   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++) 
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&qr_tmp, nbytes, "AGu");
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&tmp_vect, nbytes, "AGv");

   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&work, nbytes, "AGw");

   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */
   for (i = 0; i < aggr_count; i++) 
   {
      /* ---------------------------------------------------------- */
      /* set up the matrix we want to decompose into Q and R:       */
      /* ---------------------------------------------------------- */
      length = aggr_cnt_array[i];
      if (nullspace_vect == NULL) 
      {
         for (j = 0; j < (int) length; j++)
         {
            index = rows_in_aggs[i][j];
            for (k = 0; k < nullspace_dim; k++)
            {
	      if ( unamalg_bdry[blocks[index]] == 'T')
                qr_tmp[k*length+j] = 0.;
               else
               {
                  if (block_pde[index] == k)     
                     qr_tmp[k*length+j] = 1.0;
                  else                           
                     qr_tmp[k*length+j] = 0.0;
               }
            }
         }
      }
      else 
      {
         for (k = 0; k < nullspace_dim; k++)
         {
            for (j = 0; j < (int) length; j++)
            {
               index = rows_in_aggs[i][j];
	      if ( unamalg_bdry[blocks[index]] == 'T')
		 qr_tmp[k*length+j] = 0.;
               else {
                  if (index < Nrows) {
                     qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
                  }
                  else {
		    fprintf( stderr,
			     "*ML*ERR* in QR, index out of range\n*ML*ERR* (file %s, line %d)\n",
			     __FILE__,__LINE__ );
		    fflush(stderr);exit( EXIT_FAILURE );
                  }
               }
            }
         }
      }
      /* ---------------------------------------------------------- */
      /* for this variable block version, I think I have to put in  */
      /* a test for zero-columns.                                   */
      /*the effect of zero-columns in the nullspace is not yet clear*/
      /* to me, they might cause trouble, so better know about them */

      printflag=0;
      length = aggr_cnt_array[i];
      for (k = 0; k < nullspace_dim; k++)
      {
         index1 = -1;
         for (j = 0; j < (int) length; j++)
         {
            if (ML_dabs(qr_tmp[k*length+j])>=EPS9)
            {
               index1=j;
               break;
            }
         }
         if (index1==-1)
         {
	     fprintf( stderr,
	              "*ML*WRN* detected zero column in QR block,\n*ML*WRN* that might crash the QR-fact.\n*ML*WRN* column %d of nullspace\n*ML*WRN* (file %s, line %d)\n",
	              k,__FILE__,__LINE__ );
	     fflush(stderr);
             printflag++;
             
         }
         
      }
#if 0 /* for debugging purpose */
      if (printflag>0)
      {
         fprintf(stdout,"Number of zero columns %d\n",printflag);
         for (j = 0; j < (int) length; j++)
         {
            for (k=0; k<nullspace_dim; k++)
            {
               fprintf(stdout,"%10.6e ",qr_tmp[k*length+j]);
            }
            fprintf(stdout,"\n");
         }
         fflush(stdout);
      }
#endif      
      /* ---------------------------------------------------------- */
      /* ---------------------------------------------------------- */
      /* now calculate QR using an LAPACK routine                   */
      /* ---------------------------------------------------------- */

      if (aggr_cnt_array[i] >= nullspace_dim) {

	DGEQRF_F77(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp, 
			  &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0)
	  pr_error("%s **ERR** dgeqrf returned a non-zero %d %d\n",str,
		   aggr_cnt_array[i],i);

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGx");
	  }
	else lwork=(int) work[0];
		 
	/* ---------------------------------------------------------- */
	/* the upper triangle of qr_tmp is now R, so copy that into   */
	/* the new nullspace                                          */
	/* ---------------------------------------------------------- */

	for (j = 0; j < nullspace_dim; j++)
	  for (k = j; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
	      qr_tmp[j+aggr_cnt_array[i]*k];
		 
#if 0 /* for debugging purpose */
        if (printflag>0)
        {
           fprintf(stdout,"QR-fact: This is R:\n");
           for (j = 0; j <nullspace_dim; j++)
           {
              for (k=0; k<j; k++)
                 fprintf(stdout,"%10.6e ",0.0);
              for (k=j; k<nullspace_dim; k++)
              {
                 fprintf(stdout,"%10.6e ",qr_tmp[j+aggr_cnt_array[i]*k]);
              }
              fprintf(stdout,"\n");
           }
           fflush(stdout);
        }
#endif      
	/* ---------------------------------------------------------- */
	/* to get this block of P, need to run qr_tmp through another */
	/* LAPACK function:                                           */
	/* ---------------------------------------------------------- */

	if ( aggr_cnt_array[i] < nullspace_dim ){
	  printf("Error in dorgqr on %d row (dims are %d, %d)\n(file %s, line %d)\n",
		 nullspace_dim, i,aggr_cnt_array[i],__FILE__,__LINE__ );
	  printf("ERROR : performing QR on a MxN matrix where M<N.\n");
	  fflush(stdout);
        }
	DORGQR_F77(&(aggr_cnt_array[i]), &nullspace_dim, &nullspace_dim, 
			  qr_tmp, &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0) {
	  printf("Error in dorgqr on %d row (dims are %d, %d)\n(file %s, line %d)\n",
		 i,aggr_cnt_array[i],nullspace_dim,__FILE__,__LINE__ );
          fflush(stdout);       
	  pr_error("Error in CoarsenMIS: dorgqr returned a non-zero\n");
	}

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGy");
	  }
	else lwork=(int) work[0];

	/* ---------------------------------------------------------- */
	/* now copy Q over into the appropriate part of P:            */
	/* The rows of P get calculated out of order, so I assume the */
	/* Q is totally dense and use what I know of how big each Q   */
	/* will be to determine where in ia, ja, etc each nonzero in  */
	/* Q belongs.  If I did not assume this, I would have to keep */
	/* all of P in memory in order to determine where each entry  */
	/* should go                                                  */
	/* ---------------------------------------------------------- */

	for (j = 0; j < aggr_cnt_array[i]; j++)
	  {
	    index = rows_in_aggs[i][j];
	    
	    if ( index < Nrows )
	      {
		index3 = new_ia[index];
		for (k = 0; k < nullspace_dim; k++) 
		  {
		    new_ja [index3+k] = i * nullspace_dim + k;
		    new_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
		  }
	      }
	    else 
	      {
	     fprintf( stderr,
	              "*ML*ERR* in QR: index out of bounds (%d)\n*ML*ERR* (file %s, line %d)\n",
		      index,__FILE__,__LINE__ );
	     fflush(stderr);exit( EXIT_FAILURE );
	      }
	  }
#if 0 /* for debugging purpose */
        if (printflag>0)
        {
           fprintf(stdout,"QR-fact: This is Q:\n");
           for (j = 0; j <aggr_cnt_array[i]; j++)
           {
              for (k=0; k<nullspace_dim; k++)
              {
                 fprintf(stdout,"%12.6e ",qr_tmp[ k*aggr_cnt_array[i]+j]);
              }
              fprintf(stdout,"\n");
           }
           fflush(stdout);
           exit(0);
        }
#endif      
      }
      else {
	/* We have a small aggregate such that the QR factorization can not */
	/* be performed. Instead let us copy the null space from the fine   */
        /* into the coarse grid nullspace and put the identity for the      */
	/* prolongator????                                                  */
	fprintf(stdout,
		"*ML*WRN* in QR: aggregate (%d) smaller then nullspace\n*ML*WRN* (file %s, line %d)\n",
		 index,__FILE__,__LINE__ );
	fflush(stdout);
	
        for (j = 0; j < nullspace_dim; j++)
	  for (k = 0; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
	      qr_tmp[j+aggr_cnt_array[i]*k];
	for (j = 0; j < aggr_cnt_array[i]; j++) {
	  index = rows_in_aggs[i][j];
	  index3 = new_ia[index];
	  for (k = 0; k < nullspace_dim; k++) {
	    new_ja [index3+k] = i * nullspace_dim + k;
	    if (k == j) new_val[index3+k] = 1.;
	    else new_val[index3+k] = 0.;
	  }
	}
      }


   } /*  for (i = 0; i < aggr_count; i++) */

   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse*nullspace_dim);
   ML_memory_free( (void **) &new_null);

   /* ------------------------------------------------------------- */
   /* set up the block information for the next level.              */
   /* For the time being, we assume constant block sizes for coarse */
   /* levels.                                                       */
   /* Note that this ML_Aggregate_CoarsenVBMETIS routine works with */
   /* constant block sizes as well though slower than               */
   /* ML_Aggregate_CoarsenMETIS                                     */
   /* ------------------------------------------------------------- */
   
   if (ml_ag->cur_level+1 < ml_ag->max_levels)
   {
      new_nblocks   = Ncoarse;
      new_nrows     = new_nblocks*nullspace_dim;
      new_blocks    = (int*)ML_allocate(new_nrows*sizeof(int));
      new_block_pde = (int*)ML_allocate(new_nrows*sizeof(int));
      if (new_block_pde==NULL)
      {
      	 fprintf( stderr,
	          "*ML*ERR* ML_Aggregate_CoarsenVBMETIS: out of space\n*ML*ERR* (file %s, line %d)\n",
		  __FILE__,__LINE__ );
	 fflush(stderr);exit( EXIT_FAILURE );
      }
      index1=0;
      for (i=0; i<new_nrows; i+=nullspace_dim)
      {
         for (j=0; j<nullspace_dim; j++)
         {
            new_blocks[i+j]    = index1;
            new_block_pde[i+j] = j;
         }
         index1++;
      }
      ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS(ml_ag,ml_ag->cur_level+1,
                                                     ml_ag->max_levels,
                                                     new_nblocks,new_blocks,
                                                     new_block_pde,new_nrows);
      ML_free(new_blocks);
      ML_free(new_block_pde);
   }
      


   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

#ifdef EXTREME_DEBUGGING
   print("set up the csr_data data structure\n");
#endif
   
   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   ML_Operator_Set_ApplyFuncData( *Pmatrix, nullspace_dim*Ncoarse, Nrows, 
                                  csr_data, Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   (*Pmatrix)->getrow->pre_comm = ML_CommInfoOP_Create();
   (*Pmatrix)->max_nz_per_row = 1; /* why this? mgee */
   
   ML_Operator_Set_Getrow((*Pmatrix), Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc((*Pmatrix), CSR_matvec);
   (*Pmatrix)->max_nz_per_row = 1; /* why this? mgee */
   /* this must be set so that the hierarchy generation does not abort early
      in adaptive SA */
   (*Pmatrix)->num_PDEs = nullspace_dim;

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */
#ifdef EXTREME_DEBUGGING
   print("clean up\n");
#endif
   
   if (unamalg_bdry) ML_free(unamalg_bdry);
   ML_memory_free((void**)&aggr_index);
   if (aggr_cnt_array) ML_free(aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) 
   if (rows_in_aggs[i]) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp);
   ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);

   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->length > 0 ) ML_free( supernode->list );
      ML_free( supernode );
   }
#if defined(OUTPUT_AGGREGATES) || defined(INPUT_AGGREGATES)
   /* Print Pmatrix*v (where v is constructed using global indices) */

   dtemp = (double *) ML_allocate(sizeof(double)*(*Pmatrix)->invec_leng);
   d2temp = (double *) ML_allocate(sizeof(double)*(*Pmatrix)->outvec_leng);
   for (i = 0; i < (*Pmatrix)->outvec_leng; i++) d2temp[i] = 0.;
   for (i = 0; i < (*Pmatrix)->invec_leng; i++)
      dtemp[i] = (double) (i+agg_offset);

   sprintf(fname,"PP%d_%d",comm->ML_mypid,level_count);
   fp = fopen(fname,"w");
   ML_Operator_Apply(*Pmatrix, (*Pmatrix)->invec_leng, dtemp, 
                     (*Pmatrix)->outvec_leng, d2temp);
   for (i = 0; i < Nrows; i++) {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 1) { j = i; k = i;} 
#else
      if (level_count == 1) { j = update_index[i]; k = update[i];} 
#endif
#else
      if (level_count == 1) { j = reordered_glob_nodes[i]; k = global_node_inds[i];}
#endif /* ifndef MAXWELL */
      else                  { j = i              ; k = i+vertex_offset;}
      fprintf(fp,"PP%d(%d) (loc=%d) = %e\n",level_count,k,j, d2temp[j]);
   }
   fclose(fp);
   ML_free(dtemp);
   ML_free(d2temp);

   /*
   csr_data = (struct ML_CSR_MSRdata *) (*Pmatrix)->data;
   if (comm->ML_myp`id == 1)
   {
      printf("%d : row_ptr = %d\nrow_ptr = %d\ncol = %d\nval = %e\n\n\n",
             comm->ML_mypid,
             csr_data->rowptr[14], csr_data->rowptr[15],
             csr_data->columns[14], csr_data->values[14]);
   }
   sprintf(fname,"comm%d",comm->ML_mypid);
   ML_CommInfoOP_Print((*Pmatrix)->getrow->pre_comm, fname);
   fflush(stdout);
   */

#endif
   return Ncoarse*nullspace_dim;
   
} /* ML_Aggregate_CoarsenVBMETIS */

/* ------------------------------------------------------------------------ */
/*!
 \brief calls METIS to decompose the graph of the local matrix
 <pre>
    This function calls METIS to decompose the graph of the local matrix      
    (that is, ignoring any inter-domain connections)                          
    this function is derived from ML_DecomposeGraph_with_METIS of which it is 
    the variable block version                                                
    
    Author:  Michael W. Gee, Org. 9214, November 2004
     
</pre>
\param Amatrix              ML_Operator*       (input)     the operator of the current level
\param N_parts              int                (input)     local number of aggregates requested
\param graph_decomposition  int[]              (output)    the aggregate information
\param bdry_nodes           char[]             (output)    flag showing boundary nodes, that is rows without off-diagonal entries
\param local_or_global      int                (input)     flag to switch whether output is local or global indizes
\param offsets              int[]              (output)    ? not used here?
\param reorder_flag         int                (input)     flag indicating fill-in reducing reordering with metis
\param current_level        int                (input)     number of current level
\param total_nz             int*               (output)    total number of nonzeros

 \sa  ML_Aggregate_Set_CoarsenScheme_VBMETIS ML_Aggregate_CoarsenVBMETIS
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS ML_Aggregate_Destroy_Vblocks_CoarsenScheme_VBMETIS

 \note this function is derived from ML_DecomposeGraph_with_METIS
 
 \warning some features of the original routine ML_DecomposeGraph_with_METIS are
          in here, but they are not tested! E.g. local_or_global=ML_GLOBAL_INDICES

*/
/* ------------------------------------------------------------------------ */
static int ML_DecomposeGraph_with_VBMETIS( ML_Operator *Amatrix,
					   int N_parts,
					   int graph_decomposition[],
					   char bdry_nodes[],
					   int local_or_global,
					   int offsets[],
					   int reorder_flag,
					   int current_level,
					   int *total_nz)
{

  int i, j,jj,  count, count2;
  int Nrows, Nrows_global,NrowsMETIS, N_nonzeros, N_bdry_nodes;
  int *wgtflag=NULL, numflag, *options=NULL, edgecut;
  idxtype *xadj=NULL, *adjncy=NULL;
#ifdef HAVE_ML_METIS
  idxtype *vwgt=NULL, *adjwgt=NULL;
#endif
  idxtype *part=NULL;
  ML_Comm * comm;
  int allocated = 0;
  int * rowi_col = NULL;
  int rowi_N;
  double * rowi_val = NULL;
  int nbytes = 0, nbytes_max = 0;
  int ok = 0;
  int * nodes_per_aggre = NULL;
  double t0;
  int * perm = NULL;
  char str[80];
  char *ptrToBdry;
  
  /* ------------------- execution begins --------------------------------- */
  
  t0 = GetClock();

  sprintf( str, "VBMETIS (level %d) :", current_level );
  
  comm = Amatrix->comm;

  /* dimension of the problem (NOTE: only local matrices) */
  
  Nrows = Amatrix->getrow->Nrows;
  perm = (int *) ML_allocate( sizeof(int) * Nrows );

  /* for some Epetra_matrices, N_nonzeros is set to -1.
     In this case, get all rows to allocate memory for adjncy.
     Also, define the set of boundary nodes. NOTE: the computation of
     nonzero elements is not really needed (ML_Operator usually have
     this number already compuuted. However, I still need to
     define the boundary nodes, and to handle epetra matrices.) 
     Finally, I need to compute the number of rows to give in input to
     METIS. Those do not include Dirichlet rows. */
  
  
  /* check row with only main diagonal entry -> boundry nodes */   
  N_nonzeros = 0;
  NrowsMETIS = 0;
  ptrToBdry = ML_Operator_IdentifyDirichletRows(Amatrix);
  /* TODO  Note that N_nonzeros could be greater than what was previously
     calculated in the code that this replaces.  I don't know if this will
     be a problem.... */ 
  N_nonzeros = ML_Operator_ComputeNumNzs(Amatrix);
  for (i=0; i<Nrows; i++) {
    bdry_nodes[i] = ptrToBdry[i];
    if (ptrToBdry[i] == 'T') perm[i] = -1;
    else                     perm[i] = NrowsMETIS++;
  }
  /* N_bdry_nodes: number of boundary blocks */
  N_bdry_nodes = ML_Comm_GsumInt(comm, Nrows-NrowsMETIS);
  Nrows_global = ML_Comm_GsumInt(comm, Nrows);
  
  if( comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
    printf("%s # bdry (vblock) nodes = %d, # (vblock) nodes = %d\n",
	   str,
	   N_bdry_nodes, Nrows_global);
  }
  
  /* construct the CSR graph information of the LOCAL matrix
     using the get_row function */

  wgtflag = (idxtype *) ML_allocate (4*sizeof(idxtype));
  options = (int *)     ML_allocate (4*sizeof(int));
  
  /* set parameters */
   
  wgtflag[0] = 0;    /* no weights */
  numflag    = 0;    /* C style */
  options[0] = 0;    /* default options */
   
  xadj    = (idxtype *) ML_allocate ((NrowsMETIS+1)*sizeof(idxtype));
  adjncy  = (idxtype *) ML_allocate ((N_nonzeros)*sizeof(idxtype));
   
  if(  xadj==NULL || adjncy==NULL ) {
    fprintf( stderr,
	     "**ERR** on proc %d, not enought space for %d bytes.\nfile %s, line %d\n",
	     comm->ML_mypid, N_nonzeros,
	     __FILE__,
	     __LINE__);fflush(stderr);exit(EXIT_FAILURE);
  }
   
  count = 0; count2 = 0; xadj[0] = 0;
  
  for (i = 0; i < Nrows; i++) {

    if( bdry_nodes[i] == 'F' ) {

      xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */
    
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
		        &rowi_N, 0);

      /* need to avoid boundary nodes in METIS vectors. Skip them */
      /* (I am not pretty sure that rows with zero elements are   */
      /* well treated by METIS.) perm has been allocates of size    */
      /* Nrows, so columns corresponding to external nodes can not*/
      /* be given as input to perm                                */

      for( j=0 ; j<rowi_N ; j++ ) {
	jj = rowi_col[j];
	if( jj<Nrows ) {
	  if( jj != i && perm[jj] != -1 ) {
	    adjncy[count++] = perm[jj];
	    xadj[count2+1]++;
	  }
	}
      }
      count2++;
    }      
  }

  *total_nz = count;

#ifdef DUMP_MATLAB_FILE
      sprintf( str, "METIS_proc%d.m", comm->ML_mypid);
      fp = fopen(str,"w");
      fprintf(fp,"NrowsMETIS = %d;\n", NrowsMETIS);
      fprintf(fp,"xadj = zeros(NrowsMETIS,1);\n");
      for( i=0 ; i<NrowsMETIS+1 ; i++ ) {
	fprintf(fp,"xadj(%d) = %d;\n", i+1, xadj[i]+1);
      }
      fprintf(fp,"Nonzeros = %d;\n", count);
      fprintf(fp,"adjncy = zeros(Nonzeros,1);\n");
      for( i=0 ; i<count ; i++ ) {
	fprintf(fp,"adjncy(%d) = %d;\n", i+1, adjncy[i]+1);
      }
      fprintf(fp,"A = zeros(%d,%d)\n", NrowsMETIS,NrowsMETIS);
      for( i=0 ; i<NrowsMETIS ; i++ ) {
	for( j=xadj[i] ; j<xadj[i+1] ; j++ ) {
	  fprintf(fp,"A(%d,%d) = 1;\n",
		  i+1,adjncy[j]+1);
	}
      }
      fclose(fp);
#endif

#ifdef DUMP_WEST
      sprintf( str, "METIS_proc%d.m", comm->ML_mypid);
      fp = fopen(str,"w");
      fprintf(fp,"Nrows = %d\n", NrowsMETIS);
      for( i=0 ; i<NrowsMETIS+1 ; i++ ) {
	fprintf(fp,"%d\n", xadj[i]);
      }
      fprintf(fp,"Nonzeros = %d\n", count);
      for( i=0 ; i<count ; i++ ) {
	fprintf(fp,"%d\n", adjncy[i]);
      }
      fclose(fp);
#endif
      
  if( count > N_nonzeros || count2 != NrowsMETIS ) {
    fprintf( stderr,
	     "*ML*WRN* On proc %d, count  > N_nonzeros (%d>%d)\n"
	     "*ML*WRN* or count2 != NrowsMETIS (%d>%d)\n"
	     "a buffer overflow has probably occurred...\n",
	     comm->ML_mypid, count, N_nonzeros, count2, NrowsMETIS );
  }

  /* idxtype is by default int, but on some architectures can be
     slightly different (for instance, a short int). */
   
  part = (idxtype *) ML_allocate( sizeof(idxtype) * NrowsMETIS );
  nodes_per_aggre  = (int *) ML_allocate( sizeof(int) * N_parts );

  /* ********************************************************************** */
  /* Before calling METIS, I verify that the two extreme situations are     */
  /* handled separately.                                                    */
  /* ********************************************************************** */
  
  if( N_parts == 1 ) {

    for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = 0;
    edgecut = 0;
    
  } else if( N_parts == NrowsMETIS ) {

    fprintf( stderr,
	     "*ML*WRN*: on proc %d, N_part == N_rows_noDirichlet (%d==%d)\n",
	     comm->ML_mypid, N_parts, NrowsMETIS );
 
    for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = i;
    edgecut = 0;
  
  } else {

    ok = 0;

    while( ok == 0 ) {
      
      /* ****************************************************************** */
      /* Put -1 in part, so I can verify that METIS has filled each pos    */
      /* ****************************************************************** */

      for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = -1;
    
      /* ****************************************************************** */
      /* Estimate memory required by METIS. This memory will be dynamically */
      /* allocated inside; however this is a good estimation of how METIS   */
      /* will cost in terms of memory.                                      */
      /* Then, call METIS.                                                  */
      /* ****************************************************************** */

#ifdef HAVE_ML_METIS
      if( N_parts < 8 ) {

	i = 1; /* optype in the METIS manual */
	numflag = 0;
	METIS_EstimateMemory( &NrowsMETIS, xadj, adjncy, &numflag,
			      &i, &nbytes );
	
	METIS_PartGraphRecursive (&NrowsMETIS, xadj, adjncy, vwgt, adjwgt,
				  wgtflag, &numflag, &N_parts, options,
				  &edgecut, part);
      } else {
	
	i = 2;
	numflag = 0;

	METIS_EstimateMemory( &NrowsMETIS, xadj, adjncy, &numflag,
			      &i, &nbytes );
	
	METIS_PartGraphKway (&NrowsMETIS, xadj, adjncy, vwgt, adjwgt,
			     wgtflag, &numflag, &N_parts, options,
			     &edgecut, part);
      }
#else
      if( Amatrix->comm->ML_mypid == 0 ) {
	fprintf( stderr,
		 "*ML*WRN* This function has been compiled without the configure\n"
		 "*ML*WRN* option --with-ml_metis\n"
		 "*ML*WRN* I will put all the nodes in the same aggregate, this time...\n"
		 "*ML*WRN* (file %s, line %d)\n",
		 __FILE__,
		 __LINE__);
      }
      for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = 0;
      N_parts = 1;
#endif
      
      /* **************************************************************** */
      /* perform some checks. If aggregates with zero assigned nodes      */
      /* exist, then recall METIS, asking for a smaller number of sub     */
      /* graphs. This is the role of the `ok' variable.                   */
      /* Also, if the part vector contains some junk, recall METIS        */
      /* **************************************************************** */

      ok = 1;
      
      for( i=0 ; i<N_parts ; i++ ) nodes_per_aggre[i] = 0;
      for( i=0 ; i<NrowsMETIS ; i++ ) {
	j = part[i];
	if( j<0 || j>= N_parts ) {
	  ok = 0;
	  break;
	} 
	else nodes_per_aggre[j]++;
      }
      
      for( i=0 ; i<N_parts ; i++ ) {
	if( nodes_per_aggre[i] == 0 ) {
	  ok = 0;
	  break;
	}
      }
      
      if( ok == 0 ) {
	if( comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	  printf( "*ML*WRN* input # of (block) aggregates (%d) does not assure "
		  "non-empty aggregates.\n"
		  "*ML*WRN* Now recalling METIS with # aggregates = %d\n",
		  N_parts, N_parts/2 );
	}
	N_parts = N_parts/2;
      }
      
      if( N_parts == 0 ) {
	if( comm->ML_mypid == 0 && 9 < ML_Get_PrintLevel()) {
	  fprintf( stderr,
		   "*ML*WRN* something went **VERY** wrong in calling METIS\n"
		   "*ML*WRN* try to ask for a smaller number of subdomains\n"
		   "*ML*WRN* I will put all the nodes into one aggregate...\n"
		   "*ML*WRN* (file %s, line %d)\n",
		   __FILE__,
		   __LINE__ );
	}
	N_parts = 1;
      }
      
      /* ************************************************************** */
      /* handle the case N_parts = 1 separately. Do not recall METIS    */
      /* in this case, simply put everything to zero and continue       */
      /* ************************************************************** */
      
      if( N_parts == 1 ) {
	for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = 0;
	ok = 1;
      }
      
    } /* while( ok == 0 ) */
  
  } /* if( N_parts == 1 ) */

  /* ********************************************************************** */
  /* Some fancy output for memory usage.                                    */
  /* ********************************************************************** */

  nbytes /= 1024;
  
  nbytes_max = ML_gmax_int( nbytes, comm );
  nbytes = ML_gsum_int( nbytes, comm);

  if( Amatrix->comm->ML_mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("%s Estimated required mem for METIS = %d Kb\n"
	   "%s Max estimated mem for METIS = %d Kb\n",
	   str,
	   nbytes,
	   str,
	   nbytes_max );
  }
  
  /* ********************************************************************** */
  /* reordering using METIS to minimize the fill-in during factorization    */
  /* ********************************************************************** */

  if( reorder_flag == ML_YES ) {
    fprintf( stderr,
         "**ERR** on proc %d, Local reordering not impl. for variable blocks\n"
         "file %s, line %d\n",
         comm->ML_mypid,__FILE__,__LINE__);fflush(stderr);exit(EXIT_FAILURE);

    
  }

  /* copy back part into aggr_index, and set to -1
     the aggr_index corresponding to ghost nodes */

  for( i=0 ; i<Nrows ; i++ ) {
    j = perm[i];
    if( j != -1 ) 
      graph_decomposition[i] = (int)part[j];
    else
      graph_decomposition[i] = -1;
  }

  /* if global indices are required, modify the entries
     of graph_decomposition (only the LOCAL entries) so that
     they correspond to global indices. Also, set the array
     offsets, defined so that the global indices assigned
     to processor i are
     offsets[i] <= indices_of_i < offsets[i+1]
     Note that I do not suppose that N_parts is the same
     value among all the processors */
     
  if( local_or_global == ML_GLOBAL_INDICES ) {
    fprintf( stderr,
         "%s **WRN** on proc %d, global indizes not checked for variable blocks\n"
         "good luck! file %s, line %d\n",str,
         comm->ML_mypid,__FILE__,__LINE__);fflush(stderr);    
    ML_DecomposeGraph_BuildOffsets( N_parts, offsets, comm->ML_nprocs,
				    Amatrix->comm->USR_comm );
  }

  /* ------------------- that's all folks --------------------------------- */

  ML_free(rowi_col); ML_free(rowi_val);
  rowi_col = NULL; rowi_val = NULL;
  allocated = 0; 

  if( options != NULL ) ML_free( options );
  if( wgtflag != NULL ) ML_free( wgtflag );
  if( adjncy != NULL  ) ML_free( adjncy  );
  if( xadj != NULL    ) ML_free( xadj    );
  if( part != NULL    ) ML_free( part    );
  if( perm != NULL    ) ML_free( perm    );
  if( nodes_per_aggre != NULL ) ML_free( nodes_per_aggre );
  
  t0 = GetClock() - t0;

  if ( comm->ML_mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("%s Time to partition graph = %e (s)\n",
	   str,
	   t0 );
    
  }
  
  return N_parts;
  
} /* ML_DecomposeGraph_with_VBMETIS */
