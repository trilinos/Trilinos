/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Functions for generating grid transfer functions between 2 grids     */
/* ******************************************************************** */
/* Given 2 dependently or independently generated grids, this function  */
/* generates the grid transfer operator between the two grids.  Even    */
/* though the two grids can be input in an arbitrary manner, it is      */
/* recommended that the first argument of this function be the fine     */
/* grid for efficiency reasons.                                         */
/*                                                                      */
/* A. Inputs :                                                          */
/*                                                                      */
/*    f_grid  : the fine grid                                           */
/*    c_grid  : the coarse grid                                         */
/*                                                                      */
/* B. Output :                                                          */
/*                                                                      */
/*    xsfer_op    : the grid transfer operator                          */
/*                                                                      */
/* C. Additional Notes :                                                */
/*                                                                      */
/*    1. In order to alleviate the users of the burden to understand    */
/*       how to apply the restriction and interpolation operators, two  */
/*       other subroutines have been developed, namely,                 */
/*                                                                      */
/*         - ML_restrict(op, inlen, indata, outlen, outdata)            */
/*         - ML_interpolate(op, inlen, indata, outlen, outdata)         */
/*                                                                      */
/*       which understands the data structures for the restriction and  */
/*       interpolation operators.                                       */
/*                                                                      */
/*    2. This subroutine makes use of a local function to generate the  */
/*       interpolating coefficients, namely,                            */
/*                                                                      */
/*         USR_compute_basis_coefficients(grid,ind,coord,n,coefs,c_ptr) */
/*                                                                      */
/*       In practice, this get basis subroutine should be provided by   */
/*       users.                                                         */
/* ******************************************************************** */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : February, 1998                                       */
/* ******************************************************************** */

/* ******************************************************************** */
/* include files and global constants                                   */
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_memory.h"
#include "ml_comm.h"
#include "ml_comminfoagx.h"
#include "ml_intlist.h"
#include "ml_gridagx.h"
#include "ml_gridfunc.h"
#include "ml_operatoragx.h"
#include "ml_struct.h"
#include "ml_setup.h"

/* *********************************************************************** */
/* ML_setup : Given a fine grid f_grid and a coarse grid c_grid, set up    */
/*            the transfer operators between the 2 grids.                  */
/* *********************************************************************** */

void ML_setup_grid_xsfer_op(void *f_grid, ML_GridFunc *fgrid_fcns,
                            void *c_grid, ML_GridFunc *cgrid_fcns,
                            void **xsfer, ML_Comm *comm)
{
   int            checkpt=1, mypid;
#ifdef DEBUG
   int  status;
#endif
   ML_IntList     *lfv_list;
   ML_CommInfoAGX *comm_info;
   ML_GridAGX     *g_c_grid;
   ML_OperatorAGX *xsfer_op;

   /* -------------------------------------------------------------------- */
   /* initialize user links                                                */
   /* -------------------------------------------------------------------- */


#ifdef DEBUG
   status = ML_Comm_Check( comm );
   if ( status != 0 )
      printf("ML_setup : fail communicator check. \n");
   status = ML_GridFunc_Check( fgrid_fcns );

   if ( status != 0 )
      printf("ML_setup : fail fine grid functions check. \n");

   status = ML_GridFunc_Check( cgrid_fcns );
   if ( status != 0 )
      printf("ML_setup : fail coarse grid functions check. \n");
#endif
   gridfcns_basis = cgrid_fcns;

   /* -------------------------------------------------------------------- */
   /* Compose a complete coarse grid on each processor from one that is    */
   /* distributed. (in g_c_grid upon return)                               */
   /* -------------------------------------------------------------------- */

   mypid = comm->ML_mypid;
   /* ML_debug = -1 ; */
   /* if ( mypid == 1 ) ML_debug = 1;  turn on debug flag on processor 3 */
   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 )
   {
      printf("ML processing begins : \n");
      printf("Composing global grid ... \n");
   }
   ML_compose_global_grid( c_grid, cgrid_fcns, &g_c_grid, comm );

   /* -------------------------------------------------------------------- */
   /* With a complete global coarse grid (g_c_grid) on each processor,     */
   /* calculate the local grid transfer operator - that is, to find out    */
   /* which of my local fine nodes lie inside which of my local coarse     */
   /* elements, and compute the interpolating coefficients.                */
   /* (The partial grid transfer operator is returned in xsfer_op and the  */
   /*  outgoing coefficients are also stored temporarily in xsfer_op.)     */
   /* -------------------------------------------------------------------- */

   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 ) printf("Constructing local transfer operator ... \n");
   ML_construct_RP0( c_grid, cgrid_fcns, f_grid, fgrid_fcns, g_c_grid,
                     &xsfer_op, comm );
   (*xsfer) = (void *) xsfer_op;

   /* -------------------------------------------------------------------- */
   /* In the following, identify remote coarse elements where my local     */
   /* fine nodes possibly reside (in the form of coarse element to node    */
   /* lists).                                                              */
   /* -------------------------------------------------------------------- */

   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 ) printf("Composing grid candidates ... \n");
   ML_IntList_Create( &lfv_list, g_c_grid->Nelements, 0 );
   ML_remote_grid_candidates(f_grid, fgrid_fcns, cgrid_fcns, g_c_grid,
                             lfv_list, xsfer_op, comm);

   /* -------------------------------------------------------------------- */
   /* exchange candidates between processors                               */
   /* -------------------------------------------------------------------- */

   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 ) printf("Exchanging candidates ... \n");
   ML_CommInfoAGX_Create( &comm_info );
   ML_exchange_candidates( lfv_list, f_grid, fgrid_fcns, g_c_grid,
                           comm_info, comm );
   ML_IntList_Destroy( &lfv_list );

   /* -------------------------------------------------------------------- */
   /* calculate coefficients for remote processors                         */
   /* -------------------------------------------------------------------- */

   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 ) printf("Getting basis ... \n" );
   ML_get_basis_functions_coef( comm_info, c_grid, cgrid_fcns, xsfer_op );

   /* -------------------------------------------------------------------- */
   /* send the results of queries back to the requesting processors        */
   /* -------------------------------------------------------------------- */

   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 ) printf("Exchanging coefficients ... \n" );
   ML_exchange_coefficients( c_grid, cgrid_fcns, comm_info,
                             xsfer_op, comm );

   /* -------------------------------------------------------------------- */
   /* construct the remote part of the grid transfer operator              */
   /* -------------------------------------------------------------------- */

   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 ) printf("Constructing remote transfer operator \n");
   xsfer_op->AGX_comm = comm;
   ML_construct_RP1( f_grid, fgrid_fcns, c_grid, cgrid_fcns, g_c_grid,
                     comm_info, xsfer_op, comm);
   ML_GridAGX_Destroy( &g_c_grid );

   /* -------------------------------------------------------------------- */
   /* clean up and prepare to return                                       */
   /* -------------------------------------------------------------------- */

   ML_CommInfoAGX_Destroy( &comm_info );
   checkpt = ML_Comm_GmaxInt( comm, checkpt );
   if ( mypid == 0 ) printf( "ML processing ends. \n" );

}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ML_compose_global_grid : Given a grid distributed on all processors,    */
/*    compose in each processor the whole grid.                            */
/*   Note : local_grid can be any type of grid provided that user also     */
/*          provides a collection of grid access functions.  The global    */
/*          grid (g_c_grid), however, will be constructed in conformity    */
/*          with the grid structure in ml_grid.h.                          */
/* *********************************************************************** */

void ML_compose_global_grid( void *local_grid, ML_GridFunc *cgrid_fcns,
                             ML_GridAGX **g_c_grid, ML_Comm *comm )
{
   int     i, j, k, ndp1, proc_cnt, ivar1, ivar2, ivar3, *ibuf;
   ml_big_int *bigibuf;
   int     nproc, mypid, tot_elmnt_nvertices, leng, *tlist, tot_leng;
   int     *elmnt_proc_map, *node_proc_map, MAX_VERT_PER_ELE;
   double  *dbuf;
   ML_GridAGX *global_grid;

   /* ------------------------------------------------------------------ */
   /* Nvertices in the global grid should be the sum of Nvertices in the */
   /* local grid on all processors. (note : Nvertices does not include   */
   /* shared or external nodes)                                          */
   /* ------------------------------------------------------------------ */

   mypid = comm->ML_mypid;
   nproc = comm->ML_nprocs;
   MAX_VERT_PER_ELE = cgrid_fcns->ML_MaxElmntVert;
   if (cgrid_fcns->USR_grid_get_nvertices == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_nvertices() not found\n");
   if (cgrid_fcns->USR_grid_get_dimension == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_dimension() not found\n");
   if (cgrid_fcns->USR_grid_get_nelements == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_nelements() not found\n");
   if (cgrid_fcns->USR_grid_get_element_nvertices == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_element_nvertices() not found\n");
   if (cgrid_fcns->USR_grid_get_element_vlist == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_element_vlist() not found\n");
   if (cgrid_fcns->USR_grid_get_vertex_global_num == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_vertex_global_num() not found\n");
   if (cgrid_fcns->USR_grid_get_element_global_num == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_element_global_num() not found\n");
   if (cgrid_fcns->USR_grid_get_vertex_coordinate == NULL)
      pr_error("ML_compose_global_grid: USR_grid_get_vertex_coordinate() not found\n");


   ivar1 = cgrid_fcns->USR_grid_get_nvertices( local_grid );
   ivar2 = ML_Comm_GsumInt( comm, ivar1 );
   ML_GridAGX_Create( g_c_grid );
   global_grid = (*g_c_grid);
   global_grid->Nvertices = ivar2;
   global_grid->Nvertices_expanded = ivar2;
   global_grid->Ndim = cgrid_fcns->USR_grid_get_dimension( local_grid );

   /* ------------------------------------------------------------------ */
   /* Nelements in the global grid should be the sum of Nelements in the */
   /* local grid on all processors.                                      */
   /* ------------------------------------------------------------------ */

   ivar1 = cgrid_fcns->USR_grid_get_nelements( local_grid );
   ivar2 = ML_Comm_GsumInt( comm, ivar1 );
   global_grid->Nelements = ivar2;

   /* ------------------------------------------------------------------ */
   /* Compose the element-to-nodes list for the global grid on each      */
   /* processor                                                          */
   /* ------------------------------------------------------------------ */

   /* 1. compose the element-to-node pointer fields                      */

   ivar1++;
   ivar2  += nproc;
   ML_memory_alloc( (void**) &ibuf, ivar2 * sizeof( int ), "ibu" );
   ibuf[0] = 0;
   for ( i = 1; i < ivar1; i++ )
   {
      ibuf[i] = ibuf[i-1] +
                cgrid_fcns->USR_grid_get_element_nvertices(local_grid, i-1);
   }
   tot_elmnt_nvertices = ibuf[ivar1-1];

   ML_Comm_GappendInt( comm, ibuf, &ivar1, ivar2 );

   /* 2. allocate storage space for the coarse element to processor map  */
   /*    - this map will indicate the processors where the global        */
   /*      elements in global_grid.global_element reside                 */

   ivar1 = global_grid->Nelements;
   ML_memory_alloc( (void**) &elmnt_proc_map,ivar1*sizeof(int),"emp");

   /* 3. decode the concatencated integer vector to form a single array  */
   /*    of element-to-node pointer for the global grid (that is, the    */
   /*    variable global_grid->ele_nodes->start array). In addition, the */
   /*    element to processor map is also constructed.                   */

   k = proc_cnt = ivar1 = ivar3 = 0;
   for ( i = 1; i < ivar2; i++ )
   {
      if (ibuf[i] > ibuf[i-1])
      {
         elmnt_proc_map[ivar3++] = proc_cnt;
         ibuf[ivar3] = ibuf[i] + ivar1;
      }
      else
      {
         ivar1 += ibuf[i-1];
         proc_cnt++;
      }
   }

   ML_IntList_Create( &(global_grid->ele_nodes), 0, 0 );
   global_grid->ele_nodes->start = ibuf;
   global_grid->ele_nodes->length = global_grid->Nelements;
   global_grid->elmnt_proc_map    = elmnt_proc_map;

   /* 4. Once the element-to-node pointer is constructed, the next thing */
   /*    is to construct the global element-to-node list.                */
   /*    - concatenate the element-to-node lists across all processors   */
   /*      to give a long vector on each processor (using the special    */
   /*      special ML_GappendInt function which is different from the    */
   /*      one in Aztec).                                                */

   ivar1 = tot_elmnt_nvertices;
   ivar2 = ML_Comm_GsumInt( comm, ivar1 );
   ML_memory_alloc( (void**) &ibuf, ivar2 * sizeof( int ), "ib2" );
   tot_leng = i = 0;
   ML_memory_alloc( (void**) &tlist, MAX_VERT_PER_ELE*sizeof(int),"tl1");
   while (tot_leng < ivar1)
   {
      leng = cgrid_fcns->USR_grid_get_element_vlist(local_grid, i, tlist);
      for ( j = 0; j < leng; j++ )
      {
         ibuf[tot_leng+j] =
            cgrid_fcns->USR_grid_get_vertex_global_num(local_grid,tlist[j]);
      }
      i++;
      tot_leng += leng;
   }
   ML_memory_free( (void**) &tlist );
   ML_Comm_GappendInt( comm, ibuf, &ivar1, ivar2 );
   global_grid->ele_nodes->members = ibuf;

   /* ------------------------------------------------------------------ */
   /* Compose the global element number field of the Grid structure.     */
   /* (This again consists in first computing the data length on each    */
   /*  processor so that a buffer of appropriate size can be allocated.) */
   /* ------------------------------------------------------------------ */

   ivar1 = cgrid_fcns->USR_grid_get_nelements( local_grid );
   ivar2 = global_grid->Nelements;
   ML_memory_alloc( (void**) &bigibuf, ivar2 * sizeof( bigibuf[0] ), "ib2" );
   for ( i = 0; i < ivar1; i++ )
   {
      bigibuf[i] = cgrid_fcns->USR_grid_get_element_global_num(local_grid, i);
   }
   ML_Comm_GappendBigInt( comm, bigibuf, &ivar1, ivar2 );
   global_grid->global_element = bigibuf;

   /* ------------------------------------------------------------------ */
   /* Compose the global node  number field of the Grid structure.       */
   /* (This consists simply in setting global_vertex[i] = i since the    */
   /*  element-to-node list has been composed based on global node       */
   /*  numberings.)                                                      */
   /* ------------------------------------------------------------------ */

   ivar1 = global_grid->Nvertices;
   ML_memory_alloc((void**) &(global_grid->global_vertex),
                   ivar1 * sizeof(int), "gv1" );
   for ( i = 0; i < ivar1; i++ )
   {
      global_grid->global_vertex[i] = i;
   }

   /* ------------------------------------------------------------------ */
   /* Compose the (x,y,z) coordinate field of the Grid structure as well */
   /* as the node_proc_map field. (The global node numbers have to be    */
   /* communicated first, followed by the coordinates to ensure proper   */
   /* alignment of node to coordinates and initialization of the         */
   /* node_proc_map array).                                              */
   /* ------------------------------------------------------------------ */

   ivar1 = global_grid->Nvertices;
   ML_memory_alloc( (void**) &tlist, ivar1 * sizeof( int ), "tl2" );
   ML_memory_alloc( (void**) &node_proc_map,ivar1*sizeof(int), "nmp" );
   ivar3 = cgrid_fcns->USR_grid_get_nvertices( local_grid );
   for ( i = 0; i < ivar3; i++ )
   {
      tlist[i] = cgrid_fcns->USR_grid_get_vertex_global_num(local_grid, i);
   }
   ML_Comm_GappendInt( comm, tlist, &ivar3, ivar1 );

   ndp1    = global_grid->Ndim;
   ivar3   = cgrid_fcns->USR_grid_get_nvertices( local_grid );
   ivar1   = ivar3 * ndp1 + 1;
   ivar2   = global_grid->Nvertices * ndp1 + nproc + 1;
   ML_memory_alloc((void**) &dbuf,  ivar2 * sizeof( double ), "dbu" );
   dbuf[0] = (double) - mypid - 1000.0;
   for ( i = 0; i < ivar3; i++ )
   {
      cgrid_fcns->USR_grid_get_vertex_coordinate( local_grid,
                                        i, &(dbuf[i*ndp1+1]) );
   }
   ML_Comm_GappendDouble( comm, dbuf, &ivar1, ivar2 );
   dbuf[ivar2-1] = - 2000.0;
   k = global_grid->Nvertices;
   ML_memory_alloc((void**) &(global_grid->x), k*sizeof(double), "GGX");
   ML_memory_alloc((void**) &(global_grid->y), k*sizeof(double), "GGY");
   if ( ndp1 > 2 )
   {
      ML_memory_alloc((void**) &(global_grid->z), k*sizeof(double), "GGZ");
   }
   ivar2 = 0;
   ivar3 = 0;
   for ( i = 0; i < nproc; i++ )
   {
      ivar1 = (int) dbuf[ivar3++];
      ivar1 = - ( ivar1 + 1000 );
      while ( dbuf[ivar3] > -1000.0 )
      {
         k = tlist[ivar2++];
         global_grid->x[k] = dbuf[ivar3++];
         global_grid->y[k] = dbuf[ivar3++];
         if (ndp1 > 2) global_grid->z[k] = dbuf[ivar3++];
         if ( ivar1 != i ) printf("Error : processor no. not matched.\n");
         node_proc_map[k] = ivar1;
      }
   }
   global_grid->node_proc_map = node_proc_map;
   ML_memory_free( (void **) &dbuf );
   ML_memory_free((void **) &tlist );
}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ML_construct_RP0 : construct the local part of the grid transfer        */
/*                    operator.                                            */
/*                                                                         */
/*    INPUTS :                                                             */
/*                                                                         */
/*       c_grid      : coarse grid                                         */
/*       cgrid_fcns  : functions to access coarse grid                     */
/*       f_grid      : fine grid                                           */
/*       fgrid_fcns  : functions to access fine grid                       */
/*       g_c_grid    : global coarse grid                                  */
/*       node_flag   : indicate which fine nodes have been processed.      */
/*                     and local processor ID                              */
/*       comm        : communicator                                        */
/*    OUTPUTS :                                                            */
/*                                                                         */
/*       xsfer_op    : the grid transfer operator                          */
/*                                                                         */
/* *********************************************************************** */

void ML_construct_RP0(void *c_grid, ML_GridFunc *cgrid_fcns,
                      void *f_grid, ML_GridFunc *fgrid_fcns,
                      ML_GridAGX *g_c_grid,  ML_OperatorAGX  **xsfer_op2,
                      ML_Comm *comm)
{
   int    leng, ndim, ncvert, ncelmnts, nfvert, *elmnt_map, ngcelmnts, mypid;
   int    cur_coef_leng, nfxyz_leng, *coef_ptr, *fvlist, *vlist, gcnt, gpcnt;
   int    i, j, k, m, index, icnt, icnt2, ibegin, iend;
   ml_big_int ggelenum, gcelenum;
   int    mbegin, mend, ncnt, *fine2coarsecnts, ncand, cele_num;
   int    *coef_ptr2, cur_ptr_leng, *node_proc_map, *fnode_flag;
   int    extern_node_cnt, ext_cnt, cnodenum, cgnodenum, fnodenum;
   int    MAX_VERT_PER_ELE;
   double *coord, *coord_short, *coefs, *coefs2;
   ML_IntList     *lfv_list;
   ML_ElementAGX  *element;
   ML_OperatorAGX *xsfer_op;
   unsigned int nbytes;

   /* ----------------------------------------------------------------- */
   /* initialize the grid transfer operator (and allocate memory)       */
   /* ----------------------------------------------------------------- */

   ML_OperatorAGX_Create( xsfer_op2 );
   xsfer_op = (*xsfer_op2);

   /* ----------------------------------------------------------------- */
   /* fetch grid and processor information                              */
   /* ----------------------------------------------------------------- */
   if (fgrid_fcns->USR_grid_get_nvertices== NULL)
      pr_error("ML_construct_RP0: USR_grid_get_nvertices() not found\n");
   if (fgrid_fcns->USR_grid_get_vertex_coordinate == NULL)
      pr_error("ML_construct_RP0: USR_grid_get_vertex_coordinate() not found\n");
   if (cgrid_fcns->USR_compute_basis_coefficients == NULL)
      pr_error("ML_construct_RP0: USR_compute_basis_coefficients() not found\n");

   ndim       = cgrid_fcns->USR_grid_get_dimension( c_grid );
   ncvert     = cgrid_fcns->USR_grid_get_nvertices( c_grid );
   ncelmnts   = cgrid_fcns->USR_grid_get_nelements( c_grid );
   nfvert     = fgrid_fcns->USR_grid_get_nvertices( f_grid );
   nfxyz_leng = ndim * nfvert ;
   elmnt_map  = g_c_grid->elmnt_proc_map;
   ngcelmnts  = g_c_grid->Nelements;
   mypid      = comm->ML_mypid;
   MAX_VERT_PER_ELE = cgrid_fcns->ML_MaxElmntVert;

   /* ----------------------------------------------------------------- */
   /* allocate space to store coordinate and node information           */
   /*  fvlist - store local node numbers for the master list            */
   /*  coord  - store coordinates corresponding to fvlist               */
   /* ----------------------------------------------------------------- */

   ML_memory_alloc((void**) &fvlist, nfvert * sizeof(int), "fvl" );
   ML_memory_alloc((void**) &coord, nfxyz_leng * sizeof(double), "COO" );
   for ( i = 0; i < nfvert; i++ )
   {
      fvlist[i] = i;
      fgrid_fcns->USR_grid_get_vertex_coordinate(f_grid,i,&coord[ndim*i]);
      /* printf("%d : FVERTEX %d : %e %e \n", mypid, i, coord[ndim*i],
                                                     coord[ndim*i+1]); */
   }

   /* ----------------------------------------------------------------- */
   /* allocate space to store the local to local coefficients           */
   /* ----------------------------------------------------------------- */

   cur_ptr_leng  = nfvert + 2;
   ML_memory_alloc((void**) &coef_ptr, cur_ptr_leng*sizeof(int), "cep");
   cur_coef_leng = cur_ptr_leng * MAX_VERT_PER_ELE;
   ML_memory_alloc((void**) &coefs, cur_coef_leng*sizeof(double), "cef");
   coef_ptr[0]   = 0;

   /* ----------------------------------------------------------------- */
   /* traverse all local coarse element                                 */
   /* - at the end lfv_list stores a list of candidates for each coarse */
   /*   element in my process and coef_ptr and coefs arrays store the   */
   /*   pointers and the actual coefficients.                           */
   /* - vlist temporarily stores the candidates for each coarse element */
   /* - fnode_flag is used to keep track of which fine nodes have been  */
   /*   processed so that subsequent queries will not include them.     */
   /* ----------------------------------------------------------------- */

   ML_IntList_Create( &lfv_list, ncelmnts, 0);
   ML_ElementAGX_Create( &element, ndim, MAX_VERT_PER_ELE );
   ML_memory_alloc( (void**) &vlist, (nfvert + 2) * sizeof(int), "vlt" );
   ML_memory_alloc( (void**) &fnode_flag, nfvert * sizeof(int), "fnf" );
   xsfer_op->fnode_leng = nfvert;
   xsfer_op->fnode_flag = fnode_flag;
   gcnt  = 0;
   gpcnt = 1;
   for ( i = 0; i < nfvert; i++ ) fnode_flag[i] = -1;

   for ( i = 0; i < ngcelmnts; i++ )
   {
      if ( elmnt_map[i] == mypid )
      {
         ML_GridAGX_Get_Element( g_c_grid, i, element );

         /* search fine vertices for candidates for the given element.  */
         /* Upon return, vlist[1...] should contain the indices of the  */
         /* candidates, and ncand contains the number of such vertices) */
         /* fnode_flag is an array to tell the element to ignore those  */
         /* fine nodes with the corresponding entry not equal to -1.    */

         ML_ElementAGX_ComposeCandidates(element, nfvert, coord,
                                fvlist, fnode_flag, &ncand, &vlist[1]);

         /* if a non-empty list is returned, load the candidate info to */
         /* the int_lists structure as one row of a CSR matrix (so the  */
         /* number of rows of this list at the end is the no. of coarse */
         /* elements which have non-empty candidates). The first element*/
         /* of the list, vlist[0], is the local coarse element number.  */

         if (ncand > 0) {

            /* search for the local index in the local coarse grid      */

            index    = 0;
            ggelenum = g_c_grid->global_element[i];
            gcelenum = cgrid_fcns->USR_grid_get_element_global_num(c_grid,
                                                                   index);
            while ( gcelenum != ggelenum && index < ncelmnts )
            {
               index++;
               gcelenum =
                  cgrid_fcns->USR_grid_get_element_global_num(c_grid,index);
            }

/* #if defined(DEBUG) */
            if ( index >= ncelmnts )
            {
               printf(" Error : cannot find element in local grid.\n");
               exit(-1);
            }
/* #endif */

            leng  = ncand;
            ML_memory_alloc((void**) &coord_short,
                            ndim*leng*sizeof(double),"co2");
            icnt  = 0;
            for ( j = 1; j <= ncand; j++ )
            {
               icnt2 = vlist[j];
               fgrid_fcns->USR_grid_get_vertex_coordinate(f_grid, icnt2,
                                               &(coord_short[ndim*icnt]));
               icnt++;
            }

            /* call user function to compute coefficients        */
            /* (the length of coefs and coef_ptr are dynamically */
            /* adjusted to save memory usage.)                   */

            if ((gpcnt + leng) > cur_ptr_leng)
            {
               coef_ptr2    = coef_ptr;
               cur_ptr_leng = gpcnt + 3 * leng;
               ML_memory_alloc((void**) &coef_ptr,
                               cur_ptr_leng*sizeof(int), "cp2");
               for ( j = 0; j < gpcnt; j++ ) coef_ptr[j] = coef_ptr2[j];
               ML_memory_free( (void**) &coef_ptr2);
            }
            if ((gcnt + leng*MAX_VERT_PER_ELE) > cur_coef_leng)
            {
               coefs2  = coefs;
               cur_coef_leng = gcnt + 5 * leng * MAX_VERT_PER_ELE;
               ML_memory_alloc((void**) &coefs,
                               cur_coef_leng*sizeof(double),"ce2");
               for ( j = 0; j < gcnt; j++ ) coefs[j] = coefs2[j];
               ML_memory_free( (void **) &coefs2);
            }

            cgrid_fcns->USR_compute_basis_coefficients(c_grid, index,
                        coord_short, leng, &coefs[gcnt], &coef_ptr[gpcnt]);

            /* need to check coef_ptr to see if any one of the */
            /* candidates got in, and if so, register in the   */
            /* fnode_flag array                                */

            for ( j = 0 ; j < leng; j++ )
            {
               if (coef_ptr[gpcnt+j] > 1) fnode_flag[vlist[j+1]] = mypid;
               coef_ptr[gpcnt+j] += coef_ptr[gpcnt+j-1];
            }

            /* when debug mode is turned on, the processor that get */
            /* the mode will output basis coefficient information   */

            /*if ( ML_debug >= 0 )
            {
               printf(" ===> RP0 : element %d \n", i);
               ML_ElementAGX_Print(element);
               for ( j = 0 ; j < leng; j++ )
               {
                  if ( ndim == 2 )
                     printf("  coord %d = %e %e \n", j, coord_short[j*2],
                                                        coord_short[j*2+1]);
                  else if ( ndim == 3 )
                     printf("  coord %d = %e %e \n", j, coord_short[j*3],
                               coord_short[j*3+1], coord_short[j*3+2]);
                  if (coef_ptr[gpcnt+j] > (coef_ptr[gpcnt+j-1]+1))
                  {
                     for (k=coef_ptr[gpcnt+j-1]; k<coef_ptr[gpcnt+j]; k++)
                        printf(" coef = %e \n", coefs[k]);
                  }
               }
               printf(" <=== \n" );
            }
            */
            ML_memory_free( (void**) &coord_short);

            /* update counters and free temporary storage */

            vlist[0] = index;
            ML_IntList_Load_Sublist(lfv_list, ncand+1, vlist);
            gpcnt = gpcnt + leng;
            gcnt  = coef_ptr[gpcnt-1];

            if (gcnt > cur_coef_leng)
            {
               printf("Error : coefficient array not long enough. \n");
               exit(0);
            }
         }
      }
   }
   ML_ElementAGX_Destroy( &element );
   ML_memory_free( (void**) &vlist);

   /* ----------------------------------------------------------------- */
   /*  Calculate the storage requirements and allocate storage.         */
   /*  The storage requirement for local grid transfer are stored in    */
   /*  in the fine2coarsecnts array. Coefficients for remote processors */
   /*  are counted in extern_node_cnt (node_proc_map is used to find    */
   /*  whether a coarse node is residing within my processor or not).   */
   /* ----------------------------------------------------------------- */

   node_proc_map = g_c_grid->node_proc_map;
   ML_memory_alloc((void**) &fine2coarsecnts,ncvert * sizeof(int),"f2c");
   ML_memory_alloc((void**) &vlist,  MAX_VERT_PER_ELE*sizeof(int),"f2C");
   for ( i = 0; i < ncvert; i++ ) fine2coarsecnts[i] = 0;
   gpcnt = extern_node_cnt = 0;

   for ( i = 0; i < lfv_list->length; i++ )
   {
      k        = lfv_list->start[i];
      ibegin   = k + 1;
      iend     = lfv_list->start[i+1];
      leng     = iend - ibegin;
      cele_num = lfv_list->members[k];
      icnt     = cgrid_fcns->USR_grid_get_element_vlist(c_grid,
                                                        cele_num, vlist);

      for ( m = ibegin; m < iend; m++ )
      {
         mbegin = coef_ptr[gpcnt];
         mend   = coef_ptr[gpcnt+1];
         if ((mend - mbegin) > 1) { /* if nonempty */
            for ( j = mbegin; j < mend; j++ )
            {
               if ( coefs[j] > 1.0E-12 )
               {
                  cnodenum  = vlist[j-mbegin];
                  cgnodenum = cgrid_fcns->USR_grid_get_vertex_global_num(
                                                          c_grid,cnodenum);
                  if ( node_proc_map[cgnodenum] == mypid )
                     fine2coarsecnts[cnodenum]++;
                  else extern_node_cnt++;
               }
            }
         }
         gpcnt++;
      }
   }

   xsfer_op->Nlocal_rows = ncvert;
   nbytes = ( ncvert + 1 ) * sizeof(int);
   ML_memory_alloc((void**) &(xsfer_op->local_ia), nbytes, "XS2");
   ncnt = 0;
   for ( i = 0; i < ncvert; i++ )
   {
      xsfer_op->local_ia[i] = ncnt;
      ncnt += fine2coarsecnts[i];
   }
   xsfer_op->local_ia[ncvert] = ncnt;
   for ( i = 0; i < ncvert; i++ )
      fine2coarsecnts[i] = xsfer_op->local_ia[i];
   if ( ncnt > 0 )
   {
      nbytes = ncnt * sizeof(int);
      ML_memory_alloc((void**) &(xsfer_op->local_ja), nbytes, "XS3");
      nbytes = ncnt * sizeof(double);
      ML_memory_alloc((void**) &(xsfer_op->local_a), nbytes, "XS4");
   }
   else
   {
      xsfer_op->local_ja = NULL;
      xsfer_op->local_a  = NULL;
   }
   if ( extern_node_cnt > 0 )
   {
      icnt = extern_node_cnt;
      nbytes = icnt * sizeof(int);
      ML_memory_alloc((void**) &(xsfer_op->ext_ia), nbytes, "XS5");
      ML_memory_alloc((void**) &(xsfer_op->ext_ja), nbytes, "XS6");
      nbytes = icnt * sizeof(double);
      ML_memory_alloc((void**) &(xsfer_op->ext_a), nbytes, "XS7");
      xsfer_op->ext_cnt = icnt;
   }
   else
   {
      xsfer_op->ext_ia  = NULL;
      xsfer_op->ext_ja  = NULL;
      xsfer_op->ext_a   = NULL;
      xsfer_op->ext_cnt = 0;
   }

   /* ----------------------------------------------------------------- */
   /* fill in the local fields of the grid transfer operator            */
   /* ----------------------------------------------------------------- */

   gpcnt   = 0;
   ext_cnt = 0;
   for ( i = 0; i < lfv_list->length; i++ )
   {
      k        = lfv_list->start[i];
      ibegin   = k + 1;
      iend     = lfv_list->start[i+1];
      leng     = iend - ibegin;
      cele_num = lfv_list->members[k];
      icnt     = cgrid_fcns->USR_grid_get_element_vlist(c_grid,
                                                  cele_num, vlist);

      for ( m = ibegin; m < iend; m++ )
      {
         mbegin    = coef_ptr[gpcnt];
         mend      = coef_ptr[gpcnt+1];
         fnodenum  = lfv_list->members[m];

         if (( mend - mbegin ) > 1) { /* if nonempty */

            for ( j = mbegin; j < mend; j++ )
            {
               if ( coefs[j] > 1.0E-12 )
               {
                  cnodenum  = vlist[j-mbegin];
                  cgnodenum = cgrid_fcns->USR_grid_get_vertex_global_num(
                                                        c_grid,cnodenum);
                  if ( node_proc_map[cgnodenum] == mypid )
                  {
                     index = fine2coarsecnts[cnodenum]++;
                     xsfer_op->local_ja[index] = fnodenum;
                     xsfer_op->local_a[index]  = coefs[j];
                  }
                  else
                  {
                     xsfer_op->ext_ia[ext_cnt]  = cgnodenum;
                     xsfer_op->ext_ja[ext_cnt]  = fnodenum;
                     xsfer_op->ext_a[ext_cnt++] = coefs[j];
                  }
               }
            }
         }
         gpcnt++;
      }
   }

/* #if defined(DEBUG) */
   ncnt = xsfer_op->Nlocal_rows;
   icnt = 0;
   for ( i = 0; i < ncnt; i++ )
   {
      icnt = (fine2coarsecnts[i] - xsfer_op->local_ia[i+1]);
      if ( icnt != 0 )
         printf("Error : in local transfer operator (%d) %d.\n",i,icnt);
   }
/* #endif */
   ML_memory_free( (void **) &coef_ptr );
   ML_memory_free( (void **) &coefs);
   ML_memory_free( (void **) &vlist );
   ML_memory_free( (void **) &fvlist );
   ML_memory_free( (void **) &fine2coarsecnts);
   ML_memory_free( (void **) &coord );
   ML_IntList_Destroy( &lfv_list );
}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ML_remote_grid_candidates : Given a global coarse grid g_c_grid,        */
/*   search the local fine grid f_grid (use the given coordinates) to      */
/*   identify candidate vertices for each of the nonlocal global coarse    */
/*   grid elements. (Note : the fine grid vertices in candidates are given */
/*   local node ordering).  Upon return, the candidates structure contains */
/*   a list of sublists, each of which stores the candidates for a coarse  */
/*   element with its global (since it is referencing a nonlocal element)  */
/*   element number stored at the beginning of the sublist.                */
/* *********************************************************************** */

int ML_remote_grid_candidates(void *f_grid, ML_GridFunc *fgrid_fcns,
                           ML_GridFunc *cgrid_fcns,
                           ML_GridAGX *g_c_grid, ML_IntList *candidates,
                           ML_OperatorAGX  *xsfer_op, ML_Comm *comm)
{
   int    i, ncelmnts, nfvert, ncand, icnt, *fnode_flag;
   int    mypid, *vlist, *fvlist, ndim, nfvert_left, *element_map;
   int    nbytes, MAX_VERT_PER_ELE;
   double *coord;
   ML_ElementAGX *element;

   /* ------------------------------------------------------------------*/
   /* - fetch the number of elements in the global coarse grid          */
   /* - fetch the dimension of the fine grid                            */
   /* - fetch the number of vertices in the local fine grid             */
   /* - fetch the element to processor map                              */
   /* ----------------------------------------------------------------- */

   ncelmnts    = g_c_grid->Nelements;
   if (fgrid_fcns->USR_grid_get_nvertices == NULL)
      pr_error("ML_remote_grid_candidates: USR_grid_get_nvertices() not found\n");
   if (fgrid_fcns->USR_grid_get_dimension == NULL)
      pr_error("ML_remote_grid_candidates: USR_grid_get_dimension() not found\n");
   ndim        = fgrid_fcns->USR_grid_get_dimension(f_grid);
   nfvert      = fgrid_fcns->USR_grid_get_nvertices(f_grid);
   mypid       = comm->ML_mypid;
   element_map = g_c_grid->elmnt_proc_map;
   fnode_flag  = xsfer_op->fnode_flag;
   MAX_VERT_PER_ELE = cgrid_fcns->ML_MaxElmntVert;

   /* ----------------------------------------------------------------- */
   /* count of number of unprocessed local fine nodes (some or all of   */
   /* them may belong to some local coarse elements and thus should not */
   /* be considered any more.) If not more left, just return.           */
   /* ----------------------------------------------------------------- */

   nfvert_left = 0;
   for (i = 0; i < nfvert; i++) if ( fnode_flag[i] == -1 ) nfvert_left++;
   if ( nfvert_left == 0 ) return 0;

   /* ----------------------------------------------------------------- */
   /* - initialize the candidates structure (allocate memory)           */
   /* - initialize a simple element structure to temporarily hold the   */
   /*   information of a given coarse grid element in searching for     */
   /*   candidates                                                      */
   /* - allocate memory for temporarily storing candidate node list     */
   /* - allocate memory for fetching coordinate information             */
   /* - allocate memory for storing unprocessed fine node numbers       */
   /* ----------------------------------------------------------------- */

   ML_ElementAGX_Create(&element, ndim, MAX_VERT_PER_ELE);
   ML_memory_alloc((void**) &vlist, (nfvert_left + 1)*sizeof(int),"vl2");
   nbytes = ndim * nfvert_left * sizeof(double);
   ML_memory_alloc((void**) &coord, nbytes, "co3");
   ML_memory_alloc((void**) &fvlist, nfvert_left*sizeof(int), "fv2");

   /* ----------------------------------------------------------------- */
   /* Collect all node coordinates which have not been processed        */
   /* Note : whether a fine node has been processed or not depends on   */
   /*        whether fnode_flag[i] is -1 or not.                        */
   /* ----------------------------------------------------------------- */

   icnt = 0;
   for ( i = 0; i < nfvert; i++ )
   {
      if ( fnode_flag[i] == -1 )
      {
         fvlist[icnt] = i;
         fgrid_fcns->USR_grid_get_vertex_coordinate(f_grid, i,
                                                &coord[ndim*icnt]);
         icnt++;
      }
   }
   if ( icnt != nfvert_left )
   {
      printf("Error : in ML_remote_grid_candidates \n");
      exit(0);
   }

   /* ------------------------------------------------------------------*/
   /* process the searching for each nonlocal coarse grid element       */
   /* ------------------------------------------------------------------*/

   for ( i = 0; i < ncelmnts; i++ )
   {
      if ( element_map[i] != mypid )
      {
         /* fetch the coarse element into a simple element structure  */

         ML_GridAGX_Get_Element(g_c_grid, i, element);

         /* search fine vertices for candidates for the given element.  */
         /* Upon return, vlist[1...] should contain the indices of the  */
         /* candidates, and ncand contains the number of such vertices) */
         /* fnode_flag is an array to tell the element to ignore those  */
         /* fine nodes with the corresponding entry = -1.               */

         ML_ElementAGX_ComposeCandidates(element, nfvert_left, coord,
                                    fvlist, fnode_flag, &ncand, &vlist[1]);

         /* if a non-empty list is returned, load the candidate info.   */
         /* to the int_lists structure as one row of a CSR matrix (so   */
         /* the number of rows of this list at the end is the no. of    */
         /* coarse elements which have non-empty candidates). The first */
         /* element of the list, vlist[0], is the local coarse element  */
         /* number.                                                     */

         if (ncand > 0)
         {
            vlist[0] = i;
            ML_IntList_Load_Sublist(candidates, ncand+1, vlist);
         }
      }
   }

   /* ----------------------------------------------------------------- */
   /* clean up                                                          */
   /* ----------------------------------------------------------------- */

   ML_ElementAGX_Destroy(&element);
   ML_memory_free((void **) &coord);
   ML_memory_free((void **) &vlist);
   ML_memory_free((void **) &fvlist);
   return 0;
}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ML_exchange_candidates : The local fine grid node candidate list (for   */
/*   every nonlocal coarse grid element) in 'inlist' are exchanged between */
/*   processors (the processors where the sublist should be sent can be    */
/*   identified from the elmnt_proc_map array).  At the end the 'com'      */
/*   variable (the recv_ia, recv_list, and recv_xyz fields) contains the   */
/*   nonlocal fine grid node lists which are candidates for local coarse   */
/*   elements while the 'com' send fields contain the local fine grid node */
/*   list originally stored in lfv_list (so now lfv_list can be disposed   */
/*   of when after return from this subroutine.                            */
/*   - f_grid is passed as an argument to provide coordinate information   */
/*   - g_cgrid is passed to provide the global coarse element numbers      */
/* *********************************************************************** */

void ML_exchange_candidates(ML_IntList *inlist, void *fgrid,
                            ML_GridFunc *fgrid_fcns, ML_GridAGX *g_cgrid,
                            ML_CommInfoAGX *combuf, ML_Comm *comm)
{
   int     i, j, k1, k2, leng, nprocs, mypid, proc_id, fromproc;
   int     sendproc_cnt, recvproc_cnt, tot_recv_leng, msgtype, index= 0, length;
   int     *send_proc, *send_leng, *itmp, *proc_flag;
   ml_big_int **send_list, elenum, big_index=0;
   int     *recv_proc, *recv_leng, tot_send_leng, *trecv_proc;
   ml_big_int **recv_list;
   int     ndim, *trecv_msg, *trecv_leng, *elmnt_proc_map, nbytes, *intarray;
   double  *xyz;
   USR_REQ *Request;

   /* ----------------------------------------------------------------- */
   /* initialize an array proc_flag which is used for computing the     */
   /* communication patterns and lengths between processors             */
   /* ----------------------------------------------------------------- */

   mypid     = comm->ML_mypid;
   nprocs    = comm->ML_nprocs;
   ML_memory_alloc((void**) &proc_flag, nprocs * sizeof(int), "pf1" );
   for ( i = 0; i < nprocs; i++ ) proc_flag[i] = 0;
   combuf->proc_id = mypid;

   /* ----------------------------------------------------------------- */
   /* fetch the element to processor map                                */
   /* ----------------------------------------------------------------- */

   ndim = fgrid_fcns->USR_grid_get_dimension(fgrid);
   elmnt_proc_map = g_cgrid->elmnt_proc_map;

   /* ----------------------------------------------------------------- */
   /* now search the candidate fine node list to find the information   */
   /* about which processors requests are to be sent and length of send */
   /* (recall the first entry of each of the sublist contains the local */
   /* coarse element number and elmnt_proc_map contains the originating */
   /* processor of the corresponding coarse element.) At the end,       */
   /* proc_flag[i] contains the number of integers to be sent to        */
   /* processor i.                                                      */
   /* ----------------------------------------------------------------- */

   leng      = inlist->length;
   for ( i = 0; i < leng; i++ )
   {
      k1      = inlist->start[i];    /* get the element number pointer */
      k2      = inlist->members[k1]; /* fetch the local element index  */
      proc_id = elmnt_proc_map[k2];  /* locate where the element is */
      if (proc_id != mypid)          /* if not local, raise flag */
         proc_flag[proc_id] += (inlist->start[i+1] - k1);
   }

   /* ----------------------------------------------------------------- */
   /* Find out how many processors data are to be sent (the number of   */
   /* entries in proc_flag that are non-zero) and the total data length */
   /* to be sent to each processor.  Use this information to initialize */
   /* the communication  buffers. (Separate arrays, e.g. send_proc and  */
   /* send_list, are allocated to more conveniently set up the send     */
   /* information.  This information will be loaded to the send fields  */
   /* of the 'com' variable later.)                                     */
   /* ----------------------------------------------------------------- */

   sendproc_cnt = 0; tot_send_leng = 0;
   for ( i = 0; i < nprocs; i++)
   {
      if (proc_flag[i] != 0)
      {
         sendproc_cnt++;
         tot_send_leng += proc_flag[i];
      }
   }

   ML_CommInfoAGX_Setup_Send(combuf, sendproc_cnt, tot_send_leng);
   if ( sendproc_cnt > 0 )
   {
      nbytes = sendproc_cnt * sizeof(int);
      ML_memory_alloc((void**) &send_proc, nbytes, "sp1" );
      ML_memory_alloc((void**) &send_leng, nbytes, "sl1" );
      ML_memory_alloc((void**) &send_list, sendproc_cnt*sizeof(ml_big_int *), "st1" );
   }
   else
   {
      send_proc = NULL;
      send_leng = NULL;
      send_list = NULL;
   }
   k2 = 0;
   for ( i = 0; i < nprocs; i++ )
   {
      k1 = proc_flag[i];
      if ( k1 > 0)
      {
         send_proc[k2]   = i;
         send_leng[k2]   = 0;
         ML_memory_alloc((void**) &(send_list[k2]),k1*sizeof(ml_big_int),"sl3");
         k2++;
      }
   }

   /* ----------------------------------------------------------------- */
   /* 'sendproc_cnt' is the number of processors that data are to be    */
   /* sent for querying, 'send_proc' contains the processor numbers,    */
   /* and 'send_list' is an 2-D array such that send_list[i] has the    */
   /* coarse element to fine node candidates indices.  This 'send_list' */
   /* is to be composed next.                                           */
   /* - for each sublist in the candidate list (each sublist            */
   /*   corresponds to the node list of a coarse element)               */
   /* - encode the element number into a negative number for more       */
   /*   convenient identification later                                  */
   /* - if the coarse element does not reside locally, load the sublist  */
   /*   to the send list                                                 */
   /* ----------------------------------------------------------------- */

   leng = inlist->length;
   for ( i = 0; i < leng; i++ )
   {
      k1      = inlist->start[i];
      k2      = inlist->members[k1];
      proc_id = elmnt_proc_map[k2];
      elenum  = g_cgrid->global_element[k2];
      if ( proc_id != mypid )
      {
         for ( j = 0; j < sendproc_cnt; j++ )
            if ( proc_id == send_proc[j] ) { index = j; break; }
         k1++;
         k2 = inlist->start[i+1];

         /* convert the global element number to a negative number */

         send_list[index][send_leng[index]++] = - (elenum + 1);

         /* put the fine vertex candidates on the send list */

         for ( j = k1; j < k2; j++)
            send_list[index][send_leng[index]++] = inlist->members[j];
      }
      else
      {
         printf("Error : this should never have been reached.\n");
         exit(0);
      }
   }

   /* ----------------------------------------------------------------- */
   /* load the send information to the communication buffer             */
   /* ----------------------------------------------------------------- */

   for ( i = 0; i < sendproc_cnt; i++ )
      ML_CommInfoAGX_Load_SendList(combuf, send_proc[i], send_leng[i],
                                 send_list[i]);

   /* ----------------------------------------------------------------- */
   /* Now the send part of the communication buffer has been composed.  */
   /* Next, generate the receive part of the communication buffer.      */
   /* First, each processor sets proc_flag[i]=1 if it needs to send to  */
   /* processor i.  Then, this proc_flag is reduced (i.e. add) among    */
   /* all processors so that at each end, each processor will know how  */
   /* many processors it will receive requests from (by examining the   */
   /* entry proc_flag[mypid].)  The buffers for receive information     */
   /* will then be allocated (in recv_proc, recv_leng, and recv_list).  */
   /* Again these external (to com variable) recv variables are used    */
   /* for purposes of convenience.  These information will then be      */
   /* loaded to the recv fields of the 'com' variable.                  */
   /* ----------------------------------------------------------------- */

   for ( i = 0; i < nprocs; i++ ) if (proc_flag[i] != 0) proc_flag[i] = 1;
   ML_memory_alloc( (void **) &itmp, nprocs * sizeof(int), "itm" );
   /*ML_Comm_GsumVecInt(comm, proc_flag, itmp, nprocs );*/
   ML_gsum_vec_int(&proc_flag, &itmp, nprocs, comm );
   k1 = proc_flag[mypid];
   if ( k1 > 0 )
   {
      nbytes = k1 * sizeof(int);
      ML_memory_alloc( (void **) &recv_proc, nbytes, "rp1" );
      ML_memory_alloc( (void **) &recv_leng, nbytes, "rl1" );
      ML_memory_alloc( (void **) &recv_list, k1*sizeof(ml_big_int*), "rt1" );
   }
   else
   {
      recv_proc = NULL;
      recv_leng = NULL;
      recv_list = NULL;
   }

   /* ----------------------------------------------------------------- */
   /* Each processor sends a message to tell its destination processors */
   /* about its intention to send data over and the data length.  It    */
   /* then listens, well aware of how many messages it will receive,    */
   /* based on the data in proc_flag[mypid], to incoming messages and   */
   /* registers the incoming data sources and the integer data (which   */
   /* are the lengths of data to be received later. At the end whenever */
   /* itmp[i] is not equal to 0, processor i is a source processor and  */
   /* itmp[i] contains the length of data to be received for queries.   */
   /* ----------------------------------------------------------------- */

   msgtype = 478;
   for ( i = 0; i < nprocs; i++ ) itmp[i] = 0;
   if ( proc_flag[mypid] > 0 )
   {
      k1 = proc_flag[mypid] * sizeof(USR_REQ);
      ML_memory_alloc( (void **) &Request, k1, "ru1" );
      k1 = proc_flag[mypid] * sizeof(int);
      ML_memory_alloc( (void **) &intarray, k1, "rv1" );
   }
   else
   {
      Request  = NULL;
      intarray = NULL;
   }
   for ( i = 0; i < proc_flag[mypid]; i++)
   {
      fromproc = -1;
      comm->USR_irecvbytes((void*) &intarray[i], sizeof(int), &fromproc,
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
   }
   for ( i = 0; i < sendproc_cnt; i++ )
   {
      comm->USR_sendbytes((void*) &send_leng[i], sizeof(int), send_proc[i],
                          msgtype, comm->USR_comm);
   }
   for ( i = 0; i < proc_flag[mypid]; i++)
   {
      fromproc = -1;
      comm->USR_cheapwaitbytes((void*) &intarray[i], sizeof(int), &fromproc, 
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      itmp[fromproc] = intarray[i];
   }
   if ( proc_flag[mypid] > 0 )
   {
      ML_memory_free((void**) &Request );
      ML_memory_free((void**) &intarray );
      Request  = NULL;
      intarray = NULL;
   }

   /* ----------------------------------------------------------------- */
   /* Compress the information in itmp into the recv_leng and recv_proc */
   /* arrays.  Then set up the communication buffers (using the         */
   /* ML_CommInfoAGX_Setup_Recv function) and indexing information      */
   /* (using ML_CommInfoAGX_Load_RecvInfo) for receiving data.          */
   /* ----------------------------------------------------------------- */

   recvproc_cnt = tot_recv_leng = 0;
   for ( i = 0; i < nprocs; i++ )
   {
      if ( itmp[i] > 0 )
      {
         recv_proc[recvproc_cnt]   = i;
         recv_leng[recvproc_cnt++] = itmp[i];
         tot_recv_leng += itmp[i];
      }
   }

   ML_CommInfoAGX_Setup_Recv(combuf, recvproc_cnt, tot_recv_leng);
   for ( i = 0; i < recvproc_cnt; i++ )
      ML_CommInfoAGX_Load_RecvInfo(combuf, recv_proc[i], recv_leng[i]);

   /* ----------------------------------------------------------------- */
   /* Communicate the nonlocal fine grid node information. The incoming */
   /* data are loaded to com.recv_list.                                 */
   /* ----------------------------------------------------------------- */

   if ( recvproc_cnt > 0 )
   {
      nbytes = recvproc_cnt * sizeof(int);
      ML_memory_alloc( (void**) &trecv_proc, nbytes, "tr1" );
      ML_memory_alloc( (void**) &trecv_leng, nbytes, "tl1" );
      ML_memory_alloc( (void**) &trecv_msg,  nbytes, "tm1" );
      ML_memory_alloc( (void**) &Request, recvproc_cnt*sizeof(USR_REQ),"rv1");
   }
   else
   {
      trecv_proc = NULL;
      trecv_leng = NULL;
      trecv_msg  = NULL;
      Request    = NULL;
   }

   msgtype = 479;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      ML_CommInfoAGX_Get_RecvList(combuf,i,&proc_id,&leng,&(recv_list[i]));
      trecv_msg[i]  = 479;
      trecv_proc[i] = proc_id;
      trecv_leng[i] = leng * sizeof(recv_list[0][0]);
      comm->USR_irecvbytes((void*) (recv_list[i]), trecv_leng[i],
                           &(trecv_proc[i]), &(trecv_msg[i]),
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }
   for ( i = 0; i < sendproc_cnt; i++ )
   {
      leng = send_leng[i] * sizeof(send_list[0][0]);
      comm->USR_sendbytes((void*) (send_list[i]), leng, send_proc[i],
                          msgtype, comm->USR_comm );
   }
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      comm->USR_cheapwaitbytes((void*) (recv_list[i]), trecv_leng[i], 
                          &(trecv_proc[i]), &(trecv_msg[i]), 
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }

   /* ----------------------------------------------------------------- */
   /* In addition to the node index information, the actual coordinate  */
   /* information has also to be communicated.                          */
   /* ----------------------------------------------------------------- */

   /* ensure the xyz buffer is large enough to hold the send data */
   /* (search for maximum length among all communication) */

   k2      = 0;
   msgtype = 480;
   for ( i = 0; i < sendproc_cnt; i++ )
      if (send_leng[i] > k2) k2 = send_leng[i];
   k2  = k2 * ndim * sizeof(double);
   if ( k2 > 0 ) ML_memory_alloc( (void**) &xyz, k2, "xyz" );
   else          xyz = NULL;

   index = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      trecv_msg[i]  = 480;
      trecv_leng[i] = ndim * recv_leng[i] * sizeof(double);
      trecv_proc[i] = recv_proc[i];
      comm->USR_irecvbytes((void*)&(combuf->recv_xyz[index]), trecv_leng[i],
                           &(trecv_proc[i]), &(trecv_msg[i]),
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
      index = index + ndim * recv_leng[i];
   }

   for ( i = 0; i < sendproc_cnt; i++ )
   {
      length = send_leng[i] * sizeof(double);
      for ( j = 0; j < send_leng[i]; j++ )
      {
         big_index      = send_list[i][j];
         if (big_index >= 0)
         {
            fgrid_fcns->USR_grid_get_vertex_coordinate(fgrid, big_index,
                                                      &(xyz[j*ndim]));
         }
         else
         {
            xyz[ndim*j] = xyz[ndim*j+1] = 0.0;
            if (ndim > 2) xyz[ndim*j+2] = 0.0;
         }
      }
      comm->USR_sendbytes((void*) xyz, ndim*length, send_proc[i], msgtype,
                          comm->USR_comm );
   }

   index = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      comm->USR_cheapwaitbytes((void*)&(combuf->recv_xyz[index]), trecv_leng[i],
                           &(trecv_proc[i]), &(trecv_msg[i]), 
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
      index = index + ndim * recv_leng[i];
   }

   /* ----------------------------------------------------------------- */
   /* Finally, clean up                                                 */
   /* ----------------------------------------------------------------- */

   if ( xyz        != NULL ) ML_memory_free((void**) &xyz);
   if ( trecv_leng != NULL ) ML_memory_free((void**) &trecv_leng);
   if ( trecv_proc != NULL ) ML_memory_free((void**) &trecv_proc);
   if ( trecv_msg  != NULL ) ML_memory_free((void**) &trecv_msg);
   if ( Request    != NULL ) ML_memory_free((void**) &Request);
   if ( send_proc  != NULL ) ML_memory_free((void**) &send_proc);
   if ( send_leng  != NULL ) ML_memory_free((void**) &send_leng);
   if ( send_list  != NULL )
   {
      for ( i = 0; i < sendproc_cnt; i++ )
         ML_memory_free( (void**) &(send_list[i]) );
      ML_memory_free((void**) &send_list);
   }
   if ( recv_proc != NULL ) ML_memory_free((void **) &recv_proc);
   if ( recv_leng != NULL ) ML_memory_free((void **) &recv_leng);
   if ( recv_list != NULL ) ML_memory_free((void **) &recv_list);
   ML_memory_free((void **) &proc_flag);
   ML_memory_free((void **) &itmp);
}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ML_get_basis_functions_coef : Given a list of nodes and coordinate      */
/*    information in the 'com' variable, generate the grid transfer        */
/*    coefficients in xsfer_op ext2 variables.                             */
/* *********************************************************************** */

void ML_get_basis_functions_coef(ML_CommInfoAGX *com, void *c_grid,
                                 ML_GridFunc *cgrid_fcns,
                                 ML_OperatorAGX  *xsfer_op)
{
   int    j, m, length, mcnt, ncnt, index, nodecnt, cur_proc;
   int    cur_coef_leng, ndim, ncelmnts, *coef_ptr;
   ml_big_int celenum, gelenum, *ev_list;
   int    cur_leng, *track_array=NULL, track_leng, prev_track_leng;
   int    MAX_VERT_PER_ELE;
   double *coord, *coefs, *coefs2;

   /* ----------------------------------------------------------------- */
   /* fetch the incoming information (element to node list, coordinate  */
   /* information).                                                     */
   /* ----------------------------------------------------------------- */

   length  = com->recv_ia[com->recv_cnt];
   ev_list = com->recv_list;
   coord   = com->recv_xyz;

   /* ----------------------------------------------------------------- */
   /* allocate space for the coefficients and its index pointers        */
   /* (In order to minimize memory usage, the coefficient array is      */
   /* initially a short double array.  Its length will be incrementally */
   /* increased based on need.  'cur_coef_leng' holds the current       */
   /* length of 'coefs'.)                                               */
   /* ----------------------------------------------------------------- */

   ncelmnts  = cgrid_fcns->USR_grid_get_nelements( c_grid );
   ndim      = cgrid_fcns->USR_grid_get_dimension( c_grid );
   ML_memory_alloc((void**) &coef_ptr,(length + 1) * sizeof(int), "cp3");
   MAX_VERT_PER_ELE = cgrid_fcns->ML_MaxElmntVert;
   cur_coef_leng = length * MAX_VERT_PER_ELE + 1;
   ML_memory_alloc((void**) &coefs, cur_coef_leng*sizeof(double), "ce3");
   coef_ptr[0] = 0;

   /* ----------------------------------------------------------------- */
   /* The list consists of integer lists of node numbers.  Each sublist */
   /* is preceded by its coarse element number(with its number negated).*/
   /* ----------------------------------------------------------------- */

   mcnt     = 0;
   ncnt     = 0;
   cur_proc = 0;
   if ( com->recv_cnt > 0 )
   {
      cur_leng = com->recv_ia[1] - com->recv_ia[0];
      if ( cur_leng > 0 )
         ML_memory_alloc((void**) &track_array,cur_leng*sizeof(int),"ta1");
      else
         track_array = NULL;
      track_leng = 0;
   }
   while ( mcnt < length )
   {
      /* these array are used to track duplicates */

      if ( mcnt >= com->recv_ia[cur_proc+1] )
      {
         cur_proc++;
         cur_leng = com->recv_ia[cur_proc+1] - com->recv_ia[cur_proc];
         if ( track_array != NULL ) ML_memory_free( (void**) &track_array );
         if ( cur_leng > 0 )
            ML_memory_alloc((void**) &track_array,cur_leng*sizeof(int),"ta2");
         else
            track_array = NULL;
         track_leng = 0;
      }

      /* fetch the global element number */

      celenum = - ev_list[mcnt] - 1;

      /* find out its relative position in the local coarse grid */

      index  = 0;
      gelenum = cgrid_fcns->USR_grid_get_element_global_num(c_grid, index);
      while ( gelenum != celenum && index < ncelmnts )
      {
         index++;
         gelenum = cgrid_fcns->USR_grid_get_element_global_num(c_grid, index);
      }
/* #if defined(DEBUG) */
      if ( index >= ncelmnts ) {
         printf(" Error : cannot find element in local grid.\n");
         exit(-1);
      }
/* #endif */

      /* since this is an element number, no coefficient is assigned */

      coef_ptr[mcnt+1] = coef_ptr[mcnt];

      /* count the number of node candidates for this element */

      mcnt++;
      nodecnt = mcnt;
      while ( nodecnt < length && ev_list[nodecnt] >= 0 ) nodecnt++;
      nodecnt = nodecnt - mcnt;

      /* check to see if the coefs array is long enough to hold extra */
      /* data, and if not, reallocate a longer array.                 */
      /* (Note : MAX_VERT_PER_ELE should be changed according to the  */
      /*         types of finite element used.)                       */

      if ((ncnt + nodecnt * MAX_VERT_PER_ELE) > cur_coef_leng)
      {
         coefs2 = coefs;
         cur_coef_leng = ncnt + 5 * nodecnt * MAX_VERT_PER_ELE;
         ML_memory_alloc((void**) &coefs, cur_coef_leng*sizeof(double),"ce4");
         for ( j = 0; j < ncnt; j++ ) coefs[j] = coefs2[j];
         ML_memory_free( (void**) &coefs2);
      }

      /* call user function to generate coefficients */
      /* for the current set of node coordinates.    */

      cgrid_fcns->USR_compute_basis_coefficients(c_grid, index,
                        &coord[ndim*mcnt], nodecnt, &coefs[ncnt],
                        &coef_ptr[mcnt+1]);

      /* To interface to the user function properly,    */
      /* need to adjust the coef_ptr array accordingly. */

      for ( j = 0 ; j < nodecnt; j++ )
         coef_ptr[mcnt+j+1] += coef_ptr[mcnt+j];

      /* the processor that gets the debug flag will print out */
      /* basis coefficient information */

      /*if ( ML_debug >= 0 )
      {
         printf(" ===> (2) \n");
         printf(" Local coarse element = %d \n", index );
         for ( j = 0 ; j < nodecnt; j++ )
         {
            if ( ndim == 2 )
               printf("coord = %e %e\n",coord[(mcnt+j)*2],coord[(mcnt+j)*2+1]);
            else if ( ndim == 3 )
               printf("coord = %e %e\n",coord[(mcnt+j)*3],coord[(mcnt+j)*3+1],
                                        coord[(mcnt+j)*3+2]);
            if (coef_ptr[mcnt+j+1] > (coef_ptr[mcnt+j]+1))
            {
               for ( k = coef_ptr[mcnt+j] ; k < coef_ptr[mcnt+j+1]; k++ )
                  printf(" coef = %e \n", coefs[k]);
            }
         }
         printf(" <=== (2) \n");
      }
      */

      /* see if each fine node lies inside the coarse element */
      /* and if so, the fine node is said to be done.         */

      for ( j = 0 ; j < nodecnt; j++ )
      {
         if ((coef_ptr[mcnt+j+1] - coef_ptr[mcnt+j]) > 1 )
         {
            prev_track_leng = track_leng;
            ML_search_insert_sort( ev_list[mcnt+j], track_array,
                                   &track_leng, 0);
            if ( prev_track_leng == track_leng )
            {
               for ( m = coef_ptr[mcnt+j]; m < coef_ptr[mcnt+j+1]; m++ )
                  coefs[m] = 0.0;
            }
         }
      }

      /* finally update the counters for the next processing */

      mcnt = mcnt + nodecnt;
      ncnt = coef_ptr[mcnt];
   }

   /* store away the coefficients and their pointers to */
   /* the grid transfer object */

   xsfer_op->ext2_cnt = length;
   xsfer_op->ext2_ptr = coef_ptr;
   xsfer_op->ext2_a   = coefs;

   /* delete the coordinate information to conserve memory */

   if ( com->recv_xyz != 0 )
   {
      ML_memory_free( (void**) &com->recv_xyz );
      com->recv_xyz = 0;
   }
   if ( track_array != NULL ) ML_memory_free( (void**) &track_array );
}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ML_exchange_coefficients : shuffle back the interpolating coefficients  */
/*                            to the originating processor.                */
/*    (The processor information is stored in the communication buffer)    */
/* *********************************************************************** */

void ML_exchange_coefficients(void *c_grid, ML_GridFunc *cgrid_fcns,
                              ML_CommInfoAGX *combuf,
                              ML_OperatorAGX  *xsfer_op, ML_Comm *comm)
{
   int     i, j, k, sendproc_cnt, recvproc_cnt, tot_recv_leng;
   int     msgtype, *send_proc, *send_leng, *recv_proc, *recv_leng;
   int     leng, offset, *coefp_out, *ind_data;
   int     index, *coef_ptr, begin, end, icnt;
   int     *vlist, *index_out, ncelmnts, MAX_VERT_PER_ELE;
   ml_big_int *evlist, celenum, gelenum;
   double  *coef_out, *coefs;
   USR_REQ *Request;

   /* ----------------------------------------------------------------- */
   /* fetch the coefficients and their pointers from the grid transfer  */
   /* data structure                                                    */
   /* ----------------------------------------------------------------- */

   coef_ptr = xsfer_op->ext2_ptr;
   coefs    = xsfer_op->ext2_a;
   ncelmnts = cgrid_fcns->USR_grid_get_nelements( c_grid );
   MAX_VERT_PER_ELE = cgrid_fcns->ML_MaxElmntVert;
   ML_memory_alloc((void**) &vlist, MAX_VERT_PER_ELE*sizeof(int),"vl3");

   /* ----------------------------------------------------------------- */
   /* All the local coarse node indices to be sent back to their        */
   /* originating processors have to be converted to global coarse      */
   /* node indices first before sending.  This step is absolutely       */
   /* necessary for two reasons :                                       */
   /* 1. since the receiving processor needs some way to distinguish    */
   /*    between all incoming (remote) coarse node coefficients to be   */
   /*    able to inject the proper remote coefficients into the grid    */
   /*    transfer operator (the remote part)                            */
   /* 2. since the returned coarse node may not reside in the sending   */
   /*    processor, this is needed to identify the correct processor.   */
   /* ----------------------------------------------------------------- */

   /* the ind_data is the node list with an 1-to-1 correspondence with  */
   /* the coefs array.  This node list has to be filled before sending  */

   evlist = combuf->recv_list;

   sendproc_cnt = combuf->recv_cnt;
   leng         = combuf->recv_ia[sendproc_cnt];
   if ( leng > 0 )
   {
      icnt = coef_ptr[leng];
      if ( icnt > 0 )
         ML_memory_alloc((void**) &ind_data, icnt*sizeof(int), "id1");
   }
   else
   {
      ind_data = NULL;
   }
   icnt         = 0;

   while ( icnt < leng )
   {

      celenum = - evlist[icnt] - 1;
      index  = 0;
      gelenum = cgrid_fcns->USR_grid_get_element_global_num(c_grid, index);
      while ( gelenum != celenum && index < ncelmnts )
      {
         index++;
         gelenum = cgrid_fcns->USR_grid_get_element_global_num(c_grid,index);
      }
      k = cgrid_fcns->USR_grid_get_element_vlist( c_grid, index, vlist );
      for ( j = 0; j < k; j++ )
         vlist[j] = cgrid_fcns->USR_grid_get_vertex_global_num(c_grid,
                                                               vlist[j]);
      icnt++;

      /* return coefficient = 0 if the coarse grid */
      /* point does not reside in my processor */

      while ( evlist[icnt] >= 0 && icnt < leng )
      {
         index = coef_ptr[icnt];
         if ( ( coef_ptr[icnt+1] - index ) > 1 )
            for ( j = 0; j < k; j++ ) ind_data[index+j] = vlist[j];
         icnt++;
      }
   }

   /* ----------------------------------------------------------------- */
   /* send the coefficients to all others processors                    */
   /* ----------------------------------------------------------------- */

   /* First set up the send information                               */
   /*  - compute the length of data to be sent to other processors    */
   /*    (each fine grid node may have 1 to several coefficients      */
   /*     returned)                                                   */

   send_proc = combuf->recv_proc;
   if (sendproc_cnt > 0)
      ML_memory_alloc((void**) &send_leng, sendproc_cnt*sizeof(int), "sl4");
   else
      send_leng = NULL;

   for ( i = 0; i < sendproc_cnt; i++ )
   {
      begin = combuf->recv_ia[i];
      end   = combuf->recv_ia[i+1];
      send_leng[i] = coef_ptr[end] - coef_ptr[begin];
   }

   /* Then set up the send information                               */
   /*   - to find out the receive length, it requires communicating  */
   /*     to the sending processors.                                 */

   recvproc_cnt = combuf->send_cnt;
   recv_proc    = combuf->send_proc;
   if (recvproc_cnt > 0)
   {
      ML_memory_alloc((void**) &recv_leng,recvproc_cnt*sizeof(int),"rl3");
      ML_memory_alloc((void**) &Request,recvproc_cnt*sizeof(USR_REQ),"ro3");
   }
   else
   {
      recv_leng = NULL;
      Request   = NULL;
   }
   tot_recv_leng = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      msgtype = 13579;
      comm->USR_irecvbytes((void*) &recv_leng[i], sizeof(int), &recv_proc[i],
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
   }
   msgtype = 13579;
   for ( i = 0; i < sendproc_cnt; i++ )
   {
      comm->USR_sendbytes((void*) &send_leng[i], sizeof(int), send_proc[i],
                          msgtype, comm->USR_comm );
   }
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      comm->USR_cheapwaitbytes((void*) &recv_leng[i], sizeof(int), &recv_proc[i],
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      tot_recv_leng += recv_leng[i];
   }

   /* Finally send the coefficients and indices over */

   if (tot_recv_leng > 0)
   {
      ML_memory_alloc((void**) &coef_out,tot_recv_leng*sizeof(double),"co4");
      ML_memory_alloc((void**) &index_out, tot_recv_leng*sizeof(int), "io2");
   }
   else
   {
      coef_out  = NULL;
      index_out = NULL;
   }

   offset = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      msgtype = 13580;
      leng = recv_leng[i] * sizeof(double);
      comm->USR_irecvbytes((void*) &coef_out[offset], leng, &recv_proc[i],
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      offset = offset + recv_leng[i];
   }
   offset  = 0;
   msgtype = 13580;
   for ( i = 0; i < sendproc_cnt; i++ )
   {
      leng = send_leng[i] * sizeof(double);
      comm->USR_sendbytes((void*) &coefs[offset], leng, send_proc[i],
                          msgtype, comm->USR_comm );
      offset = offset + send_leng[i];
   }
   offset = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      msgtype = 13580;
      leng = recv_leng[i] * sizeof(double); 
      comm->USR_cheapwaitbytes((void*) &coef_out[offset], leng, &recv_proc[i], 
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      offset = offset + recv_leng[i];
   }

   offset = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      msgtype = 13581;
      leng = recv_leng[i] * sizeof(int);
      comm->USR_irecvbytes((void*) &index_out[offset], leng, &recv_proc[i],
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      offset = offset + recv_leng[i];
   }
   offset  = 0;
   msgtype = 13581;
   for ( i = 0; i < sendproc_cnt; i++ )
   {
      leng = send_leng[i] * sizeof(int);
      comm->USR_sendbytes((void*) &ind_data[offset], leng, send_proc[i],
                          msgtype, comm->USR_comm );
      offset = offset + send_leng[i];
   }
   offset = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      msgtype = 13581;
      leng = recv_leng[i] * sizeof(int); 
      comm->USR_cheapwaitbytes((void*) &index_out[offset], leng, &recv_proc[i], 
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      offset = offset + recv_leng[i];
   }
   if ( recvproc_cnt > 0 )
   {
      ML_memory_free( (void**) &Request);
   }

   /* ----------------------------------------------------------------- */
   /* send the coefficient pointers to all others processors            */
   /* (coefficient pointers are needed to group coefficients with       */
   /*  vertices)                                                        */
   /* ----------------------------------------------------------------- */

   /* First set up the send information and receive information         */

   sendproc_cnt  = combuf->recv_cnt;
   send_proc     = combuf->recv_proc;
   for ( i = 0; i < sendproc_cnt; i++ )
   {
      send_leng[i] = combuf->recv_ia[i+1] - combuf->recv_ia[i] + 1;
   }
   recvproc_cnt  = combuf->send_cnt;
   recv_proc     = combuf->send_proc;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      recv_leng[i] = combuf->send_ia[i+1] - combuf->send_ia[i] + 1;
   }
   tot_recv_leng = 0;
   for ( i = 0; i < recvproc_cnt; i++ ) tot_recv_leng += recv_leng[i];
   if (tot_recv_leng > 0)
   {
      ML_memory_alloc((void**) &coefp_out,tot_recv_leng*sizeof(int),"co5");
   }
   else
   {
      coefp_out = NULL;
   }

   /* Then send the coefficient pointers */

   if ( recvproc_cnt > 0 )
   {
      ML_memory_alloc((void**) &Request,recvproc_cnt*sizeof(USR_REQ),"cv5");
   }
   else
   {
      Request = NULL;
   }

   offset = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      msgtype = 13582;
      leng = recv_leng[i] * sizeof(int);
      comm->USR_irecvbytes((char*) &coefp_out[offset], leng, &recv_proc[i],
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      offset = offset + recv_leng[i];
   }
   msgtype = 13582;
   offset  = 0;
   for ( i = 0; i < sendproc_cnt; i++ )
   {
      leng = send_leng[i] * sizeof(int);
      comm->USR_sendbytes((void*) &coef_ptr[offset], leng, send_proc[i],
                          msgtype, comm->USR_comm );
      offset = offset + send_leng[i] - 1;
   }
   offset = 0;
   for ( i = 0; i < recvproc_cnt; i++ )
   {
      msgtype = 13582;
      leng = recv_leng[i] * sizeof(int); 
      comm->USR_cheapwaitbytes((char*) &coefp_out[offset], leng, &recv_proc[i], 
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      offset = offset + recv_leng[i];
   }

   /* ----------------------------------------------------------------- */
   /* cleaning up (the setting of pcoefs pointer must be done after the */
   /* rearranging the communication buffer)                             */
   /* ----------------------------------------------------------------- */

   if ( Request != NULL ) ML_memory_free( (void**) &Request);
   ML_memory_free( (void**) &combuf->recv_list);
   ML_memory_free( (void**) &send_leng);
   ML_memory_free( (void**) &recv_leng);
   ML_memory_free( (void**) &ind_data);
   ML_memory_free( (void**) &coefs);
   ML_memory_free( (void**) &coef_ptr);
   ML_memory_free( (void**) &vlist);

   /* ----------------------------------------------------------------- */
   /* finally store the coefficients and pointer to xsfer               */
   /* ----------------------------------------------------------------- */

   xsfer_op->ext2_cnt   = tot_recv_leng;
   xsfer_op->ext2_ptr   = coefp_out;
   xsfer_op->ext2_a     = coef_out;
   xsfer_op->ext2_index = index_out;
}

/* *********************************************************************** */
/* *********************************************************************** */
/* *********************************************************************** */
/* ML_construct_RP1 : use the coefficients stored in temporary arrays in   */
/*    xsfer_op to construct the send and receive fields of the grid        */
/*    transfer operator.                                                   */
/* *********************************************************************** */

void ML_construct_RP1(void *f_grid, ML_GridFunc *fgrid_fcns,
                      void *c_grid, ML_GridFunc *cgrid_fcns,
                      ML_GridAGX *g_c_grid, ML_CommInfoAGX *combuf,
                      ML_OperatorAGX  *xsfer_op, ML_Comm *comm)
{
   int     i, j, k, m, ext_cnt, *ext_ia, *ext_ja, *ext2_ptr;
   int     *ext2_ind, mypid, nprocs, fnvert, cnvert, gcnvert, *node_proc_map;
   int     *cv_leng, *cv_cur_leng, **cv_list, **cv_cnt, *cv_ia, colind;
   ml_big_int fvnum;
   int     com_cnt, cvnum, pnum, *tcv_list, begin, end, rownum, lcnt;
   int     index, proc_cnt, colbase, colend, colcnt, *tmp_ia, *proc_flag;
   int     *recv_list32;
   int     int_size=sizeof(int), msgtype, fromsize, fromproc;
   int     *fnode_flag, err_flag, ndim, nbytes, *intarray;
   double  *ext2_a, *ext_a, coord[3];
   USR_REQ *Request;
   ML_Operator *oper;

   /* ------------------------------------------------------------*/
   /* fetch the data to construct the grid transfer operator      */
   /* ------------------------------------------------------------*/

   ext_cnt   = xsfer_op->ext_cnt;
   ext_a     = xsfer_op->ext_a;
   ext_ia    = xsfer_op->ext_ia;
   ext_ja    = xsfer_op->ext_ja;
   ext2_a    = xsfer_op->ext2_a;
   ext2_ptr  = xsfer_op->ext2_ptr;
   ext2_ind  = xsfer_op->ext2_index;

   /* ------------------------------------------------------------*/
   /* fetch processor, grid, and other information                */
   /* ------------------------------------------------------------*/

   mypid   = comm->ML_mypid;
   nprocs  = comm->ML_nprocs;
   xsfer_op->proc_id = mypid;
   xsfer_op->num_procs = nprocs;
   fnvert  = fgrid_fcns->USR_grid_get_nvertices(f_grid);
   cnvert  = cgrid_fcns->USR_grid_get_nvertices(c_grid);
   gcnvert = ML_Comm_GsumInt( comm, cnvert );
   node_proc_map = g_c_grid->node_proc_map;
   fnode_flag    = xsfer_op->fnode_flag;

   /* ------------------------------------------------------------*/
   /* allocate storage for bookkeeping                            */
   /* - cv_leng : data length for each processor                  */
   /* - cv_cur_leng : current allocated length for cv_list        */
   /* - cv_list : coarse indices for each processor               */
   /* - cv_cnt : no. of times coarse indices are accessed         */
   /* ------------------------------------------------------------*/

   ML_memory_alloc((void**) &cv_leng, nprocs * int_size, "cv1" );
   for ( i = 0; i < nprocs; i++ ) cv_leng[i] = 0;
   ML_memory_alloc((void**) &cv_cur_leng, nprocs * int_size, "cl1" );
   m = gcnvert * 3 / nprocs + 1;
   for ( i = 0; i < nprocs; i++ ) cv_cur_leng[i] = m;
   ML_memory_alloc((void**) &cv_list, nprocs * sizeof(int*), "cv5" );
   ML_memory_alloc((void**) &cv_cnt, nprocs * sizeof(int*), "cc5" );
   for ( i = 0; i < nprocs; i++ )
   {
      m = cv_cur_leng[i];
      ML_memory_alloc((void**) &(cv_list[i]), m * int_size, "cl6" );
      ML_memory_alloc((void**) &(cv_cnt[i]), m * int_size, "cc6" );
      for ( j = 0; j < m; j++ ) cv_cnt[i][j] = 0;
   }
   ML_memory_alloc((void**) &cv_ia, (nprocs + 1) * int_size, "cia" );

   /* ------------------------------------------------------------*/
   /* process the local coefficient array first                   */
   /* ------------------------------------------------------------*/

   for ( i = 0; i < ext_cnt; i++ )
   {
      cvnum = ext_ia[i];
      pnum  = node_proc_map[cvnum];
      if ( cv_leng[pnum] >= cv_cur_leng[pnum] )
      {
         tcv_list = cv_list[pnum];
         cv_cur_leng[pnum] += 20;
         m = cv_cur_leng[pnum];
         ML_memory_alloc((void**) &(cv_list[pnum]), m * int_size, "cl7" );
         for ( j = 0; j < cv_leng[pnum]; j++ )
            cv_list[pnum][j] = tcv_list[j];
         ML_memory_free( (void **) &tcv_list );
         tcv_list = cv_cnt[pnum];
         ML_memory_alloc((void**) &(cv_cnt[pnum]), m * int_size, "cc4" );
         for ( j = 0; j < cv_leng[pnum]; j++ )
            cv_cnt[pnum][j] = tcv_list[j];
         ML_memory_free( (void **) &tcv_list );
      }
      ML_search_insert_sort( cvnum, cv_list[pnum], &(cv_leng[pnum]),
                            cv_cnt[pnum] );
   }

   /* ------------------------------------------------------------*/
   /* now search the incoming coefficient array and its pointer   */
   /* array (from queries) to be added to cv_leng and cv_list.    */
   /* ------------------------------------------------------------*/

   index = proc_cnt = colbase = colend = 0;

   com_cnt = combuf->send_cnt;
   for ( i = 0; i < com_cnt; i++)
   {

      /* the following integer variables are used to set up the */
      /* offset to the coefficient array, due to the fact that  */
      /* each source processor delivers its subarray according  */
      /* to its own offset - so that first thing is to fetch    */
      /* this offset from its first location and compute its    */
      /* true offset                                            */

      colbase = colbase + colend;
      colcnt  = colbase - ext2_ptr[index] ;
      colend  = 0;

      /* now for each node sent for query to processor */
      /* com.send_proc[i]                              */

      for ( j = combuf->send_ia[i]; j < combuf->send_ia[i+1]; j++ )
      {
         fvnum = combuf->send_list[j];

         /* if it is an element number, ignore; */
         /* otherwise it is a node number       */

         if ( fvnum >= 0 )
         {

            /* fetch the beginning and end of */
            /* the coefficient pointers       */

            begin = ext2_ptr[index] + colcnt;
            end   = ext2_ptr[index+1] + colcnt;

            /* if the coefficient feedback has more than one */
            /* element, it is possible that one or more of   */
            /* the coarse grid vertices for the inquired     */
            /* element returns a nonzero, so register it     */

            if ((end - begin) > 1)
            {
               for ( k = begin; k < end; k++ )
               {
                  if ( ext2_a[k] > 1.0E-12 )
                  {
                     cvnum = ext2_ind[k];  /* global coarse node number */
                     pnum  = node_proc_map[cvnum];
                     if ( cv_leng[pnum] >= cv_cur_leng[pnum] )
                     {
                        tcv_list = cv_list[pnum];
                        cv_cur_leng[pnum] += 20;
                        m = cv_cur_leng[pnum];
                        ML_memory_alloc((void**) &(cv_list[pnum]),
                                        m*int_size, "cl7");
                        for ( j = 0; j < cv_leng[pnum]; j++ )
                           cv_list[pnum][j] = tcv_list[j];
                        ML_memory_free( (void**) &tcv_list );
                        tcv_list = cv_cnt[pnum];
                        ML_memory_alloc((void**) &(cv_cnt[pnum]),
                                        m*int_size, "cc8");
                        for ( j = 0; j < cv_leng[pnum]; j++ )
                           cv_cnt[pnum][j] = tcv_list[j];
                        ML_memory_free( (void**) &tcv_list );
                     }
                     ML_search_insert_sort(cvnum, cv_list[pnum],
                                       &(cv_leng[pnum]), cv_cnt[pnum]);
                  }
               }
            }
            colend = colend + end - begin;
         }

         /* increment to the next set of coefficients */

         index++;
      }

      /* increment to the next element */

      index++;
   }

   /* ------------------------------------------------------------*/
   /* see if there is anything to be sent at all, and also        */
   /* compose the cv_ia array                                     */
   /* ------------------------------------------------------------*/

   proc_cnt = 0;
   cv_ia[0] = 0;
   for ( i = 0; i < nprocs; i++ )
   {
      if ( cv_leng[i] > 0 ) proc_cnt++;
      cv_ia[i+1] = cv_ia[i] + cv_leng[i];
   }

   /* ------------------------------------------------------------*/
   /* now let the remote processor of my intention to send during */
   /* restriction and to receive during prolongation.             */
   /* (i.e. to set up the recv_list in xsfer_op)                  */
   /* ------------------------------------------------------------*/

   /* find out how many processors to receive from and register*/

   ML_memory_alloc((void**) &proc_flag, nprocs * int_size, "pf3" );
   ML_memory_alloc((void**) &tmp_ia,    nprocs * int_size, "tia" );
   for ( i = 0; i < nprocs; i++ )
   {
      if ( cv_leng[i] > 0 ) proc_flag[i] = 1;
      else                  proc_flag[i] = 0;
   }
   /*ML_Comm_GsumVecInt( comm, proc_flag, tmp_ia, nprocs );*/
   ML_gsum_vec_int( &proc_flag, &tmp_ia, nprocs, comm );
   m = proc_flag[mypid];
   xsfer_op->com->recv_cnt = m;
   if ( m > 0 )
   {
      ML_memory_alloc((void**)&(xsfer_op->com->recv_proc),m*int_size,"XS8");
      ML_memory_alloc((void**)&(xsfer_op->com->recv_ia),(m+1)*int_size,"XS9");
   }
   else
   {
      xsfer_op->com->recv_proc = NULL;
      xsfer_op->com->recv_ia   = NULL;
      xsfer_op->com->recv_list = NULL;
   }

   /* get the length of data to be received from each remote processor */

   for ( i = 0; i < nprocs; i++ ) tmp_ia[i] = 0;
   if ( proc_flag[mypid] > 0 )
   {
      k = proc_flag[mypid] * sizeof(USR_REQ);
      ML_memory_alloc((void**)&Request, k, "XT9");
      k = proc_flag[mypid] * sizeof(int);
      ML_memory_alloc((void**)&intarray, k, "XU9");
   }
   else
   {
      Request  = NULL;
      intarray = NULL;
   }

   msgtype = 5346;
   for ( i = 0; i < proc_flag[mypid]; i++ )
   {
      fromproc = -1;
      comm->USR_irecvbytes((void*) &intarray[i], int_size, &fromproc,
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
   }
   for ( i = 0; i < nprocs; i++ )
   {
      if ( cv_leng[i] > 0 )
         comm->USR_sendbytes((void*) &cv_leng[i], int_size, i, msgtype,
                             comm->USR_comm );
   }
   for ( i = 0; i < proc_flag[mypid]; i++ )
   {
      fromproc = -1;
      comm->USR_cheapwaitbytes((void*) &intarray[i], int_size, &fromproc, 
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i]);
#else
                           &msgtype, comm->USR_comm, (void *)&Request[i]);
#endif
      tmp_ia[fromproc] = intarray[i];
   }
   if ( proc_flag[mypid] > 0 )
   {
      ML_memory_free( (void**) &Request );
      ML_memory_free( (void**) &intarray );
   }

   k = 0;
   for ( i = 0; i < nprocs; i++ )
   {
      if ( tmp_ia[i] > 0 ) {
         xsfer_op->com->recv_proc[k] = i;
         xsfer_op->com->recv_ia[k]   = tmp_ia[i];;
         k++;
      }
   }
   if ( xsfer_op->com->recv_cnt > 0 )
   {
      for ( i = xsfer_op->com->recv_cnt; i >= 1; i-- )
         xsfer_op->com->recv_ia[i] = xsfer_op->com->recv_ia[i-1];
      xsfer_op->com->recv_ia[0] = 0;
      for ( i = 1; i <= xsfer_op->com->recv_cnt; i++ )
         xsfer_op->com->recv_ia[i] += xsfer_op->com->recv_ia[i-1];
      m = xsfer_op->com->recv_ia[xsfer_op->com->recv_cnt];
      ML_memory_alloc((void**)&(xsfer_op->com->recv_list),m*sizeof(xsfer_op->com->recv_list[0]),"XSa");
      ML_memory_alloc((void**) &recv_list32, m*int_size, "tia" );
   }
   ML_memory_free( (void**) &tmp_ia );
   ML_memory_free( (void**) &proc_flag );

   /* send and receive node lists (coarse nodes in global ordering) */

   msgtype = 5347;
   k = xsfer_op->com->recv_cnt;
   if ( k > 0 )
   {
      ML_memory_alloc( (void**) &Request, k*sizeof(USR_REQ), "XW0" );
   }
   else
   {
      Request = NULL;
   }
   for ( i = 0; i < k; i++ )
   {
      fromproc = xsfer_op->com->recv_proc[i];
      fromsize = (xsfer_op->com->recv_ia[i+1] - xsfer_op->com->recv_ia[i]) *
                  int_size;
      tmp_ia = &(recv_list32[xsfer_op->com->recv_ia[i]]);
      comm->USR_irecvbytes((void*) tmp_ia, fromsize, &fromproc, &msgtype,
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }
   for ( i = 0; i < nprocs; i++ )
   {
      if ( cv_leng[i] > 0 )
      {
         m = cv_leng[i] * int_size;
         comm->USR_sendbytes((void*) cv_list[i], m, i, msgtype,
                             comm->USR_comm);
      }
   }
   for ( i = 0; i < xsfer_op->com->recv_cnt; i++ )
   {
      fromproc = xsfer_op->com->recv_proc[i];
      fromsize = (xsfer_op->com->recv_ia[i+1] - xsfer_op->com->recv_ia[i]) *
                  int_size;
      tmp_ia = &(recv_list32[xsfer_op->com->recv_ia[i]]);
      comm->USR_cheapwaitbytes((void*) tmp_ia, fromsize, &fromproc, &msgtype, 
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }
   k = xsfer_op->com->recv_cnt;
   if ( k > 0 )
   {
      ML_memory_free( (void**) &Request );
   }

   /* convert the incoming indices into local indices */

   if ( xsfer_op->com->recv_cnt > 0 )
   {
      m      = xsfer_op->com->recv_ia[xsfer_op->com->recv_cnt];
      tmp_ia = recv_list32;
      for ( i = 0; i < m; i++ )
      {
         index = tmp_ia[i];
         j = 0;
         k = cgrid_fcns->USR_grid_get_vertex_global_num( c_grid, j );
         while ( index != k && j < cnvert )
         {
            j++;
            k = cgrid_fcns->USR_grid_get_vertex_global_num( c_grid, j );
         }
         if ( j >= cnvert )
         {
            printf("RP1 %d : something wrong. \n", mypid );
            exit(0);
         }
         tmp_ia[i] = j;
      }

      for (i=0; i<m; ++i)
        xsfer_op->com->recv_list[i]=recv_list32[i];
   }


   if ( xsfer_op->com->recv_cnt > 0 )
     ML_memory_free( (void**) &recv_list32 );

   /* ------------------------------------------------------------*/
   /* use the cv_leng, cv_cnt, and cv_list arrays to compute the  */
   /* amount of memory needed to store the remote part of the     */
   /* grid transfer operator.                                     */
   /* ------------------------------------------------------------*/

   if ( proc_cnt > 0 )
   {
      xsfer_op->com->send_cnt  = proc_cnt;
      nbytes = proc_cnt * sizeof(int);
      ML_memory_alloc((void**)&(xsfer_op->com->send_proc), nbytes, "XSb");
      nbytes = (proc_cnt + 1) * sizeof(int);
      ML_memory_alloc((void**)&(xsfer_op->com->send_ia), nbytes, "XSc");
      xsfer_op->com->send_ia[0] = 0;
   }
   else
   {
      xsfer_op->com->send_cnt  = 0;
      xsfer_op->com->send_proc = NULL;
      xsfer_op->com->send_ia   = NULL;
   }

   /* set up the row index field */

   k = 0;
   for ( i = 0; i < nprocs; i++ )
   {
      if (cv_leng[i] > 0)
      {
         xsfer_op->com->send_proc[k++] = i;
         xsfer_op->com->send_ia[k] = xsfer_op->com->send_ia[k-1] + cv_leng[i];
      }
   }
   if ( proc_cnt > 0 ) lcnt = xsfer_op->com->send_ia[k];
   else                lcnt = 0;
   xsfer_op->Nremote_rows   = lcnt;
   xsfer_op->com->send_list  = NULL;

   /* set up the matrix operator for restriction */

   if ( lcnt > 0 )
   {
      nbytes = ( lcnt + 1 ) * sizeof(int);
      ML_memory_alloc((void**) &(xsfer_op->remote_ia), nbytes, "XSd");
      xsfer_op->remote_ia[0] = 0;
      k = 0;
      for ( i = 0; i < nprocs; i++ )
      {
         for ( j = 0; j < cv_leng[i]; j++ )
         {
            xsfer_op->remote_ia[k+1] = xsfer_op->remote_ia[k] + cv_cnt[i][j];
            k++;
         }
      }
      m = xsfer_op->remote_ia[k];
      nbytes = m * sizeof(int);
      ML_memory_alloc((void**) &(xsfer_op->remote_ja), nbytes, "XSe");
      nbytes = m * sizeof(double);
      ML_memory_alloc((void**) &(xsfer_op->remote_a), nbytes, "XSf");
   }
   else
   {
      xsfer_op->remote_ia = NULL;
      xsfer_op->remote_ja = NULL;
      xsfer_op->remote_a  = NULL;
   }

   /* ------------------------------------------------------------*/
   /* insert the operator coefficients into the operator          */
   /* ------------------------------------------------------------*/

   tmp_ia = NULL;
   if ( lcnt > 0 )
   {

      ML_memory_alloc((void**) &tmp_ia,  lcnt * sizeof(int), "tia" );
      for ( i = 0; i < lcnt; i++ ) tmp_ia[i] = xsfer_op->remote_ia[i];

      /* insert the coefficients coming from local queries */

      for ( i = 0; i < ext_cnt; i++ )
      {
         cvnum  = ext_ia[i];
         pnum   = node_proc_map[cvnum];
         index  = ML_sorted_search( cvnum, cv_leng[pnum], cv_list[pnum]);
         rownum = cv_ia[pnum] + index;
         colind = tmp_ia[rownum]++;
         xsfer_op->remote_ja[colind] = ext_ja[i];
         xsfer_op->remote_a[colind]  = ext_a[i];
      }

      /* insert the coefficients coming back from queries */

      index = colbase = colend = 0;
      for ( i = 0; i < com_cnt; i++ )
      {

         colbase = colbase + colend;
         colcnt  = colbase - ext2_ptr[index] ;
         colend  = 0;

         for ( j = combuf->send_ia[i]; j < combuf->send_ia[i+1]; j++ )
         {
            fvnum = combuf->send_list[j];

            if (fvnum >= 0)
            {
               begin = ext2_ptr[index] + colcnt;
               end   = ext2_ptr[index+1] + colcnt;

               if ((end - begin) > 1 && (fnode_flag[fvnum] == -1 ) )
               {
                  fnode_flag[fvnum] = 1;
                  for ( k = begin; k < end; k++ )
                  {
                     if ( ext2_a[k] > 1.0E-12 )
                     {
                        cvnum = ext2_ind[k];
                        pnum  = node_proc_map[cvnum];
                        m = ML_sorted_search(cvnum, cv_leng[pnum],
                                             cv_list[pnum]);
                        if ( m < 0 )
                        {
                           printf("%d : impossible : searching for %d \n",
                                  mypid,  cvnum );
                           for ( j = 0; j < cv_leng[pnum]; j++ )
                              printf("    %d : available = %d \n", mypid,
                                     cv_list[pnum][j]);
                        }
                        rownum = cv_ia[pnum] + m;
                        colind = tmp_ia[rownum]++;
                        xsfer_op->remote_ja[colind] = fvnum;
                        xsfer_op->remote_a[colind]  = ext2_a[k];
                     }
                  }
               }
               colend = colend + end - begin;
            }
            index++;
         }
         index++;
      }
   }

   /* ------------------------------------------------------------*/
   /* compress the remote part of the operator, if needed.        */
   /* ------------------------------------------------------------*/

   k = xsfer_op->Nremote_rows;
   m = 0;
   for ( i = 0; i < k; i++ )
      m += (tmp_ia[i] - xsfer_op->remote_ia[i]);
   if ( m > 0 )
   {
      ML_memory_alloc((void**) &ext_ja,  m * sizeof(int), "XSh" );
      ML_memory_alloc((void**) &ext_a,  m * sizeof(double), "XSi" );
      m = 0;
      for ( i = 0; i < k; i++ )
      {
         for ( j = xsfer_op->remote_ia[i]; j < tmp_ia[i]; j++)
         {
            ext_ja[m]  = xsfer_op->remote_ja[j];
            ext_a[m++] = xsfer_op->remote_a[j];
         }
         tmp_ia[i] = m;
      }
      xsfer_op->remote_ia[0] = 0;
      for ( i = 1; i <= k; i++ )
         xsfer_op->remote_ia[i] = tmp_ia[i-1];
      ML_memory_free( (void**) &(xsfer_op->remote_ja) );
      ML_memory_free( (void**) &(xsfer_op->remote_a) );
      xsfer_op->remote_ja = ext_ja;
      xsfer_op->remote_a  = ext_a;
   }
   if ( tmp_ia != NULL ) ML_memory_free((void **) &tmp_ia);

   /* ------------------------------------------------------------*/
   /* clean up before moving on                                   */
   /* ------------------------------------------------------------*/

   ML_memory_free((void **) &cv_leng);
   ML_memory_free((void **) &cv_cur_leng);
   for ( i = 0; i < nprocs; i++ ) ML_memory_free((void **) &(cv_list[i]));
   ML_memory_free( (void **) &cv_list );
   for ( i = 0; i < nprocs; i++ ) ML_memory_free((void **) &(cv_cnt[i]));
   ML_memory_free( (void **) &cv_cnt );
   ML_memory_free( (void **) &cv_ia );

   /* ------------------------------------------------------------*/
   /* find normalization factor for each row for restriction      */
   /* ------------------------------------------------------------*/

   ML_memory_alloc((void**) &ext_a, fnvert * sizeof(double), "eta" );
   ML_memory_alloc((void**) &ext2_a, cnvert * sizeof(double), "e2a" );
   for ( i = 0; i < fnvert; i++ ) ext_a[i] = 1.0;
   nbytes = cnvert * sizeof(double);
   ML_memory_alloc((void**) &(xsfer_op->restrict_wgts), nbytes, "XSg");
   for ( i = 0; i < cnvert; i++ ) xsfer_op->restrict_wgts[i] = 1.0;
   oper = (ML_Operator *) ML_allocate( sizeof(ML_Operator) );
   oper->data = (void *) xsfer_op;
   oper->invec_leng = fnvert;
   oper->outvec_leng = cnvert;
   xsfer_op->AGX_stride = 1;
   ML_OperatorAGX_Restrict( oper, fnvert, ext_a, cnvert, ext2_a);
   ML_free( oper );

   /* ------------------------------------------------------------*/
   /* check normalization factors                                 */
   /* ------------------------------------------------------------*/

   err_flag = 0;
   for ( i = 0; i < cnvert; i++ )
   {
      if ( ext2_a[i] == 0.0 ) err_flag++;
   }
   err_flag = ML_Comm_GmaxInt( comm, err_flag );
   if ( err_flag > 0 )
   {
      ndim = cgrid_fcns->USR_grid_get_dimension( c_grid );
      for ( i = 0; i < cnvert; i++ )
      {
         cgrid_fcns->USR_grid_get_vertex_coordinate( c_grid, i, coord);
         j = cgrid_fcns->USR_grid_get_vertex_global_num( c_grid, i );
         if ( ndim == 2 )
            printf("%3d : coord, wgt %3d = (%12.6e %12.6e) %12.6e \n",
                   mypid, j, coord[0], coord[1], ext2_a[i] );
         else
            printf("%3d : coord, wgt %3d = (%12.6e %12.6e %12.6e) %12.6e \n",
                   mypid, j, coord[0], coord[1], coord[2], ext2_a[i] );
      }
      exit(1);
   }

   /* ------------------------------------------------------------*/
   /* if normalization factors are nonzero, get their inverses    */
   /* ------------------------------------------------------------*/

   for ( i = 0; i < cnvert; i++ )
      xsfer_op->restrict_wgts[i] = 1.0 / ext2_a[i];

   /* ------------------------------------------------------------*/
   /* final clean up                                              */
   /* ------------------------------------------------------------*/

   ML_memory_free( (void**) &ext_a );
   ML_memory_free( (void**) &ext2_a );
   xsfer_op->ext_cnt = 0;
   if ( xsfer_op->ext_ia != 0 )
   {
      ML_memory_free( (void**) &(xsfer_op->ext_ia) );
   }
   if ( xsfer_op->ext_ja != 0 )
   {
      ML_memory_free( (void**) &(xsfer_op->ext_ja) );
   }
   if ( xsfer_op->ext_a != 0 )
   {
      ML_memory_free( (void**) &(xsfer_op->ext_a) );
   }
   xsfer_op->ext2_cnt = 0;
   if ( xsfer_op->ext2_a != 0 )
   {
      ML_memory_free( (void**) &(xsfer_op->ext2_a) );
   }
   if ( xsfer_op->ext2_index != 0 )
   {
      ML_memory_free( (void**) &(xsfer_op->ext2_index) );
   }
   if ( xsfer_op->ext2_ptr != 0 )
   {
      ML_memory_free( (void**) &(xsfer_op->ext2_ptr) );
   }
   xsfer_op->fnode_leng = 0;
   if ( xsfer_op->fnode_flag != 0 )
   {
      ML_memory_free( (void**) &(xsfer_op->fnode_flag) );
   }
}

