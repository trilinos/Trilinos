/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for manipulating ML_OperatorAGX objects                    */
/* ******************************************************************** */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : August, 1998                                         */
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include "ml_operatoragx.h"

/* ******************************************************************** */
/* initialization function                                              */
/* -------------------------------------------------------------------- */

int ML_OperatorAGX_Create(ML_OperatorAGX **opp)
{
   ML_OperatorAGX *op;

   ML_memory_alloc( (void **) opp, sizeof(ML_OperatorAGX), "OP1");
   op =              (*opp);
   op->ML_id         = ML_ID_OPAGX;
   op->proc_id       = -1;
   op->step          = 0;
   op->num_procs     = 0;
   op->Nlocal_rows   = 0;
   op->local_ia      = 0;
   op->local_ja      = 0;
   op->local_a       = 0;
   op->Nremote_rows  = 0;
   op->remote_ia     = 0;
   op->remote_ja     = 0; 
   op->remote_a      = 0; 
   op->restrict_wgts = 0;
   op->ext_cnt       = 0;
   op->ext_ia        = 0;
   op->ext_ja        = 0;
   op->ext_a         = 0;
   op->ext2_cnt      = 0;
   op->ext2_ptr      = 0; 
   op->ext2_index    = 0;
   op->ext2_a        = 0;
   op->fnode_flag    = 0;
   op->fnode_leng    = 0;
   op->remote_restrict_wgts = 0;
   op->coarse_bdry_list = 0;
   op->coarse_bdry_leng = -1;
   op->AGX_comm      = 0;
   op->AGX_stride    = 0;

   ML_CommInfoAGX_Create(&(op->com));
   return 0;
}

/* ******************************************************************** */
/* destructor function                                                  */
/* -------------------------------------------------------------------- */
void ML_Operator2AGX_Destroy(void *data)
{
   ML_OperatorAGX *opagx;
   opagx = (ML_OperatorAGX *) data;
   ML_OperatorAGX_Destroy( &opagx);
}
int ML_OperatorAGX_Destroy(ML_OperatorAGX **opp)
{
   ML_OperatorAGX *op;

   op = (*opp);
   if ( op->ML_id != ML_ID_OPAGX )
   {
      printf("ML_OperatorAGX_Destroy : Wrong object. \n");
      exit(1);
   }
   op->Nlocal_rows  = 0;
   op->Nremote_rows = 0;
   op->ext_cnt      = 0;
   op->ext2_cnt     = 0;
   op->fnode_leng   = 0;
   if (op->local_ia   != 0) ML_memory_free( (void **) &(op->local_ia) ); 
   if (op->local_ja   != 0) ML_memory_free( (void **) &(op->local_ja) );
   if (op->local_a    != 0) ML_memory_free( (void **) &(op->local_a) );
   if (op->remote_ia  != 0) ML_memory_free( (void **) &(op->remote_ia) );
   if (op->remote_ja  != 0) ML_memory_free( (void **) &(op->remote_ja) );
   if (op->remote_a   != 0) ML_memory_free( (void **) &(op->remote_a) ); 
   if (op->ext_ia     != 0) ML_memory_free( (void **) &(op->ext_ia) );
   if (op->ext_ja     != 0) ML_memory_free( (void **) &(op->ext_ja) ); 
   if (op->ext_a      != 0) ML_memory_free( (void **) &(op->ext_a) ); 
   if (op->ext2_ptr   != 0) ML_memory_free( (void **) &(op->ext2_ptr) ); 
   if (op->ext2_index != 0) ML_memory_free( (void **) &(op->ext2_index) ); 
   if (op->ext2_a     != 0) ML_memory_free( (void **) &(op->ext2_a) ); 
   if (op->fnode_flag != 0) ML_memory_free( (void **) &(op->fnode_flag) ); 
   if (op->restrict_wgts != 0) ML_memory_free((void **) &(op->restrict_wgts)); 
   if (op->remote_restrict_wgts != 0) ML_free(op->remote_restrict_wgts);
   if (op->coarse_bdry_list!= 0) ML_free(op->coarse_bdry_list);
   op->remote_restrict_wgts = 0;
   op->coarse_bdry_list     = 0;
   op->coarse_bdry_leng     = -1;
   ML_CommInfoAGX_Destroy(&(op->com));
   op->ML_id = -1;
   ML_memory_free( (void **) opp);
   return 0;
}

/* ******************************************************************** */
/* print operator data function                                         */
/* -------------------------------------------------------------------- */

int ML_OperatorAGX_Print(ML_OperatorAGX *op)
{
   int i, j;

   if ( op->ML_id != ML_ID_OPAGX )
   {
      printf("ML_OperatorAGX_Print : Wrong object. \n");
      exit(1);
   }
   printf("**** Local operator : numrows = %d \n",op->Nlocal_rows);
   for ( i = 0; i < op->Nlocal_rows; i++ ) 
   {
      for (j=op->local_ia[i]; j<op->local_ia[i+1]; j++) 
         printf("   (%4d,%4d) = %e\n",i,op->local_ja[j],op->local_a[j]);
   }
   printf("**** Remote operator : numrows = %d \n",op->Nremote_rows);
   for ( i = 0; i < op->Nremote_rows; i++ ) 
   {
      for (j=op->remote_ia[i]; j<op->remote_ia[i+1]; j++) 
         printf("   (%4d,%4d) = %e\n",i,op->remote_ja[j],op->remote_a[j]);
   }
   printf("**** Remote data : no. entries = %d \n",op->ext_cnt);
   for ( i = 0; i < op->ext_cnt; i++ ) 
   {
      printf("   (%4d,%4d) = %e\n", op->ext_ia[i], op->ext_ja[i], 
                                       op->ext_a[i]);
   }
   ML_CommInfoAGX_Print( op->com );
   return 0;
}

/* **********************************************************************/
/* ML_OperatorAGX_Restrict : Given an operator and a data array, this   */
/* subroutine uses the given operator to transform it to another data   */
/* array                                                                */
/* -------------------------------------------------------------------- */

int ML_OperatorAGX_Restrict(ML_Operator *vop, int inlen, double *din, int outlen,
                            double *dout)
{
   int     i, j, k, index, icnt, icnt2, cnum, leng; 
   int     *send_proc = NULL, *recv_proc = NULL, *send_leng, *recv_leng = NULL;
   int ibegin, iend;
   int     pid, send_cnt, recv_cnt, mtype, step, istep, jstep;
   double  *send_buf, *recv_buf, mult, dtmp, *darray;
   ML_Operator    *op;
   ML_OperatorAGX *oplocal;
   USR_REQ    *Request;
   ML_Comm *comm;

   /* ----------------------------------------------------------------- */
   /* extract the operator data object                                  */
   /* ----------------------------------------------------------------- */
   op = (ML_Operator *) vop;
   oplocal = (ML_OperatorAGX *) op->data;
   if ( oplocal->ML_id != ML_ID_OPAGX )
   {
      printf("ML_OperatorAGX_Restrict : Wrong object. \n");
      ML_avoid_unused_param((void *) &inlen);
      ML_avoid_unused_param((void *) &outlen);
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* other checks                                                      */
   /* ----------------------------------------------------------------- */
   comm = oplocal->AGX_comm;
   icnt = oplocal->Nlocal_rows;
   step = oplocal->AGX_stride;

   /* ----------------------------------------------------------------- */
   /* allocate receive buffer                                           */
   /* ----------------------------------------------------------------- */

   recv_cnt  = oplocal->com->recv_cnt;
   if (recv_cnt > 0)
   {
      recv_proc = oplocal->com->recv_proc;
      recv_leng = oplocal->com->recv_ia;
      icnt = oplocal->com->recv_ia[recv_cnt];
      icnt2 = icnt * step;
      ML_memory_alloc((void **) &recv_buf, icnt2 * sizeof(double), "OP2");
   }

   /* ----------------------------------------------------------------- */
   /* perform matrix vector product to be sent to remote processors     */
   /* ----------------------------------------------------------------- */

   icnt = oplocal->Nremote_rows;
   icnt2 = icnt * step;
   if (icnt > 0)
   {
      ML_memory_alloc((void **) &send_buf, icnt2 * sizeof(double), "OP3");
      for ( i = 0; i < icnt; i++ )
      {
         istep  = i * step;
         ibegin = oplocal->remote_ia[i];
         iend   = oplocal->remote_ia[i+1];
         for ( k = 0; k < step; k++ ) send_buf[istep+k] = 0.0;
         for ( j = ibegin; j < iend; j++ )
         { 
            cnum  = oplocal->remote_ja[j];
            mult  = oplocal->remote_a[j];
            jstep = cnum * step;
            for ( k = 0; k < step; k++ ) 
              send_buf[istep+k] += ( mult * din[jstep+k] );
         }
      }
   }

   /* ----------------------------------------------------------------- */
   /* set up send information                                           */
   /* ----------------------------------------------------------------- */

   send_cnt = oplocal->com->send_cnt;
   if (send_cnt > 0)
   {
      send_proc = oplocal->com->send_proc;
      ML_memory_alloc((void **) &send_leng,send_cnt*sizeof(int),"OP4");
      for ( i = 0; i < send_cnt; i++ )
         send_leng[i] = oplocal->com->send_ia[i+1] - oplocal->com->send_ia[i];
   }   

   /* ----------------------------------------------------------------- */
   /* post a receive for incoming data                                  */
   /* ----------------------------------------------------------------- */

   if ( recv_cnt > 0 )
   {
      ML_memory_alloc((void **) &Request, recv_cnt*sizeof(USR_REQ),"OP9");
   }
   for ( i = 0; i < recv_cnt; i++ )
   {
      mtype = 767;
      leng  = (recv_leng[i+1] - recv_leng[i]) * sizeof(double) * step;
      pid   = recv_proc[i];
      index = recv_leng[i] * step;
      comm->USR_irecvbytes((void*) &(recv_buf[index]), leng, &pid, &mtype, 
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }

   /* ----------------------------------------------------------------- */
   /* send border data                                                  */
   /* ----------------------------------------------------------------- */

   mtype = 767;
   icnt  = 0;
   for ( i = 0; i < send_cnt; i++ )
   {
      leng = send_leng[i] * sizeof(double) * step;
      pid  = send_proc[i];
      comm->USR_sendbytes((void*)&(send_buf[icnt]), leng, pid, mtype, 
                          comm->USR_comm );
      icnt += (send_leng[i] * step);
   }

   /* ----------------------------------------------------------------- */
   /* local multiplication                                              */
   /* ----------------------------------------------------------------- */

   icnt   = oplocal->Nlocal_rows;
   icnt2  = oplocal->Nlocal_rows * step;
   ML_memory_alloc((void **) &darray, step * sizeof(double),"OP5");

   for ( i = 0; i < icnt; i++ )
   {
      for ( k = 0; k < step; k++ ) darray[k] = 0.0;
      ibegin = oplocal->local_ia[i];
      iend   = oplocal->local_ia[i+1];
      for ( j = ibegin; j < iend; j++ )
      {
         cnum  = oplocal->local_ja[j];
         mult  = oplocal->local_a[j];
         jstep = cnum * step;
         for ( k = 0; k < step; k++ ) 
            darray[k] += ( mult * din[jstep+k] );
      }
      istep = i * step;
      for ( k = 0; k < step; k++ ) dout[istep+k] = darray[k];
   }
   ML_memory_free( (void **) &darray );

   /* ----------------------------------------------------------------- */
   /* wait for incoming data                                            */
   /* ----------------------------------------------------------------- */

   for ( i = 0; i < recv_cnt; i++ )
   {
      mtype = 767;
      leng  = (recv_leng[i+1] - recv_leng[i]) * sizeof(double) * step;
      pid   = recv_proc[i];
      index = recv_leng[i] * step;
      comm->USR_cheapwaitbytes((void*) &(recv_buf[index]), leng, &pid, &mtype, 
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }

   /* ----------------------------------------------------------------- */
   /* multiplication for the external data                              */
   /* ----------------------------------------------------------------- */

   if (recv_cnt > 0)
   {
      icnt = recv_leng[recv_cnt];
      for ( i = 0; i < icnt; i++ )
      {
         index = oplocal->com->recv_list[i];
         istep = i * step;
         jstep = index * step;
         for ( k = 0; k < step; k++ ) 
            dout[jstep+k] += recv_buf[istep+k];
      }
   }

   /* ----------------------------------------------------------------- */
   /* scaling                                                           */
   /* ----------------------------------------------------------------- */

   icnt = oplocal->Nlocal_rows;
   for ( i = 0; i < icnt; i++ )
   { 
      dtmp  = oplocal->restrict_wgts[i];
      istep = i * step;
      for ( k = 0; k < step; k++ ) dout[istep+k] *= dtmp;
   }

   /* ----------------------------------------------------------------- */
   /* clean up                                                          */
   /* ----------------------------------------------------------------- */

   icnt = oplocal->Nremote_rows;
   if (icnt > 0)       ML_memory_free( (void **) &send_buf);
   if ( send_cnt > 0 ) ML_memory_free( (void **) &send_leng);
   if (recv_cnt > 0)   ML_memory_free( (void **) &recv_buf);
   if (recv_cnt > 0)   ML_memory_free( (void **) &Request);
   return 0;
}

/* **********************************************************************/
/* ML_OperatorAGX_Prolongate : Given an operator and a data array, this */
/*  subroutine uses the given operator to transposely transform it to   */
/*  another data array                                                  */
/* -------------------------------------------------------------------- */

int ML_OperatorAGX_Prolongate(ML_Operator *vop, int inlen, double *din, 
                              int outlen, double *dout)
{
   int    i, j, k, index, icnt, icnt2, rnum, leng; 
   int    *send_proc = NULL, *recv_proc = NULL, *send_leng, *recv_leng = NULL;
   int    ibegin, iend;
   int    pid, send_cnt, recv_cnt, mtype, step, istep, jstep;
   double *send_buf, *recv_buf, mult;
   ML_Operator    *op;
   ML_OperatorAGX *oplocal;
   USR_REQ    *Request;
   ML_Comm    *comm;

   /* ----------------------------------------------------------------- */
   /* extract the operator data object                                  */
   /* ----------------------------------------------------------------- */

   op = (ML_Operator *) vop;
   oplocal = (ML_OperatorAGX *) op->data;
   if ( oplocal->ML_id != ML_ID_OPAGX )
   {
      printf("ML_OperatorAGX_Prolongate : Wrong object. \n");
      ML_avoid_unused_param((void *) &inlen);
      exit(1);
   }
   comm = oplocal->AGX_comm;
   icnt = oplocal->Nlocal_rows;
   step = oplocal->AGX_stride;

   /* ----------------------------------------------------------------- */
   /* allocate receive buffer                                           */
   /* ----------------------------------------------------------------- */

   recv_cnt  = oplocal->com->send_cnt;
   if (recv_cnt > 0)
   {
      recv_proc = oplocal->com->send_proc;
      recv_leng = oplocal->com->send_ia;
      icnt = oplocal->com->send_ia[recv_cnt];
      icnt2 = icnt * step;
      ML_memory_alloc( (void **) &recv_buf, icnt2*sizeof(double), "OP6");
   }

   /* ----------------------------------------------------------------- */
   /* set up send information                                           */
   /* ----------------------------------------------------------------- */

   send_cnt = oplocal->com->recv_cnt;
   if (send_cnt > 0)
   {
      send_proc = oplocal->com->recv_proc;
      ML_memory_alloc( (void **) &send_leng,send_cnt*sizeof(int), "OP7");
      for ( i = 0; i < send_cnt; i++ )
         send_leng[i] = oplocal->com->recv_ia[i+1] - oplocal->com->recv_ia[i];
   }   

   /* ----------------------------------------------------------------- */
   /* collect send data                                                 */
   /* ----------------------------------------------------------------- */

   if (send_cnt > 0)
   {
      icnt  = oplocal->com->recv_ia[send_cnt];
      icnt2 = icnt * step;
      ML_memory_alloc( (void **) &send_buf,icnt2*sizeof(double), "OP8");
      for ( i = 0; i < icnt; i++ )
      {
         istep = i * step;
         index = oplocal->com->recv_list[i];
         jstep = index * step;
         for ( k = 0; k < step; k++ ) send_buf[istep+k] = din[jstep+k];
      }
   }

   /* ----------------------------------------------------------------- */
   /* post a receive for incoming data                                  */
   /* ----------------------------------------------------------------- */

   if ( recv_cnt > 0 )
   {
      ML_memory_alloc((void **) &Request, recv_cnt*sizeof(USR_REQ),"OQ9");
   }
   for ( i = 0; i < recv_cnt; i++ )
   {
      mtype = 769;
      leng  = (recv_leng[i+1] - recv_leng[i]) * sizeof(double) * step;
      pid   = recv_proc[i];
      index = recv_leng[i] * step;
      comm->USR_irecvbytes((void*) &(recv_buf[index]), leng, &pid, &mtype, 
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }

   /* ----------------------------------------------------------------- */
   /* send border data                                                  */
   /* ----------------------------------------------------------------- */

   mtype = 769;
   icnt  = 0;
   for ( i = 0; i < send_cnt; i++ )
   {
      leng = send_leng[i] * sizeof(double) * step;
      pid  = send_proc[i];
      comm->USR_sendbytes((void*)&(send_buf[icnt]), leng, pid, mtype, 
                          comm->USR_comm );
      icnt += (send_leng[i] * step);
   }

   /* ----------------------------------------------------------------- */
   /* local multiplication                                              */
   /* ----------------------------------------------------------------- */

   icnt  = oplocal->Nlocal_rows;
   icnt2 = icnt * step;

   for ( i = 0; i < outlen; i++ ) dout[i] = 0.0;
   for ( i = 0; i < icnt; i++ )
   {
      ibegin = oplocal->local_ia[i];
      iend   = oplocal->local_ia[i+1];
      istep  = i * step;
      for ( j = ibegin; j < iend; j++ )
      {
         rnum  = oplocal->local_ja[j];
         mult  = oplocal->local_a[j];
         jstep = rnum * step;
         for ( k = 0; k < step; k++ ) 
            dout[jstep+k] += ( mult * din[istep+k] );
      }
   }

   /* ----------------------------------------------------------------- */
   /* wait for incoming data                                            */
   /* ----------------------------------------------------------------- */

   for ( i = 0; i < recv_cnt; i++ )
   {
      mtype = 769;
      leng  = (recv_leng[i+1] - recv_leng[i]) * sizeof(double) * step;
      pid   = recv_proc[i];
      index = recv_leng[i] * step;
      comm->USR_cheapwaitbytes((void*) &(recv_buf[index]), leng, &pid, &mtype, 
#ifdef ML_CPP
                           comm->USR_comm, &Request[i]);
#else
                           comm->USR_comm, (void *)&Request[i]);
#endif
   }

   /* ----------------------------------------------------------------- */
   /* multiplication for external data                                  */
   /* ----------------------------------------------------------------- */

   if (recv_cnt > 0)
   {
      icnt = oplocal->Nremote_rows;
      for ( i = 0; i < icnt; i++ )
      {
         ibegin = oplocal->remote_ia[i];
         iend   = oplocal->remote_ia[i+1];
         istep  = i * step;
         for ( j = ibegin; j < iend; j++ )
         {
            rnum  = oplocal->remote_ja[j];
            mult  = oplocal->remote_a[j];
            jstep = rnum * step;
            for ( k = 0; k < step; k++ ) 
               dout[jstep+k] += ( mult * recv_buf[istep+k] );
         }
      }
   }

   /* ----------------------------------------------------------------- */
   /* clean up                                                          */
   /* ----------------------------------------------------------------- */

   if (send_cnt > 0)
   {
      ML_memory_free( (void **) &send_buf);
      ML_memory_free( (void **) &send_leng);
   }
   if (recv_cnt > 0)
   {
      ML_memory_free( (void **) &recv_buf);
      ML_memory_free( (void **) &Request);
   }
   return 0;
}

/* **********************************************************************/
/* Ray's version of getrows (This is the one to use for now.)           */
/* ---------------------------------------------------------------------*/

int ML_OperatorAGX_Getrows(ML_Operator *data, int N_requested_rows, 
                int requested_rows[], int allocated_space, int columns[], 
                double values[], int row_lengths[])
{
   int            i, j, *ia, *ja, ncount;
   int            blk_ind, step, offset;
   double         *a;
   ML_OperatorAGX *local_op;
   int            ii, dlength, *dlist;
   int            coarse_bdry_leng, *coarse_bdry_list;
   double         *remote_restrict_wgts;
   int            *temp_local_ia, *temp_local_ja, rowind;
   int            *temp_remote_ia, *temp_remote_ja;
   int            Nlocal, *new_local_ia = NULL, *new_local_ja;
   int            Nremote, *new_remote_ia = NULL, *new_remote_ja = NULL;
   int            ncnt, k, istep;
   double         *invec, *outvec, *new_remote_a = NULL, *temp_remote_a, dtmp;
   ML_1Level      *fine, *coarse;
   ML_Operator    *Rmat;

   /* --------------------------------------------- */
   /* generate coarse boundary list                 */
   /* --------------------------------------------- */

   Rmat     = (ML_Operator *) data;
   local_op = (ML_OperatorAGX *) Rmat->data;
   coarse   = Rmat->to;
   fine     = Rmat->from;
   step     = local_op->AGX_stride;

   coarse_bdry_leng = local_op->coarse_bdry_leng;
   coarse_bdry_list = local_op->coarse_bdry_list;
   remote_restrict_wgts = local_op->remote_restrict_wgts;

   if ( coarse_bdry_leng == -1 )
   {
      temp_local_ia = local_op->local_ia;
      temp_local_ja = local_op->local_ja;
      temp_remote_ia = local_op->remote_ia;
      temp_remote_ja = local_op->remote_ja;
      temp_remote_a = local_op->remote_a;
      Nlocal = local_op->Nlocal_rows;
      if (local_op->com->send_cnt == 0) Nremote = 0;
      else Nremote = local_op->com->send_ia[local_op->com->send_cnt];
      if ( Nlocal > 0 ) {
         new_local_ia = (int *) ML_allocate( (Nlocal + 1) * sizeof(int) );
         for ( i = 0; i <= Nlocal; i++ ) new_local_ia[i] = 0;
         new_local_ja = 0;
         local_op->local_ia = new_local_ia;
         local_op->local_ja = new_local_ja;
      }
      if ( Nremote > 0 ) {
         new_remote_ia = (int *) ML_allocate( (Nremote + 1) * sizeof(int) );
         for ( i = 0; i <= Nremote; i++ ) new_remote_ia[i] = i;
         new_remote_ja = (int *) ML_allocate( Nremote * sizeof(int) );
         for ( i = 0; i < Nremote; i++ ) new_remote_ja[i] = i;
         new_remote_a = (double *) ML_allocate( Nremote * sizeof(double) );
         for ( i = 0; i < Nremote; i++ ) new_remote_a[i] = 1.0;
         local_op->remote_ia = new_remote_ia;
         local_op->remote_ja = new_remote_ja;
         local_op->remote_a = new_remote_a;
      }
      if ( Nlocal > 0 ) {
         invec = (double *) ML_allocate( Nlocal * step * sizeof(double) );
         for ( i = 0; i < Nlocal; i++ ) {
            dtmp  = local_op->restrict_wgts[i];
            istep = i * step;
            for ( k = 0; k < step; k++ ) invec[istep+k] = dtmp;
         }
      } else invec = 0;
      if ( Nremote > 0 ) {
         outvec = (double *) ML_allocate( Nremote * step * sizeof(double) );
         for ( i = 0; i < Nremote*step; i++ ) outvec[i] = 0.0;
      } else outvec = 0;
      ncnt = 0;
      ML_BdryPts_Get_Dirichlet_Grid_Info(coarse->BCs, &dlength, &dlist);
      for ( i = 0; i < dlength; i++ ) {
         invec[dlist[i]] = -invec[dlist[i]];
         ncnt++;
      }

      ML_OperatorAGX_Prolongate( (void *) coarse->Pmat,
/* (void *) local_op, */ Nlocal*step, invec, 
                           Nremote*step, outvec);
      ML_BdryPts_Get_Dirichlet_Grid_Info(coarse->BCs, &j, &dlist);
      for ( i = 0; i < Nremote*step; i++ )
         if ( outvec[i] < 0.0 ) j++;
      if ( j > 0 ) coarse_bdry_list = (int *) ML_allocate( j * sizeof(int) );
      else         coarse_bdry_list = 0;
      remote_restrict_wgts = (double *) ML_allocate((1+Nremote)*sizeof(double));
      if (remote_restrict_wgts == NULL) { printf("Out of space\n"); exit(1); }
      
      coarse_bdry_leng = 0;
      ML_BdryPts_Get_Dirichlet_Grid_Info(coarse->BCs, &dlength, &dlist);
      if ( dlist != NULL) {
         for (j = 0; j < dlength; j++)
           coarse_bdry_list[coarse_bdry_leng++] = dlist[j];
      }
      for (i = 0; i < Nremote; i++) {
         remote_restrict_wgts[i] = outvec[i*step];
         if (remote_restrict_wgts[i] < 0.0) 
            remote_restrict_wgts[i] = -remote_restrict_wgts[i];
      }
      for ( i = 0; i < Nremote*step; i++ )
         if ( outvec[i] < 0.0 ) 
	    coarse_bdry_list[coarse_bdry_leng++]= i+Nlocal*step;
      if ( Nlocal > 0 ) {
         ML_free( new_local_ia );
         ML_free( invec );
      }
      if ( Nremote > 0 ) {
         ML_free( new_remote_ia );
         ML_free( new_remote_ja );
         ML_free( new_remote_a );
         ML_free( outvec );
      }
      local_op->local_ia = temp_local_ia;
      local_op->local_ja = temp_local_ja;
      local_op->remote_ia = temp_remote_ia;
      local_op->remote_ja = temp_remote_ja;
      local_op->remote_a = temp_remote_a;
      local_op->coarse_bdry_leng = coarse_bdry_leng;
      local_op->coarse_bdry_list = coarse_bdry_list;
      local_op->remote_restrict_wgts = remote_restrict_wgts;
   }

   /* --------------------------------------------- */
   /* the rest                                      */
   /* --------------------------------------------- */

   ncount = 0;
   for (k = 0; k < N_requested_rows; k++) {
      rowind = requested_rows[k];
      blk_ind = rowind/step;
      offset  = rowind - blk_ind*step;

      if ( blk_ind < local_op->Nlocal_rows) {
         ia = local_op->local_ia;
         ja = local_op->local_ja;
         a  = local_op->local_a;
         dtmp = local_op->restrict_wgts[blk_ind];
      }
      else {
         blk_ind -= local_op->Nlocal_rows;
         ia = local_op->remote_ia;
         ja = local_op->remote_ja;
         a  = local_op->remote_a;
         dtmp = remote_restrict_wgts[blk_ind];
      }

      row_lengths[k] = ia[blk_ind+1] - ia[blk_ind];
      if (ncount + row_lengths[k] > allocated_space) return(0);
      for ( j = ia[blk_ind]; j < ia[blk_ind+1]; j++ ) {
         columns[ncount]  = ja[j]*step + offset;
         values[ncount++] = a[j]*dtmp;
      }

      for ( ii = 0; ii < coarse_bdry_leng; ii++ )
         if ( rowind == coarse_bdry_list[ii] ) break;

      if (ii < coarse_bdry_leng) {
         ML_BdryPts_Get_Dirichlet_Grid_Info(fine->BCs, &dlength, &dlist);
         for (j = ncount - row_lengths[k]; j < ncount; j++) {
            if ( ML_find_index(columns[j], dlist, dlength) == -1) {
               values[j] = 0.0;
            }
         }
      }
   }
   return(1);
}

/* **********************************************************************/
/* Ray's version of getcols (This is the one to use for now.)           */
/* ---------------------------------------------------------------------*/

int ML_OperatorAGX_Getcols(ML_Operator *data, int N_requested_rows, 
                    int requested_rows[], int allocated_space, 
                    int columns[], double values[], int row_lengths[])

{
   int          i, j, ncount = 0, kk, colind;
   int          blk_ind, step, offset, fnodes, *patemp;
   int          Nlocal, *local_ia, *local_ja;
   int          Nremote, *remote_ia, *remote_ja;
   int           *ext_ia, *ext_ja, nbytes, col, row;
   double          *local_a, *remote_a, *ext_a;
   ML_OperatorAGX  *local_op;
   ML_Operator     *Pmat;

   /* --------------------------------------------------------- */
   /* extract the operator                                      */
   /* --------------------------------------------------------- */
   Pmat     = (ML_Operator *) data;
   local_op = (ML_OperatorAGX *) Pmat->data;
   step = local_op->AGX_stride;

   /* --------------------------------------------------------- */
   /* if column ordering has been done, return column           */
   /* --------------------------------------------------------- */
   if ( local_op->ext_cnt == 0   || local_op->ext_ia == NULL || 
        local_op->ext_ja == NULL || local_op->ext_a  == NULL )
   {
      Nlocal   = local_op->Nlocal_rows;
      local_ia = local_op->local_ia;
      local_ja = local_op->local_ja;
      local_a  = local_op->local_a;
      Nremote  = local_op->Nremote_rows;
      remote_ia = local_op->remote_ia;
      remote_ja = local_op->remote_ja;
      remote_a  = local_op->remote_a;
      fnodes = 0;
      for ( i = 0; i < local_ia[Nlocal]; i++ )
      {
         if ( local_ja[i] >= fnodes ) fnodes = local_ja[i];
      }
      if (Nremote > 0)
      {
         for ( i = 0; i < remote_ia[Nremote]; i++ )
         {
            if ( remote_ja[i] >= fnodes ) fnodes = remote_ja[i];
         }
      }
      fnodes++;
      ML_memory_alloc( (void**) &patemp, fnodes*sizeof(int), "ot1" );
      for ( i = 0; i < fnodes; i++) patemp[i] = 0;
      for ( i = 0; i < local_ia[Nlocal]; i++ ) patemp[local_ja[i]]++;
      if ( Nremote > 0 )
      {
         for (i = 0; i < remote_ia[Nremote]; i++) patemp[remote_ja[i]]++;
      }

      /* --------------------------------------------------------- */
      /* allocate memory for column matrix                         */
      /* --------------------------------------------------------- */
      nbytes = ( fnodes + 1 ) * sizeof(int);
      ML_memory_alloc( (void**) &(local_op->ext_ja), nbytes, "ot2" );
      if ( Nremote > 0 ) ncount = local_ia[Nlocal] + remote_ia[Nremote];
      else               ncount = local_ia[Nlocal];
      nbytes = ncount * sizeof(int);
      ML_memory_alloc( (void**) &(local_op->ext_ia), nbytes, "ot3" );
      nbytes = ncount * sizeof(double);
      ML_memory_alloc( (void**) &(local_op->ext_a), nbytes, "ot4" );

      /* --------------------------------------------------------- */
      /* fill in the matrix                                        */
      /* --------------------------------------------------------- */
      local_op->ext_cnt = fnodes;
      ext_ja      = local_op->ext_ja;
      ext_ia      = local_op->ext_ia;
      ext_a       = local_op->ext_a;
      ext_ja[0]   = 0; 
      for ( i = 1; i <= fnodes; i++ )
         ext_ja[i] = patemp[i-1] + ext_ja[i-1]; 
      for ( i = 0; i < Nlocal; i++ )
      {
         for ( j = local_ia[i]; j < local_ia[i+1]; j++ )
         {
            col = i;
            row = local_ja[j];
            ext_ia[ext_ja[row]]  = col;
            ext_a[ext_ja[row]++] = local_a[j];
         }
      }
      for ( i = 0; i < Nremote; i++ )
      {
         for ( j = remote_ia[i]; j < remote_ia[i+1]; j++ )
         {
            col = Nlocal + i;
            row = remote_ja[j];
            ext_ia[ext_ja[row]]  = col;
            ext_a[ext_ja[row]++] = remote_a[j];
         }
      }
      for ( i = fnodes; i > 0; i-- ) ext_ja[i] = ext_ja[i-1];
      ext_ja[0] = 0;
      ML_memory_free( (void**) &patemp );
   }

   /* --------------------------------------------------------- */
   /* if column ordering has been done, return column           */
   /* --------------------------------------------------------- */

   ext_ja = local_op->ext_ja;
   ext_ia = local_op->ext_ia;
   ext_a  = local_op->ext_a;
   ncount = 0;

   for (kk = 0; kk < N_requested_rows; kk++) {
      colind = requested_rows[kk];
      blk_ind = colind/step;
      offset  = colind - blk_ind*step;

      row_lengths[kk] = ext_ja[blk_ind+1] - ext_ja[blk_ind];
      if (ncount + row_lengths[kk] > allocated_space) return(0);

      for (i = ext_ja[blk_ind]; i < ext_ja[blk_ind+1]; i++) {
         columns[ncount] = ext_ia[i]*step + offset;
         values[ncount++] = ext_a[i];
      }
   }

   return(1);
}

/* **********************************************************************/
/* Ray's version (This is the one to use for now.)                      */
/* ---------------------------------------------------------------------*/

int ML_OperatorAGX_Clean_Getrows(ML_Operator **opp)
{
/* Clean up stuff that is only needed for RAP */

   ML_OperatorAGX *op;

   op = (ML_OperatorAGX *) (*opp)->data;
   if ( op->ML_id != ML_ID_OPAGX )
   {
      printf("ML_OperatorAGX_Clean_Getrows: Wrong object. \n");
      exit(1);
   }
   if (op->remote_restrict_wgts != 0) ML_free(op->remote_restrict_wgts);
   if (op->coarse_bdry_list!= 0) ML_free(op->coarse_bdry_list);
   op->remote_restrict_wgts = 0;
   op->coarse_bdry_list     = 0;
   op->coarse_bdry_leng     = -1;
   op->ext_cnt      = 0;
   if (op->ext_ja     != 0) ML_memory_free( (void **) &(op->ext_ja) ); 
   if (op->ext_ia     != 0) ML_memory_free( (void **) &(op->ext_ia) );
   if (op->ext_a      != 0) ML_memory_free( (void **) &(op->ext_a) ); 
   return 0;
}


int ML_OperatorAGX_Gen_ComminfoOp(ML_OperatorAGX *vop, ML_Operator *Rmat,
	ML_Operator *Pmat)
{
   ML_CommInfoAGX *com;
   int            i, j, k, kk, flag, step;
   int            send_count, recv_count, *send_data = NULL, *recv_data = NULL, 
                  *neighbors = NULL, Nneighbors, *remap;

   step       = vop->AGX_stride;
   com        = vop->com;
 
   /* convert between Charles communication information and RST's */

   send_count = 1;
   for (i = 0; i < com->send_cnt; i++) {
     kk = com->send_ia[i+1] - com->send_ia[i];
     if (kk > send_count) send_count = kk;
   }
   send_count *= step;

   recv_count = 1;
   for (i = 0; i < com->recv_cnt; i++) {
     kk = com->recv_ia[i+1] - com->recv_ia[i];
     if (kk > recv_count) recv_count = kk;
   }
   recv_count *= step;

   send_data = (int *) ML_allocate( send_count*sizeof(int));
   recv_data = (int *) ML_allocate( recv_count*sizeof(int));
   neighbors = (int *) ML_allocate((com->send_cnt+com->recv_cnt)*sizeof(int));

   /* copy the send neighbors into the neighbor array */

   Nneighbors = com->send_cnt;
   for (i = 0; i < Nneighbors; i++) neighbors[i] = com->send_proc[i];

   /* Add any receive neighbors who are not send neighbors */

   for (i = 0; i < com->recv_cnt; i++) {
      flag = 0;
      for (j = 0; j < com->send_cnt; j++)
         if (neighbors[j] == com->recv_proc[i]) {flag = 1; break;}
      if (flag == 0) neighbors[Nneighbors++] = com->recv_proc[i];
   }
   ML_CommInfoOP_Set_neighbors(&(Pmat->getrow->pre_comm),Nneighbors,
                               neighbors, ML_OVERWRITE,NULL,0);
   kk = (vop->Nlocal_rows + vop->Nremote_rows)*step;
   remap = (int *) ML_allocate(kk*sizeof(int));
   for (i = 0; i < Rmat->outvec_leng; i++) remap[i    ] =  i;
   for (i = 0; i < kk-Rmat->outvec_leng; i++) remap[i+Rmat->outvec_leng] = -1;
   ML_CommInfoOP_Set_neighbors(&(Rmat->getrow->post_comm), Nneighbors,
                               neighbors, ML_ADD, remap, kk);

   for (i = 0; i < Nneighbors; i++ ) {
      /* send list for ith neighbor */

      send_count = 0;
      if (i < com->send_cnt) {
         j = i;
         for (k = com->send_ia[j]; k < com->send_ia[j+1]; k++) {
            for (kk = 0; kk < step; kk++)
               send_data[send_count++]  = kk + (vop->Nlocal_rows +
                                                k) *step;
         }
      }

      /* recv list for ith neighbor */

      recv_count = 0;
      for (j = 0; j < com->recv_cnt; j++) {
         if ( com->recv_proc[j] == neighbors[i] ) {
            for (k = com->recv_ia[j]; k < com->recv_ia[j+1]; k++) {
               for (kk = 0; kk < step; kk++)
                  recv_data[recv_count++]  = kk + com->recv_list[k]*step;
            }
            break;
         }
      }

      ML_CommInfoOP_Set_exch_info(Pmat->getrow->pre_comm, neighbors[i],
                                  send_count, send_data, recv_count, recv_data);
      ML_CommInfoOP_Set_exch_info(Rmat->getrow->post_comm,neighbors[i],
                                  recv_count, recv_data,send_count, send_data);
   }
   if (neighbors != NULL) ML_free(neighbors);
   if (recv_data != NULL) ML_free(recv_data);
   if (send_data != NULL) ML_free(send_data);

   return 0;
}
