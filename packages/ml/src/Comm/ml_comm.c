/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions to communicate between processors                          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : Sept, 1998                                           */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_comm.h"
#include "ml_utils.h"

ML_Comm *global_comm = NULL; /* should not be used to avoid side effect */

/* ******************************************************************** */
/* Create a communicator for ML                                         */
/* -------------------------------------------------------------------- */

int ML_Comm_Create( ML_Comm ** com )
{
   ML_Comm *com_ptr;

   ML_memory_alloc( (void **) com, sizeof (ML_Comm), "CM1" );
   com_ptr              = (*com);
   com_ptr->ML_id       = ML_ID_COMM;
   com_ptr->ML_mypid    = 0;
   com_ptr->ML_nprocs   = 1;
   com_ptr->USR_comm    = 0;
   com_ptr->USR_sendbytes  = ML_Comm_Send;
   com_ptr->USR_irecvbytes = ML_Comm_Irecv;
   com_ptr->USR_waitbytes  = ML_Comm_Wait;
   com_ptr->USR_cheapwaitbytes  = ML_Comm_CheapWait;
   com_ptr->USR_errhandler = NULL;

#ifdef ML_MPI
   MPI_Comm_size(MPI_COMM_WORLD, &(com_ptr->ML_nprocs));
   MPI_Comm_rank(MPI_COMM_WORLD, &(com_ptr->ML_mypid));
   com_ptr->USR_sendbytes  = ML_Comm_Send;
   com_ptr->USR_irecvbytes = ML_Comm_Irecv;
   com_ptr->USR_waitbytes  = ML_Comm_Wait;
   com_ptr->USR_cheapwaitbytes  = ML_Comm_CheapWait;
   com_ptr->USR_comm       = MPI_COMM_WORLD;
#ifdef ML_CATCH_MPI_ERRORS_IN_DEBUGGER
   /* register the error handling function */
   ML_Comm_ErrorHandlerCreate((USR_ERRHANDLER_FUNCTION *) ML_Comm_ErrorHandler,
                              com_ptr->USR_errhandler);
   /* associate the error handling function with the communicator */
   ML_Comm_ErrorHandlerSet(com_ptr->USR_comm, com_ptr->USR_errhandler);
#endif
#endif /*ifdef ML_MPI*/

   return 0;
}

/* ******************************************************************** */
/* destroy a data structure for grid communication functions            */
/* -------------------------------------------------------------------- */

int ML_Comm_Destroy( ML_Comm ** com )
{
   if ( (*com) != NULL )
   {
      if ( (*com)->ML_id != ML_ID_COMM )
      {
         printf("ML_Comm_Destroy : Wrong Comm object to destroy. \n");
         return -1;
      }
      (*com)->ML_id = -1;
#ifdef ML_CATCH_MPI_ERRORS_IN_DEBUGGER
      ML_Comm_ErrorHandlerDestroy(&((*com)->USR_errhandler));
#endif
      ML_memory_free( (void **) com );
   }
   return 0;
}

/* ******************************************************************** */
/* Check that all communicator functions and variables have been set    */
/* (return -1 if not)                                                   */
/* -------------------------------------------------------------------- */

int ML_Comm_Check( ML_Comm *com_ptr )
{
   int ready_flag = 1;

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_Check : Wrong Comm object to check. \n");
      return -1;
   }
   if ( com_ptr->USR_irecvbytes == NULL ) ready_flag = 0;
   if ( com_ptr->USR_sendbytes  == NULL ) ready_flag = 0;
   if ( com_ptr->USR_waitbytes  == NULL ) ready_flag = 0;
   if ( com_ptr->USR_cheapwaitbytes  == NULL ) ready_flag = 0;
   if ( com_ptr->USR_errhandler == NULL ) ready_flag = 0;
   if ( com_ptr->ML_mypid  < 0 )          ready_flag = 0;
   if ( com_ptr->ML_nprocs < 0 )          ready_flag = 0;
   if ( com_ptr->USR_comm == 0 )          ready_flag = 0;
   if ( ready_flag == 1 ) return 0;
   else                   return -1;
}

/* ******************************************************************** */
/* Functions to set communicator parameters                             */
/* -------------------------------------------------------------------- */

int ML_Comm_Set_UsrComm( ML_Comm *com_ptr, USR_COMM com )
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_Set_UsrComm : Wrong object (%d).\n",com_ptr->ML_id);
      exit(1);
   }
   com_ptr->USR_comm = com;
   return 0;
}

/* -------------------------------------------------------------------- */

int ML_Comm_Set_Mypid( ML_Comm *com_ptr, int mypid )
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_Set_Mypid : Wrong object (%d).\n",com_ptr->ML_id);
      exit(1);
   }
   com_ptr->ML_mypid = mypid;
   return 0;
}

/* -------------------------------------------------------------------- */

int ML_Comm_Set_Nprocs( ML_Comm *com_ptr, int nprocs )
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_Set_Nprocs : Wrong object (%d).\n",com_ptr->ML_id);
      exit(1);
   }
   com_ptr->ML_nprocs = nprocs;
   return 0;
}

/* -------------------------------------------------------------------- */

int ML_Comm_Set_SendFcn( ML_Comm *com_ptr, int (*func)(void*,unsigned int,int,int,USR_COMM))
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_Set_SendFcn : Wrong object (%d).\n",com_ptr->ML_id);
      exit(1);
   }
   com_ptr->USR_sendbytes = func;
   return 0;
}

/* -------------------------------------------------------------------- */

int ML_Comm_Set_RecvFcn( ML_Comm *com_ptr, int (*func)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*))
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_Set_RecvFcn : Wrong object (%d).\n",com_ptr->ML_id);
      exit(1);
   }
   com_ptr->USR_irecvbytes = func;
   return 0;
}

/* -------------------------------------------------------------------- */

int ML_Comm_Set_WaitFcn( ML_Comm *com_ptr, int (*func)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*))
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_Set_WaitFcn : Wrong object (%d).\n",com_ptr->ML_id);
      exit(1);
   }
   com_ptr->USR_waitbytes = func;
   return 0;
}

/************************************************************************/
/* find the maximum integer among the ones in each processor            */
/* -------------------------------------------------------------------- */

int ML_Comm_GmaxInt(ML_Comm *com_ptr, int idata)

{
   int     mask, partner, hbit, msgtype, msgbase=245, nprocs, mypid;
   int     i, k, indata, outdata;
   USR_REQ Request;

#ifdef ML_USEMPIFUNCTIONS
   return ML_gmax_int(idata, com_ptr);
#endif

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GmaxInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   outdata = idata;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = sizeof(int);
            com_ptr->USR_irecvbytes((void*) &indata, k, &partner,
                                    &msgtype, com_ptr->USR_comm,
#ifdef ML_CPP
                                    &Request );
#else
                                    (void *)&Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &indata, k, &partner, 
                                   &msgtype, com_ptr->USR_comm, 
#ifdef ML_CPP
                                   &Request );
#else
                                   (void *)&Request );
#endif
            if ( indata > outdata ) outdata = indata;
         }
         else if (partner < nprocs)
         {
            k = sizeof(int);
            com_ptr->USR_sendbytes((void*) &outdata, k, partner,
                                   msgtype, com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 538;
   mask    = 32767;
   k       = sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
         }
      }
   }
   return outdata;
}

/************************************************************************/
/* find the sum of integers residing in each processor                  */
/* -------------------------------------------------------------------- */

int ML_Comm_GsumInt(ML_Comm *com_ptr, int idata)
{
   int     mask, partner, hbit, msgtype, msgbase=246;
   int     i, k, indata, outdata, nprocs, mypid;
   USR_REQ Request;

#ifdef ML_USEMPIFUNCTIONS
  MPI_Allreduce((void *) &idata,(void *) &i, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  return i;
#endif

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GsumInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   outdata = idata;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = sizeof(int);
            com_ptr->USR_irecvbytes((void*) &indata, k, &partner, &msgtype,
#ifdef ML_CPP
                                    com_ptr->USR_comm, &Request );
#else
                                    com_ptr->USR_comm, (void *)&Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &indata, k, &partner, &msgtype, 
#ifdef ML_CPP
                                    com_ptr->USR_comm, &Request );
#else
                                    com_ptr->USR_comm, (void *)&Request );
#endif
            outdata = outdata + indata;
         }
         else if (partner < nprocs)
         {
            k = sizeof(int);
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 539;
   mask    = 32767;
   k       = sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                    com_ptr->USR_comm, &Request );
#else
                                    com_ptr->USR_comm, (void *)&Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                    com_ptr->USR_comm, &Request );
#else
                                    com_ptr->USR_comm, (void *)&Request );
#endif
         }
      }
   }
   return outdata;
}

/************************************************************************/
/* find the sum of doubles  residing in each processor                  */
/* -------------------------------------------------------------------- */

double ML_Comm_GsumDouble(ML_Comm *com_ptr, double ddata)
{
   int     mask, partner, hbit, msgtype, msgbase=247;
   int     i, k, nprocs, mypid;
   double  indata, outdata;
   USR_REQ Request;
#ifdef ML_USEMPIFUNCTIONS
   return ML_gsum_double(ddata, com_ptr);
#endif

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GsumDouble : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapse (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   outdata = ddata;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = sizeof(double);
            com_ptr->USR_irecvbytes((void*) &indata, k, &partner, &msgtype,
#ifdef ML_CPP
                                    com_ptr->USR_comm, &Request );
#else
                                    com_ptr->USR_comm, (void *)&Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &indata, k, &partner, &msgtype, 
#ifdef ML_CPP
                                    com_ptr->USR_comm, &Request );
#else
                                    com_ptr->USR_comm, (void *)&Request );
#endif
            outdata = outdata + indata;
         }
         else if (partner < nprocs)
         {
            k = sizeof(double);
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 539;
   mask    = 32767;
   k       = sizeof(double);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
         }
      }
   }
   return outdata;
}

/************************************************************************/
/* find the max of doubles  residing in each processor                  */
/* -------------------------------------------------------------------- */

double ML_Comm_GmaxDouble(ML_Comm *com_ptr, double ddata)
{
   int     mask, partner, hbit, msgtype, msgbase=247;
   int     i, k, nprocs, mypid;
   double  indata, outdata;
   USR_REQ Request;

#ifdef ML_USEMPIFUNCTIONS
   return ML_gmax_double(ddata, com_ptr);
#endif

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GmaxDouble : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapse (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   outdata = ddata;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = sizeof(double);
            com_ptr->USR_irecvbytes((void*) &indata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &indata, k, &partner, &msgtype, 
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            outdata = (outdata > indata) ? outdata : indata;
         }
         else if (partner < nprocs)
         {
            k = sizeof(double);
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 539;
   mask    = 32767;
   k       = sizeof(double);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) &outdata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
         }
      }
   }
   return outdata;
}

/************************************************************************/
/* sum an integer vector across all processors                          */
/*----------------------------------------------------------------------*/

int ML_Comm_GsumVecInt(ML_Comm *com_ptr, int *idata, int *itmp, int leng)
{
   int     mask, partner, hbit, msgtype, msgbase=247;
   int     i, j, k, nprocs, mypid;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GsumVecInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;

   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = leng * sizeof(int);
            com_ptr->USR_irecvbytes((void*) itmp, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) itmp, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            for ( j = 0; j < leng; j++ ) idata[j] += itmp[j];
         }
         else if (partner < nprocs)
         {
            k = leng * sizeof(int);
            com_ptr->USR_sendbytes((void*) idata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 540;
   mask    = 32767;
   k       = leng * sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) idata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) idata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) idata, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
         }
      }
   }
   return 0;
}

/************************************************************************/
/* This is a modification from Tuminaro's AZ_gappend_int subroutine.    */
/* The modification is done so that the data are arranged according to  */
/* processor number.                                                    */
/*----------------------------------------------------------------------*/

int ML_Comm_GappendInt(ML_Comm *com_ptr, int *vals, int *cur_length,
                    int total_length)
{
   int     mask, partner, hbit, msgtype, msgbase=145;
   int     i, k, nbytes, mypid, nprocs;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GappendInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = (total_length - (*cur_length)) * sizeof(int);
            com_ptr->USR_irecvbytes((void*)&(vals[*cur_length]), k,
                                   &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            nbytes = com_ptr->USR_waitbytes((void*)&(vals[*cur_length]), k,
                                   &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            (*cur_length) += (nbytes / sizeof(int));
         }
         else if (partner < nprocs)
         {
            k = (*cur_length) * sizeof(int);
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 438;
   mask    = 32767;
   k       = total_length * sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) vals, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_waitbytes((void*) vals, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
         }
      }
   }
   (*cur_length) = total_length;
   return 0;
}



/************************************************************************/
/* This is a modification from Tuminaro's AZ_gappend_int subroutine.    */
/* The modification is done so that the data are arranged according to  */
/* processor number.                                                    */
/* This is the same function as ML_Comm_GappendInt, except the          */
/* parameter vals is of type ml_big_int                                 */
/*----------------------------------------------------------------------*/

int ML_Comm_GappendBigInt(ML_Comm *com_ptr, ml_big_int *vals, int *cur_length,
                    int total_length)
{
   int     mask, partner, hbit, msgtype, msgbase=145;
   int     i, k, nbytes, mypid, nprocs;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GappendInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = (total_length - (*cur_length)) * sizeof(vals[0]);
            com_ptr->USR_irecvbytes((void*)&(vals[*cur_length]), k,
                                   &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            nbytes = com_ptr->USR_waitbytes((void*)&(vals[*cur_length]), k,
                                   &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            (*cur_length) += (nbytes / sizeof(vals[0]));
         }
         else if (partner < nprocs)
         {
            k = (*cur_length) * sizeof(vals[0]);
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 438;
   mask    = 32767;
   k       = total_length * sizeof(vals[0]);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) vals, k, &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            com_ptr->USR_cheapwaitbytes((void*) vals, k, &partner, &msgtype, 
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
         }
      }
   }
   (*cur_length) = total_length;
   return 0;
}




/************************************************************************/
/* This is a modification from Tuminaro's AZ_gappend_double subroutine. */
/* The modification is done so that the data are arranged according to  */
/* processor number.                                                    */
/*----------------------------------------------------------------------*/

int ML_Comm_GappendDouble(ML_Comm *com_ptr, double *vals, int *cur_length,
                            int total_length)
{
   int     mask, partner, hbit, msgtype, msgbase=245;
   int     i, k, nbytes, nprocs, mypid;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GappendDouble : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = (total_length - (*cur_length)) * sizeof(double);
            com_ptr->USR_irecvbytes((void*)&(vals[*cur_length]), k,
                                    &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            nbytes = com_ptr->USR_waitbytes((void*)&(vals[*cur_length]), k,
                                    &partner, &msgtype,
#ifdef ML_CPP
                                     com_ptr->USR_comm, &Request );
#else
                                     com_ptr->USR_comm, (void *) &Request );
#endif
            (*cur_length) += (nbytes / sizeof(double));
         }
         else if (partner < nprocs)
         {
            k = (*cur_length) * sizeof(double);
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 538;
   mask    = 32767;
   k       = total_length * sizeof(double);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) vals, k, &partner, &msgtype,
                                    com_ptr->USR_comm, (void *) &Request );
            com_ptr->USR_cheapwaitbytes((void*) vals, k, &partner, &msgtype, 
                                    com_ptr->USR_comm, (void *) &Request );
         }
      }
   }
   (*cur_length) = total_length;
   return 0;
}

/**************************************************************************/
/* communication subroutines                                              */
/*------------------------------------------------------------------------*/

int ML_Comm_Irecv(void* buf, unsigned int count, int *src,
                  int *mid, USR_COMM comm, USR_REQ *request )
{
   int err = 0;
#ifdef ML_MPI
   if (*mid == -1) {
     *mid = MPI_ANY_TAG;
     /* bogus code to avoid warnings */
     if (*src == -59)  ML_avoid_unused_param((void *) &comm);
   }
   if (*src == -1) *src = MPI_ANY_SOURCE;
   err = MPI_Irecv(buf,count,MPI_BYTE,*src,*mid,MPI_COMM_WORLD,request);
#else
   /* bogus code to avoid warnings */
   if (*mid == -59) {
     ML_avoid_unused_param(buf);
     ML_avoid_unused_param((void *) &count);
     ML_avoid_unused_param((void *) src);
     ML_avoid_unused_param((void *) &comm);
     ML_avoid_unused_param((void *) request);
   }
#endif
   return err;
}

/*------------------------------------------------------------------------*/

int ML_Comm_Wait (void* buf, unsigned int count, int *src,
                  int *mid, USR_COMM comm, USR_REQ *request )
{
   int        return_cnt = 0;
#ifdef ML_MPI
   MPI_Status status;
   MPI_Wait(request, &status);
   MPI_Get_count(&status, MPI_BYTE, &return_cnt);
   *src = status.MPI_SOURCE;
   *mid = status.MPI_TAG;
   /* bogus code to avoid warnings */
   if (*mid == -59) {
     ML_avoid_unused_param(buf);
     ML_avoid_unused_param((void *) &count);
     ML_avoid_unused_param((void *) &comm);
   }
#else
   /* bogus code to avoid warnings */
   if (*mid == -59) {
     ML_avoid_unused_param(buf);
     ML_avoid_unused_param((void *) &count);
     ML_avoid_unused_param((void *) src);
     ML_avoid_unused_param((void *) &comm);
     ML_avoid_unused_param((void *) request);
   }
#endif
   return return_cnt;
}

/*------------------------------------------------------------------------*/

/* Identical to ML_Comm_Wait, but  the message length is not calculated. */

void ML_Comm_CheapWait (void* buf, unsigned int count, int *src,
                  int *mid, USR_COMM comm, USR_REQ *request )
{
#ifdef ML_MPI
   MPI_Status status;
   MPI_Wait(request, &status);
   *src = status.MPI_SOURCE;
   *mid = status.MPI_TAG;
   /* bogus code to avoid warnings */
   if (*mid == -59) {
     ML_avoid_unused_param(buf);
     ML_avoid_unused_param((void *) &count);
     ML_avoid_unused_param((void *) &comm);
   }
#else
   /* bogus code to avoid warnings */
   if (*mid == -59) {
     ML_avoid_unused_param(buf);
     ML_avoid_unused_param((void *) &count);
     ML_avoid_unused_param((void *) src);
     ML_avoid_unused_param((void *) &comm);
     ML_avoid_unused_param((void *) request);
   }
#endif
}

/*------------------------------------------------------------------------*/

int ML_Comm_Send(void* buf, unsigned int count, int dest, int mid,
                 USR_COMM comm )
{
   int err = 0;
#ifdef ML_MPI
   err = MPI_Send(buf, count, MPI_BYTE, dest, mid, MPI_COMM_WORLD);
   /* bogus code to avoid warnings */
   if (mid == -59) {
     ML_avoid_unused_param((void *) &comm);
   }
#else
   /* bogus code to avoid warnings */
   if (mid == -59) {
     ML_avoid_unused_param(buf);
     ML_avoid_unused_param((void *) &count);
     ML_avoid_unused_param((void *) &dest);
     ML_avoid_unused_param((void *) &comm);
   }
#endif

   return err;
}

#if defined(ML_MPI) && defined(ML_CATCH_MPI_ERRORS_IN_DEBUGGER)

/*------------------------------------------------------------------------------

 Briefly, this function associates *errhandler with the communicator comm.

 Sets the behavior of message-passing error trapping.  For MPI, the default
 behavior is to exit if an error is detected.   This can be changed
 to pass errors back to the calling function via return code, or to use a
 user-defined error handling function.

 input:
    comm            communicator (for MPI, type MPI_Comm)
    errhandler      specifies how errors should be handled
                    (for MPI, type MPI_Errhandler)

 output:
    err             error code


 MPI-specific options for errhandler:

    MPI_ERRORS_ARE_FATAL        (default)
    MPI_ERRORS_RETURN

    MPE_Errors_call_dbx_in_xterm (mpich only)
    MPE_Signals_call_debugger    (mpich only)

 Note!  There may very well be other implementation-specific options for
 errhandler.

------------------------------------------------------------------------------*/

int ML_Comm_ErrorHandlerSet(USR_COMM comm, USR_ERRHANDLER *errhandler)
{
   int err = 0;
#ifdef ML_MPI
   /*err = MPI_Errhandler_set(comm, *errhandler);*/
   err = USR_ERRHANDLER_SET(comm, *errhandler);
#endif
   return err;
}

/*------------------------------------------------------------------------------
 Wrapper for registration of MPI error-handling function.
------------------------------------------------------------------------------*/

int ML_Comm_ErrorHandlerCreate(void (*fcn)(USR_COMM*,int*),
/*int ML_Comm_ErrorHandlerCreate(USR_ERRHANDLER_FUNCTION *fcn,*/
                               USR_ERRHANDLER *errhandler)
{
   int err = 0;
#ifdef ML_MPI
   /*err = USR_ERRHANDLER_CREATE(fcn, errhandler);*/
   err = USR_ERRHANDLER_CREATE(fcn, errhandler);
#endif
   return err;
}

/*------------------------------------------------------------------------------
 Wrapper for destruction of handle to MPI error-handling function.
------------------------------------------------------------------------------*/

int ML_Comm_ErrorHandlerDestroy(USR_ERRHANDLER **errhandler)
{
  int err = 0;
#ifdef ML_MPI
  err = MPI_Errhandler_free(*errhandler);
#endif
  *errhandler = NULL;
  return err; }

/*------------------------------------------------------------------------------
 ML's very own error-handling routine for MPI.
 Useful if trying to trap error with a debugger.
------------------------------------------------------------------------------*/

void ML_Comm_ErrorHandler(USR_COMM *comm, int *error_code)
{
#ifdef ML_MPI
   int message_size;
   char message[MPI_MAX_ERROR_STRING];

   MPI_Error_string(*error_code,message,&message_size);
   fprintf(stderr,"(In ML_Comm_ErrorHandler) MPI error: %s\n",message);
#endif
   /* we call abort so we can trap the error with a debugger */
   abort();
}

#endif /* if defined(ML_MPI) && defined(ML_CATCH_MPI_ERRORS_IN_DEBUGGER) */
