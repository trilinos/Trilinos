/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for manipulating the ML_CommInfo_AGX object                */
/* (used in the automatic grid transfer structure)                      */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : February, 1999                                       */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_comminfoagx.h"

/* ******************************************************************** */
/* Create a local communication buffer object                           */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Create(ML_CommInfoAGX **combuf)
{
   ML_CommInfoAGX *com;

   ML_memory_alloc((void**) combuf, sizeof(ML_CommInfoAGX),"CI1");
   com            = (*combuf);
   com->ML_id     = ML_ID_COMMINFOAGX;
   com->proc_id   = 0;
   com->send_cur  = 0;
   com->send_cnt  = 0;
   com->recv_cur  = 0;
   com->recv_cnt  = 0;
   com->send_proc = 0;
   com->send_ia   = 0;
   com->send_list = 0;
   com->recv_proc = 0;
   com->recv_ia   = 0;
   com->recv_list = 0;
   com->recv_xyz  = 0;
   return 0;
}

/* ******************************************************************** */
/* Deallocate the storage and destroy this object                       */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Destroy(ML_CommInfoAGX **combuf)
{
   ML_CommInfoAGX *com;

   com = (*combuf);
   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Destroy : wrong object. \n");
      exit(1);
   }
   if (com->send_cnt > 0)
   {
      if (com->send_proc != 0) ML_memory_free((void**)&(com->send_proc));
      if (com->send_list != 0) ML_memory_free((void**)&(com->send_list));
   }
   if (com->recv_cnt > 0)
   {
      if (com->recv_proc != 0) ML_memory_free((void**)&(com->recv_proc));
      if (com->recv_list != 0) ML_memory_free((void**)&(com->recv_list));
   }
   if (com->send_ia  != 0) ML_memory_free((void**)&(com->send_ia));
   if (com->recv_ia  != 0) ML_memory_free((void**)&(com->recv_ia));
   if (com->recv_xyz != 0) ML_memory_free((void**)&(com->recv_xyz));
   com->ML_id = -1;
   ML_memory_free( (void **) combuf );
   return 0;
}

/* ******************************************************************** */
/* Sets up buffer space to store the send information                   */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Setup_Send(ML_CommInfoAGX *com, int count, int count2)
{
   int nbytes;

   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Setup_Send : wrong object. \n");
      exit(1);
   }
   com->send_cur = 0;
   com->send_cnt = count;
   nbytes        = ( count + 1 ) * sizeof(int);
   ML_memory_alloc( (void**) &(com->send_ia), (unsigned int) nbytes, "CB1" );
   if ( count > 0 )
   {
      nbytes = count * sizeof(int);
      ML_memory_alloc( (void**) &(com->send_proc), (unsigned int) nbytes, "CB2" );
   }
   else
   {
      com->send_proc  = 0;
   }
   if ( count2 > 0 )
   {
      nbytes = count2 * sizeof(ml_big_int);
      ML_memory_alloc( (void**) &(com->send_list), (unsigned int) nbytes, "CB3" );
   }
   else
   {
      com->send_list  = 0;
   }
   com->send_ia[0] = 0;
   return 0;
}

/* ******************************************************************** */
/*    Load a destination ID and a list (which are indices to a data     */
/*    array) to the structure.                                          */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Load_SendList(ML_CommInfoAGX *com,int proc,int leng,
                                 ml_big_int *list)
{
   int  i, k, offset;

   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Load_SendList : wrong object. \n");
      exit(1);
   }
   k = (com->send_cur)++;
   com->send_proc[k] = proc;
   offset = com->send_ia[k];
   for ( i = 0; i < leng; i++ ) com->send_list[offset+i] = list[i];
   com->send_ia[k+1] = offset + leng;
   return 0;
}

/* ******************************************************************** */
/*    Given an index to the array, return the send information.         */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Get_SendList(ML_CommInfoAGX *com, int index, int *proc,
                              int *leng, ml_big_int **list)
{
   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Get_SendList : wrong object. \n");
      exit(1);
   }
   (*proc) = com->send_proc[index];
   (*leng) = com->send_ia[index+1] - com->send_ia[index];
   (*list) = &(com->send_list[com->send_ia[index]]);
   return 0;
}

/* ******************************************************************** */
/*  Sets up buffer space for receive information                        */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Setup_Recv(ML_CommInfoAGX *com, int count, int count2)
{
   int nbytes;

   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Setup_Recv : wrong object. \n");
      exit(1);
   }
   com->recv_cur = 0;
   com->recv_cnt = count;
   nbytes        = ( count + 1 ) * sizeof(int);
   ML_memory_alloc( (void **) &(com->recv_ia), (unsigned int) nbytes, "CB4" );
   if ( count > 0 )
   {
      nbytes = count * sizeof(int);
      ML_memory_alloc( (void **) &(com->recv_proc), (unsigned int) nbytes, "CB5" );
   }
   else
   {
      com->recv_proc = 0;
   }
   if ( count2 > 0 )
   {
      nbytes = count2 * sizeof(ml_big_int);
      ML_memory_alloc( (void **) &(com->recv_list), (unsigned int) nbytes, "CB6" );
      nbytes = 3 * count2 * sizeof(double);
      ML_memory_alloc( (void **) &(com->recv_xyz), (unsigned int) nbytes, "CB7" );
   }
   else
   {
      com->recv_list = 0;
      com->recv_xyz  = 0;
   }
   com->recv_ia[0] = 0;
   return 0;
}

/* ******************************************************************** */
/*    Load a source ID and the corresponding receive length to the      */
/*    structure.                                                        */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Load_RecvInfo( ML_CommInfoAGX *com, int proc, int leng )
{
   int  k;

   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Load_RecvInfo : wrong object. \n");
      exit(1);
   }
   k = (com->recv_cur)++;
   com->recv_proc[k] = proc;
   com->recv_ia[k+1] = com->recv_ia[k] + leng;
   return 0;
}

/* ******************************************************************** */
/*    Load the receive integer and double data (coordinates) to the     */
/*    structure.                                                        */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Load_RecvData(ML_CommInfoAGX *com, int proc, ml_big_int *list,
                                double *x, double *y, double *z)
{
   int  i, k = 0, flag = 0, j = 0;

   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Load_RecvData : wrong object. \n");
      exit(1);
   }
   while (k < com->recv_cnt && flag == 0)
   {
      if (com->recv_proc[k] == proc) flag = 1;
      else                           k++;
   }
   for ( i =com->recv_ia[k]; i < com->recv_ia[k+1]; i++ )
   {
      com->recv_list[i]    = list[j];
      com->recv_xyz[3*i]   = x[j];
      com->recv_xyz[3*i+1] = y[j];
      com->recv_xyz[3*i+2] = z[j++];
   }
   return 0;
}

/* ******************************************************************** */
/*    Given an index to the array, return the recv information.         */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Get_RecvList(ML_CommInfoAGX *com, int index, int *proc,
                            int *leng, ml_big_int **list)
{
   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Get_RecvList : wrong object. \n");
      exit(1);
   }
   (*proc) = com->recv_proc[index];
   (*leng) = com->recv_ia[index+1] - com->recv_ia[index];
   (*list) = &(com->recv_list[com->recv_ia[index]]);
   return 0;
}

/* ******************************************************************** */
/*    This function print all information about the communication       */
/*    buffer standard output device.                                    */
/* -------------------------------------------------------------------- */

int ML_CommInfoAGX_Print(ML_CommInfoAGX *com)
{
   int  i, j, leng;

   if ( com->ML_id != ML_ID_COMMINFOAGX )
   {
      printf("ML_CommInfoAGX_Print : wrong object. \n");
      exit(1);
   }
   printf("ML_CommInfoAGX : number of destinations = %d \n",com->send_cnt);
   for ( i = 0; i < com->send_cnt; i++ )
   {
      leng = com->send_ia[i+1] - com->send_ia[i];
      printf("   To : %d , leng = %d \n", com->send_proc[i], leng);
      if (com->send_list != 0)
      {
         for (j=com->send_ia[i]; j<com->send_ia[i+1]; j++)
            printf("    index = %d \n", com->send_list[j]);
      }
   }
   printf("ML_CommInfoAGX : number of sources = %d \n",com->recv_cnt);
   for ( i = 0; i < com->recv_cnt; i++ )
   {
      leng = com->recv_ia[i+1] - com->recv_ia[i];
      printf("   From : %d , leng = %d \n", com->recv_proc[i], leng);
      if (com->recv_list != 0)
      {
         for (j=com->recv_ia[i]; j<com->recv_ia[i+1]; j++)
            printf("    index = %d \n", com->recv_list[j]);
      }
   }
   return 0;
}

