/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "ExternalRows_dh.h"
#include "Factor_dh.h"
#include "Euclid_dh.h"
#include "SubdomainGraph_dh.h"
#include "Mem_dh.h"
#include "Parser_dh.h"
#include "Hash_dh.h"

 /* tags for MPI comms */
enum{ ROW_CT_TAG, NZ_CT_TAG, ROW_LENGTH_TAG, ROW_NUMBER_TAG,
      CVAL_TAG, FILL_TAG, AVAL_TAG };

#undef __FUNC__
#define __FUNC__ "ExternalRows_dhCreate"
void ExternalRows_dhCreate(ExternalRows_dh *er)
{
  START_FUNC_DH
  struct _extrows_dh* tmp = (struct _extrows_dh*)MALLOC_DH(sizeof(struct _extrows_dh)); CHECK_V_ERROR;
  *er = tmp;

  if (MAX_MPI_TASKS < np_dh) {
    SET_V_ERROR("MAX_MPI_TASKS is too small; change, then recompile!");
  }

  { int i;
    for (i=0; i<MAX_MPI_TASKS; ++i) {
      tmp->rcv_row_lengths[i] = NULL;
      tmp->rcv_row_numbers[i] = NULL;
    }
  }

  tmp->cvalExt = NULL;
  tmp->fillExt = NULL;
  tmp->avalExt = NULL;
  tmp->my_row_counts = NULL;
  tmp->my_row_numbers = NULL;
  tmp->cvalSend = NULL;
  tmp->fillSend = NULL;
  tmp->avalSend = NULL;
  tmp->rowLookup = NULL;
  tmp->sg = NULL;
  tmp->F = NULL;
  tmp->debug = Parser_dhHasSwitch(parser_dh, "-debug_ExtRows");
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "ExternalRows_dhDestroy"
void ExternalRows_dhDestroy(ExternalRows_dh er)
{
  START_FUNC_DH
  int i;

  for (i=0; i<MAX_MPI_TASKS; ++i) {
    if (er->rcv_row_lengths[i] != NULL) {
      FREE_DH(er->rcv_row_lengths[i]); CHECK_V_ERROR;
    }
    if (er->rcv_row_numbers[i] != NULL) {
      FREE_DH(er->rcv_row_numbers[i]); CHECK_V_ERROR;
    }
  }

  if (er->cvalExt != NULL) { FREE_DH(er->cvalExt); CHECK_V_ERROR; }
  if (er->fillExt != NULL) { FREE_DH(er->fillExt); CHECK_V_ERROR; }
  if (er->avalExt != NULL) { FREE_DH(er->avalExt); CHECK_V_ERROR; }

  if (er->my_row_counts != NULL) { FREE_DH(er->my_row_counts); CHECK_V_ERROR; }
  if (er->my_row_numbers != NULL) { FREE_DH(er->my_row_numbers); CHECK_V_ERROR; }

  if (er->cvalSend != NULL) { FREE_DH(er->cvalSend); CHECK_V_ERROR; }
  if (er->fillSend != NULL) { FREE_DH(er->fillSend); CHECK_V_ERROR; }
  if (er->avalSend != NULL) { FREE_DH(er->avalSend); CHECK_V_ERROR; }

  if (er->rowLookup != NULL) { Hash_dhDestroy(er->rowLookup); CHECK_V_ERROR; }
  FREE_DH(er); CHECK_V_ERROR;
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "ExternalRows_dhInit"
void ExternalRows_dhInit(ExternalRows_dh er, Euclid_dh ctx)
{
  START_FUNC_DH
  er->sg = ctx->sg;
  er->F = ctx->F;
  END_FUNC_DH
}

/*=====================================================================
 * method for accessing external rows
 *=====================================================================*/

#undef __FUNC__
#define __FUNC__ "ExternalRows_dhGetRow"
void ExternalRows_dhGetRow(ExternalRows_dh er, int globalRow,
                            int *len, int **cval, int **fill, REAL_DH **aval)
{
  START_FUNC_DH
  if (er->rowLookup == NULL) {
    *len = 0;
  } 

  else {
    HashData *r = NULL;
    r = Hash_dhLookup(er->rowLookup, globalRow); CHECK_V_ERROR;
    if (r != NULL) {
      *len = r->iData;
      if (cval != NULL) *cval = r->iDataPtr;
      if (fill != NULL) *fill = r->iDataPtr2;
      if (aval != NULL) *aval = r->fDataPtr;
    } else {
      *len = 0;
    }
  }
  END_FUNC_DH
}

/*=====================================================================
 * methods for receiving  external rows from lower-ordered subdomains
 *=====================================================================*/
static void rcv_ext_storage_private(ExternalRows_dh extRows); 
static void build_hash_table_private(ExternalRows_dh er); 
static void rcv_external_rows_private(ExternalRows_dh er);
static void allocate_ext_row_storage_private(ExternalRows_dh er);
static void print_received_rows_private(ExternalRows_dh er);

#undef __FUNC__
#define __FUNC__ "ExternalRows_dhRecvRows"
void ExternalRows_dhRecvRows(ExternalRows_dh er)
{
  START_FUNC_DH
  bool debug = false;
  if (logFile != NULL && er->debug) debug = true;

  if (er->sg->loCount > 0) {
    /* get number of rows and length of each row to be received
       from each lower ordered nabor.
       (allocates: *rcv_row_lengths[], *rcv_row_numbers[])
    */
    rcv_ext_storage_private(er); CHECK_V_ERROR;


    /* allocate data structures for receiving the rows (no comms) 
       (allocates: cvalExt, fillExt, avalExt)
       (no communications)
     */
    allocate_ext_row_storage_private(er); CHECK_V_ERROR;


    /* construct hash table for  external row lookup (no comms) 
       (Creates/allocates: rowLookup)
       (no communications)
     */
    build_hash_table_private(er); CHECK_V_ERROR;

    /* receive the actual row structures and values 
       from lower ordered neighbors 
     */
    rcv_external_rows_private(er); CHECK_V_ERROR;

    if (debug) {
      print_received_rows_private(er); CHECK_V_ERROR;
    }
  } 
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "rcv_ext_storage_private"
void rcv_ext_storage_private(ExternalRows_dh er)
{
  START_FUNC_DH
  int i; 
  int loCount = er->sg->loCount, *loNabors = er->sg->loNabors;
  int *rcv_row_counts = er->rcv_row_counts;
  int *rcv_nz_counts = er->rcv_nz_counts;
  int **lengths = er->rcv_row_lengths, **numbers = er->rcv_row_numbers;
  bool debug = false;

  if (logFile != NULL && er->debug) debug = true;

  /* get number of rows, and total nonzeros, that each lo-nabor will send */
  for (i=0; i<loCount; ++i) {
    int nabor = loNabors[i];
    MPI_Irecv(rcv_row_counts+i, 1, MPI_INT, nabor, ROW_CT_TAG, comm_dh, er->req1+i);
    MPI_Irecv(rcv_nz_counts+i,  1, MPI_INT, nabor, NZ_CT_TAG,  comm_dh, er->req2+i);
  }    
  MPI_Waitall(loCount, er->req1, er->status);
  MPI_Waitall(loCount, er->req2, er->status);

  if (debug) {
    fprintf(logFile, "\nEXR rcv_ext_storage_private:: <nabor,rowCount,nzCount>\nEXR ");
    for (i=0; i<loCount; ++i) {
      fprintf(logFile, "<%i,%i,%i> ", loNabors[i], rcv_row_counts[i], rcv_nz_counts[i]);
    }
  }

  /* get lengths and global number of each row to be received */
  for (i=0; i<loCount; ++i) {
    int nz = rcv_nz_counts[i];
    int nabor = loNabors[i];
    lengths[i] =  (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
    numbers[i] =  (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
    MPI_Irecv(lengths[i], nz, MPI_INT, nabor, ROW_LENGTH_TAG, comm_dh, er->req1+i);
    MPI_Irecv(numbers[i], nz, MPI_INT, nabor, ROW_NUMBER_TAG, comm_dh, er->req2+i);
  }
  MPI_Waitall(loCount, er->req1, er->status);
  MPI_Waitall(loCount, er->req2, er->status);

  if (debug) {
    int j, nz;
    for (i=0; i<loCount; ++i) { 
      fprintf(logFile, "\nEXR rows <number,length> to be received from P_%i\nEXR ", loNabors[i]);
      nz = rcv_row_counts[i];
      for (j=0; j<nz; ++j) fprintf(logFile, "<%i,%i> ", numbers[i][j], lengths[i][j]);
      fprintf(logFile, "\n");
    }
  }

  END_FUNC_DH
}

/* allocates: cvalExt, fillExt, avalExt */
#undef __FUNC__
#define __FUNC__ "allocate_ext_row_storage_private"
void allocate_ext_row_storage_private(ExternalRows_dh er)
{
  START_FUNC_DH
  int i, nz = 0;
  int loCount = er->sg->loCount;
  int *rcv_nz_counts = er->rcv_nz_counts;

  /* count total number of nonzeros to be received */
  for (i=0; i<loCount; ++i) nz += rcv_nz_counts[i];

  /* allocate buffers */
  er->cvalExt = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  er->fillExt = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  er->avalExt = (REAL_DH*)MALLOC_DH(nz*sizeof(REAL_DH)); CHECK_V_ERROR;
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "build_hash_table_private"
void build_hash_table_private(ExternalRows_dh er)
{
  START_FUNC_DH
  int loCount = er->sg->loCount;
  int i, j, offset, rowCt = 0;
  Hash_dh table;
  HashData record;
  int *extRowCval = er->cvalExt, *extRowFill = er->fillExt;
  REAL_DH *extRowAval = er->avalExt;
  int *rcv_row_counts = er->rcv_row_counts;
  int **rcv_row_numbers = er->rcv_row_numbers;
  int **rcv_row_lengths = er->rcv_row_lengths;

  /* count total number of rows to be received */
  for (i=0; i<loCount; ++i) rowCt += rcv_row_counts[i];

  /* build table for looking up external rows */
  Hash_dhCreate(&table, rowCt); CHECK_V_ERROR;
  er->rowLookup = table;
  offset = 0;

  /* loop over lower ordered nabors in subdomain graph */
  for (i=0; i<loCount; ++i) { 

    /* number of rows to be received from nabor(i) */
    int rowCount = rcv_row_counts[i]; 

    /* loop over rows to be received from nabor(i) */
    for (j=0; j<rowCount; ++j) { 

      /* insert a record to locate row(j) in the hash table */
      int row = rcv_row_numbers[i][j];
      int rowLength = rcv_row_lengths[i][j];
      record.iData     = rowLength;
      record.iDataPtr  = extRowCval + offset;
      record.iDataPtr2 = extRowFill + offset;
      record.fDataPtr  = extRowAval + offset;
      Hash_dhInsert(table, row, &record); CHECK_V_ERROR;
      offset += rowLength;
    }
  }
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "rcv_external_rows_private"
void rcv_external_rows_private(ExternalRows_dh er)
{
  START_FUNC_DH
  int *rcv_nz_counts = er->rcv_nz_counts;
  int i, loCount = er->sg->loCount, *loNabors = er->sg->loNabors;
  int nabor, nz = 0, offset = 0;
  int *extRowCval = er->cvalExt, *extRowFill = er->fillExt;
  double *extRowAval = er->avalExt;

  /* start receives of external rows */
  nz = 0;
  for (i=0; i<loCount; ++i) {
    nabor = loNabors[i];
    nz = rcv_nz_counts[i];
    MPI_Irecv(extRowCval+offset, nz, MPI_INT,    nabor, CVAL_TAG, comm_dh, er->req1+i);
    MPI_Irecv(extRowFill+offset, nz, MPI_INT,    nabor, FILL_TAG, comm_dh, er->req2+i);
    MPI_Irecv(extRowAval+offset, nz, MPI_DOUBLE, nabor, AVAL_TAG, comm_dh, er->req3+i);
    offset += nz;
  }

  /* wait for external rows to arrive */
  MPI_Waitall(loCount, er->req1, er->status);
  MPI_Waitall(loCount, er->req2, er->status);
  MPI_Waitall(loCount, er->req3, er->status);
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "print_received_rows_private"
void print_received_rows_private(ExternalRows_dh er)
{
  START_FUNC_DH
  bool noValues = (Parser_dhHasSwitch(parser_dh, "-noValues"));
  int i, j, k, rwCt, idx = 0, nabor;
  int loCount = er->sg->loCount, *loNabors = er->sg->loNabors;
  int n = er->F->n;

  fprintf(logFile, "\nEXR ================= received rows, printed from buffers ==============\n");

  /* loop over nabors from whom we received rows */
  for (i=0; i<loCount; ++i) {
    rwCt = er->rcv_row_counts[i];
    nabor = loNabors[i];
    fprintf(logFile, "\nEXR Rows received from P_%i:\n", nabor);

    /* loop over each row to be received from this nabor */
    for (j=0; j<rwCt; ++j) {  
      int rowNum = er->rcv_row_numbers[i][j];  
      int rowLen  = er->rcv_row_lengths[i][j];
      fprintf(logFile, "EXR %i :: ", 1+rowNum);
      for (k=0; k<rowLen; ++k) {
        if (noValues) {
          fprintf(logFile, "%i,%i ; ", er->cvalExt[idx], er->fillExt[idx]);
        } else {
          fprintf(logFile, "%i,%i,%g ; ", er->cvalExt[idx], er->fillExt[idx], er->avalExt[idx]);
        }
        ++idx;
      }
      fprintf(logFile, "\n");
    }
  }

  fprintf(logFile, "\nEXR =============== received rows, printed from hash table =============\n");
  for (i=0; i<n; ++i) {
    int len, *cval, *fill;
    REAL_DH *aval;
    ExternalRows_dhGetRow(er, i, &len, &cval, &fill, &aval); CHECK_V_ERROR;
    if (len > 0) {
      fprintf(logFile, "EXR %i :: ", i+1);
      for (j=0; j<len; ++j) {
        if (noValues) {
          fprintf(logFile, "%i,%i ; ", cval[j], fill[j]);
        } else {
          fprintf(logFile, "%i,%i,%g ; ", cval[j], fill[j], aval[j]);
        }
      }
      fprintf(logFile, "\n");
    }
  }

  END_FUNC_DH
}

/*=====================================================================
 * methods for sending rows to higher ordered nabors in subdomain graph
 *=====================================================================*/

static void send_ext_storage_private(ExternalRows_dh er);
static void send_external_rows_private(ExternalRows_dh er);
static void waitfor_sends_private(ExternalRows_dh er);

#undef __FUNC__
#define __FUNC__ "ExternalRows_dhSendRows"
void ExternalRows_dhSendRows(ExternalRows_dh er)
{
  START_FUNC_DH
  if (er->sg->hiCount > 0) {
    /* send number of rows and length of each row to be sent
       to each higher ordered nabor.
    */
    send_ext_storage_private(er); CHECK_V_ERROR;

    /* send the row's colum indices, fill levels, and values */
    send_external_rows_private(er); CHECK_V_ERROR;

    waitfor_sends_private(er); CHECK_V_ERROR;
  }
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "send_ext_storage_private"
void send_ext_storage_private(ExternalRows_dh er)
{
  START_FUNC_DH
  int nz, i, j;
  int *nzCounts, *nzNumbers;
  int hiCount = er->sg->hiCount, *hiNabors = er->sg->hiNabors;
  int *rp = er->F->rp, *diag = er->F->diag;
  int m = er->F->m;
  int beg_row = er->F->beg_row;
  int rowCount = er->F->bdry_count;  /* number of boundary rows */
  int first_bdry = er->F->first_bdry; 
  bool debug = false;

  if (logFile != NULL && er->debug) debug = true;

  /* allocate storage to hold nz counts for each row */
  nzCounts =  er->my_row_counts = (int*)MALLOC_DH(rowCount*sizeof(int)); CHECK_V_ERROR;
  nzNumbers =  er->my_row_numbers = (int*)MALLOC_DH(rowCount*sizeof(int)); CHECK_V_ERROR;

  /* count nonzeros in upper triangular portion of each boundary row */
  nz = 0;
  for (i=first_bdry, j=0; i<m; ++i, ++j) {
    int tmp = (rp[i+1] - diag[i]);
    nz += tmp;
    nzCounts[j] = tmp;
  }
  er->nzSend = nz;

  if (debug) {
    fprintf(logFile, "EXR send_ext_storage_private:: rowCount = %i\n", rowCount);
    fprintf(logFile, "EXR send_ext_storage_private:: nz Count = %i\n", nz);
  }

  /* send  number of rows, and total nonzeros, to higher ordered nabors */
  for (i=0; i<hiCount; ++i) {
    int nabor = hiNabors[i];
    MPI_Isend(&rowCount, 1, MPI_INT, nabor, ROW_CT_TAG, comm_dh, er->req1+i);
    MPI_Isend(&nz,       1, MPI_INT, nabor, NZ_CT_TAG,  comm_dh, er->req2+i);
  }    

  /* set up array for global row numbers */
  for (i=0, j=first_bdry; j<m; ++i, ++j) {
    nzNumbers[i] = j+beg_row;
  }

  /* start sends of length and global number of each of this processor's
     boundary row to higher ordered nabors; the receiving processor will
     use this information to allocate storage buffers for the actual
     row structures and values.
   */
  for (i=0; i<hiCount; ++i) {
    int nabor = hiNabors[i];
    MPI_Isend(nzNumbers, rowCount, MPI_INT, nabor, ROW_NUMBER_TAG, comm_dh, er->req3+i);
    MPI_Isend(nzCounts,  rowCount, MPI_INT, nabor, ROW_LENGTH_TAG, comm_dh, er->req4+i);
  }

  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "send_external_rows_private"
void send_external_rows_private(ExternalRows_dh er)
{
  START_FUNC_DH
  int i, j, hiCount = er->sg->hiCount, *hiNabors = er->sg->hiNabors;
  int offset, nz = er->nzSend;
  int *cvalSend, *fillSend; 
  REAL_DH *avalSend; 
  int *cval = er->F->cval, *fill = er->F->fill;
  int m = er->F->m;
  int *rp = er->F->rp, *diag = er->F->diag;
  int first_bdry = er->F->first_bdry;
  REAL_DH *aval = er->F->aval;
  bool debug = false;

  if (logFile != NULL && er->debug) debug = true;

  /* allocate buffers to hold upper triangular portion of boundary rows */
  cvalSend = er->cvalSend = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  fillSend = er->fillSend = (int*)MALLOC_DH(nz*sizeof(int)); CHECK_V_ERROR;
  avalSend = er->avalSend = (double*)MALLOC_DH(nz*sizeof(double)); CHECK_V_ERROR;

  /* copy upper triangular portion of boundary rows int send buffers */
  offset = 0;
  for (i=first_bdry, j=0; i<m; ++i, ++j) {
    int tmp = (rp[i+1] - diag[i]);

    memcpy(cvalSend+offset, cval+diag[i], tmp*sizeof(int));
    memcpy(fillSend+offset, fill+diag[i], tmp*sizeof(int));
    memcpy(avalSend+offset, aval+diag[i], tmp*sizeof(double));
    offset += tmp;
  }

  if (debug) {
    int beg_row = er->F->beg_row;
    int idx = 0;
    bool noValues = (Parser_dhHasSwitch(parser_dh, "-noValues"));

    fprintf(logFile, "\nEXR ======================= send buffers ======================\n");

    for (i=first_bdry, j=0; i<m; ++i, ++j) {
      int tmp = (rp[i+1] - diag[i]);
      fprintf(logFile, "EXR %i :: ", i+beg_row);

      for (j=0; j<tmp; ++j) {
        if (noValues) {
          fprintf(logFile, "%i,%i ; ", cvalSend[idx], fillSend[idx]);
        } else {
          fprintf(logFile, "%i,%i,%g ; ", cvalSend[idx], fillSend[idx], avalSend[idx]);
        }
        ++idx;
      }
      fprintf(logFile, "\n");
    }
  }

  /* start sends to higher-ordred nabors */
  for (i=0; i<hiCount; ++i) {
    int nabor = hiNabors[i];
    MPI_Isend(cvalSend, nz, MPI_INT,    nabor, CVAL_TAG, comm_dh, er->cval_req+i);
    MPI_Isend(fillSend, nz, MPI_INT,    nabor, FILL_TAG, comm_dh, er->fill_req+i); 
    MPI_Isend(avalSend, nz, MPI_DOUBLE, nabor, AVAL_TAG, comm_dh, er->aval_req+i);
  }
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "waitfor_sends_private"
void waitfor_sends_private(ExternalRows_dh er)
{
  START_FUNC_DH
  MPI_Status *status = er->status;
  int hiCount = er->sg->hiCount;

  if (hiCount) {
    MPI_Waitall(hiCount, er->req1, status);
    MPI_Waitall(hiCount, er->req2, status);
    MPI_Waitall(hiCount, er->req3, status);
    MPI_Waitall(hiCount, er->req4, status);
    MPI_Waitall(hiCount, er->cval_req, status);
    MPI_Waitall(hiCount, er->fill_req, status);
    MPI_Waitall(hiCount, er->aval_req, status);
  }
  END_FUNC_DH
}
