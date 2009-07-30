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

/* to do: re-integrate fix-smalll-pivots */

#include "ilu_dh.h"
#include "Mem_dh.h"
#include "Parser_dh.h"
#include "Euclid_dh.h"
#include "getRow_dh.h"
#include "Factor_dh.h"
#include "SubdomainGraph_dh.h"

static bool check_constraint_private(Euclid_dh ctx, int b, int j);

static int symbolic_row_private(int localRow, 
                 int *list, int *marker, int *tmpFill,
                 int len, int *CVAL, double *AVAL,
                 int *o2n_col, Euclid_dh ctx, bool debug);

static int numeric_row_private(int localRow, 
                        int len, int *CVAL, double *AVAL,
                        REAL_DH *work, int *o2n_col, Euclid_dh ctx, bool debug);


#undef __FUNC__
#define __FUNC__ "compute_scaling_private"
void compute_scaling_private(int row, int len, double *AVAL, Euclid_dh ctx)
{
  START_FUNC_DH
  double tmp = 0.0;
  int j;

  for (j=0; j<len; ++j) tmp = MAX( tmp, fabs(AVAL[j]) );
  if (tmp) {
    ctx->scale[row] = 1.0/tmp;
  }
  END_FUNC_DH
}

#if 0

/* not used ? */
#undef __FUNC__
#define __FUNC__ "fixPivot_private"
double fixPivot_private(int row, int len, float *vals)
{
  START_FUNC_DH
  int i;
  float max = 0.0;
  bool debug = false;

  for (i=0; i<len; ++i) {
    float tmp = fabs(vals[i]);
    max = MAX(max, tmp);
  }
  END_FUNC_VAL(max* ctxPrivate->pivotFix)
}

#endif




#undef __FUNC__
#define __FUNC__ "iluk_seq"
void iluk_seq(Euclid_dh ctx)
{
  START_FUNC_DH
  int      *rp, *cval, *diag;
  int      *CVAL;
  int      i, j, len, count, col, idx = 0;
  int      *list, *marker, *fill, *tmpFill;
  int      temp, m, from = ctx->from, to = ctx->to;
  int      *n2o_row, *o2n_col, beg_row, beg_rowP;
  double   *AVAL;
  REAL_DH  *work, *aval;
  Factor_dh F = ctx->F;
  SubdomainGraph_dh sg = ctx->sg;
  bool debug = false;

  if (logFile != NULL  &&  Parser_dhHasSwitch(parser_dh, "-debug_ilu")) debug = true;

  m = F->m;
  rp = F->rp;
  cval = F->cval;
  fill = F->fill;
  diag = F->diag;
  aval = F->aval;
  work = ctx->work;
  count = rp[from];

  if (sg == NULL) {
    SET_V_ERROR("subdomain graph is NULL");
  }

  n2o_row = ctx->sg->n2o_row;
  o2n_col = ctx->sg->o2n_col;
  beg_row  = ctx->sg->beg_row[myid_dh];
  beg_rowP  = ctx->sg->beg_rowP[myid_dh];

  /* allocate and initialize working space */
  list   = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  marker = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  tmpFill = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) marker[i] = -1;

  /* working space for values */
  for (i=0; i<m; ++i) work[i] = 0.0; 

/*    printf_dh("====================== starting iluk_seq; level= %i\n\n", ctx->level); 
*/


  /*---------- main loop ----------*/

  for (i=from; i<to; ++i) {
    int row = n2o_row[i];             /* local row number */
    int globalRow = row+beg_row;      /* global row number */

/*fprintf(logFile, "--------------------------------- localRow= %i\n", 1+i);
*/

    if (debug) {
	fprintf(logFile, "ILU_seq ================================= starting local row: %i, (global= %i) level= %i\n", i+1, i+1+sg->beg_rowP[myid_dh], ctx->level);
    }

    EuclidGetRow(ctx->A, globalRow, &len, &CVAL, &AVAL); CHECK_V_ERROR;

    /* compute scaling value for row(i) */
    if (ctx->isScaled) { 
      compute_scaling_private(i, len, AVAL, ctx); CHECK_V_ERROR; 
    }

    /* Compute symbolic factor for row(i);
       this also performs sparsification
     */
    count = symbolic_row_private(i, list, marker, tmpFill, 
                                 len, CVAL, AVAL,
                                 o2n_col, ctx, debug); CHECK_V_ERROR;

    /* Ensure adequate storage; reallocate, if necessary. */
    if (idx + count > F->alloc) {
      Factor_dhReallocate(F, idx, count); CHECK_V_ERROR;
      SET_INFO("REALLOCATED from ilu_seq");
      cval = F->cval;
      fill = F->fill;
      aval = F->aval;
    }

    /* Copy factored symbolic row to permanent storage */
    col = list[m];
    while (count--) {
      cval[idx] = col;  
      fill[idx] = tmpFill[col];
      ++idx;
/*fprintf(logFile, "  col= %i\n", 1+col);
*/
      col = list[col];
    }

    /* add row-pointer to start of next row. */
    rp[i+1] = idx;

    /* Insert pointer to diagonal */
    temp = rp[i]; 
    while (cval[temp] != i) ++temp; 
    diag[i] = temp;

/*fprintf(logFile, "  diag[i]= %i\n", diag);
*/

    /* compute numeric factor for current row */
     numeric_row_private(i, len, CVAL, AVAL, 
                          work, o2n_col, ctx, debug); CHECK_V_ERROR
    EuclidRestoreRow(ctx->A, globalRow, &len, &CVAL, &AVAL); CHECK_V_ERROR;

    /* Copy factored numeric row to permanent storage,
       and re-zero work vector
     */
    if (debug) {
      fprintf(logFile, "ILU_seq:  ");
      for (j=rp[i]; j<rp[i+1]; ++j) {
        col = cval[j];  
        aval[j] = work[col];  
        work[col] = 0.0;
        fprintf(logFile, "%i,%i,%g ; ", 1+cval[j], fill[j], aval[j]);
        fflush(logFile);
      }
      fprintf(logFile, "\n");
    } else {
      for (j=rp[i]; j<rp[i+1]; ++j) {
        col = cval[j];  
        aval[j] = work[col];  
        work[col] = 0.0;
      }
    }

    /* check for zero diagonal */
    if (! aval[diag[i]]) {
      sprintf(msgBuf_dh, "zero diagonal in local row %i", i+1);
      SET_V_ERROR(msgBuf_dh);
    }
  }

  FREE_DH(list); CHECK_V_ERROR;
  FREE_DH(tmpFill); CHECK_V_ERROR;
  FREE_DH(marker); CHECK_V_ERROR;

  /* adjust column indices back to global */
  if (beg_rowP) {
    int start = rp[from];
    int stop = rp[to];
    for (i=start; i<stop; ++i) cval[i] += beg_rowP;
  }

  /* for debugging: this is so the Print methods will work, even if
     F hasn't been fully factored
  */
  for (i=to+1; i<m; ++i) rp[i] = 0;

  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "iluk_seq_block"
void iluk_seq_block(Euclid_dh ctx)
{
  START_FUNC_DH
  int      *rp, *cval, *diag;
  int      *CVAL;
  int      h, i, j, len, count, col, idx = 0;
  int      *list, *marker, *fill, *tmpFill;
  int      temp, m;
  int      *n2o_row, *o2n_col, *beg_rowP, *n2o_sub, blocks;
  int      *row_count, *dummy = NULL, dummy2[1];
  double   *AVAL;
  REAL_DH  *work, *aval;
  Factor_dh F = ctx->F;
  SubdomainGraph_dh sg = ctx->sg;
  bool bj = false, constrained = false;
  int discard = 0;
  int gr = -1;  /* globalRow */
  bool debug = false;

  if (logFile != NULL  &&  Parser_dhHasSwitch(parser_dh, "-debug_ilu")) debug = true;

/*fprintf(stderr, "====================== starting iluk_seq_block; level= %i\n\n", ctx->level); 
*/

  if (!strcmp(ctx->algo_par, "bj")) bj = true;
  constrained = ! Parser_dhHasSwitch(parser_dh, "-unconstrained");

  m    = F->m;
  rp   = F->rp;
  cval = F->cval;
  fill = F->fill;
  diag = F->diag;
  aval = F->aval;
  work = ctx->work;

  if (sg != NULL) {
    n2o_row   = sg->n2o_row;
    o2n_col   = sg->o2n_col;
    row_count = sg->row_count;
    /* beg_row   = sg->beg_row ; */
    beg_rowP  = sg->beg_rowP;
    n2o_sub   = sg->n2o_sub;
    blocks    = sg->blocks;
  }

  else {
    dummy = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
    for (i=0; i<m; ++i) dummy[i] = i;
    n2o_row   = dummy;
    o2n_col   = dummy;
    dummy2[0] = m; row_count = dummy2;
    /* beg_row   = 0; */
    beg_rowP  = dummy;
    n2o_sub   = dummy;
    blocks    = 1;
  }

  /* allocate and initialize working space */
  list   = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  marker = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  tmpFill = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) marker[i] = -1;

  /* working space for values */
  for (i=0; i<m; ++i) work[i] = 0.0; 

  /*---------- main loop ----------*/

 for (h=0; h<blocks; ++h) {
  /* 1st and last row in current block, with respect to A */
  int curBlock = n2o_sub[h];
  int first_row = beg_rowP[curBlock];
  int end_row   = first_row + row_count[curBlock];

    if (debug) {
        fprintf(logFile, "\n\nILU_seq BLOCK: %i @@@@@@@@@@@@@@@ \n", curBlock);
    }

  for (i=first_row; i<end_row; ++i) {
    int row = n2o_row[i];  
    ++gr;

    if (debug) {
      fprintf(logFile, "ILU_seq  global: %i  local: %i =================================\n", 1+gr, 1+i-first_row);
    }

/*prinft("first_row= %i  end_row= %i\n", first_row, end_row);
*/

    EuclidGetRow(ctx->A, row, &len, &CVAL, &AVAL); CHECK_V_ERROR;

    /* compute scaling value for row(i) */
    if (ctx->isScaled) { 
      compute_scaling_private(i, len, AVAL, ctx); CHECK_V_ERROR; 
    }

    /* Compute symbolic factor for row(i);
       this also performs sparsification
     */
    count = symbolic_row_private(i, list, marker, tmpFill, 
                                 len, CVAL, AVAL,
                                 o2n_col, ctx, debug); CHECK_V_ERROR;

    /* Ensure adequate storage; reallocate, if necessary. */
    if (idx + count > F->alloc) {
      Factor_dhReallocate(F, idx, count); CHECK_V_ERROR;
      SET_INFO("REALLOCATED from ilu_seq");
      cval = F->cval;
      fill = F->fill;
      aval = F->aval;
    }

    /* Copy factored symbolic row to permanent storage */
    col = list[m];
    while (count--) {

      /* constrained pilu */
      if (constrained && !bj) {
        if (col >= first_row && col < end_row) {
          cval[idx] = col;  
          fill[idx] = tmpFill[col];
          ++idx;
        } else { 
          if (check_constraint_private(ctx, curBlock, col)) {
            cval[idx] = col;  
            fill[idx] = tmpFill[col];
            ++idx;
          } else {
            ++discard; 
          }
        }
        col = list[col];
      }

      /* block jacobi case */
      else if (bj) {
        if (col >= first_row && col < end_row) {
          cval[idx] = col;  
          fill[idx] = tmpFill[col];
          ++idx;
        } else { 
          ++discard; 
        }
        col = list[col];
      }

      /* general case */
      else {
        cval[idx] = col;  
        fill[idx] = tmpFill[col];
        ++idx;
        col = list[col];
      }
    }

    /* add row-pointer to start of next row. */
    rp[i+1] = idx;

    /* Insert pointer to diagonal */
    temp = rp[i]; 
    while (cval[temp] != i) ++temp; 
    diag[i] = temp;

    /* compute numeric factor for current row */
    numeric_row_private(i, len, CVAL, AVAL, 
                          work, o2n_col, ctx, debug); CHECK_V_ERROR
    EuclidRestoreRow(ctx->A, row, &len, &CVAL, &AVAL); CHECK_V_ERROR;

    /* Copy factored numeric row to permanent storage,
       and re-zero work vector
     */
    if (debug) {
      fprintf(logFile, "ILU_seq: ");
      for (j=rp[i]; j<rp[i+1]; ++j) {
        col = cval[j];  
        aval[j] = work[col];  
        work[col] = 0.0;
        fprintf(logFile, "%i,%i,%g ; ", 1+cval[j], fill[j], aval[j]);
      }
      fprintf(logFile, "\n");
     } 

     /* normal operation */
     else {
      for (j=rp[i]; j<rp[i+1]; ++j) {
        col = cval[j];  
        aval[j] = work[col];  
        work[col] = 0.0;
      } 
    }

    /* check for zero diagonal */
    if (! aval[diag[i]]) {
      sprintf(msgBuf_dh, "zero diagonal in local row %i", i+1);
      SET_V_ERROR(msgBuf_dh);
    }
  }
 }

/*  printf("bj= %i  constrained= %i  discarded= %i\n", bj, constrained, discard); */

  if (dummy != NULL) { FREE_DH(dummy); CHECK_V_ERROR; }
  FREE_DH(list); CHECK_V_ERROR;
  FREE_DH(tmpFill); CHECK_V_ERROR;
  FREE_DH(marker); CHECK_V_ERROR;

  END_FUNC_DH
}



/* Computes ILU(K) factor of a single row; returns fill 
   count for the row.  Explicitly inserts diag if not already 
   present.  On return, all column indices are local 
   (i.e, referenced to 0).
*/
#undef __FUNC__
#define __FUNC__ "symbolic_row_private"
int symbolic_row_private(int localRow, 
                 int *list, int *marker, int *tmpFill,
                 int len, int *CVAL, double *AVAL,
                 int *o2n_col, Euclid_dh ctx, bool debug)
{
  START_FUNC_DH
  int level = ctx->level, m = ctx->F->m;
  int *cval = ctx->F->cval, *diag = ctx->F->diag, *rp = ctx->F->rp; 
  int *fill = ctx->F->fill;
  int count = 0;
  int j, node, tmp, col, head;
  int fill1, fill2, beg_row;
  double val;
  double thresh = ctx->sparseTolA;
  REAL_DH scale;

  scale = ctx->scale[localRow]; 
  ctx->stats[NZA_STATS] += (double)len;
  beg_row  = ctx->sg->beg_row[myid_dh];

  /* Insert col indices in linked list, and values in work vector.
   * List[m] points to the first (smallest) col in the linked list.
   * Column values are adjusted from global to local numbering.
   */
  list[m] = m;
  for (j=0; j<len; ++j) {
    tmp = m;
    col = *CVAL++;
    col -= beg_row;     /* adjust to zero based */
    col = o2n_col[col]; /* permute the column */
    val = *AVAL++;
    val *= scale;       /* scale the value */

    if (fabs(val) > thresh || col == localRow) {  /* sparsification */
      ++count;
      while (col > list[tmp]) tmp = list[tmp];
      list[col]   = list[tmp];
      list[tmp]   = col;
      tmpFill[col] = 0;
      marker[col] = localRow;
    }
  }

  /* insert diag if not already present */
  if (marker[localRow] != localRow) {
    tmp = m;
    while (localRow > list[tmp]) tmp = list[tmp];
    list[localRow]    = list[tmp];
    list[tmp]    = localRow;
    tmpFill[localRow] = 0;
    marker[localRow]  = localRow;
    ++count;
  }
  ctx->stats[NZA_USED_STATS] += (double)count;

  /* update row from previously factored rows */
  head = m;
  if (level > 0) {
    while (list[head] < localRow) {
      node = list[head];
      fill1 = tmpFill[node];

      if (debug) {
        fprintf(logFile, "ILU_seq   sf updating from row: %i\n", 1+node);
      }

      if (fill1 < level) {
        for (j = diag[node]+1; j<rp[node+1]; ++j) {
          col = cval[j];
          fill2 = fill1 + fill[j] + 1;

          if (fill2 <= level) {
            /* if newly discovered fill entry, mark it as discovered;
             * if entry has level <= K, add it to the linked-list.
             */
            if (marker[col] < localRow) {
              tmp = head;
              marker[col] = localRow;
              tmpFill[col] = fill2;
              while (col > list[tmp]) tmp = list[tmp];
              list[col] = list[tmp];
              list[tmp]    = col;
              ++count; /* increment fill count */
            }

          /* if previously-discovered fill, update the entry's level. */
            else {
              tmpFill[col] = (fill2 < tmpFill[col]) ? fill2 : tmpFill[col];
            }
          }
        }
      } /* fill1 < level  */
      head = list[head];  /* advance to next item in linked list */
    }
  }
  END_FUNC_VAL(count)
}


#undef __FUNC__
#define __FUNC__ "numeric_row_private"
int numeric_row_private(int localRow, 
                        int len, int *CVAL, double *AVAL,
                        REAL_DH *work, int *o2n_col, Euclid_dh ctx, bool debug)
{
  START_FUNC_DH
  double  pc, pv, multiplier;
  int     j, k, col, row;
  int     *rp = ctx->F->rp, *cval = ctx->F->cval;
  int     *diag = ctx->F->diag;
  int     beg_row;
  double  val;
  REAL_DH *aval = ctx->F->aval, scale;

  scale = ctx->scale[localRow]; 
  beg_row  = ctx->sg->beg_row[myid_dh];

  /* zero work vector */
  /* note: indices in col[] are already permuted. */
  for (j=rp[localRow]; j<rp[localRow+1]; ++j) { 
    col = cval[j];  
    work[col] = 0.0; 
  }

  /* init work vector with values from A */
  /* (note: some values may be na due to sparsification; this is O.K.) */
  for (j=0; j<len; ++j) {
    col = *CVAL++;
    col -= beg_row;
    val = *AVAL++;
    col = o2n_col[col];  /* note: we permute the indices from A */
    work[col] = val*scale;
  }



/*fprintf(stderr, "local row= %i\n", 1+localRow);
*/


  for (j=rp[localRow]; j<diag[localRow]; ++j) {
    row = cval[j];     /* previously factored row */
    pc = work[row];


      pv = aval[diag[row]]; /* diagonal of previously factored row */

/*
if (pc == 0.0 || pv == 0.0) {
fprintf(stderr, "pv= %g; pc= %g\n", pv, pc);
}
*/

    if (pc != 0.0 && pv != 0.0) {
      multiplier = pc / pv;
      work[row] = multiplier;

      if (debug) {
        fprintf(logFile, "ILU_seq   nf updating from row: %i; multiplier= %g\n", 1+row, multiplier);
      }

      for (k=diag[row]+1; k<rp[row+1]; ++k) {
        col = cval[k];
        work[col] -= (multiplier * aval[k]);
      }
    } else  {
      if (debug) {
        fprintf(logFile, "ILU_seq   nf NO UPDATE from row %i; pc = %g; pv = %g\n", 1+row, pc, pv);
      }
    }
  }

  /* check for zero or too small of a pivot */
#if 0
  if (fabs(work[i]) <= pivotTol) {
    /* yuck! assume row scaling, and just stick in a value */
    aval[diag[i]] = pivotFix;
  }
#endif

  END_FUNC_VAL(0)
}


/*-----------------------------------------------------------------------*
 * ILUT starts here
 *-----------------------------------------------------------------------*/
int ilut_row_private(int localRow, int *list, int *o2n_col, int *marker,
                     int len, int *CVAL, double *AVAL,
                     REAL_DH *work, Euclid_dh ctx, bool debug);

#undef __FUNC__
#define __FUNC__ "ilut_seq"
void ilut_seq(Euclid_dh ctx)
{
  START_FUNC_DH
  int      *rp, *cval, *diag, *CVAL;
  int      i, len, count, col, idx = 0;
  int      *list, *marker;
  int      temp, m, from, to;
  int      *n2o_row, *o2n_col, beg_row, beg_rowP;
  double   *AVAL, droptol; 
  REAL_DH *work, *aval, val;
  Factor_dh F = ctx->F;
  SubdomainGraph_dh sg = ctx->sg;
  bool debug = false;

  if (logFile != NULL  &&  Parser_dhHasSwitch(parser_dh, "-debug_ilu")) debug = true;

  m = F->m;
  rp = F->rp;
  cval = F->cval;
  diag = F->diag;
  aval = F->aval;
  work = ctx->work;
  from = ctx->from;
  to = ctx->to;
  count = rp[from];
  droptol = ctx->droptol;

  if (sg == NULL) {
    SET_V_ERROR("subdomain graph is NULL");
  }

  n2o_row = ctx->sg->n2o_row;
  o2n_col = ctx->sg->o2n_col;
  beg_row  = ctx->sg->beg_row[myid_dh];
  beg_rowP  = ctx->sg->beg_rowP[myid_dh];


  /* allocate and initialize working space */
  list   = (int*)MALLOC_DH((m+1)*sizeof(int)); CHECK_V_ERROR;
  marker = (int*)MALLOC_DH(m*sizeof(int)); CHECK_V_ERROR;
  for (i=0; i<m; ++i) marker[i] = -1;
  rp[0] = 0;

  /* working space for values */
  for (i=0; i<m; ++i) work[i] = 0.0; 

  /* ----- main loop start ----- */
  for (i=from; i<to; ++i) {
    int row = n2o_row[i];             /* local row number */
    int globalRow = row + beg_row;    /* global row number */
    EuclidGetRow(ctx->A, globalRow, &len, &CVAL, &AVAL); CHECK_V_ERROR;

    /* compute scaling value for row(i) */
    compute_scaling_private(i, len, AVAL, ctx); CHECK_V_ERROR; 

    /* compute factor for row i */
    count = ilut_row_private(i, list, o2n_col, marker,
                         len, CVAL, AVAL, work, ctx, debug); CHECK_V_ERROR;

    EuclidRestoreRow(ctx->A, globalRow, &len, &CVAL, &AVAL); CHECK_V_ERROR;

    /* Ensure adequate storage; reallocate, if necessary. */
    if (idx + count > F->alloc) {
      Factor_dhReallocate(F, idx, count); CHECK_V_ERROR;
      SET_INFO("REALLOCATED from ilu_seq");
      cval = F->cval;
      aval = F->aval;
    }

    /* Copy factored row to permanent storage,
       apply 2nd drop test,
       and re-zero work vector
     */
    col = list[m];
    while (count--) {
      val = work[col];
      if (col == i || fabs(val) > droptol) {
        cval[idx] = col;   
        aval[idx++] = val;
        work[col] = 0.0;
      }
      col = list[col];
    }

    /* add row-pointer to start of next row. */
    rp[i+1] = idx;

    /* Insert pointer to diagonal */
    temp = rp[i]; 
    while (cval[temp] != i) ++temp;
    diag[i] = temp;

    /* check for zero diagonal */
    if (! aval[diag[i]]) {
      sprintf(msgBuf_dh, "zero diagonal in local row %i", i+1);
      SET_V_ERROR(msgBuf_dh);
    }
  } /* --------- main loop end --------- */

  /* adjust column indices back to global */
  if (beg_rowP) {
    int start = rp[from];
    int stop = rp[to];
    for (i=start; i<stop; ++i) cval[i] += beg_rowP;
  }

  FREE_DH(list);
  FREE_DH(marker);
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "ilut_row_private"
int ilut_row_private(int localRow, int *list, int *o2n_col, int *marker,
                     int len, int *CVAL, double *AVAL,
                     REAL_DH *work, Euclid_dh ctx, bool debug)
{
  START_FUNC_DH
  Factor_dh F = ctx->F;
  int     j, col, m = ctx->m, *rp = F->rp, *cval = F->cval;
  int     tmp, *diag = F->diag;
  int     head;
  int     count = 0, beg_row;
  double  val;
  double  mult, *aval = F->aval;
  double  scale, pv, pc;
  double  droptol = ctx->droptol;
  double thresh = ctx->sparseTolA;

  scale = ctx->scale[localRow];
  ctx->stats[NZA_STATS] += (double)len;
  beg_row  = ctx->sg->beg_row[myid_dh];


  /* Insert col indices in linked list, and values in work vector.
   * List[m] points to the first (smallest) col in the linked list.
   * Column values are adjusted from global to local numbering.
   */
  list[m] = m;
  for (j=0; j<len; ++j) {
    tmp = m;
    col = *CVAL++;
    col -= beg_row;     /* adjust to zero based */
    col = o2n_col[col]; /* permute the column */
    val = *AVAL++;
    val *= scale;       /* scale the value */

    if (fabs(val) > thresh || col == localRow) {  /* sparsification */
      ++count;
      while (col > list[tmp]) tmp = list[tmp];
      list[col]   = list[tmp];
      list[tmp]   = col;
      work[col] = val;
      marker[col] = localRow;
    }
  }

  /* insert diag if not already present */
  if (marker[localRow] != localRow) {
    tmp = m;
    while (localRow > list[tmp]) tmp = list[tmp];
    list[localRow]    = list[tmp];
    list[tmp]    = localRow;
    marker[localRow]  = localRow;
    ++count;
  }

  /* update current row from previously factored rows */
  head = m;
  while (list[head] < localRow) {
    int row = list[head];

    /* get the multiplier, and apply 1st drop tolerance test */
    pc = work[row];
    if (pc != 0.0) {
      pv = aval[diag[row]];  /* diagonal (pivot) of previously factored row */
      mult = pc / pv;

      /* update localRow from previously factored "row" */
      if (fabs(mult) > droptol) {
        work[row] = mult;

        for (j=diag[row]+1; j<rp[row+1]; ++j) {
          col = cval[j];
          work[col] -= (mult * aval[j]);

          /* if col isn't already present in the linked-list, insert it.  */
          if (marker[col] < localRow) {
            marker[col] = localRow;     /* mark the column as known fill */
            tmp = head;             /* insert in list [this and next 3 lines] */
            while (col > list[tmp]) tmp = list[tmp];
            list[col] = list[tmp];
            list[tmp] = col;
            ++count;                /* increment fill count */
          }
        } 
      }
    }
    head = list[head];  /* advance to next item in linked list */
  }

  END_FUNC_VAL(count)
}


#undef __FUNC__
#define __FUNC__ "check_constraint_private"
bool check_constraint_private(Euclid_dh ctx, int p1, int j) 
{
  START_FUNC_DH
  bool retval = false;
  int i, p2;
  int *nabors, count;
  SubdomainGraph_dh sg = ctx->sg;

  if (sg == NULL) {
    SET_ERROR(-1, "ctx->sg == NULL");
  }
  
  p2 = SubdomainGraph_dhFindOwner(ctx->sg, j, true);


  nabors = sg->adj + sg->ptrs[p1];
  count = sg->ptrs[p1+1]  - sg->ptrs[p1];

/*
printf("p1= %i, p2= %i;  p1's nabors: ", p1, p2);
for (i=0; i<count; ++i) printf("%i ", nabors[i]);
printf("\n");
*/

  for (i=0; i<count; ++i) {
/* printf("  @@@ next nabor= %i\n", nabors[i]);
*/
    if (nabors[i] == p2) {
      retval = true;
      break;
    }
  }

  END_FUNC_VAL(retval)
}
