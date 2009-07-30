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

#include "Euclid_dh.h"
#include "Mem_dh.h"
#include "Mat_dh.h"
#include "Vec_dh.h"
#include "Factor_dh.h"
#include "getRow_dh.h"
#include "ilu_dh.h"
#include "Parser_dh.h"
#include "SortedList_dh.h"
#include "SubdomainGraph_dh.h"
#include "ExternalRows_dh.h"
#include "krylov_dh.h"

static void get_runtime_params_private(Euclid_dh ctx);
static void invert_diagonals_private(Euclid_dh ctx);
static void compute_rho_private(Euclid_dh ctx); 
static void factor_private(Euclid_dh ctx); 
/* static void discard_indices_private(Euclid_dh ctx); */
static void reduce_timings_private(Euclid_dh ctx);

#undef __FUNC__
#define __FUNC__ "Euclid_dhCreate"
void Euclid_dhCreate(Euclid_dh *ctxOUT)
{
  START_FUNC_DH
  struct _mpi_interface_dh * ctx =
     (struct _mpi_interface_dh*)MALLOC_DH(sizeof(struct _mpi_interface_dh)); CHECK_V_ERROR;
  *ctxOUT = ctx;

  ctx->isSetup = false;

  ctx->rho_init = 2.0;
  ctx->rho_final = 0.0;

  ctx->m = 0;
  ctx->n = 0;
  ctx->rhs = NULL;
  ctx->A = NULL;
  ctx->F = NULL;
  ctx->sg = NULL;

  ctx->scale = NULL;
  ctx->isScaled = false;
  ctx->work = NULL;
  ctx->work2 = NULL;
  ctx->from = 0;
  ctx->to = 0;

  strcpy(ctx->algo_par, "pilu");
  strcpy(ctx->algo_ilu, "iluk");
  ctx->level = 1;
  ctx->droptol = DEFAULT_DROP_TOL;
  ctx->sparseTolA = 0.0;
  ctx->sparseTolF = 0.0;
  ctx->pivotMin = 0.0;
  ctx->pivotFix = PIVOT_FIX_DEFAULT;
  ctx->maxVal = 0.0;

  ctx->slist = NULL;
  ctx->extRows = NULL;

  strcpy(ctx->krylovMethod, "bicgstab");
  ctx->maxIts = 200;
  ctx->rtol = 1e-5;
  ctx->atol = 1e-50;
  ctx->its = 0;
  ctx->itsTotal = 0;
  ctx->setupCount = 0;
  ctx->logging = 0;
  ctx->printStats = (Parser_dhHasSwitch(parser_dh, "-printStats"));

  { int i;
    for (i=0; i<TIMING_BINS; ++i) ctx->timing[i] = 0.0;
    for (i=0; i<STATS_BINS; ++i) ctx->stats[i] = 0.0;
  }
  ctx->timingsWereReduced = false;

  ++ref_counter;
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "Euclid_dhDestroy"
void Euclid_dhDestroy(Euclid_dh ctx)
{
  START_FUNC_DH

  if (Parser_dhHasSwitch(parser_dh, "-eu_stats")
      || ctx->logging) {
    /* insert switch so memory report will also be printed */
    Parser_dhInsert(parser_dh, "-eu_mem", "1"); CHECK_V_ERROR;
    Euclid_dhPrintHypreReport(ctx, stdout); CHECK_V_ERROR;
  }

  if (ctx->setupCount > 1 && ctx->printStats) {
    Euclid_dhPrintStatsShorter(ctx, stdout); CHECK_V_ERROR;
  }

  if (ctx->F != NULL) { Factor_dhDestroy(ctx->F); CHECK_V_ERROR; }
  if (ctx->sg != NULL) { SubdomainGraph_dhDestroy(ctx->sg); CHECK_V_ERROR; }
  if (ctx->scale != NULL) { FREE_DH(ctx->scale); CHECK_V_ERROR; }
  if (ctx->work != NULL) { FREE_DH(ctx->work); CHECK_V_ERROR; }
  if (ctx->work2 != NULL) { FREE_DH(ctx->work2); CHECK_V_ERROR; }
  if (ctx->slist != NULL) { SortedList_dhDestroy(ctx->slist); CHECK_V_ERROR; }
  if (ctx->extRows != NULL) { ExternalRows_dhDestroy(ctx->extRows); CHECK_V_ERROR; }
  FREE_DH(ctx); CHECK_V_ERROR; 

  --ref_counter;
  END_FUNC_DH
}


/* on entry, "A" must have been set.  If this context is being
   reused, user must ensure ???
*/
#undef __FUNC__
#define __FUNC__ "Euclid_dhSetup"
void Euclid_dhSetup(Euclid_dh ctx)
{
  START_FUNC_DH
  int m, n, beg_row;
  double t1;
  bool isSetup = ctx->isSetup;
  bool bj = false;

  /*----------------------------------------------------
   * If Euclid was previously setup, print summary of
   * what happened during previous setup/solve
   *----------------------------------------------------*/
  if (ctx->setupCount && ctx->printStats) {
    Euclid_dhPrintStatsShorter(ctx, stdout); CHECK_V_ERROR;
    ctx->its = 0;
  }

  /*----------------------------------------------------
   * zero array for statistical reporting
   *----------------------------------------------------*/
  { int i;
    for (i=0; i<STATS_BINS; ++i) ctx->stats[i] = 0.0;
  }

  /*----------------------------------------------------
   * internal timing
   *----------------------------------------------------*/
  ctx->timing[SOLVE_START_T] = MPI_Wtime();
  /* sum timing from last linear solve cycle, if any */
  ctx->timing[TOTAL_SOLVE_T] += ctx->timing[TOTAL_SOLVE_TEMP_T];
  ctx->timing[TOTAL_SOLVE_TEMP_T] = 0.0;

  if (ctx->F != NULL) {
    Factor_dhDestroy(ctx->F); CHECK_V_ERROR;
    ctx->F = NULL;
  }

  if (ctx->A == NULL) {
    SET_V_ERROR("must set ctx->A before calling init");
  }
  EuclidGetDimensions(ctx->A, &beg_row, &m, &n); CHECK_V_ERROR;
  ctx->m = m;
  ctx->n = n;

  if (Parser_dhHasSwitch(parser_dh, "-print_size")) {
    printf_dh("setting up linear system; global rows: %i  local rows: %i (on P_0)\n", n,m);
  }

  sprintf(msgBuf_dh, "localRow= %i;  globalRows= %i;  beg_row= %i", m, n, beg_row);
  SET_INFO(msgBuf_dh);

  bj = Parser_dhHasSwitch(parser_dh, "-bj");

  /*------------------------------------------------------------------------
   * Setup the SubdomainGraph, which contains connectivity and
   * and permutation information.  If this context is being reused,
   * this may already have been done; if being resused, the underlying
   * subdomain graph cannot change (user's responsibility?)
   *------------------------------------------------------------------------*/
  if (ctx->sg == NULL) {
    int blocks = np_dh;
    t1 = MPI_Wtime();
    if (np_dh == 1) {
      Parser_dhReadInt(parser_dh, "-blocks", &blocks); CHECK_V_ERROR;
      SubdomainGraph_dhCreate(&(ctx->sg)); CHECK_V_ERROR;
      SubdomainGraph_dhInit(ctx->sg, blocks, bj, ctx->A); CHECK_V_ERROR;
    } else {
      SubdomainGraph_dhCreate(&(ctx->sg)); CHECK_V_ERROR;
      SubdomainGraph_dhInit(ctx->sg, -1, bj, ctx->A); CHECK_V_ERROR;
    }
    ctx->timing[SUB_GRAPH_T] += (MPI_Wtime() - t1);
  }


/* SubdomainGraph_dhDump(ctx->sg, "SG.dump"); CHECK_V_ERROR; */

  /*----------------------------------------------------
   * for debugging
   *----------------------------------------------------*/
  if (Parser_dhHasSwitch(parser_dh, "-doNotFactor")) {
    goto END_OF_FUNCTION;
  }


  /*----------------------------------------------------
   * query parser for runtime parameters
   *----------------------------------------------------*/
  if (! isSetup) {
    get_runtime_params_private(ctx); CHECK_V_ERROR;
  }
  if (! strcmp(ctx->algo_par, "bj")) bj = false;

  /*--------------------------------------------------------- 
   * allocate and initialize storage for row-scaling
   * (ctx->isScaled is set in get_runtime_params_private(); )
   *---------------------------------------------------------*/
  if (ctx->scale == NULL) {
    ctx->scale = (REAL_DH*)MALLOC_DH(m*sizeof(REAL_DH)); CHECK_V_ERROR;
  }
  { int i; for (i=0; i<m; ++i) ctx->scale[i] = 1.0; }

  /*------------------------------------------------------------------ 
   * allocate work vectors; used in factorization and triangular solves;
   *------------------------------------------------------------------*/
  if ( ctx->work == NULL) { 
    ctx->work = (REAL_DH*)MALLOC_DH(m*sizeof(REAL_DH)); CHECK_V_ERROR;
  }
  if ( ctx->work2 == NULL) { 
    ctx->work2 = (REAL_DH*)MALLOC_DH(m*sizeof(REAL_DH)); CHECK_V_ERROR;
  }

  /*-----------------------------------------------------------------
   * perform the incomplete factorization (this should be, at least
   * for higher level ILUK, the most time-intensive portion of setup)
   *-----------------------------------------------------------------*/
  t1 = MPI_Wtime();    
  factor_private(ctx); CHECK_V_ERROR;
  ctx->timing[FACTOR_T] += (MPI_Wtime() - t1);

  /*-------------------------------------------------------------- 
   * invert diagonals, for faster triangular solves
   *--------------------------------------------------------------*/
  if (strcmp(ctx->algo_par, "none")) {
    invert_diagonals_private(ctx); CHECK_V_ERROR; 
  }

  /*-------------------------------------------------------------- 
   * compute rho_final: global ratio of nzF/nzA
   * also, if -sparseA > 0,  compute ratio of nzA 
   * used in factorization
   *--------------------------------------------------------------*/
  /* for some reason compute_rho_private() was expensive, so now it's
     an option, unless there's only one mpi task.
   */
  if (Parser_dhHasSwitch(parser_dh, "-computeRho") || np_dh == 1) {
   if (strcmp(ctx->algo_par, "none")) { 
     t1 = MPI_Wtime();
     compute_rho_private(ctx); CHECK_V_ERROR;  
     ctx->timing[COMPUTE_RHO_T] += (MPI_Wtime() - t1);
   }
 }

  /*-------------------------------------------------------------- 
   * if using PILU, set up persistent comms and global-to-local
   * number scheme, for efficient triangular solves.
   * (Thanks to Edmond Chow for these algorithmic ideas.)
   *--------------------------------------------------------------*/

  if (! strcmp(ctx->algo_par, "pilu")  &&  np_dh > 1) { 
    t1 = MPI_Wtime();
    Factor_dhSolveSetup(ctx->F, ctx->sg); CHECK_V_ERROR;
    ctx->timing[SOLVE_SETUP_T] += (MPI_Wtime() - t1);
  }

END_OF_FUNCTION: ;

  /*-------------------------------------------------------
   * internal timing
   *-------------------------------------------------------*/
  ctx->timing[SETUP_T] += (MPI_Wtime() - ctx->timing[SOLVE_START_T]);
  ctx->setupCount += 1;  

  ctx->isSetup = true;

  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "get_runtime_params_private"
void get_runtime_params_private(Euclid_dh ctx)
{
  START_FUNC_DH
  char *tmp;

  /* params for use of internal solvers */
  Parser_dhReadInt(parser_dh,    "-maxIts",&(ctx->maxIts)); 
  Parser_dhReadDouble(parser_dh, "-rtol", &(ctx->rtol));
  Parser_dhReadDouble(parser_dh, "-atol", &(ctx->atol));

  /* parallelization strategy (bj, pilu, none) */
  tmp = NULL;
  Parser_dhReadString(parser_dh, "-par", &tmp);
  if (tmp != NULL) {
    strcpy(ctx->algo_par, tmp);
  }
  if (Parser_dhHasSwitch(parser_dh, "-bj")) {
    strcpy(ctx->algo_par, "bj");
  }


  /* factorization parameters */
  Parser_dhReadDouble(parser_dh, "-rho", &(ctx->rho_init)); 
                                  /* inital storage allocation for factor */
  Parser_dhReadInt(parser_dh, "-level", &ctx->level);
  Parser_dhReadInt(parser_dh, "-pc_ilu_levels", &ctx->level);

  if (Parser_dhHasSwitch(parser_dh, "-ilut")) {
    Parser_dhReadDouble(parser_dh, "-ilut", &ctx->droptol);
    ctx->isScaled = true;
    strcpy(ctx->algo_ilu, "ilut");
  }

  /* make sure both algo_par and algo_ilu are set to "none,"
     if at least one is.
  */
  if (! strcmp(ctx->algo_par, "none")) { 
    strcpy(ctx->algo_ilu, "none");
  }
  else if (! strcmp(ctx->algo_ilu, "none")) { 
    strcpy(ctx->algo_par, "none");
  }


  Parser_dhReadDouble(parser_dh, "-sparseA",&(ctx->sparseTolA));  
                                        /* sparsify A before factoring */
  Parser_dhReadDouble(parser_dh, "-sparseF",&(ctx->sparseTolF));  
                                        /* sparsify after factoring */
  Parser_dhReadDouble(parser_dh, "-pivotMin", &(ctx->pivotMin));    
                                        /* adjust pivots if smaller than this */
  Parser_dhReadDouble(parser_dh, "-pivotFix", &(ctx->pivotFix));    
                                        /* how to adjust pivots */

  /* set row scaling for mandatory cases */
  if (ctx->sparseTolA || ! strcmp(ctx->algo_ilu, "ilut")) {
    ctx->isScaled = true;
  }

  /* solve method */
  tmp = NULL;
  Parser_dhReadString(parser_dh, "-ksp_type", &tmp);
  if (tmp != NULL) {
    strcpy(ctx->krylovMethod, tmp);

    if (! strcmp(ctx->krylovMethod, "bcgs")) {
      strcpy(ctx->krylovMethod, "bicgstab");
    }
  }
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "invert_diagonals_private"
void invert_diagonals_private(Euclid_dh ctx)
{
  START_FUNC_DH
  REAL_DH *aval = ctx->F->aval;
  int *diag = ctx->F->diag;
  if (aval == NULL || diag == NULL) {
    SET_INFO("can't invert diags; either F->aval or F->diag is NULL");
  } else {
    int i, m = ctx->F->m;
    for (i=0; i<m; ++i) {
        aval[diag[i]] = 1.0/aval[diag[i]]; 
    }
  }
  END_FUNC_DH
}


/* records rho_final (max for all processors) */
#undef __FUNC__
#define __FUNC__ "compute_rho_private"
void compute_rho_private(Euclid_dh ctx)
{
  START_FUNC_DH
  if (ctx->F != NULL) {
    double bufLocal[3], bufGlobal[3];
    int m = ctx->m;

    ctx->stats[NZF_STATS] = (double)ctx->F->rp[m];
    bufLocal[0] = ctx->stats[NZA_STATS];      /* nzA */
    bufLocal[1] = ctx->stats[NZF_STATS];      /* nzF */
    bufLocal[2] = ctx->stats[NZA_USED_STATS]; /* nzA used */

    if (np_dh == 1) {
      bufGlobal[0] = bufLocal[0];
      bufGlobal[1] = bufLocal[1];
      bufGlobal[2] = bufLocal[2];
    } else {
      MPI_Reduce(bufLocal, bufGlobal, 3, MPI_DOUBLE, MPI_SUM, 0, comm_dh);
    }

    if (myid_dh == 0) {

      /* compute rho */
      if (bufGlobal[0] && bufGlobal[1]) {
        ctx->rho_final = bufGlobal[1]/bufGlobal[0];
      } else {
        ctx->rho_final = -1;
      }

      /* compute ratio of nonzeros in A that were used */
      if (bufGlobal[0] && bufGlobal[2]) {
        ctx->stats[NZA_RATIO_STATS] = 100.0*bufGlobal[2]/bufGlobal[0];
      } else {
        ctx->stats[NZA_RATIO_STATS] = 100.0;
      }
    }
  }
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "factor_private"
void factor_private(Euclid_dh ctx) 
{
  START_FUNC_DH
  /*-------------------------------------------------------------
   * special case, for testing/debugging: no preconditioning
   *-------------------------------------------------------------*/
  if (! strcmp(ctx->algo_par, "none")) {
    goto DO_NOTHING;
  }

  /*-------------------------------------------------------------
   * Initialize object to hold factor.
   *-------------------------------------------------------------*/
  { int br = 0;
    int id = np_dh;
    if (ctx->sg != NULL) {
      br = ctx->sg->beg_rowP[myid_dh];
      id = ctx->sg->o2n_sub[myid_dh];
    }
    Factor_dhInit(ctx->A, true, true, ctx->rho_init, id, br, &(ctx->F)); CHECK_V_ERROR;
    ctx->F->bdry_count = ctx->sg->bdry_count[myid_dh];
    ctx->F->first_bdry = ctx->F->m - ctx->F->bdry_count;
    if (! strcmp(ctx->algo_par, "bj")) ctx->F->blockJacobi = true;
    if (Parser_dhHasSwitch(parser_dh, "-bj")) ctx->F->blockJacobi = true;
  }

  /*-------------------------------------------------------------
   * single mpi task with single or multiple subdomains
   *-------------------------------------------------------------*/
  if (np_dh == 1) {

    /* ILU(k) factorization */
    if (! strcmp(ctx->algo_ilu, "iluk")) {
      ctx->from = 0;
      ctx->to = ctx->m;

      /* only for debugging: use ilu_mpi_pilu */
      if (Parser_dhHasSwitch(parser_dh, "-mpi")) {
        if (ctx->sg != NULL && ctx->sg->blocks > 1) {
          SET_V_ERROR("only use -mpi, which invokes ilu_mpi_pilu(), for np = 1 and -blocks 1");
        }
        iluk_mpi_pilu(ctx); CHECK_V_ERROR;
      }

      /* "normal" operation */
      else  {
        iluk_seq_block(ctx); CHECK_V_ERROR;
        /* note: iluk_seq_block() performs block jacobi iluk if ctx->algo_par == bj.  */
      } 
    } 

    /* ILUT factorization */
    else if (! strcmp(ctx->algo_ilu, "ilut")) {
      ctx->from = 0;
      ctx->to = ctx->m;
      ilut_seq(ctx); CHECK_V_ERROR;
    }

    /* all other factorization methods */
    else {
        sprintf(msgBuf_dh, "factorization method: %s is not implemented", 
                                                                ctx->algo_ilu);
        SET_V_ERROR(msgBuf_dh);
    }
  }

  /*-------------------------------------------------------------
   * multiple mpi tasks with multiple subdomains
   *-------------------------------------------------------------*/
  else {
    /* block jacobi */
    if (! strcmp(ctx->algo_par, "bj")) {
      ctx->from = 0;
      ctx->to = ctx->m;
      iluk_mpi_bj(ctx); CHECK_V_ERROR;
    }

    /* iluk */
    else if (! strcmp(ctx->algo_ilu, "iluk")) {
      bool bj = ctx->F->blockJacobi;  /* for debugging */

      /* printf_dh("\n@@@ starting ilu_mpi_pilu @@@\n"); */

      SortedList_dhCreate(&(ctx->slist)); CHECK_V_ERROR;
      SortedList_dhInit(ctx->slist, ctx->sg); CHECK_V_ERROR;
      ExternalRows_dhCreate(&(ctx->extRows)); CHECK_V_ERROR;
      ExternalRows_dhInit(ctx->extRows, ctx); CHECK_V_ERROR;

      /* factor interior rows */
      ctx->from = 0;
      ctx->to = ctx->F->first_bdry;

/*
if (Parser_dhHasSwitch(parser_dh, "-test")) {
       printf("[%i] Euclid_dh :: TESTING ilu_seq\n", myid_dh);
       iluk_seq(ctx); CHECK_V_ERROR; 
} else {
       iluk_mpi_pilu(ctx); CHECK_V_ERROR; 
}
*/

       iluk_seq(ctx); CHECK_V_ERROR; 

      /* get external rows from lower ordered neighbors in the
         subdomain graph; these rows are needed for factoring
         this subdomain's boundary rows.
      */
      if (! bj) {
        ExternalRows_dhRecvRows(ctx->extRows); CHECK_V_ERROR;
      }

      /* factor boundary rows */
      ctx->from = ctx->F->first_bdry;
      ctx->to = ctx->F->m;
      iluk_mpi_pilu(ctx); CHECK_V_ERROR;

      /* send this processor's boundary rows to higher ordered
         neighbors in the subdomain graph.
      */
      if (! bj) {
        ExternalRows_dhSendRows(ctx->extRows); CHECK_V_ERROR;
      }

      /* discard column indices in factor if they would alter
         the subdomain graph (any such elements are in upper
         triangular portion of the row)
       */


      SortedList_dhDestroy(ctx->slist); CHECK_V_ERROR;
      ctx->slist = NULL;
      ExternalRows_dhDestroy(ctx->extRows); CHECK_V_ERROR;
      ctx->extRows = NULL;
    } 

    /* all other factorization methods */
    else {
        sprintf(msgBuf_dh, "factorization method: %s is not implemented", 
                                                                ctx->algo_ilu);
        SET_V_ERROR(msgBuf_dh);
    }
  }

DO_NOTHING: ;

  END_FUNC_DH
}

#if 0

#undef __FUNC__
#define __FUNC__ "discard_indices_private"
void discard_indices_private(Euclid_dh ctx)
{
  START_FUNC_DH
#if 0
  int *rp = ctx->F->rp, *cval = ctx->F->cval;
  double *aval = ctx->F->aval;
  int m = F->m, *nabors = ctx->nabors, nc = ctx->naborCount;
  int i, j, k, idx, count = 0, start_of_row;
  int beg_row = ctx->beg_row, end_row = beg_row + m;
  int *diag = ctx->F->diag;

  /* if col is not locally owned, and doesn't belong to a
   * nabor in the (original) subdomain graph, we need to discard
   * the column index and associated value.  First, we'll flag all
   * such indices for deletion.
   */
  for (i=0; i<m; ++i) {
    for (j=rp[i]; j<rp[i+1]; ++j) {
      int col = cval[j];
      if (col < beg_row  || col >= end_row) {
        bool flag = true;
        int owner = find_owner_private_mpi(ctx, col); CHECK_V_ERROR;

        for (k=0; k<nc; ++k) {
          if (nabors[k] == owner) {
            flag = false;
            break;
          }
        }

        if (flag) {
          cval[j] = -1; 
          ++count;
        }
      }
    }
  }

  sprintf(msgBuf_dh, "deleting %i indices that would alter the subdomain graph", count);
  SET_INFO(msgBuf_dh);

  /* Second, perform the actual deletion */
  idx = 0;
  start_of_row = 0;
  for (i=0; i<m; ++i) {
    for (j=start_of_row; j<rp[i+1]; ++j) {
      int    col = cval[j];
      double val = aval[j];
      if (col != -1) { 
        cval[idx] = col;
        aval[idx] = val;
        ++idx;
      } 
    }
    start_of_row = rp[i+1];
    rp[i+1] = idx;
  }

  /* rebuild diagonal pointers */
  for (i=0; i<m; ++i) {
    for (j=rp[i]; j<rp[i+1]; ++j) {
      if (cval[j] == i+beg_row) {
        diag[i] = j;
        break;
      }
    }
  }
#endif
  END_FUNC_DH
}
#endif

#undef __FUNC__
#define __FUNC__ "Euclid_dhSolve"
void Euclid_dhSolve(Euclid_dh ctx, Vec_dh x, Vec_dh b, int *its)
{
  START_FUNC_DH
  int itsOUT;
  Mat_dh A = (Mat_dh)ctx->A;

  if (! strcmp(ctx->krylovMethod, "cg")) {
    cg_euclid(A, ctx, x->vals, b->vals, &itsOUT); ERRCHKA;
  } else if (! strcmp(ctx->krylovMethod, "bicgstab")) {
    bicgstab_euclid(A, ctx, x->vals, b->vals, &itsOUT); ERRCHKA;
  } else {
    sprintf(msgBuf_dh, "unknown krylov solver: %s", ctx->krylovMethod);
    SET_V_ERROR(msgBuf_dh);
  }
  *its = itsOUT;
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "Euclid_dhPrintStats"
void Euclid_dhPrintStats(Euclid_dh ctx, FILE *fp)
{
  START_FUNC_DH
  double *timing;
  int nz;

  nz = Factor_dhReadNz(ctx->F); CHECK_V_ERROR;
  timing = ctx->timing;

  /* add in timing from lasst setup (if any) */
  ctx->timing[TOTAL_SOLVE_T] += ctx->timing[TOTAL_SOLVE_TEMP_T];
  ctx->timing[TOTAL_SOLVE_TEMP_T] = 0.0;

  reduce_timings_private(ctx); CHECK_V_ERROR;

  fprintf_dh(fp, "\n==================== Euclid report (start) ====================\n");
  fprintf_dh(fp, "\nruntime parameters\n");
  fprintf_dh(fp, "------------------\n");
  fprintf_dh(fp, "   setups:                 %i\n", ctx->setupCount);
  fprintf_dh(fp, "   tri solves:             %i\n", ctx->itsTotal);
  fprintf_dh(fp, "   parallelization method: %s\n", ctx->algo_par);
  fprintf_dh(fp, "   factorization method:   %s\n", ctx->algo_ilu);
  fprintf_dh(fp, "   matrix was row scaled:  %i\n", ctx->isScaled);

  fprintf_dh(fp, "   matrix row count:       %i\n", ctx->n);
  fprintf_dh(fp, "   nzF:                    %i\n", nz);
  fprintf_dh(fp, "   rho:                    %g\n", ctx->rho_final);
  fprintf_dh(fp, "   level:                  %i\n", ctx->level);
  fprintf_dh(fp, "   sparseA:                %g\n", ctx->sparseTolA);

  fprintf_dh(fp, "\nEuclid timing report\n");
  fprintf_dh(fp, "--------------------\n");
  fprintf_dh(fp, "   solves total:  %0.2f (see docs)\n", timing[TOTAL_SOLVE_T]);
  fprintf_dh(fp, "   tri solves:    %0.2f\n", timing[TRI_SOLVE_T]);
  fprintf_dh(fp, "   setups:        %0.2f\n", timing[SETUP_T]);
  fprintf_dh(fp, "      subdomain graph setup:  %0.2f\n", timing[SUB_GRAPH_T]);
  fprintf_dh(fp, "      factorization:          %0.2f\n", timing[FACTOR_T]);
  fprintf_dh(fp, "      solve setup:            %0.2f\n", timing[SOLVE_SETUP_T]);
  fprintf_dh(fp, "      rho:                    %0.2f\n", ctx->timing[COMPUTE_RHO_T]);
  fprintf_dh(fp, "      misc (should be small): %0.2f\n", 
                timing[SETUP_T] - 
               (timing[SUB_GRAPH_T]+timing[FACTOR_T]+
                timing[SOLVE_SETUP_T]+timing[COMPUTE_RHO_T]));

  if (ctx->sg != NULL) {
    SubdomainGraph_dhPrintStats(ctx->sg, fp); CHECK_V_ERROR;
    SubdomainGraph_dhPrintRatios(ctx->sg, fp); CHECK_V_ERROR;
  }


  fprintf_dh(fp, "\nApplicable if Euclid's internal solvers were used:\n");
  fprintf_dh(fp, "---------------------------------------------------\n");
  fprintf_dh(fp, "   solve method: %s\n", ctx->krylovMethod);
  fprintf_dh(fp, "   maxIts:       %i\n", ctx->maxIts);
  fprintf_dh(fp, "   rtol:         %g\n", ctx->rtol);
  fprintf_dh(fp, "   atol:         %g\n", ctx->atol);
  fprintf_dh(fp, "\n==================== Euclid report (end) ======================\n");
  END_FUNC_DH
}


/* nzA ratio and rho refer to most recent solve, if more than
   one solve (call to Setup) was performed.  Other stats
   are cumulative.
*/
#undef __FUNC__
#define __FUNC__ "Euclid_dhPrintStatsShort"
void Euclid_dhPrintStatsShort(Euclid_dh ctx, double setup, double solve, FILE *fp)
{
  START_FUNC_DH
  double *timing = ctx->timing;
  /* double *stats = ctx->stats; */
  /* double setup_factor; */
  /* double setup_other; */
  double apply_total;
  double apply_per_it;
  /* double nzUsedRatio; */
  double perIt;
  int blocks = np_dh;

  if (np_dh == 1) blocks = ctx->sg->blocks;

  reduce_timings_private(ctx); CHECK_V_ERROR;

  /* setup_factor  = timing[FACTOR_T]; */
  /* setup_other   = timing[SETUP_T] - setup_factor; */
  apply_total   = timing[TRI_SOLVE_T];
  apply_per_it  = apply_total/(double)ctx->its;
  /* nzUsedRatio   = stats[NZA_RATIO_STATS]; */
  perIt         = solve/(double)ctx->its;

  fprintf_dh(fp, "\n");
  fprintf_dh(fp, "%6s %6s %6s %6s %6s %6s %6s %6s %6s %6s XX\n",
                 "method", "subdms", "level", "its", "setup", "solve", "total", "perIt", "perIt", "rows");
  fprintf_dh(fp, "------  -----  -----  -----  -----  -----  -----  -----  -----  -----  XX\n");
  fprintf_dh(fp, "%6s %6i %6i %6i %6.2f %6.2f %6.2f %6.4f %6.5f %6g  XXX\n",
                 ctx->algo_par,     /* parallelization strategy [pilu, bj] */
                 blocks,            /* number of subdomains */
                 ctx->level,        /* level, for ILU(k) */
                 ctx->its,          /* iterations */
                 setup,             /* total setup time, from caller */
                 solve,             /* total setup time, from caller */
                 setup+solve,       /* total time, from caller */
                 perIt,              /* time per iteration, solver+precond. */
                 apply_per_it,     /* time per iteration, solver+precond. */
                 (double)ctx->n    /* global unknnowns */
             );




#if 0
  fprintf_dh(fp, "\n");
  fprintf_dh(fp, "%6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s XX\n",
             "", "","","","","setup","setup","","","","","","");

  fprintf_dh(fp, "%6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s XX\n",
                 "method", "subdms", "level", "its", "total", "factor",
                 "other", "apply", "perIt", "rho", "A_tol", "A_%", "rows");
  fprintf_dh(fp, "------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  ----- XX\n");


  fprintf_dh(fp, "%6s %6i %6i %6i %6.2f %6.2f %6.2f %6.2f %6.4f %6.1f %6g %6.2f %6g  XXX\n",
                 ctx->algo_par,     /* parallelization strategy [pilu, bj] */
                 blocks,            /* number of subdomains */
                 ctx->level,        /* level, for ILU(k) */
                 ctx->its,          /* iterations */
                 setup,             /* total setup time, from caller */
                 solve,             /* total setup time, from caller */
                 setup_factor,      /* pc solve: factorization */
                 setup_other,       /* pc setup: other */
                 apply_total,       /* triangular solve time */
                 apply_per_it,      /* time for one triangular solve */
                 ctx->rho_final,    /* rho */
                 ctx->sparseTolA,   /* sparseA tolerance */
                 nzUsedRatio,       /* percent of A that was used */
                 (double)ctx->n     /* global unknnowns */
            );
#endif

#if 0
  /* special: for scalability studies */
  fprintf_dh(fp, "\n%6s %6s %6s %6s %6s %6s WW\n", "method",  "level", "subGph", "factor", "solveS", "perIt");
  fprintf_dh(fp, "------  -----  -----  -----  -----  -----  WW\n");
  fprintf_dh(fp, "%6s %6i %6.2f %6.2f %6.2f %6.4f  WWW\n", 
               ctx->algo_par,
               ctx->level,
               timing[SUB_GRAPH_T], 
               timing[FACTOR_T], 
               timing[SOLVE_SETUP_T], 
               apply_per_it);
#endif
  END_FUNC_DH
}


/* its during last solve; rho; nzaUsed */
#undef __FUNC__
#define __FUNC__ "Euclid_dhPrintStatsShorter"
void Euclid_dhPrintStatsShorter(Euclid_dh ctx, FILE *fp)
{
  START_FUNC_DH
  double *stats = ctx->stats;

  int    its           = ctx->its;
  double rho           = ctx->rho_final;
  double nzUsedRatio   = stats[NZA_RATIO_STATS];


  fprintf_dh(fp, "\nStats from last linear solve: YY\n");
  fprintf_dh(fp, "%6s %6s %6s     YY\n", "its", "rho","A_%");
  fprintf_dh(fp, " -----  -----  -----     YY\n");
  fprintf_dh(fp, "%6i %6.2f %6.2f     YYY\n", its, rho, nzUsedRatio);

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "Euclid_dhPrintScaling"
void Euclid_dhPrintScaling(Euclid_dh ctx, FILE *fp)
{
  START_FUNC_DH
  int i, m = ctx->m;

  if (m > 10) m = 10;

  if (ctx->scale == NULL) {
    SET_V_ERROR("ctx->scale is NULL; was Euclid_dhSetup() called?");
  }

  fprintf(fp, "\n---------- 1st %i row scaling values:\n", m);
  for (i=0; i<m; ++i) {
    fprintf(fp, "   %i  %g  \n", i+1, ctx->scale[i]);
  }
  END_FUNC_DH
}


#undef __FUNC__
#define __FUNC__ "reduce_timings_private"
void reduce_timings_private(Euclid_dh ctx)
{
  START_FUNC_DH
  if (np_dh > 1) {
    double bufOUT[TIMING_BINS];

    memcpy(bufOUT, ctx->timing, TIMING_BINS*sizeof(double));
    MPI_Reduce(bufOUT, ctx->timing, TIMING_BINS, MPI_DOUBLE, MPI_MAX, 0, comm_dh);
  }

  ctx->timingsWereReduced = true;
  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "Euclid_dhPrintHypreReport"
void Euclid_dhPrintHypreReport(Euclid_dh ctx, FILE *fp)
{
  START_FUNC_DH
  double *timing;
  int nz;

  nz = Factor_dhReadNz(ctx->F); CHECK_V_ERROR;
  timing = ctx->timing;

  /* add in timing from lasst setup (if any) */
  ctx->timing[TOTAL_SOLVE_T] += ctx->timing[TOTAL_SOLVE_TEMP_T];
  ctx->timing[TOTAL_SOLVE_TEMP_T] = 0.0;

  reduce_timings_private(ctx); CHECK_V_ERROR;

 if (myid_dh == 0) {

  fprintf(fp, "@@@@@@@@@@@@@@@@@@@@@@ Euclid statistical report (start)\n");
  fprintf_dh(fp, "\nruntime parameters\n");
  fprintf_dh(fp, "------------------\n");
  fprintf_dh(fp, "   setups:                 %i\n", ctx->setupCount);
  fprintf_dh(fp, "   tri solves:             %i\n", ctx->itsTotal);
  fprintf_dh(fp, "   parallelization method: %s\n", ctx->algo_par);
  fprintf_dh(fp, "   factorization method:   %s\n", ctx->algo_ilu);
  if (! strcmp(ctx->algo_ilu, "iluk")) {
    fprintf_dh(fp, "      level:               %i\n", ctx->level);
  }

  if (ctx->isScaled) {
    fprintf_dh(fp, "   matrix was row scaled\n");
  }

  fprintf_dh(fp, "   global matrix row count: %i\n", ctx->n);
  fprintf_dh(fp, "   nzF:                     %i\n", nz);
  fprintf_dh(fp, "   rho:                     %g\n", ctx->rho_final);
  fprintf_dh(fp, "   sparseA:                 %g\n", ctx->sparseTolA);

  fprintf_dh(fp, "\nEuclid timing report\n");
  fprintf_dh(fp, "--------------------\n");
  fprintf_dh(fp, "   solves total:  %0.2f (see docs)\n", timing[TOTAL_SOLVE_T]);
  fprintf_dh(fp, "   tri solves:    %0.2f\n", timing[TRI_SOLVE_T]);
  fprintf_dh(fp, "   setups:        %0.2f\n", timing[SETUP_T]);
  fprintf_dh(fp, "      subdomain graph setup:  %0.2f\n", timing[SUB_GRAPH_T]);
  fprintf_dh(fp, "      factorization:          %0.2f\n", timing[FACTOR_T]);
  fprintf_dh(fp, "      solve setup:            %0.2f\n", timing[SOLVE_SETUP_T]);
  fprintf_dh(fp, "      rho:                    %0.2f\n", ctx->timing[COMPUTE_RHO_T]);
  fprintf_dh(fp, "      misc (should be small): %0.2f\n", 
                timing[SETUP_T] - 
               (timing[SUB_GRAPH_T]+timing[FACTOR_T]+
                timing[SOLVE_SETUP_T]+timing[COMPUTE_RHO_T]));

  if (ctx->sg != NULL) {
    SubdomainGraph_dhPrintStats(ctx->sg, fp); CHECK_V_ERROR;
    SubdomainGraph_dhPrintRatios(ctx->sg, fp); CHECK_V_ERROR;
  }

  fprintf(fp, "@@@@@@@@@@@@@@@@@@@@@@ Euclid statistical report (end)\n");

 }

  END_FUNC_DH
}

#undef __FUNC__
#define __FUNC__ "Euclid_dhPrintTestData"
void Euclid_dhPrintTestData(Euclid_dh ctx, FILE *fp)
{
  START_FUNC_DH
  /* Print data that should remain that will hopefully 
     remain the same for any platform.
     Possibly "tri solves" may change . . .
  */
  if (myid_dh == 0) {
    fprintf(fp, "   setups:                 %i\n", ctx->setupCount);
    fprintf(fp, "   tri solves:             %i\n", ctx->its);
    fprintf(fp, "   parallelization method: %s\n", ctx->algo_par);
    fprintf(fp, "   factorization method:   %s\n", ctx->algo_ilu);
    fprintf(fp, "   level:                  %i\n", ctx->level);
    fprintf(fp, "   row scaling:            %i\n", ctx->isScaled);
  }
  SubdomainGraph_dhPrintRatios(ctx->sg, fp); CHECK_V_ERROR;
  END_FUNC_DH
}
