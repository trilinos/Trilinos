/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * ParaSails.h header file.
 *
 *****************************************************************************/

#include "Matrix.h"
#include "Numbering.h"
#include "PrunedRows.h"
#include "StoredRows.h"
#include "RowPatt.h"
#include "LoadBal.h"

#ifndef _PARASAILS_H
#define _PARASAILS_H

typedef struct
{
    int        symmetric;
    double     thresh;
    int        num_levels;
    double     filter;
    double     loadbal_beta;

    double     cost;          /* cost for this processor */
    double     setup_pattern_time;
    double     setup_values_time;

    Numbering *numb;
    Matrix    *M;             /* preconditioner */

    MPI_Comm   comm;
    int        beg_row;
    int        end_row;
    int       *beg_rows;
    int       *end_rows;
}
ParaSails;

ParaSails *ParaSailsCreate(MPI_Comm comm, int beg_row, int end_row, int sym);
void ParaSailsDestroy(ParaSails *ps);
void ParaSailsSetupPattern(ParaSails *ps, Matrix *A, 
  double thresh, int num_levels);
void ParaSailsSetupPatternExt(ParaSails *ps, Matrix *A, 
  double thresh_global, double thresh_local, int num_levels);
void ParaSailsSetupValues(ParaSails *ps, Matrix *A, double filter);
void ParaSailsApply(ParaSails *ps, double *u, double *v);
double ParaSailsStatsPattern(ParaSails *ps, Matrix *A);
void ParaSailsStatsValues(ParaSails *ps, Matrix *A);

#endif /* _PARASAILS_H */
