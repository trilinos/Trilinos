#ifndef _NAME_SPACE_SINGULAR_COARSE
#define _NAME_SPACE_SINGULAR_COARSE

#include "AztecOO.h"

#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_LAPACK.h"

namespace singularCoarse {
  double *Qcoarse;
  double *QcoarseTQcoarse;
  int rowQcoarse;
  int colQcoarse;
  const Epetra_Comm *commCoarse;
  //---------------------
  void setNullSpace(double *V, int row, int col, double *VtV, const Epetra_Comm *_Comm);
  void projection(double *z, int *options, int *proc_config, double *params,
                  AZ_MATRIX_STRUCT *Amat, AZ_PREC_STRUCT *prec);
}

#endif
