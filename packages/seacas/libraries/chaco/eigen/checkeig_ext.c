/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"
#include <stdio.h> // for fprintf, printf, NULL, FILE

/* Check an extended eigenpair of A by direct multiplication. Uses
   the Ay = extval*y + Dg form of the problem for convenience. */

double checkeig_ext(double *err, double *work, /* work vector of length n */
                    struct vtx_data **A, double *y, int n, double extval, double *vwsqrt,
                    double *gvec, double eigtol,
                    int warnings /* don't want to see warning messages in one of the
                                    contexts this is called */
)
{
  extern FILE *Output_File;   /* output file or null */
  extern int   DEBUG_EVECS;   /* print debugging output? */
  extern int   WARNING_EVECS; /* print warning messages? */
  double       resid;         /* the extended eigen residual */
  double       ch_norm();     /* vector norm */
  void         splarax();     /* sparse matrix vector mult */
  void         scadd();       /* scaled vector add */
  void         scale_diag();  /* scale vector by another's elements */
  void         cpvec();       /* vector copy */

  splarax(err, A, n, y, vwsqrt, work);
  scadd(err, 1, n, -extval, y);
  cpvec(work, 1, n, gvec); /* only need if going to re-use gvec */
  scale_diag(work, 1, n, vwsqrt);
  scadd(err, 1, n, -1.0, work);
  resid = ch_norm(err, 1, n);

  if (DEBUG_EVECS > 0) {
    printf("  extended residual: %g\n", resid);
    if (Output_File != NULL) {
      fprintf(Output_File, "  extended residual: %g\n", resid);
    }
  }
  if (warnings && WARNING_EVECS > 0 && resid > eigtol) {
    printf("WARNING: Extended residual (%g) greater than tolerance (%g).\n", resid, eigtol);
    if (Output_File != NULL) {
      fprintf(Output_File, "WARNING: Extended residual (%g) greater than tolerance (%g).\n", resid,
              eigtol);
    }
  }

  return (resid);
}
