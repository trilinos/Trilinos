/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"

/* Check an eigenpair of A by direct multiplication.  */
double checkeig(double *err, struct vtx_data **A, double *y, int n, double lambda, double *vwsqrt,
                double *work)
{
  double resid;
  double normy;
  double ch_norm();
  void   splarax(), scadd();

  splarax(err, A, n, y, vwsqrt, work);
  scadd(err, 1, n, -lambda, y);
  normy = ch_norm(y, 1, n);
  resid = ch_norm(err, 1, n) / normy;
  return (resid);
}
