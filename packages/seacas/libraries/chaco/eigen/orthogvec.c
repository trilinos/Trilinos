/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

void orthogvec(double *vec1,     /* vector to be orthogonalized */
               int beg, int end, /* start and stop range for vector */
               double *vec2      /* vector to be orthogonalized against */
)
{
  double alpha;
  double dot();
  void   scadd();

  alpha = -dot(vec1, beg, end, vec2) / dot(vec2, beg, end, vec2);
  scadd(vec1, beg, end, alpha, vec2);
}

void orthogvec_float(float *vec1,      /* vector to be orthogonalized */
                     int beg, int end, /* start and stop range for vector */
                     float *vec2       /* vector to be orthogonalized against */
)
{
  float  alpha;
  double dot_float();
  void   scadd_float();

  alpha = -dot_float(vec1, beg, end, vec2) / dot_float(vec2, beg, end, vec2);
  scadd_float(vec1, beg, end, alpha, vec2);
}
