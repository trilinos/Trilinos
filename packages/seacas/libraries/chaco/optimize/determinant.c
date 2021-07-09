/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

double determinant(double M[3][3], int ndims)
{
  if (ndims == 1) {
    return (M[0][0]);
  }
  if (ndims == 2) {
    return (M[0][0] * M[1][1] - M[0][1] * M[1][0]);
  }

  else if (ndims == 3) {
    return (M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) -
            M[1][0] * (M[0][1] * M[2][2] - M[2][1] * M[0][2]) +
            M[2][0] * (M[0][1] * M[1][2] - M[1][1] * M[0][2]));
  }

  else {
    return (0.0);
  }
}
