/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef SEAMS_ARRAY_H
#define SEAMS_ARRAY_H

namespace SEAMS {
  struct array;

  double array_value(array *arr, double row, double col);
  array *array_add(const array *a, const array *b);
  array *array_sub(const array *a, const array *b);
  array *array_scale(const array *a, double s);
  array *array_mult(const array *a, const array *b);
} // namespace SEAMS
#endif
