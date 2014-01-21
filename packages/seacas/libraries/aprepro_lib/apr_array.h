#ifndef SEAMS_ARRAY_H
#define SEAMS_ARRAY_H

/* Array structure */
struct array {
  double *data;
  int rows;
  int cols;
};

namespace SEAMS {
  array *array_add(const array *a, const array *b);
  array *array_sub(const array *a, const array *b);
  array *array_scale(const array *a, double s);
  array *array_mult(const array *a, const array *b);
}
#endif
