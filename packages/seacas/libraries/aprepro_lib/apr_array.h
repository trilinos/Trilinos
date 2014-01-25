#ifndef SEAMS_ARRAY_H
#define SEAMS_ARRAY_H

namespace SEAMS {
  struct array;
  
  array *array_add(const array *a, const array *b);
  array *array_sub(const array *a, const array *b);
  array *array_scale(const array *a, double s);
  array *array_mult(const array *a, const array *b);
}
#endif
