#include "aprepro.h"

namespace SEAMS {
  array *array_add(const array *a, const array *b)
  {
    array *array_data = new array(a->rows, a->cols);
    for (int i=0; i < a->rows*a->cols; i++) {
      array_data->data[i] = a->data[i] + b->data[i];
    }
    return array_data;
  }

  array *array_sub(const array *a, const array *b)
  {
    array *array_data = new array(a->rows, a->cols);

    for (int i=0; i < a->rows*a->cols; i++) {
      array_data->data[i] = a->data[i] - b->data[i];
    }
    return array_data;
  }

  array *array_scale(const array *a, double s)
  {
    array *array_data = new array(a->rows, a->cols);

    for (int i=0; i < a->rows*a->cols; i++) {
      array_data->data[i] = a->data[i] * s;
    }

    return array_data;
  }

  array *array_mult(const array *a, const array *b)
  {
    int ac = a->cols;
    int bc = b->cols;

    array *array_data = new array(a->rows, b->cols);

    for (int i=0; i<b->cols; i++) {
      for (int j=0; j<a->rows; j++) {
	double sum = 0.0;
	for (int k=0; k<a->cols; k++) {
	  sum += a->data[j*ac+k] * b->data[k*bc+i];
	}
	array_data->data[j*bc+i] = sum;
      }
    }
    return array_data;
  }
}
