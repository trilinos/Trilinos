#include "aprepro.h"

namespace SEAMS {
  array *array_add(const array *a, const array *b)
  {
    int i;
    array *array_data = (array*) malloc(sizeof(array));
    array_data->rows = a->rows;
    array_data->cols = a->cols;

    /* Allocate space to store data... */
    array_data->data = (double*) calloc(a->rows*a->cols,sizeof(double));

    for (i=0; i < a->rows*a->cols; i++) {
      array_data->data[i] = a->data[i] + b->data[i];
    }
    return array_data;
  }

  array *array_sub(const array *a, const array *b)
  {
    int i;
    array *array_data = (array*) malloc(sizeof(array));
    array_data->rows = a->rows;
    array_data->cols = a->cols;

    /* Allocate space to store data... */
    array_data->data = (double*) calloc(a->rows*a->cols,sizeof(double));

    for (i=0; i < a->rows*a->cols; i++) {
      array_data->data[i] = a->data[i] - b->data[i];
    }
    return array_data;
  }

  array *array_scale(const array *a, double s)
  {
    int i;
    array *array_data = (array*) malloc(sizeof(array));
    array_data->rows = a->rows;
    array_data->cols = a->cols;

    /* Allocate space to store data... */
    array_data->data = (double*) calloc(a->rows*a->cols,sizeof(double));

    for (i=0; i < a->rows*a->cols; i++) {
      array_data->data[i] = a->data[i] * s;
    }

    return array_data;
  }

  array *array_mult(const array *a, const array *b)
  {
    int i,j, k;
    int ac = a->cols;
    int bc = b->cols;

    array *array_data = (array*) malloc(sizeof(array));
    array_data->rows = a->rows;
    array_data->cols = b->cols;

    /* Allocate space to store data... */
    array_data->data = (double*) calloc(a->rows*b->cols,sizeof(double));

    for (i=0; i<b->cols; i++) {
      for (j=0; j<a->rows; j++) {
	double sum = 0.0;
	for (k=0; k<a->cols; k++) {
	  sum += a->data[j*ac+k] * b->data[k*bc+i];
	}
	array_data->data[j*bc+i] = sum;
      }
    }
    return array_data;
  }
}
