#include "my_aprepro.h"
#include <stdlib.h>

array *array_construct(int rows, int cols)
{
  array *array_data = (array*) malloc(sizeof(array));
  if (array_data == NULL) {
    yyerror("Error allocating memory.");
  }
  array_data->rows = rows;
  array_data->cols = cols;

  /* Allocate space to store data... */
  array_data->data = (double*) calloc(rows*cols,sizeof(double));
  if (array_data->data == NULL) {
    yyerror("Error allocating memory.");
  }

  return array_data;
}

array *array_add(array *a, array *b)
{
  int i;
  array *array_data = array_construct(a->rows, a->cols);

  for (i=0; i < a->rows*a->cols; i++) {
    array_data->data[i] = a->data[i] + b->data[i];
  }
  return array_data;
}

array *array_sub(array *a, array *b)
{
  int i;
  array *array_data = array_construct(a->rows, a->cols);

  for (i=0; i < a->rows*a->cols; i++) {
    array_data->data[i] = a->data[i] - b->data[i];
  }
  return array_data;
}

array *array_scale(array *a, double s)
{
  int i;
  array *array_data = array_construct(a->rows, a->cols);

  for (i=0; i < a->rows*a->cols; i++) {
    array_data->data[i] = a->data[i] * s;
  }

  return array_data;
}

array *array_mult(array *a, array *b)
{
  int i,j, k;
  int ac = a->cols;
  int bc = b->cols;

  array *array_data = array_construct(a->rows, b->cols);

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
