#include "my_aprepro.h"
#include <stdlib.h>

extern aprepro_options ap_options;

double array_interpolate(array *arr, double row, double col)
{
  /*
   * Bilinear interpolation.
   * Assumes equal grid spacing over the range
   * (0.0 -> rows-1) (0.0 -> cols-1)
   */
  
  if (ap_options.one_based_index == True) {
    row--;
    col--;
  }

  int irl = row;
  int irh = irl+1;
  int icl = col;
  int ich = icl+1;

  int cols = arr->cols;
  int rows = arr->rows;

  double value = 0.0;
  
  if (irh < rows && ich < cols) {
    double v11 = arr->data[irl*cols+icl];
    double v21 = arr->data[irh*cols+icl];
    double v12 = arr->data[irl*cols+ich];
    double v22 = arr->data[irh*cols+ich];
    value =
      v11 * (irh - row) * (ich - col) + v21 * (row - irl) * (ich - col) +
      v12 * (irh - row) * (col - icl) + v22 * (row - irl) * (col - icl);
  }
  else {
    yyerror("Row or Column index out of range"); 
  }
  return value;
}

double array_value(array *arr, double row, double col)
{
  if (ap_options.one_based_index == True) {
    row--;
    col--;
  }
  double value = 0.0;
  int cols = arr->cols;
  int rows = arr->rows;
  if (row >= 0 && row < rows && col >= 0 && col < cols) {
    if (row != (int)row || col != (int)col) {
      value = array_interpolate(arr, row, col);
    }
    else {
      int irow = row;
      int icol = col;
      int offset = irow*cols+icol;
      value = arr->data[offset];
    }
  }
  else {
    yyerror("Row or Column index out of range"); 
  }
  return value;
}

array *array_construct(int rows, int cols)
{
  array *array_data = (array*) malloc(sizeof(array));
  if (array_data == NULL) {
    yyerror("Error allocating memory.");
    exit(EXIT_FAILURE);
  }
  array_data->rows = rows;
  array_data->cols = cols;

  /* Allocate space to store data... */
  array_data->data = (double*) calloc(rows*cols,sizeof(double));
  if (array_data->data == NULL) {
    yyerror("Error allocating memory.");
    exit(EXIT_FAILURE);
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
