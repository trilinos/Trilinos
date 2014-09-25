// Copyright (c) 2014, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <vector>                       // for vector
#include "aprepro.h"                    // for array, Aprepro, etc

namespace SEAMS {
  extern SEAMS::Aprepro *aprepro;

  double array_interpolate(const array *arr, double row, double col)
  {
    /*
     * Bilinear interpolation.
     * Assumes equal grid spacing over the range
     * (0.0 -> rows-1) (0.0 -> cols-1)
     */
  
    if (aprepro->ap_options.one_based_index) {
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
      aprepro->error("Row or Column index out of range"); 
    }
    return value;
  }

  double array_value(array *arr, double row, double col)
  {
    if (aprepro->ap_options.one_based_index) {
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
      aprepro->error("Row or Column index out of range"); 
    }
    return value;
  }

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
