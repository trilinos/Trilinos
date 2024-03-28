// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "aprepro.h" // for array, Aprepro, etc

#include <vector> // for vector

namespace SEAMS {
  extern SEAMS::Aprepro *aprepro;

  double array_interpolate(const array *arr, double row, double col)
  {
    /*
     * Bilinear interpolation.
     * Assumes equal grid spacing over the range
     * (0.0 -> rows-1) (0.0 -> cols-1)
     */

    int cols = arr->cols;
    int rows = arr->rows;

    int irl = row;
    int irh = rows > 1 ? irl + 1 : irl;
    int icl = col;
    int ich = cols > 1 ? icl + 1 : icl;

    double value = 0.0;

    if (irh < rows && ich < cols) {
      double v11 = arr->data[irl * cols + icl];
      double v21 = arr->data[irh * cols + icl];
      double v12 = arr->data[irl * cols + ich];
      double v22 = arr->data[irh * cols + ich];
      if (rows > 1 && cols > 1) {
        value = (v11 * (irh - row) + v21 * (row - irl)) * (ich - col) +
                (v12 * (irh - row) * v22 * (row - irl)) * (col - icl);
      }
      else if (rows > 1 && cols == 1) {
        value = v11 * (irh - row) + v21 * (row - irl);
      }
      else if (cols > 1 && rows == 1) {
        value = v11 * (ich - col) + v12 * (col - icl);
      }
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
    int    cols  = arr->cols;
    int    rows  = arr->rows;
    if (row >= 0 && row < rows && col >= 0 && col < cols) {
      if (row != static_cast<int>(row) || col != static_cast<int>(col)) {
        value = array_interpolate(arr, row, col);
      }
      else {
        int irow   = row;
        int icol   = col;
        int offset = irow * cols + icol;
        value      = arr->data[offset];
      }
    }
    else {
      aprepro->error("Row or Column index out of range");
    }
    return value;
  }

  array *array_add(const array *a, const array *b)
  {
    auto array_data = aprepro->make_array(a->rows, a->cols);
    for (int i = 0; i < a->rows * a->cols; i++) {
      array_data->data[i] = a->data[i] + b->data[i];
    }
    return array_data;
  }

  array *array_sub(const array *a, const array *b)
  {
    auto array_data = aprepro->make_array(a->rows, a->cols);

    for (int i = 0; i < a->rows * a->cols; i++) {
      array_data->data[i] = a->data[i] - b->data[i];
    }
    return array_data;
  }

  array *array_scale(const array *a, double s)
  {
    auto array_data = aprepro->make_array(a->rows, a->cols);

    for (int i = 0; i < a->rows * a->cols; i++) {
      array_data->data[i] = a->data[i] * s;
    }

    return array_data;
  }

  array *array_mult(const array *a, const array *b)
  {
    int ac = a->cols;
    int bc = b->cols;

    auto array_data = aprepro->make_array(a->rows, b->cols);

    for (int i = 0; i < b->cols; i++) {
      for (int j = 0; j < a->rows; j++) {
        double sum = 0.0;
        for (int k = 0; k < a->cols; k++) {
          sum += a->data[j * ac + k] * b->data[k * bc + i];
        }
        array_data->data[j * bc + i] = sum;
      }
    }
    return array_data;
  }
} // namespace SEAMS
