// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MasterElementCalc_h
#define Akri_MasterElementCalc_h

namespace krino {

class MasterElementCalc
{
public:
  static void scalar_gradient(
      const int nint,       //:  number of intg points
      const int nelem,      //:  number of elements to process
      const int ndims,
      const int nnodes,
      const double* gradop, //: (nvec,npe,nelem,nint)
      const double* det_J,  //: (nelem,nint)
      const double* sfield, //: (npe,nelem)
      double* vector );     //: (nvec,nelem,nint)

  static void determinant(
      const int num_elem_dims,
      const int num_coord_dims,
      const int  nint,
      const int npe_g,
      const double* deriv_g,  // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,   // (nvec,npe,nelem)
      double* det_J,          // (nelem,nint)
      double* error );        // (nelem)

  static void gradient_operator(
      const int num_elem_dims,
      const int nint,
      const int npe_g,
      const double* deriv_g,  // (nvec,npe_g,nint)
      const int npe_f,
      const double* deriv_f,  // (nvec,npe_f,nint)
      const int   nelem,
      const double* coords,   // (nvec,npe,nelem)
      double* gradop,         // (nvec,npe,nelem,nint)
      double* det_J,          // (nelem,nint)
      double* error);         // (nelem)

  static void determinant_2d(
      const int nint,
      const int npe,
      const double* deriv,   // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,  // (nvec,npe,nelem)
      double* det_J,         // (nelem,nint)
      double* error);        // (nelem)
  static void determinant_3d(
      const int nint,
      const int npe,
      const double* deriv,   // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,  // (nvec,npe,nelem)
      double* det_J,         // (nelem,nint)
      double* error);        // (nelem)

  static void determinant_element2d_in_3d(
      const int nint,
      const int npe,
      const double* deriv,   // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,  // (nvec,npe,nelem)
      double* det_J,         // (nelem,nint)
      double* error);        // (nelem)

  static void determinant_element1d_in_3d(
      const int nint,
      const int npe,
      const double* deriv,   // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,  // (nvec,npe,nelem)
      double* det_J,         // (nelem,nint)
      double* error);        // (nelem)

  static void determinant_element1d_in_2d(
      const int nint,
      const int npe,
      const double* deriv,   // (nvec,npe_g,nint)
      const int   nelem,
      const double* coords,  // (nvec,npe,nelem)
      double* det_J,         // (nelem,nint)
      double* error);        // (nelem)

  static void gradient_operator_2d(
      const int nint,
      const int npe_g,
      const double* deriv_g, // (nvec,npe_g,nint)
      const int npe_f,
      const double* deriv_f, // (nvec,npe_f,nint)
      const int   nelem,
      const double* coords,  // (nvec,npe,nelem)
      double* gradop,        // (nvec,npe,nelem,nint)
      double* det_J,         // (nelem,nint)
      double* error);        // (nelem)

  static void gradient_operator_3d(
      const int nint,
      const int npe_g,
      const double* deriv_g, // (nvec,npe_g,nint)
      const int npe_f,
      const double* deriv_f, // (nvec,npe_f,nint)
      const int   nelem,
      const double* coords,  // (nvec,npe,nelem)
      double* gradop,        // (nvec,npe,nelem,nint)
      double* det_J,         // (nelem,nint)
      double* error);        // (nelem)
};

}  // end namespace krino

#endif // Akri_MasterElementCalc_h
