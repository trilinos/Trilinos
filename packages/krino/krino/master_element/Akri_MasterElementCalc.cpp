// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MasterElementCalc.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <cmath>

namespace krino {

void
MasterElementCalc::scalar_gradient(
        const int nint,       //:  number of intg points
        const int nelem,      //:  number of elements to process
        const int ndims,
        const int nnodes,
        const double* gradop, //: (nvec,npe,nelem,nint)
        const double* /*det_J*/,  //: (nelem,nint)
        const double* sfield, //: (npe,nelem)
        double* vector )      //: (nvec,nelem,nint)
{
  for ( int ip(0); ip < nint; ++ip ) {
    for ( int elem(0); elem < nelem; ++elem) {
      for ( int dim(0); dim < ndims; ++dim ) {
        double & val = vector[ip*ndims + dim];
        val = 0.0;
        for ( int node(0); node < nnodes; ++node ) {
          val += gradop[( (ip*nelem + elem)*nnodes + node)*ndims + dim] * sfield[ node ];
        }
      }
    }
  }
}

void
MasterElementCalc::determinant(
      const int num_elem_dims,
      const int num_coord_dims,
      const int  nint,
      const int npe_g,
      const double* deriv_g, // (num_elem_dims,npe_g,nint)
      const int   nelem,
      const double* coords,  // (num_coord_dims,npe,nelem)
      double* det_J,         // (nelem,nint)
      double* error )        // (nelem)
{
  if (num_elem_dims != num_coord_dims)
  {
    if (2 == num_elem_dims)
    {
      STK_ThrowAssert(3 == num_coord_dims);
      MasterElementCalc::determinant_element2d_in_3d(nint, npe_g, deriv_g, nelem, coords, det_J, error);
    }
    else
    {
      STK_ThrowAssert(1 == num_elem_dims);
      if (2 == num_coord_dims)
      {
        MasterElementCalc::determinant_element1d_in_2d(nint, npe_g, deriv_g, nelem, coords, det_J, error);
      }
      else
      {
        MasterElementCalc::determinant_element1d_in_3d(nint, npe_g, deriv_g, nelem, coords, det_J, error);
      }
    }
  }
  else
  {
    STK_ThrowAssert(num_elem_dims >= 2 && num_elem_dims <= 3);
    if (2 == num_elem_dims) MasterElementCalc::determinant_2d(nint, npe_g, deriv_g, nelem, coords, det_J, error);
    else MasterElementCalc::determinant_3d(nint, npe_g, deriv_g, nelem, coords, det_J, error);
  }
}

void
MasterElementCalc::gradient_operator(
        const int num_elem_dims,
        const int nint,
        const int npe_g,
        const double* deriv_g,   // (nvec,npe_g,nint)
        const int npe_f,
        const double* deriv_f,   // (nvec,npe_f,nint)
        const int   nelem,
        const double* coords,    // (nvec,npe,nelem)
        double* gradop,          // (nvec,npe,nelem,nint)
        double* det_J,           // (nelem,nint)
        double* error)
{
  STK_ThrowAssert(num_elem_dims >= 2 && num_elem_dims <= 3);
  if (2 == num_elem_dims) gradient_operator_2d(nint, npe_g, deriv_g, npe_f, deriv_f, nelem, coords, gradop, det_J, error);
  else gradient_operator_3d(nint, npe_g, deriv_g, npe_f, deriv_f, nelem, coords, gradop, det_J, error);
}

void
MasterElementCalc::determinant_2d(
        const int nint,
        const int npe,
        const double* deriv,   // (nvec,npe_g,nint)
        const int   nelem,
        const double* coords,    // (nvec,npe,nelem)
        double* det_J,           // (nelem,nint)
        double* error)
{
  auto c2d = [coords,npe](int d, int i, int e) { return coords[d + 2*(i + npe*e)]; };
  auto d2d = [deriv,npe](int d, int n, int q) { return deriv[d + 2*(n + npe*q)]; };
  auto detj = [det_J,nelem](int e, int q) -> double& { return det_J[e+nelem*q]; };

  for (int elem(0); elem < nelem; ++elem) error[elem] = 0.;

  for (int ke(0); ke < nelem; ++ke) {
    for ( int ki(0); ki < nint; ++ki ) {
      double dx_ds0 = 0.;
      double dx_ds1 = 0.;
      double dy_ds0 = 0.;
      double dy_ds1 = 0.;

      for ( int kn(0); kn < npe; ++kn ) {
        dx_ds0 += d2d(0,kn,ki)*c2d(0,kn,ke);
        dx_ds1 += d2d(1,kn,ki)*c2d(0,kn,ke);

        dy_ds0 += d2d(0,kn,ki)*c2d(1,kn,ke);
        dy_ds1 += d2d(1,kn,ki)*c2d(1,kn,ke);
      }

      detj(ke,ki) = dx_ds0*dy_ds1 - dy_ds0*dx_ds1;

      if ( detj(ke,ki) <= 0. )
      {
        error[ke] = 1.;
      }
    }
  }
}

void
MasterElementCalc::determinant_3d(
        const int nint,
        const int npe,
        const double* deriv,   // (nvec,npe_g,nint)
        const int   nelem,
        const double* coords,    // (nvec,npe,nelem)
        double* det_J,           // (nelem,nint)
        double* error)
{
  auto c3d = [coords,npe](int d, int i, int e) { return coords[d + 3*(i + npe*e)]; };
  auto d3d = [deriv,npe](int d, int n, int q) { return deriv[d + 3*(n + npe*q)]; };
  auto detj = [det_J,nelem](int e, int q) -> double& { return det_J[e+nelem*q]; };

  for (int elem(0); elem < nelem; ++elem) error[elem] = 0.;

  for (int ke(0); ke < nelem; ++ke) {
    for ( int ki(0); ki < nint; ++ki ) {
      double dx_ds0 = 0.;
      double dx_ds1 = 0.;
      double dx_ds2 = 0.;
      double dy_ds0 = 0.;
      double dy_ds1 = 0.;
      double dy_ds2 = 0.;
      double dz_ds0 = 0.;
      double dz_ds1 = 0.;
      double dz_ds2 = 0.;

      for ( int kn(0); kn < npe; ++kn ) {
        dx_ds0 += d3d(0,kn,ki)*c3d(0,kn,ke);
        dx_ds1 += d3d(1,kn,ki)*c3d(0,kn,ke);
        dx_ds2 += d3d(2,kn,ki)*c3d(0,kn,ke);

        dy_ds0 += d3d(0,kn,ki)*c3d(1,kn,ke);
        dy_ds1 += d3d(1,kn,ki)*c3d(1,kn,ke);
        dy_ds2 += d3d(2,kn,ki)*c3d(1,kn,ke);

        dz_ds0 += d3d(0,kn,ki)*c3d(2,kn,ke);
        dz_ds1 += d3d(1,kn,ki)*c3d(2,kn,ke);
        dz_ds2 += d3d(2,kn,ki)*c3d(2,kn,ke);
      }

      detj(ke,ki) = dx_ds0*( dy_ds1*dz_ds2 - dz_ds1*dy_ds2 )
                  + dy_ds0*( dz_ds1*dx_ds2 - dx_ds1*dz_ds2 )
                  + dz_ds0*( dx_ds1*dy_ds2 - dy_ds1*dx_ds2 );

      if ( detj(ke,ki) <= 0. )
      {
        error[ke] = 1.;
      }
    }
  }
}

void
MasterElementCalc::determinant_element2d_in_3d(
  const int nint,
  const int npe,
  const double* deriv,   // (nvec,npe_g,nint)
  const int   nelem,
  const double* coords,    // (nvec,npe,nelem)
  double* det_J,           // (nelem,nint)
  double* error)
{
  auto c3d = [coords,npe](int d, int i, int e) { return coords[d + 3*(i + npe*e)]; };
  auto d2d = [deriv,npe](int d, int n, int q) { return deriv[d + 2*(n + npe*q)]; };
  auto detj = [det_J,nelem](int e, int q) -> double& { return det_J[e+nelem*q]; };

  for (int elem(0); elem < nelem; ++elem) error[elem] = 0.;

  for (int ke(0); ke < nelem; ++ke) {
    for ( int ki(0); ki < nint; ++ki ) {
      double dx_ds0 = 0.;
      double dx_ds1 = 0.;
      double dy_ds0 = 0.;
      double dy_ds1 = 0.;
      double dz_ds0 = 0.;
      double dz_ds1 = 0.;

      for ( int kn(0); kn < npe; ++kn ) {
        dx_ds0 += d2d(0,kn,ki)*c3d(0,kn,ke);
        dx_ds1 += d2d(1,kn,ki)*c3d(0,kn,ke);

        dy_ds0 += d2d(0,kn,ki)*c3d(1,kn,ke);
        dy_ds1 += d2d(1,kn,ki)*c3d(1,kn,ke);

        dz_ds0 += d2d(0,kn,ki)*c3d(2,kn,ke);
        dz_ds1 += d2d(1,kn,ki)*c3d(2,kn,ke);
      }

      const double detXY = dx_ds0*dy_ds1 - dx_ds1*dy_ds0;
      const double detYZ = dy_ds0*dz_ds1 - dy_ds1*dz_ds0;
      const double detXZ =-dx_ds0*dz_ds1 + dx_ds1*dz_ds0;

      detj(ke,ki) = std::sqrt(detXY*detXY + detYZ*detYZ + detXZ*detXZ);

      if ( detj(ke,ki) <= 0. )
      {
        error[ke] = 1.;
      }
    }
  }
}

void
MasterElementCalc::determinant_element1d_in_3d(
  const int nint,
  const int npe,
  const double* deriv,   // (nvec,npe_g,nint)
  const int   nelem,
  const double* coords,    // (nvec,npe,nelem)
  double* det_J,           // (nelem,nint)
  double* error)
{
  auto c3d = [coords,npe](int d, int i, int e) { return coords[d + 3*(i + npe*e)]; };
  auto d1d = [deriv,npe](int d, int n, int q) { return deriv[d + n + npe*q]; };
  auto detj = [det_J,nelem](int e, int q) -> double& { return det_J[e+nelem*q]; };

  for (int elem(0); elem < nelem; ++elem) error[elem] = 0.;

  for (int ke(0); ke < nelem; ++ke) {
    for ( int ki(0); ki < nint; ++ki ) {
      double dx_ds = 0.;
      double dy_ds = 0.;
      double dz_ds = 0.;

      for ( int kn(0); kn < npe; ++kn ) {
        dx_ds += d1d(0,kn,ki)*c3d(0,kn,ke);
        dy_ds += d1d(0,kn,ki)*c3d(1,kn,ke);
        dz_ds += d1d(0,kn,ki)*c3d(2,kn,ke);
      }

      detj(ke,ki) = std::sqrt(dx_ds*dx_ds + dy_ds*dy_ds + dz_ds*dz_ds);

      if ( detj(ke,ki) <= 0. )
      {
        error[ke] = 1.;
      }
    }
  }
}

void
MasterElementCalc::determinant_element1d_in_2d(
  const int nint,
  const int npe,
  const double* deriv,   // (nvec,npe_g,nint)
  const int   nelem,
  const double* coords,    // (nvec,npe,nelem)
  double* det_J,           // (nelem,nint)
  double* error)
{
  auto c2d = [coords,npe](int d, int i, int e) { return coords[d + 2*(i + npe*e)]; };
  auto d1d = [deriv,npe](int d, int n, int q) { return deriv[d + n + npe*q]; };
  auto detj = [det_J,nelem](int e, int q) -> double& { return det_J[e+nelem*q]; };

  for (int elem(0); elem < nelem; ++elem) error[elem] = 0.;

  for (int ke(0); ke < nelem; ++ke) {
    for ( int ki(0); ki < nint; ++ki ) {
      double dx_ds = 0.;
      double dy_ds = 0.;

      for ( int kn(0); kn < npe; ++kn ) {
        dx_ds += d1d(0,kn,ki)*c2d(0,kn,ke);
        dy_ds += d1d(0,kn,ki)*c2d(1,kn,ke);
      }

      detj(ke,ki) = std::sqrt(dx_ds*dx_ds + dy_ds*dy_ds);

      if ( detj(ke,ki) <= 0. )
      {
        error[ke] = 1.;
      }
    }
  }
}

void
MasterElementCalc::gradient_operator_2d(
        const int nint,
        const int npe_g,
        const double* deriv_g,   // (nvec,npe_g,nint)
        const int npe_f,
        const double* deriv_f,   // (nvec,npe_f,nint)
        const int   nelem,
        const double* coords,    // (nvec,npe,nelem)
        double* gradop,          // (nvec,npe,nelem,nint)
        double* det_J,           // (nelem,nint)
        double* error)
{
  auto c2d = [coords,npe_g](int d, int i, int e) { return coords[d + 2*(i + npe_g*e)]; };
  auto d2d = [deriv_g,npe_g](int d, int n, int q) { return deriv_g[d + 2*(n + npe_g*q)]; };
  auto d2df = [deriv_f,npe_f](int d, int n, int q) { return deriv_f[d + 2*(n + npe_f*q)]; };
  auto detj = [det_J,nelem](int e, int q) -> double& { return det_J[e+nelem*q]; };
  auto g2d = [gradop,npe_f,nelem](int d, int n, int e, int q) -> double& { return gradop[d + 2*(n + npe_f*(e+nelem*q))]; };

  for (int elem(0); elem < nelem; ++elem) error[elem] = 0.;

  for (int ke(0); ke < nelem; ++ke) {
    for ( int ki(0); ki < nint; ++ki ) {
      double dx_ds0 = 0.;
      double dx_ds1 = 0.;
      double dy_ds0 = 0.;
      double dy_ds1 = 0.;

      for ( int kn(0); kn < npe_g; ++kn ) {
        dx_ds0 = dx_ds0+d2d(0,kn,ki)*c2d(0,kn,ke);
        dx_ds1 = dx_ds1+d2d(1,kn,ki)*c2d(0,kn,ke);

        dy_ds0 = dy_ds0+d2d(0,kn,ki)*c2d(1,kn,ke);
        dy_ds1 = dy_ds1+d2d(1,kn,ki)*c2d(1,kn,ke);
      }

      detj(ke,ki) = dx_ds0*dy_ds1 - dy_ds0*dx_ds1;

      double denom = 0.0;
      if ( detj(ke,ki) <= 0. )
      {
        error[ke] = 1.;
      }
      else
      {
        denom = 1./detj(ke,ki);
      }

      //  compute the gradient operators at the integration station -
      for ( int kn(0); kn < npe_f; ++kn ) {
        g2d(0,kn,ki,ke) = denom * ( d2df(0,kn,ki)*dy_ds1 - d2df(1,kn,ki)*dy_ds0 );
        g2d(1,kn,ki,ke) = denom * ( d2df(1,kn,ki)*dx_ds0 - d2df(0,kn,ki)*dx_ds1 );
      }
    }
  }
}

void
MasterElementCalc::gradient_operator_3d(
        const int nint,
        const int npe_g,
        const double* deriv_g,   // (nvec,npe_g,nint)
        const int npe_f,
        const double* deriv_f,   // (nvec,npe_f,nint)
        const int   nelem,
        const double* coords,    // (nvec,npe,nelem)
        double* gradop,          // (nvec,npe,nelem,nint)
        double* det_J,           // (nelem,nint)
        double* error)
{
  auto c3d = [coords,npe_g](int d, int i, int e) { return coords[d + 3*(i + npe_g*e)]; };
  auto d3d = [deriv_g,npe_g](int d, int n, int q) { return deriv_g[d + 3*(n + npe_g*q)]; };
  auto d3df = [deriv_f,npe_f](int d, int n, int q) { return deriv_f[d + 3*(n + npe_f*q)]; };
  auto detj = [det_J,nelem](int e, int q) -> double& { return det_J[e+nelem*q]; };
  auto g3d = [gradop,npe_f,nelem](int d, int n, int e, int q) -> double& { return gradop[d + 3*(n + npe_f*(e+nelem*q))]; };

  for (int elem(0); elem < nelem; ++elem) error[elem] = 0.;

  for (int ke(0); ke < nelem; ++ke) {
    for ( int ki(0); ki < nint; ++ki ) {
      double dx_ds0 = 0.;
      double dx_ds1 = 0.;
      double dx_ds2 = 0.;
      double dy_ds0 = 0.;
      double dy_ds1 = 0.;
      double dy_ds2 = 0.;
      double dz_ds0 = 0.;
      double dz_ds1 = 0.;
      double dz_ds2 = 0.;

      for ( int kn(0); kn < npe_g; ++kn ) {
        dx_ds0 = dx_ds0+d3d(0,kn,ki)*c3d(0,kn,ke);
        dx_ds1 = dx_ds1+d3d(1,kn,ki)*c3d(0,kn,ke);
        dx_ds2 = dx_ds2+d3d(2,kn,ki)*c3d(0,kn,ke);

        dy_ds0 = dy_ds0+d3d(0,kn,ki)*c3d(1,kn,ke);
        dy_ds1 = dy_ds1+d3d(1,kn,ki)*c3d(1,kn,ke);
        dy_ds2 = dy_ds2+d3d(2,kn,ki)*c3d(1,kn,ke);

        dz_ds0 = dz_ds0+d3d(0,kn,ki)*c3d(2,kn,ke);
        dz_ds1 = dz_ds1+d3d(1,kn,ki)*c3d(2,kn,ke);
        dz_ds2 = dz_ds2+d3d(2,kn,ki)*c3d(2,kn,ke);
      }

      detj(ke,ki) = dx_ds0*( dy_ds1*dz_ds2 - dz_ds1*dy_ds2 )
                  + dy_ds0*( dz_ds1*dx_ds2 - dx_ds1*dz_ds2 )
                  + dz_ds0*( dx_ds1*dy_ds2 - dy_ds1*dx_ds2 );

      double denom = 0.0;
      if ( detj(ke,ki) <= 0. )
      {
        error[ke] = 1.;
      }
      else
      {
        denom = 1./detj(ke,ki);
      }

      //  compute the gradient operators at the integration station -
      for ( int kn(0); kn < npe_f; ++kn ) {
        g3d(0,kn,ki,ke) = denom *
          ( d3df(0,kn,ki)*(dy_ds1*dz_ds2 - dz_ds1*dy_ds2)
          + d3df(1,kn,ki)*(dz_ds0*dy_ds2 - dy_ds0*dz_ds2)
          + d3df(2,kn,ki)*(dy_ds0*dz_ds1 - dz_ds0*dy_ds1) );

        g3d(1,kn,ki,ke) = denom *
          ( d3df(0,kn,ki)*(dz_ds1*dx_ds2 - dx_ds1*dz_ds2)
          + d3df(1,kn,ki)*(dx_ds0*dz_ds2 - dz_ds0*dx_ds2)
          + d3df(2,kn,ki)*(dz_ds0*dx_ds1 - dx_ds0*dz_ds1) );

        g3d(2,kn,ki,ke) = denom *
          ( d3df(0,kn,ki)*(dx_ds1*dy_ds2 - dy_ds1*dx_ds2)
          + d3df(1,kn,ki)*(dy_ds0*dx_ds2 - dx_ds0*dy_ds2)
          + d3df(2,kn,ki)*(dx_ds0*dy_ds1 - dy_ds0*dx_ds1) );
      }
    }
  }
}

} // namespace krino
