// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef HAVE_QUADRATURE_QUAD_H
#define HAVE_QUADRATURE_QUAD_H

#include "phx_core_Constants.h"
#include "phx_quadrature_Element.h"

namespace phx {

namespace quadrature {

class Quad : public Element
{
public:

  Quad(const int numQuadrNodes)
  {
    numQuadrNodes_ = numQuadrNodes;
    if (numQuadrNodes_ == phx::core::Constants::MIN) numQuadrNodes_ = 1;
    if (numQuadrNodes_ == phx::core::Constants::MAX) numQuadrNodes_ = 9;

    numLocalNodes_ = 4;
    numBasisFunctions_ = 4;

    J_.Reshape(3,3);
    basis_rs_.Reshape(numLocalNodes_, numQuadrNodes_);
    basis_dr_.Reshape(numLocalNodes_, numQuadrNodes_);
    basis_ds_.Reshape(numLocalNodes_, numQuadrNodes_);
    basis_dt_.Reshape(numLocalNodes_, numQuadrNodes_);

    basis_xy_.Reshape(numLocalNodes_, 1);
    basis_dx_.Reshape(numLocalNodes_, 1);
    basis_dy_.Reshape(numLocalNodes_, 1);
    basis_dz_.Reshape(numLocalNodes_, 1);

    basis_rs_temp_.Reshape(numLocalNodes_, 1);
    basis_dr_temp_.Reshape(numLocalNodes_, 1);
    basis_ds_temp_.Reshape(numLocalNodes_, 1);
    basis_dt_temp_.Reshape(numLocalNodes_, 1);

    weight_.Reshape(numQuadrNodes_, 1);

    coord_.Reshape(numLocalNodes_, 3);
    for (int i = 0; i < numLocalNodes_; ++i)
      for (int j = 0; j < 3; ++j)
        coord_(i, j) = 0.0;

    qr_.Reshape(numQuadrNodes_, 1);
    qs_.Reshape(numQuadrNodes_, 1);

    switch (numQuadrNodes_) {
    case 1:      
      qr_[0]     = 0.0;
      qs_[0]     = 0.0;
      weight_[0] = 4.0;
      break;

    case 4:

      qr_[0]     =  -sqrt(3.0) / 3.0;
      qs_[0]     =  -sqrt(3.0) / 3.0;
      weight_[0] =   1.0;

      qr_[1]     =   sqrt(3.0) / 3.0;
      qs_[1]     =  -sqrt(3.0) / 3.0;
      weight_[1] =   1.0;

      qr_[2]     =   sqrt(3.0) / 3.0;
      qs_[2]     =   sqrt(3.0) / 3.0;
      weight_[2] =   1.0;

      qr_[3]     =  -sqrt(3.0) / 3.0;
      qs_[3]     =   sqrt(3.0) / 3.0;
      weight_[3] =   1.0;

      break;

    case 9:

      qr_[0]     =  -sqrt(3.0) / 5.0;
      qs_[0]     =  -sqrt(3.0) / 5.0;
      weight_[0] =   +25.0/81;

      qr_[1]     =  +0.0;
      qs_[1]     =  -sqrt(3.0) / 5.0;
      weight_[1] =  +40.0/81;

      qr_[2]     =  +sqrt(3.0) / 5.0;
      qs_[2]     =  -sqrt(3.0) / 5.0;
      weight_[2] =  +25.0/81;

      qr_[3]     =  -sqrt(3.0) / 5.0;
      qs_[3]     =  +0.0;
      weight_[3] =  +40.0/81;

      qr_[4]     =  +0.0;
      qs_[4]     =  +0.0;
      weight_[4] =  +64.0/81;

      qr_[5]     =  +sqrt(3.0) / 5.0;
      qs_[5]     =  +0.0;
      weight_[5] =  +40.0/81;

      qr_[6]     =  -sqrt(3.0) / 5.0;
      qs_[6]     =  +sqrt(3.0) / 5.0;
      weight_[6] =  +25.0/81;

      qr_[7]     =  +0.0;
      qs_[7]     =  +sqrt(3.0) / 5.0;
      weight_[7] =  +40.0/81;

      qr_[8]     =  +sqrt(3.0) / 5.0;
      qs_[8]     =  +sqrt(3.0) / 5.0;
      weight_[8] =  +25.0/81;

      break;

    default:
      TEST_FOR_EXCEPTION(true, std::out_of_range,
                         "Selected quadrature nodes, " << numQuadrNodes_ <<
                         ", not defined. Available choices are 1, 4, 9");
    }

    double x[4], y[4];

    x[0] = -1.0;  x[1] =  1.0;  x[2] =  1.0;  x[3] = -1.0;
    y[0] = -1.0;  y[1] = -1.0;  y[2] =  1.0;  y[3] =  1.0;

    for (int k = 0; k < numQuadrNodes_; k++) {
      for (int i = 0; i < 4; i++) {
        basis_rs_(i,k) = 0.25*(1+x[i] * qr_[k])*(1 + y[i] * qs_[k]);
        basis_dr_(i,k) = 0.25*   x[i]          *(1 + y[i] * qs_[k]);
        basis_ds_(i,k) = 0.25*(1+x[i] * qr_[k])*     y[i];
      }
    }
  }

  //! Destructor.
  ~Quad()
  {}

  virtual void computeJacobian(const int quadrNode) const
  {
    const double& x_0 = coord_(0, 0);
    const double& x_1 = coord_(1, 0);
    const double& x_2 = coord_(2, 0);
    const double& x_3 = coord_(3, 0);

    const double& y_0 = coord_(0, 1);
    const double& y_1 = coord_(1, 1);
    const double& y_2 = coord_(2, 1);
    const double& y_3 = coord_(3, 1);

    double divide_by;
    double ijacobian[2][2];

    double qr = qr_[quadrNode];
    double qs = qs_[quadrNode];

    /* transformation from the reference square to the actual one */
    ijacobian[0][0] = 0.25 * (-x_0 * (1-qr) + x_1 * (1-qr) + x_2 * (1+qr) - x_3 * (1+qr));
    ijacobian[0][1] = 0.25 * (-x_0 * (1-qs) - x_1 * (1+qs) + x_2 * (1+qs) + x_3 * (1-qs));
    ijacobian[1][0] = 0.25 * (-y_0 * (1-qr) + y_1 * (1-qr) + y_2 * (1+qr) - y_3 * (1+qr));
    ijacobian[1][1] = 0.25 * (-y_0 * (1-qs) - y_1 * (1+qs) + y_2 * (1+qs) + y_3 * (1-qs));

    det_J_ = ijacobian[0][0] * ijacobian[1][1] - ijacobian[0][1] * ijacobian[1][0];

    TEST_FOR_EXCEPTION(det_J_ == 0, std::logic_error,
                       "element has zero determinant, " << endl <<
                       "x = (" << x_0 << ", " << x_1 << ", " << x_2 << ", " << x_3 << "); "
                       "y = (" << y_0 << ", " << y_1 << ", " << y_2 << ", " << y_3 << "); ");

    divide_by = 1.0 / (det_J_);

    /* transformation from the actual to the reference */
    J_(0,0) =   divide_by * ijacobian[0][0];
    J_(1,0) = - divide_by * ijacobian[1][0];
    J_(0,1) = - divide_by * ijacobian[0][1];
    J_(1,1) =   divide_by * ijacobian[1][1];
  }

}; // class Triangle

} // namespace quadrature

} // namespace phx

#endif

