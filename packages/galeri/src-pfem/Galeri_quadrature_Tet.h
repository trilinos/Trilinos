// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HAVE_QUADRATURE_TET_H
#define HAVE_QUADRATURE_TET_H

#include "Galeri_quadrature_Element.h"

namespace Galeri {

namespace quadrature {

class Tet : public Element
{
public:

  Tet(const int numQuadrNodes)
  {
    numQuadrNodes_ = numQuadrNodes;
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
      qr_[0]    = 1.0/4;
      qs_[0]    = 1.0/4;
      qt_[0]    = 1.0/4;
      weight_[0] = 1.0/6;
      break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::out_of_range,
                         "Selected quadrature nodes, " << numQuadrNodes_ <<
                         ", not defined. Available choices are: 1");
    }

    for (int k = 0 ; k < numQuadrNodes_ ; k++) 
    {
      basis_rs_(0,k) = 1.0 - qr_[k] -  qs_[k] - qt_[k];
      basis_rs_(1,k) = qr_[k];
      basis_rs_(2,k) = qs_[k];
      basis_rs_(3,k) = qt_[k];

      basis_dr_(0,k) = -1.0;
      basis_dr_(1,k) = 1.0;
      basis_dr_(2,k) = 0.0;
      basis_dr_(3,k) = 0.0;

      basis_ds_(0,k) = -1.0;
      basis_ds_(1,k) = 0.0;
      basis_ds_(2,k) = 1.0;
      basis_ds_(3,k) = 0.0;

      basis_dt_(0,k) = -1.0;
      basis_dt_(1,k) = 0.0;
      basis_dt_(2,k) = 0.0;
      basis_dt_(3,k) = 1.0;
    }
  }

  //! Destructor.
  ~Tet()
  {}

  virtual void computeJacobian(const int quadrNode) const
  {
    double a, b, c, d, e, f, g, h, l;
    double divide_by;

    /* jacobian^{-1} is the matrix

                   | a b c |
       jacobian =  | d e f |
                   | g h l |

       which transforms from the reference tet to the "real" tet.
       inv_jacobian is the Jacobian of the transformation from the
       "real" tetrahedron to the reference one.
       Finally, det_jacobian is the determinant of this latter transformation.
       */
    const double& x_0 = coord_(0, 0);
    const double& x_1 = coord_(1, 0);
    const double& x_2 = coord_(2, 0);
    const double& x_3 = coord_(3, 0);

    const double& y_0 = coord_(0, 1);
    const double& y_1 = coord_(1, 1);
    const double& y_2 = coord_(2, 1);
    const double& y_3 = coord_(3, 1);

    const double& z_0 = coord_(0, 2);
    const double& z_1 = coord_(1, 2);
    const double& z_2 = coord_(2, 2);
    const double& z_3 = coord_(3, 2);

    a = x_1 - x_0;
    b = x_2 - x_0;
    c = x_3 - x_0;

    d = y_1 - y_0;
    e = y_2 - y_0;
    f = y_3 - y_0;

    g = z_1 - z_0;
    h = z_2 - z_0;
    l = z_3 - z_0;

    det_J_ = (a * e * l - a * f * h - d * b * l + d * c * h +
              g * b * f - g * c * e);

    if (det_J_ < 0) det_J_ = - det_J_;

    TEUCHOS_TEST_FOR_EXCEPTION(det_J_ == 0, std::logic_error,
                       "element has zero determinant, " << endl <<
                       "x = (" << x_0 << ", " << x_1 << ", " << x_2 << ", " << x_3 << "); "
                       "y = (" << y_0 << ", " << y_1 << ", " << y_2 << ", " << y_3 << "); "
                       "z = (" << z_0 << ", " << z_1 << ", " << z_2 << ", " << z_3 << "); ");

    divide_by = - 1.0/(det_J_);

    J_(0,0) = divide_by * (-e * l + f * h);
    J_(1,0) = divide_by * ( b * l - c * h); 
    J_(2,0) = divide_by * (-b * f + c * e);

    J_(0,1) = divide_by * ( d * l - f * g);
    J_(1,1) = divide_by * (-a * l + c * g);
    J_(2,1) = divide_by * ( a * f - c * d);

    J_(0,2) = divide_by * (-d * h + e * g);
    J_(1,2) = divide_by * ( a * h - b * g);
    J_(2,2) = divide_by * (-a * e + b * d);
  }

}; // class Tet

} // namespace quadrature

} // namespace Galeri

#endif
