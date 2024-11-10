// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HAVE_QUADRATURE_SEGMENT_H
#define HAVE_QUADRATURE_SEGMENT_H

#include "Galeri_quadrature_Element.h"

namespace Galeri {

namespace quadrature {

class Segment : public Element
{
public:

  Segment(const int numQuadrNodes)
  {
    numQuadrNodes_ = numQuadrNodes;
    numLocalNodes_ = 2;
    numBasisFunctions_ = 2;

    J_.Reshape(3, 3);
    basis_rs_.Reshape(numLocalNodes_,numQuadrNodes_);
    basis_dr_.Reshape(numLocalNodes_,numQuadrNodes_);
    basis_ds_.Reshape(numLocalNodes_,numQuadrNodes_);
    basis_dt_.Reshape(numLocalNodes_,numQuadrNodes_);

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

    // I got the following from:
    // http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
    
    switch (numQuadrNodes_) {
    case 1:
      // up to order 1
      qr_[0]     = 0.0;
      weight_[0] = 2.0;
      break;

    case 2:
      qr_[0]     = -(1.0 / 3) * sqrt(3.0);
      weight_[0] = 1.0;

      qr_[1]     = +(1.0 / 3) * sqrt(3.0);
      weight_[1] = 1.0;

      break;

    case 3:
      qr_[0]     = 0.0;
      weight_[0] = 8.0/ 9;

      qr_[1]     = -(1.0 / 5.0) * sqrt(15.0);
      weight_[1] = 5.0/ 9;

      qr_[2]     = +(1.0 / 5.0) * sqrt(15.0);
      weight_[2] = 5.0/ 9;

      break;

    case 4:
      qr_[0]     = -(1.0 / 35) * sqrt(525 - 70 * sqrt(30.0));
      weight_[0] =  (1.0 / 36) * (18 + sqrt(30.0));

      qr_[1]     = +(1.0 / 35) * sqrt(525 - 70 * sqrt(30.0));
      weight_[1] =  (1.0 / 36) * (18 + sqrt(30.0));

      qr_[2]     = -(1.0 / 35) * sqrt(525 + 70 * sqrt(30.0));
      weight_[2] =  (1.0 / 36) * (18 - sqrt(30.0));

      qr_[3]     = +(1.0 / 35) * sqrt(525 + 70 * sqrt(30.0));
      weight_[3] =  (1.0 / 36) * (18 - sqrt(30.0));

      break;

    case 5:
      qr_[0]     = 0.0;
      weight_[0] = 128.0 / 225.0;

      qr_[1]     = +(1.0 / 21 ) * sqrt(245 - 14 * sqrt(70.0));
      weight_[1] =  (1.0 / 900) * (322 + 13 * sqrt(70.0));

      qr_[2]     = -(1.0 / 21 ) * sqrt(245 - 14 * sqrt(70.0));
      weight_[2] =  (1.0 / 900) * (322 + 13 * sqrt(70.0));

      qr_[3]     = +(1.0 / 21 ) * sqrt(245 + 14 * sqrt(70.0));
      weight_[3] =  (1.0 / 900) * (322 - 13 * sqrt(70.0));

      qr_[4]     = -(1.0 / 21 ) * sqrt(245 + 14 * sqrt(70.0));
      weight_[4] =  (1.0 / 900) * (322 - 13 * sqrt(70.0));

      break;

    default:
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::out_of_range,
                         "Selected quadrature nodes, " << numQuadrNodes_ <<
                         ", not defined. Available choices are: 1, 2, 3, 4, 5");
    }

    for (int k = 0; k < numQuadrNodes_; k++) 
    {
      basis_rs_(0,k) = 0.5 - 0.5 * qr_[k];
      basis_rs_(1,k) = 0.5 + 0.5 * qr_[k];

      basis_dr_(0,k) = -0.5;
      basis_dr_(1,k) = 0.5;
    }
  }

  //! Destructor.
  ~Segment()
  {}

  virtual void computeJacobian(const int quadrNode) const
  {
    const double& x_0 = coord_(0, 0);
    const double& x_1 = coord_(1, 0);

    det_J_ = fabs(x_1 - x_0) / 2.0;

    TEUCHOS_TEST_FOR_EXCEPTION(det_J_ == 0, std::logic_error,
                       "element has zero determinant, " 
                       "x_0 = " << x_0 << ", x_1 = " << x_1);

    double divide_by = 1.0 / det_J_;

    J_(0,0) = divide_by;
  }

}; // class Segment

} // namespace quadrature

} // namespace Galeri

#endif
