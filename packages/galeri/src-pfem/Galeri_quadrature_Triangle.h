// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef HAVE_QUADRATURE_TRIANGLE_H
#define HAVE_QUADRATURE_TRIANGLE_H

#include "Galeri_core_Workspace.h"
#include "Galeri_quadrature_Element.h"

namespace Galeri {

namespace quadrature {

class Triangle : public Element
{
public:

  Triangle(const int numQuadrNodes)
  {
    numQuadrNodes_ = numQuadrNodes;
    if (numQuadrNodes_ == Galeri::core::Workspace::MIN) numQuadrNodes_ = 1;
    if (numQuadrNodes_ == Galeri::core::Workspace::MAX) numQuadrNodes_ = 7;

    numLocalNodes_ = 3;
    numBasisFunctions_ = 3;

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
    qs_.Reshape(numQuadrNodes_, 1);

    switch (numQuadrNodes_) {
    case 1:
      // up to order 1
      qs_[0]    = 1.0/3;
      qr_[0]    = 1.0/3;
      weight_[0] = 0.5;
      break;

    case 3:
      // up to order 2
      qr_[0] =  0.5;
      qr_[1] =  0.5;
      qr_[2] =  0.0;

      qs_[0] =  0.0;  
      qs_[1] =  0.5;
      qs_[2] =  0.5;

      weight_[0] = 1.0/6.0;
      weight_[1] = 1.0/6.0;    
      weight_[2] = 1.0/6.0;
      break;

    case 4:
      // up to order 3
      qr_[0] =  1.0/3;
      qr_[1] =  1.0/5;
      qr_[2] =  3.0/5;
      qr_[3] =  1.0/5;

      qs_[0] =  1.0/3;
      qs_[1] =  1.0/5;
      qs_[2] =  3.0/5;
      qs_[3] =  1.0/5;

      weight_[0] = -27.0/96.0;
      weight_[1] =  25.0/96.0;    
      weight_[2] =  25.0/96.0;
      weight_[2] =  25.0/96.0;
      break;

    case 7:
      // up to order 5
      qr_[0] =  0.0;
      qr_[1] =  0.5;
      qr_[2] =  1.0;
      qr_[3] =  0.5;
      qr_[4] =  0.0;
      qr_[5] =  0.0;
      qr_[6] =  1.0/3.0;

      qs_[0] =  0.0;  
      qs_[1] =  0.0;
      qs_[2] =  0.0;
      qs_[3] =  0.5;
      qs_[4] =  1.0;
      qs_[5] =  0.5;
      qs_[6] =  1.0/3.0;

      weight_[0] = 1.0/40.0; 
      weight_[1] = 1.0/15.0;    
      weight_[2] = 1.0/40.0;
      weight_[3] = 1.0/15.0;
      weight_[4] = 1.0/40.0;
      weight_[5] = 1.0/15.0;
      weight_[6] = 9.0/40.0;
      break;

    default:
      TEST_FOR_EXCEPTION(false, std::out_of_range,
                         "Selected quadrature nodes, " << numQuadrNodes_ <<
                         ", not defined. Available choices are 1, 3, 4, 7");
    }

    for (int k = 0; k < numQuadrNodes_; k++) 
    {
      basis_rs_(0,k) = 1 - qr_[k] - qs_[k];
      basis_rs_(1,k) = qr_[k];
      basis_rs_(2,k) = qs_[k];

      basis_dr_(0,k) = -1;
      basis_dr_(1,k) = 1;
      basis_dr_(2,k) = 0;

      basis_ds_(0,k) = -1;
      basis_ds_(1,k) = 0;
      basis_ds_(2,k) = 1;
    }
  }

  //! Destructor.
  ~Triangle()
  {}

  virtual void computeJacobian(const int quadrNode) const
  {
    const double& x_triangle_0 = coord_(0, 0);
    const double& x_triangle_1 = coord_(1, 0);
    const double& x_triangle_2 = coord_(2, 0);

    const double& y_triangle_0 = coord_(0, 1);
    const double& y_triangle_1 = coord_(1, 1);
    const double& y_triangle_2 = coord_(2, 1);

    det_J_ = (- y_triangle_0 + y_triangle_2) *
             (- x_triangle_0 + x_triangle_1) -
             (  y_triangle_0 - y_triangle_1) *
             (  x_triangle_0 - x_triangle_2);

    TEST_FOR_EXCEPTION(det_J_ == 0, std::logic_error,
                       "element has zero determinant, " << 
                       "x = (" << x_triangle_0 << ", " << x_triangle_1 << ", " << x_triangle_2 << "); " << 
                       "y = (" << y_triangle_0 << ", " << y_triangle_1 << ", " << y_triangle_2 << "); ");

    double divide_by = 1.0 / det_J_;

    J_(0,0) = divide_by * (- y_triangle_0 + y_triangle_2);
    J_(1,0) = divide_by * (  x_triangle_0 - x_triangle_2);
    J_(0,1) = divide_by * (  y_triangle_0 - y_triangle_1);
    J_(1,1) = divide_by * (- x_triangle_0 + x_triangle_1);
  }

}; // class Triangle

} // namespace quadrature

} // namespace Galeri

#endif
