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

#ifndef HAVE_QUADRATURE_HEX_H
#define HAVE_QUADRATURE_HEX_H

#include "Galeri_core_Workspace.h"
#include "Galeri_quadrature_Element.h"

namespace Galeri {

namespace quadrature {

class Hex : public Element
{
public:

  Hex(const int numQuadrNodes)
  {
    numQuadrNodes_ = numQuadrNodes;
    if (numQuadrNodes_ == Galeri::core::Workspace::MIN) numQuadrNodes_ = 1;
    if (numQuadrNodes_ == Galeri::core::Workspace::MAX) numQuadrNodes_ = 8;

    numLocalNodes_ = 8;
    numBasisFunctions_ = 8;

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

    qr_.Reshape(numQuadrNodes_, 1);
    qs_.Reshape(numQuadrNodes_, 1);
    qt_.Reshape(numQuadrNodes_, 1);

    switch (numQuadrNodes_) {
    case 1:      
      qs_[0]    = 0.0;
      qr_[0]    = 0.0;
      qt_[0]    = 0.0;
      weight_[0] = 8.0;
      break;

    case 8:

      qs_[0]     =  -0.57735026918963;
      qr_[0]     =  -0.57735026918963;
      qt_[0]     =  -0.57735026918963;
      weight_[0] = 1;

      qs_[1]     =   0.57735026918963;
      qr_[1]     =  -0.57735026918963;
      qt_[1]     =  -0.57735026918963;
      weight_[1] = 1;

      qs_[2]     =   0.57735026918963;
      qr_[2]     =   0.57735026918963;
      qt_[2]     =  -0.57735026918963;    
      weight_[2] = 1;

      qs_[3]     =  -0.57735026918963;
      qr_[3]     =   0.57735026918963;
      qt_[3]     =  -0.57735026918963;
      weight_[3] = 1;

      qs_[4]     =  -0.57735026918963;
      qr_[4]     =  -0.57735026918963;
      qt_[4]     =   0.57735026918963;
      weight_[4] = 1;

      qs_[5]     =   0.57735026918963;
      qr_[5]     =  -0.57735026918963;
      qt_[5]     =   0.57735026918963;
      weight_[5] = 1;
      
      qs_[6]     =   0.57735026918963;
      qr_[6]     =   0.57735026918963;
      qt_[6]     =   0.57735026918963;    
      weight_[6] = 1;

      qs_[7]     =  -0.57735026918963;
      qr_[7]     =   0.57735026918963;
      qt_[7]     =   0.57735026918963;
      weight_[7] = 1;

      break;

    default:
      TEST_FOR_EXCEPTION(false, std::out_of_range,
                         "Selected quadrature nodes, " << numQuadrNodes_ <<
                         ", not defined. Available choices are: 1, 8");
    }

    double x[8], y[8], z[8];

    x[0] = -1.0; x[1] =  1.0; x[2] =  1.0; x[3] = -1.0; x[4] = -1.0; x[5] =  1.0; x[6] =  1.0; x[7] = -1.0;
    y[0] = -1.0; y[1] = -1.0; y[2] =  1.0; y[3] =  1.0; y[4] = -1.0; y[5] = -1.0; y[6] =  1.0; y[7] =  1.0;
    z[0] = -1.0; z[1] = -1.0; z[2] = -1.0; z[3] = -1.0; z[4] =  1.0; z[5] =  1.0; z[6] =  1.0; z[7] =  1.0;

    for (int k = 0 ; k < numQuadrNodes_ ; k++) 
    {
      for (int i = 0 ; i < 8 ; i++) 
      { 
        basis_rs_(i,k) = 0.125 * (1+x[i] * qr_[k]) * (1+y[i] * qs_[k]) * (1+z[i] * qt_[k]);
        basis_dr_(i,k) = 0.125*     x[i]           * (1+y[i] * qs_[k]) * (1+z[i] * qt_[k]);
        basis_ds_(i,k) = 0.125 * (1+x[i] * qr_[k]) *   y[i]            * (1+z[i] * qt_[k]);
        basis_dt_(i,k) = 0.125 * (1+x[i] * qr_[k]) * (1+y[i] * qs_[k]) *    z[i]          ;
      }
    }
  }

  //! Destructor.
  ~Hex()
  {}

  virtual void computeJacobian(const int quadrNode) const
  {
    double a = 0.0, b = 0.0, c = 0.0;
    double d = 0.0, e = 0.0, f = 0.0;
    double g = 0.0, h = 0.0, l = 0.0;
    double divide_by;

    /* jacobian^{-1} is the matrix

                   [ a b c ]
       jacobian =  [ d e f ]
                   [ g h l ]
     */

    // FIXME: MAKE STATIC VARIABLES
    double x[8], y[8], z[8];

    x[0] = -1.0; x[1] =  1.0; x[2] =  1.0; x[3] = -1.0; x[4] = -1.0; x[5] =  1.0; x[6] =  1.0; x[7] = -1.0;
    y[0] = -1.0; y[1] = -1.0; y[2] =  1.0; y[3] =  1.0; y[4] = -1.0; y[5] = -1.0; y[6] =  1.0; y[7] =  1.0;
    z[0] = -1.0; z[1] = -1.0; z[2] = -1.0; z[3] = -1.0; z[4] =  1.0; z[5] =  1.0; z[6] =  1.0; z[7] =  1.0;

    double qr = qr_[quadrNode];
    double qs = qs_[quadrNode];
    double qt = qt_[quadrNode];

    for (int i = 0 ; i < 8 ; i++) 
    {
      a += 0.125 * coord_(i, 0) *    x[i]     * (1+y[i]*qs) * (1+z[i]*qt);
      b += 0.125 * coord_(i, 0) * (1+x[i]*qr) *    y[i]     * (1+z[i]*qt);
      c += 0.125 * coord_(i, 0) * (1+x[i]*qr) * (1+y[i]*qs) *    z[i]    ;

      d += 0.125 * coord_(i, 1) *    x[i]     * (1+y[i]*qs) * (1+z[i]*qt);
      e += 0.125 * coord_(i, 1) * (1+x[i]*qr) *    y[i]     * (1+z[i]*qt);
      f += 0.125 * coord_(i, 1) * (1+x[i]*qr) * (1+y[i]*qs) *    z[i]    ;

      g += 0.125 * coord_(i, 2) *    x[i]     * (1+y[i]*qs) * (1+z[i]*qt);
      h += 0.125 * coord_(i, 2) * (1+x[i]*qr) *    y[i]     * (1+z[i]*qt);
      l += 0.125 * coord_(i, 2) * (1+x[i]*qr) * (1+y[i]*qs) *    z[i]    ;
    }

    det_J_ = ( a * e * l - a * f * h - d * b * l + d * c * h + g * b * f - g * c * e );

    if (det_J_ < 0) det_J_ = - det_J_;

    TEST_FOR_EXCEPTION(det_J_ == 0, std::logic_error,
                       "element has zero determinant");

    divide_by = - 1.0 / (det_J_);

    J_(0,0) = divide_by * (-e * l + f * h);
    J_(0,1) = divide_by * ( b * l - c * h); 
    J_(0,2) = divide_by * (-b * f + c * e);

    J_(1,0) = divide_by * ( d * l - f * g);
    J_(1,1) = divide_by * (-a * l + c * g);
    J_(1,2) = divide_by * ( a * f - c * d);

    J_(2,0) = divide_by * (-d * h + e * g);
    J_(2,1) = divide_by * ( a * h - b * g);
    J_(2,2) = divide_by * (-a * e + b * d);
  }

}; // class Tet

} // namespace quadrature

} // namespace Galeri

#endif


