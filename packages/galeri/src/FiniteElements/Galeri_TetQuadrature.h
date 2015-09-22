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

#ifndef GALERI_TETQUADRATURE_H
#define GALERI_TETQUADRATURE_H

/*!
 * \file Galeri_TetQuadrature.h
 */

#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Galeri_AbstractQuadrature.h"

namespace Galeri {
namespace FiniteElements {

/*!
 * \class TetQuadrature
 *
 * \brief Quadrature formula on tetrahedra.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \author Last updated on Apr-05.
 *
 */
class TetQuadrature : public AbstractQuadrature
{
    
public:

  //! Constructor.
  /*! 
   * \param NumQuadrNodes - (In) Number of quadrature nodes per element.
   *                             Valid choices are: 1.
   */
  TetQuadrature(const int NumQuadrNodes)
  {
    NumQuadrNodes_ = NumQuadrNodes;
    NumDimensions_ = 3;
    NumLocalNodes_ = 4;

    J_.Reshape(NumDimensions_,NumDimensions_);
    basis_rs_.Reshape(NumLocalNodes_,NumQuadrNodes_);
    basis_dr_.Reshape(NumLocalNodes_,NumQuadrNodes_);
    basis_ds_.Reshape(NumLocalNodes_,NumQuadrNodes_);
    basis_dt_.Reshape(NumLocalNodes_,NumQuadrNodes_);

    basis_xy_.Reshape(NumLocalNodes_, 1);
    basis_dx_.Reshape(NumLocalNodes_, 1);
    basis_dy_.Reshape(NumLocalNodes_, 1);
    basis_dz_.Reshape(NumLocalNodes_, 1);

    basis_rs_temp_.Reshape(NumLocalNodes_, 1);
    basis_dr_temp_.Reshape(NumLocalNodes_, 1);
    basis_ds_temp_.Reshape(NumLocalNodes_, 1);
    basis_dt_temp_.Reshape(NumLocalNodes_, 1);

    Weight_.Reshape(NumQuadrNodes_, 1);

    qr_.Reshape(NumQuadrNodes_, 1);
    qs_.Reshape(NumQuadrNodes_, 1);
    qt_.Reshape(NumQuadrNodes_, 1);

    switch (NumQuadrNodes_) {
    case 1:      
      qr_[0]    = 1.0/4;
      qs_[0]    = 1.0/4;
      qt_[0]    = 1.0/4;
      Weight_[0] = 1.0/6;
      break;

    default:
      cerr << "The selected number of quadrature nodes ("
           << NumQuadrNodes_ << " is not available" << endl;
      cerr << "Valid choices are: 1." << endl;
      throw(-1);
    }

    for (int k = 0 ; k < NumQuadrNodes_ ; k++) 
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

  ~TetQuadrature()
  {}

  void ComputeJacobian(const int QuadrNode, 
                       const double* x, 
                       const double* y,
                       const double* z) const
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

    a = x[1]-x[0];
    b = x[2]-x[0];
    c = x[3]-x[0];

    d = y[1]-y[0];
    e = y[2]-y[0];
    f = y[3]-y[0];

    g = z[1]-z[0];
    h = z[2]-z[0];
    l = z[3]-z[0];

    det_J_ = (a * e * l - a * f * h - d * b * l + d * c * h +
              g * b * f - g * c * e);

    if (det_J_ < 0) det_J_ = - det_J_;

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

  void ComputeQuadrNodes(const int ii,
                         const double* x, const double* y, const double* z,
                         double& xq, double& yq, double& zq) const
  {
    xq = 0.0;
    yq = 0.0;
    zq = 0.0;

    for (int k = 0 ; k < NumLocalNodes_ ; k++) 
    {
      xq += basis_rs_(k,ii) * x[k];
      yq += basis_rs_(k,ii) * y[k];
      zq += basis_rs_(k,ii) * z[k];
      basis_dr_temp_[k] = basis_dr_(k,ii);
      basis_ds_temp_[k] = basis_ds_(k,ii);
      basis_dt_temp_[k] = basis_dt_(k,ii);
      basis_xy_[k] = basis_rs_(k,ii);
    }
  }

  void ComputeDerivatives(const int QuadrNode) const
  {
    for (int i = 0 ; i < NumLocalNodes_ ; i++) 
    {
      basis_dx_(i) = basis_dr_(i,QuadrNode) * J_(0,0) +
                     basis_ds_(i,QuadrNode) * J_(0,1) +
                     basis_dt_(i,QuadrNode) * J_(0,2);
      basis_dy_(i) = basis_dr_(i,QuadrNode) * J_(1,0) +
                     basis_ds_(i,QuadrNode) * J_(1,1) +
                     basis_dt_(i,QuadrNode) * J_(1,2);
      basis_dz_(i) = basis_dr_(i,QuadrNode) * J_(2,0) +
                     basis_ds_(i,QuadrNode) * J_(2,1) +
                     basis_dt_(i,QuadrNode) * J_(2,2);
    }
  }

  inline double QuadrWeight(const int QuadrNode) const
  {
    return(Weight_[QuadrNode]);
  }

  inline double DetJacobian(const int QuadrNode) const
  {
    return(det_J_);
  }

  inline double Phi(const int i) const
  {
    return(basis_xy_(i));
  }

  inline double PhiX(const int i) const
  {
    return(basis_dx_(i));
  }

  inline double PhiY(const int i) const
  {
    return(basis_dy_(i));
  }

  inline double PhiZ(const int i) const
  {
    return(basis_dz_(i));
  }

  inline double Psi(const int i) const
  {
    return(basis_xy_(i));
  }

  inline double PsiX(const int i) const
  {
    return(basis_dx_(i));
  }

  inline double PsiY(const int i) const
  {
    return(basis_dy_(i));
  }

  inline double PsiZ(const int i) const
  {
    return(basis_dz_(i));
  }

  inline int NumQuadrNodes() const
  {
    return NumQuadrNodes_;
  }

  inline int NumPhiFunctions() const 
  {
    return(4);
  }

  inline int NumPsiFunctions() const 
  {
    return(4);
  }

protected:

  int NumQuadrNodes_;
  int NumDimensions_;
  int NumLocalNodes_;

  mutable double det_J_;

  mutable Epetra_SerialDenseMatrix J_;

  mutable Epetra_SerialDenseMatrix basis_rs_;
  mutable Epetra_SerialDenseMatrix basis_dr_;
  mutable Epetra_SerialDenseMatrix basis_ds_;
  mutable Epetra_SerialDenseMatrix basis_dt_;

  mutable Epetra_SerialDenseVector basis_xy_;
  mutable Epetra_SerialDenseVector basis_dx_;
  mutable Epetra_SerialDenseVector basis_dy_;
  mutable Epetra_SerialDenseVector basis_dz_;

  mutable Epetra_SerialDenseVector basis_rs_temp_;
  mutable Epetra_SerialDenseVector basis_dr_temp_;
  mutable Epetra_SerialDenseVector basis_ds_temp_;
  mutable Epetra_SerialDenseVector basis_dt_temp_;

  mutable Epetra_SerialDenseVector Weight_;

  mutable Epetra_SerialDenseVector qr_;
  mutable Epetra_SerialDenseVector qs_;
  mutable Epetra_SerialDenseVector qt_;

};

} // namespace FiniteElements
} // namespace Galeri
#endif
