// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_HEXQUADRATURE_H
#define GALERI_HEXQUADRATURE_H

/*!
 * \file Galeri_HexQuadrature.h
 */

#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Galeri_AbstractQuadrature.h"

namespace Galeri {
namespace FiniteElements {

/*!
 * \class HexQuadrature
 *
 * \brief Quadrature formula on hexahedra.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \author Last updated on 02-Apr-05.
 *
 */
class HexQuadrature : public AbstractQuadrature
{
    
public:

  //! Constructor.
  /*! 
   * \param NumQuadrNodes - (In) Number of quadrature nodes per element.
   *                             Valid choices are: 1.
   */
  HexQuadrature(const int NumQuadrNodes)
  {
    NumQuadrNodes_ = NumQuadrNodes;
    NumDimensions_ = 3;
    NumLocalNodes_ = 8;

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
      qs_[0]    = 0.0;
      qr_[0]    = 0.0;
      qt_[0]    = 0.0;
      Weight_[0] = 8.0;
      break;

    case 8:

      qs_[0]     =  -0.57735026918963;
      qr_[0]     =  -0.57735026918963;
      qt_[0]     =  -0.57735026918963;
      Weight_[0] = 1;
      qs_[1]     =   0.57735026918963;
      qr_[1]     =  -0.57735026918963;
      qt_[1]     =  -0.57735026918963;
      Weight_[1] = 1;
      qs_[2]     =   0.57735026918963;
      qr_[2]     =   0.57735026918963;
      qt_[2]     =  -0.57735026918963;    
      Weight_[2] = 1;
      qs_[3]     =   0.57735026918963;
      qr_[3]     =  -0.57735026918963;
      qt_[3]     =  -0.57735026918963;
      Weight_[3] = 1;
      qs_[4]     =  -0.57735026918963;
      qr_[4]     =  -0.57735026918963;
      qt_[4]     =   0.57735026918963;
      Weight_[4] = 1;
      qs_[5]     =   0.57735026918963;
      qr_[5]     =  -0.57735026918963;
      qt_[5]     =   0.57735026918963;
      Weight_[5] = 1;
      qs_[6]     =   0.57735026918963;
      qr_[6]     =   0.57735026918963;
      qt_[6]     =   0.57735026918963;    
      Weight_[6] = 1;
      qs_[7]     =   0.57735026918963;
      qr_[7]     =  -0.57735026918963;
      qt_[7]     =   0.57735026918963;
      Weight_[7] = 1;
      break;

    default:
      cerr << "The selected number of quadrature nodes ("
           << NumQuadrNodes_ << " is not available" << endl;
      cerr << "Valid choices are: 1, 8." << endl;
      throw(-1);
    }

    double x[8], y[8], z[8];

    x[0] = -1.0;
    x[1] =  1.0;
    x[2] =  1.0;
    x[3] = -1.0;
    x[4] = -1.0;
    x[5] =  1.0;
    x[6] =  1.0;
    x[7] = -1.0;

    y[0] = -1.0;
    y[1] = -1.0;
    y[2] =  1.0;
    y[3] =  1.0;
    y[4] = -1.0;
    y[5] = -1.0;
    y[6] =  1.0;
    y[7] =  1.0;

    z[0] = -1.0;
    z[1] = -1.0;
    z[2] = -1.0;
    z[3] = -1.0;
    z[4] =  1.0;
    z[5] =  1.0;
    z[6] =  1.0;
    z[7] =  1.0;

    for (int k = 0 ; k < NumQuadrNodes_ ; k++) 
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

  ~HexQuadrature()
  {}

  void ComputeJacobian(const int QuadrNode, 
                       const double* x_hex, 
                       const double* y_hex,
                       const double* z_hex) const
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

    double x[8], y[8], z[8];

    x[0] = -1.0;
    x[1] =  1.0;
    x[2] =  1.0;
    x[3] = -1.0;
    x[4] = -1.0;
    x[5] =  1.0;
    x[6] =  1.0;
    x[7] = -1.0;

    y[0] = -1.0;
    y[1] = -1.0;
    y[2] =  1.0;
    y[3] =  1.0;
    y[4] = -1.0;
    y[5] = -1.0;
    y[6] =  1.0;
    y[7] =  1.0;

    z[0] = -1.0;
    z[1] = -1.0;
    z[2] = -1.0;
    z[3] = -1.0;
    z[4] =  1.0;
    z[5] =  1.0;
    z[6] =  1.0;
    z[7] =  1.0;

    double qr = qr_[QuadrNode];
    double qs = qs_[QuadrNode];
    double qt = qt_[QuadrNode];

    for (int i = 0 ; i < 8 ; i++) 
    {
      a += 0.125 * x_hex[i] *    x[i]     * (1+y[i]*qs) * (1+z[i]*qt);
      b += 0.125 * x_hex[i] * (1+x[i]*qr) *    y[i]     * (1+z[i]*qt);
      c += 0.125 * x_hex[i] * (1+x[i]*qr) * (1+y[i]*qs) *    z[i]    ;

      d += 0.125 * y_hex[i] *    x[i]     * (1+y[i]*qs) * (1+z[i]*qt);
      e += 0.125 * y_hex[i] * (1+x[i]*qr) *    y[i]     * (1+z[i]*qt);
      f += 0.125 * y_hex[i] * (1+x[i]*qr) * (1+y[i]*qs) *    z[i]    ;

      g += 0.125 * z_hex[i] *    x[i]     * (1+y[i]*qs) * (1+z[i]*qt);
      h += 0.125 * z_hex[i] * (1+x[i]*qr) *    y[i]     * (1+z[i]*qt);
      l += 0.125 * z_hex[i] * (1+x[i]*qr) * (1+y[i]*qs) *    z[i]    ;
    }

    det_J_ = ( a * e * l - a * f * h - d * b * l + d * c * h +
              g * b * f - g * c * e );

    if (det_J_ < 0) det_J_ = - det_J_;

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
    return(NumQuadrNodes_);
  }

  inline int NumPhiFunctions() const 
  {
    return(8);
  }

  inline int NumPsiFunctions() const 
  {
    return(8);
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
