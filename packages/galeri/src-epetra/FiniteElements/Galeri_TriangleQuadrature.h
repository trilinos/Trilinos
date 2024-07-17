// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_TRIANGLEQUADRATURE_H
#define GALERI_TRIANGLEQUADRATURE_H

/*!
 * \file Galeri_TriangleQuadrature.h
 */

#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Galeri_AbstractQuadrature.h"

namespace Galeri {
namespace FiniteElements {

/*!
 * \class TriangleQuadrature
 *
 * \brief Quadrature formula on triangles.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \author Last updated on Apr-05.
 *
 */
class TriangleQuadrature : public AbstractQuadrature
{
    
public:

  //! Constructor.
  /*! 
   * \param NumQuadrNodes - (In) Number of quadrature nodes per element.
   *                             Valid choices are: 1, 3, 4, 7.
   */
  TriangleQuadrature(const int NumQuadrNodes)
  {
    NumQuadrNodes_ = NumQuadrNodes;
    NumDimensions_ = 2;
    NumLocalNodes_ = 3;

    J_.Reshape(NumDimensions_,NumDimensions_);
    basis_rs_.Reshape(NumLocalNodes_,NumQuadrNodes_);
    basis_dr_.Reshape(NumLocalNodes_,NumQuadrNodes_);
    basis_ds_.Reshape(NumLocalNodes_,NumQuadrNodes_);

    basis_xy_.Reshape(NumLocalNodes_, 1);
    basis_dx_.Reshape(NumLocalNodes_, 1);
    basis_dy_.Reshape(NumLocalNodes_, 1);

    basis_rs_temp_.Reshape(NumLocalNodes_, 1);
    basis_dr_temp_.Reshape(NumLocalNodes_, 1);
    basis_ds_temp_.Reshape(NumLocalNodes_, 1);

    Weight_.Reshape(NumQuadrNodes_, 1);

    qr_.Reshape(NumQuadrNodes_, 1);
    qs_.Reshape(NumQuadrNodes_, 1);

    switch (NumQuadrNodes_) {
    case 1:
      // up to order 1
      qs_[0]    = 1.0/3;
      qr_[0]    = 1.0/3;
      Weight_[0] = 0.5;
      break;

    case 3:
      // up to order 2
      qr_[0] =  0.5;
      qr_[1] =  0.5;
      qr_[2] =  0.0;

      qs_[0] =  0.0;  
      qs_[1] =  0.5;
      qs_[2] =  0.5;

      Weight_[0] = 1.0/6.0;
      Weight_[1] = 1.0/6.0;    
      Weight_[2] = 1.0/6.0;
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

      Weight_[0] = -27.0/96.0;
      Weight_[1] =  25.0/96.0;    
      Weight_[2] =  25.0/96.0;
      Weight_[2] =  25.0/96.0;
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

      Weight_[0] = 1.0/40.0; 
      Weight_[1] = 1.0/15.0;    
      Weight_[2] = 1.0/40.0;
      Weight_[3] = 1.0/15.0;
      Weight_[4] = 1.0/40.0;
      Weight_[5] = 1.0/15.0;
      Weight_[6] = 9.0/40.0;
      break;

    default:
      cerr << "Specified number of quadrature points ("
           << NumQuadrNodes << ") is not available." << endl;
      cerr << "Valid choices are: 1, 3, 4, 7." << endl;
      throw(-1);
    }

    for (int k = 0 ; k < NumQuadrNodes_ ; k++) 
    {
      basis_rs_(0,k) = 1 - qr_(k) - qs_(k);
      basis_rs_(1,k) = qr_(k);
      basis_rs_(2,k) = qs_(k);

      basis_dr_(0,k) = -1;
      basis_dr_(1,k) = 1;
      basis_dr_(2,k) = 0;

      basis_ds_(0,k) = -1;
      basis_ds_(1,k) = 0;
      basis_ds_(2,k) = 1;
    }
  }

  //! Deastructor.
  ~TriangleQuadrature()
  {}

  void ComputeJacobian(const int QuadrNode, 
                       const double* x_triangle, 
                       const double* y_triangle,
                       const double* z_triangle) const
  {

    det_J_ = (- y_triangle[0] + y_triangle[2]) *
             (- x_triangle[0] + x_triangle[1]) -
             (  y_triangle[0] - y_triangle[1]) *
             (  x_triangle[0] - x_triangle[2]);

    if (det_J_ == 0) 
    {
      cerr << "TriangleQuadrature: element has zero determinant" << endl;
      cerr << "Coordinates are: " << endl;
      cerr << "x = " << x_triangle[0] << " " << x_triangle[1] << " " << x_triangle[2] << endl;
      cerr << "y = " << y_triangle[0] << " " << y_triangle[1] << " " << y_triangle[2] << endl;
      cerr << "z = " << z_triangle[0] << " " << z_triangle[1] << " " << z_triangle[2] << endl;
      throw(-1);
    }

    double divide_by = 1.0 / det_J_;

    J_(0,0) = divide_by * (- y_triangle[0] + y_triangle[2]);
    J_(1,0) = divide_by * (  x_triangle[0] - x_triangle[2]);
    J_(0,1) = divide_by * (  y_triangle[0] - y_triangle[1]);
    J_(1,1) = divide_by * (- x_triangle[0] + x_triangle[1]);
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
      basis_dr_temp_[k] = basis_dr_(k,ii);
      basis_ds_temp_[k] = basis_ds_(k,ii);
      basis_xy_[k] = basis_rs_(k,ii);
    }
  }

  void ComputeDerivatives(const int QuadrNode) const
  {
    for (int i = 0 ;  i <NumLocalNodes_ ; i++) 
    {
      basis_dx_(i) = basis_dr_(i,QuadrNode) * J_(0,0) +
                     basis_ds_(i,QuadrNode) * J_(0,1);
      basis_dy_(i) = basis_dr_(i,QuadrNode) * J_(1,0) +
                     basis_ds_(i,QuadrNode) * J_(1,1);
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
    return(0.0);
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
    return(0.0);
  }

  inline int NumQuadrNodes() const
  {
    return(NumQuadrNodes_);
  }

  inline int NumPhiFunctions() const 
  {
    return(3);
  }

  inline int NumPsiFunctions() const 
  {
    return(3);
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
