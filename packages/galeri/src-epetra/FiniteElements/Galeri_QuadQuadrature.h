// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_QUADQUADRATURE_H
#define GALERI_QUADQUADRATURE_H

/*!
 * \file Galeri_QuadQuadrature.h
 */

#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Galeri_AbstractQuadrature.h"

namespace Galeri {
namespace FiniteElements {

/*!
 * \class QuadQuadrature
 *
 * \brief Quadrature formula on quadrilaterals.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \author Last updated on 31-Mar-05.
 *
 */
class QuadQuadrature : public AbstractQuadrature
{
    
public:

  //! Constructor.
  /*! 
   * \param NumQuadrNodes - (In) Number of quadrature nodes per element.
   *                             Valid choices are: 1, 4, 9.
   */
  QuadQuadrature(const int NumQuadrNodes)
  {
    NumQuadrNodes_ = NumQuadrNodes;
    NumDimensions_ = 2;
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

    switch (NumQuadrNodes_) {
    case 1:      
      qs_[0]     = 0.0;
      qr_[0]     = 0.0;
      Weight_[0] = 4.0;
      break;

    case 4:

      qs_[0]     =  -0.57735026918963;
      qr_[0]     =  -0.57735026918963;
      Weight_[0] =   1.0;

      qs_[1]     =   0.57735026918963;
      qr_[1]     =  -0.57735026918963;
      Weight_[1] =   1.0;
      qs_[2]     =   0.57735026918963;
      qr_[2]     =   0.57735026918963;

      Weight_[2] =   1.0;
      qs_[3]     =   0.57735026918963;
      qr_[3]     =  -0.57735026918963;
      Weight_[3] =   1.0;

      break;

    case 9:

      qs_[0]     =  -0.34641016151378; /* -sqrt(3)/5 */
      qr_[0]     =  -0.34641016151378;
      Weight_[0] =  25.0/81;

      qs_[1]     =  0.0;
      qr_[1]     =  -0.34641016151378;
      Weight_[1] = 40.0/81;

      qs_[2]     =  +0.34641016151378;
      qr_[2]     =  -0.34641016151378;
      Weight_[2] = 25.0/81;

      qs_[3]     =  -0.34641016151378;
      qr_[3]     =  0.0;
      Weight_[3] = 40.0/81;

      qs_[4]     =  0.0;
      qr_[4]     =  0.0;
      Weight_[4] = 64.0/81;

      qs_[5]     =  +0.34641016151378;
      qr_[5]     =  0.0;
      Weight_[5] = 40.0/81;

      qs_[6]     =  -0.34641016151378;
      qr_[6]     =  +0.34641016151378;
      Weight_[6] = 25.0/81;

      qs_[7]     =  0.0;
      qr_[7]     =  +0.34641016151378;
      Weight_[7] = 40.0/81;

      qs_[8]     =  +0.34641016151378;
      qr_[8]     =  +0.34641016151378;
      Weight_[8] = 25.0/81;

      break;

    default:
      cerr << "The selected number of quadrature nodes ("
           << NumQuadrNodes_ << " is not available" << endl;
      cerr << "Valid choices are: 1, 4, 9." << endl;
      throw(-1);
    }

    double x[4], y[4];

    x[0] = -1.0;  x[1] =  1.0;  x[2] =  1.0;  x[3] = -1.0;
    y[0] = -1.0;  y[1] = -1.0;  y[2] =  1.0;  y[3] =  1.0;

    for (int k = 0 ; k < NumQuadrNodes_ ; k++) {
      for (int i = 0 ; i < 4 ; i++) {
        basis_rs_(i,k) = 0.25*(1+x[i] * qr_[k])*(1 + y[i] * qs_[k]);
        basis_dr_(i,k) = 0.25*   x[i]          *(1 + y[i] * qs_[k]);
        basis_ds_(i,k) = 0.25*(1+x[i] * qr_[k])*     y[i];
      }
    }
  }

  //! Destructor.
  ~QuadQuadrature()
  {}

  void ComputeJacobian(const int QuadrNode, 
                       const double* x, 
                       const double* y,
                       const double* z) const
  {
    double divide_by;
    double ijacobian[2][2];

    double qr = qr_[QuadrNode];
    double qs = qs_[QuadrNode];

    /* transformation from the reference square to the actual one */
    ijacobian[0][0] = 0.25 * (-x[0] * (1-qs) + x[1] * (1-qs) + x[2] * (1+qs) - x[3] * (1+qs));
    ijacobian[0][1] = 0.25 * (-x[0] * (1-qr) - x[1] * (1+qr) + x[2] * (1+qr) + x[3] * (1-qr));
    ijacobian[1][0] = 0.25 * (-y[0] * (1-qs) + y[1] * (1-qs) + y[2] * (1+qs) - y[3] * (1+qs));
    ijacobian[1][1] = 0.25 * (-y[0] * (1-qr) - y[1] * (1+qr) + y[2] * (1+qr) + y[3] * (1-qr));

    det_J_ = ijacobian[0][0] * ijacobian[1][1] - ijacobian[0][1] * ijacobian[1][0];

    divide_by = 1.0 / (det_J_);

    /* transformation from the actual to the reference */
    J_(1,1) = divide_by * ijacobian[0][0];
    J_(0,1) = - divide_by * ijacobian[0][1];
    J_(1,0) = - divide_by * ijacobian[1][0];
    J_(0,0) = divide_by * ijacobian[1][1];
  }

  void ComputeQuadrNodes(const int QuadrNode,
                        const double* x, const double* y, const double* z,
                        double & xq, double & yq, double & zq ) const
  {
    xq = 0.0;
    yq = 0.0;
    zq = 0.0;

    for (int k = 0 ; k < NumLocalNodes_ ; k++) 
    {
      xq += basis_rs_(k,QuadrNode) * x[k];
      yq += basis_rs_(k,QuadrNode) * y[k];
      basis_dr_temp_[k] = basis_dr_(k,QuadrNode);
      basis_ds_temp_[k] = basis_ds_(k,QuadrNode);
      basis_xy_[k] = basis_rs_(k,QuadrNode);
    }
  }

  void ComputeDerivatives(const int QuadrNode) const
  {
    for (int i = 0 ; i < NumLocalNodes_ ; i++) 
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

  int NumPhiFunctions() const 
  {
    return(4);
  }

  int NumPsiFunctions() const 
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
