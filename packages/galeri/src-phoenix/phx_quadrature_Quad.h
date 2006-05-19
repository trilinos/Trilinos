#ifndef HAVE_QUADRATURE_QUAD_H
#define HAVE_QUADRATURE_QUAD_H

#include "phx_quadrature_Element.h"

namespace phx {

namespace quadrature {

class Quad : public Element
{
public:

  Quad(const int NumQuadrNodes)
  {
    NumQuadrNodes_ = NumQuadrNodes;
    NumDimensions_ = 2;
    NumLocalNodes_ = 4;
    NumBasisFunctions_ = 4;

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

    coord_.Reshape(NumLocalNodes_, 3);
    for (int i = 0; i < NumLocalNodes_; ++i)
      for (int j = 0; j < 3; ++j)
        coord_(i, j) = 0.0;

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
  ~Quad()
  {}

  virtual void ComputeJacobian(const int QuadrNode) const
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

    double qr = qr_[QuadrNode];
    double qs = qs_[QuadrNode];

    /* transformation from the reference square to the actual one */
    ijacobian[0][0] = 0.25 * (-x_0 * (1-qs) + x_1 * (1-qs) + x_2 * (1+qs) - x_3 * (1+qs));
    ijacobian[0][1] = 0.25 * (-x_0 * (1-qr) - x_1 * (1+qr) + x_2 * (1+qr) + x_3 * (1-qr));
    ijacobian[1][0] = 0.25 * (-y_0 * (1-qs) + y_1 * (1-qs) + y_2 * (1+qs) - y_3 * (1+qs));
    ijacobian[1][1] = 0.25 * (-y_0 * (1-qr) - y_1 * (1+qr) + y_2 * (1+qr) + y_3 * (1-qr));

    det_J_ = ijacobian[0][0] * ijacobian[1][1] - ijacobian[0][1] * ijacobian[1][0];

    assert (det_J_ != 0.0);

    divide_by = 1.0 / (det_J_);

    /* transformation from the actual to the reference */
    J_(0,0) =   divide_by * ijacobian[0][0];
    J_(1,0) = - divide_by * ijacobian[0][1];
    J_(0,1) = - divide_by * ijacobian[1][0];
    J_(1,1) =   divide_by * ijacobian[1][1];
  }

}; // class Triangle

} // namespace quadrature

} // namespace phx

#endif

