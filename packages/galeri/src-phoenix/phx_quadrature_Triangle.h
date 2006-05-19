#ifndef HAVE_QUADRATURE_TRIANGLE_H
#define HAVE_QUADRATURE_TRIANGLE_H

#include "phx_quadrature_Element.h"

namespace phx {

namespace quadrature {

class Triangle : public Element
{
public:

  Triangle(const int NumQuadrNodes)
  {
    NumQuadrNodes_ = NumQuadrNodes;
    NumDimensions_ = 2;
    NumLocalNodes_ = 3;
    NumBasisFunctions_ = 3;

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

    for (int k = 0; k < NumQuadrNodes_; k++) 
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

  //! Deastructor.
  ~Triangle()
  {}

  virtual void ComputeJacobian(const int QuadrNode) const
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

    if (det_J_ == 0) 
    {
      cerr << "TriangleQuadrature: element has zero determinant" << endl;
      cerr << "Coordinates are: " << endl;
      cerr << "x = " << x_triangle_0 << " " << x_triangle_1 << " " << x_triangle_2 << endl;
      cerr << "y = " << y_triangle_0 << " " << y_triangle_1 << " " << y_triangle_2 << endl;
      throw(-1);
    }

    double divide_by = 1.0 / det_J_;

    J_(0,0) = divide_by * (- y_triangle_0 + y_triangle_2);
    J_(1,0) = divide_by * (  x_triangle_0 - x_triangle_2);
    J_(0,1) = divide_by * (  y_triangle_0 - y_triangle_1);
    J_(1,1) = divide_by * (- x_triangle_0 + x_triangle_1);
  }

}; // class Triangle

} // namespace quadrature

} // namespace phx

#endif
