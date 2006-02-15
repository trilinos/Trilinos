#ifndef FNX_ELEMENTS_H
#define FNX_ELEMENTS_H

#include "Epetra_Object.h"

namespace fnx
{
  class QuadrElement : public Epetra_Object
  {
    public:
      virtual ~QuadrElement() {}

      virtual void SetQuadrNodes(const int NumQuadrNodes) = 0;
      virtual int NumQuadrNodes() const = 0;

      virtual double Phi(const int QuadrNode) const = 0;

      virtual double PhiDer(const int dim, const int QuadrNode) const = 0;

      //! Computes the Jacobian at the specified quadrature node.
      virtual void ComputeJacobian(const int QuadrNode,
                                   const double* x) = 0;
  
      //! Maps the quadrature nodes from the reference element to the actual one.
      virtual void ComputeQuadrNodes(const int QuadrNode, 
                                     const double* x,
                                     const double* xq) const = 0;
    
      //! Computes the derivatives at the specified quadrature node.
      virtual void ComputeDerivatives(const int QuadrNode) const = 0;

      //! Computes the weight at the specified quadrature node.
      virtual double QuadrWeight(const int QuadrNode) const = 0;

      //! Computes the determinant of the Jacobian matrix at the quadrature node.
      virtual double DetJacobian(const int QuadrNode) const = 0;
  };
  
  class TriangleQuadrElement : public QuadrElement
  {
    public:
      virtual TriangleQuadrElement() {}

      virtual ~TriangleQuadrElement() {}

      virtual void SetQuadrNodes(const int NumQuadrNodes) const
      {
        NumQuadrNodes_ = NumQuadrNodes;

        NumDimensions_ = 2;
        NumLocalNodes_ = 3;

        Weight_.resize(NumQuadrNodes_);

        qr_.Reshape(NumQuadrNodes_);
        qs_.Reshape(NumQuadrNodes_);

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
          FNX_THROW("Specified number of quadrature points not correct");
        }

        J_.Reshape(NumDimensions_ * NumDimensions_);
        phi_.Reshape(NumLocalNodes_, NumQuadrNodes_);
        phi_dr_.Reshape(NumLocalNodes_, NumQuadrNodes_);
        phi_ds_.Reshape(NumLocalNodes_, NumQuadrNodes_);

        phi_dx_.Reshape(NumLocalNodes_, 1);
        phi_dy_.Reshape(NumLocalNodes_, 1);

        for (int k = 0 ; k < NumQuadrNodes_ ; k++) 
        {
          phi_(0,k) = 1 - qr_(k) - qs_(k);
          phi_(1,k) = qr_(k);
          phi_(2,k) = qs_(k);

          phi_dr_(0,k) = -1;
          phi_dr_(1,k) = 1;
          phi_dr_(2,k) = 0;

          phi_ds_(0,k) = -1;
          phi_ds_(1,k) = 0;
          phi_ds_(2,k) = 1;
        }
      }

      inline double Phi(const int i) const
      {
        return(phi_(i));
      }

      inline double PhiDer(const int dim, const int i) const
      {
        if (dim == 0) return(phi_dx_(i));
        if (dim == 1) return(phi_dy_(i));
      }

      // coord(i,j) = j-th coord of i-th node
      void ComputeJacobian(const int QuadrNode, 
                           Epetra_SerialDenseMatrix& x) const
      {
        det_J_ = (- x(1,0) + x(1,2)) * (- x(0,0) + x(0,1)) -
                 (  x(1,0) - x(1,1)) * (  x(0,0) - x(0,2));

        if (det_J_ == 0) FNX_THROW("element with zero determinant!");

        double divide_by = 1.0 / det_J_;

        J_(0,0) = divide_by * (- x(1,0) + x(1,2));
        J_(1,0) = divide_by * (  x(0,0) - x(0,2));
        J_(0,1) = divide_by * (  x(1,0) - x(1,1));
        J_(1,1) = divide_by * (- x(0,0) + x(0,1));
      }

      void ComputeQuadrNodes(const int ii,
                             const Epetra_SerialDenseMatrix& x,
                             Epetra_SerialDenseMatrix& xq) const
      {
        for (int k = 0 ; k < NumLocalNodes_ ; k++)  xq(0,k) = xq(1,k) = 0.0;

        for (int k = 0 ; k < NumLocalNodes_ ; k++) 
        {
          xq(0, k) += phi_(k,ii) * x(0, k);
          xq(1, k) += phi_(k,ii) * x(1, k);
        }
      }

      void ComputeDerivatives(const int QuadrNode) const
      {
        for (int i = 0 ;  i <NumLocalNodes_ ; i++) 
        {
          phi_dx_(i) = phi_dr_(i,QuadrNode) * J_(0,0) + phi_ds_(i,QuadrNode) * J_(0,1);
          phi_dy_(i) = phi_dr_(i,QuadrNode) * J_(1,0) + phi_ds_(i,QuadrNode) * J_(1,1);
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

    private:
      int NumQuadrNodes_;
      int NumDimensions_;
      int NumLocalNodes_;

      double det_J_;

      Epetra_SerialDenseMatrix J_;
      Epetra_SerialDenseMatrix phi_;
      Epetra_SerialDenseMatrix phi_dr_;
      Epetra_SerialDenseMatrix phi_ds_;
      Epetra_SerialDenseVector Weight_, qr_, qs_;
      Epetra_SerialDenseVector phi_dx_;
      Epetra_SerialDenseVector phi_dy_;
  }
      
}; // namespace fnx

#endif
