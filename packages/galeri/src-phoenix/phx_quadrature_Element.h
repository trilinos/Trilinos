#ifndef PHX_QUADRATURE_ELEMENT_H
#define PHX_QUADRATURE_ELEMENT_H

#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "phx_core_Object.h"

namespace phx {
namespace quadrature {

class Element : public core::Object
{
public:
  //! Deastructor.
  virtual ~Element()
  {}

  inline int& operator()(const int i)
  {
    return(ID_);
  }

  inline double& operator()(const int i, const int j)
  {
    return(coord_(i, j));
  }

  inline const double& operator()(const int i, const int j) const
  {
    return(coord_(i, j));
  }

  virtual void ComputeJacobian(const int QuadrNode) const = 0;

  void ComputeQuadrNodes(const int ii, double& xq, double& yq, double& zq) const
  {
    xq = 0.0;
    yq = 0.0;
    zq = 0.0;

    for (int k = 0; k < NumLocalNodes_; k++) 
    {
      xq += basis_rs_(k,ii) * coord_(k, 0);
      yq += basis_rs_(k,ii) * coord_(k, 1);
      zq += basis_rs_(k,ii) * coord_(k, 2);
      basis_dr_temp_[k] = basis_dr_(k,ii);
      basis_ds_temp_[k] = basis_ds_(k,ii);
      basis_dt_temp_[k] = basis_dt_(k,ii);
      basis_xy_[k] = basis_rs_(k,ii);
    }
  }

  void ComputeDerivatives(const int QuadrNode) const
  {
    for (int i = 0; i < NumLocalNodes_; i++) 
    {
      basis_dx_[i] = basis_dr_(i,QuadrNode) * J_(0,0) +
                     basis_ds_(i,QuadrNode) * J_(0,1);
      basis_dy_[i] = basis_dr_(i,QuadrNode) * J_(1,0) +
                     basis_ds_(i,QuadrNode) * J_(1,1);
    }
  }

  inline double getQuadrWeight(const int QuadrNode) const
  {
    return(Weight_[QuadrNode]);
  }

  inline double getDetJacobian(const int QuadrNode) const
  {
    return(det_J_);
  }

  inline double getPhi(const int i) const
  {
    return(basis_xy_[i]);
  }

  inline double getPhiX(const int i) const
  {
    return(basis_dx_[i]);
  }

  inline double getPhiY(const int i) const
  {
    return(basis_dy_[i]);
  }

  inline double getPhiZ(const int i) const
  {
    return(basis_dy_[i]);
  }

  inline int getNumQuadrNodes() const
  {
    return(NumQuadrNodes_);
  }

  inline int getNumBasisFunctions() const 
  {
    return(NumBasisFunctions_);
  }

  virtual void Print(ostream & os) const
  {
    cout << "Number of quadrature nodes = " << NumQuadrNodes_ << endl;
    cout << "Number of dimensions = " << NumDimensions_ << endl;
    cout << coord_;

  }

protected:

  int NumQuadrNodes_;
  int NumDimensions_;
  int NumLocalNodes_;
  int NumBasisFunctions_;

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

  Epetra_SerialDenseMatrix coord_; 

  int ID_;
};

} // namespace quadrature
} // namespace phx
#endif
