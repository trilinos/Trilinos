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

#ifndef GALERI_QUADRATURE_ELEMENT_H
#define GALERI_QUADRATURE_ELEMENT_H

#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Galeri_core_Object.h"

namespace Galeri {
namespace quadrature {

class Element : public core::Object
{
public:
  //! Deastructor.
  virtual ~Element()
  {}

  inline double& operator()(const int i, const int j)
  {
    return(coord_(i, j));
  }

  inline const double& operator()(const int i, const int j) const
  {
    return(coord_(i, j));
  }

  virtual void computeJacobian(const int quadrNode) const = 0;

  void computeQuadrNodes(const int ii, double& xq, double& yq, double& zq) const
  {
    xq = 0.0;
    yq = 0.0;
    zq = 0.0;

    for (int k = 0; k < numLocalNodes_; k++) 
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

  void computeDerivatives(const int quadrNode) const
  {
    for (int i = 0; i < numLocalNodes_; i++) 
    {
      basis_dx_[i] = basis_dr_(i, quadrNode) * J_(0,0) +
                     basis_ds_(i, quadrNode) * J_(0,1) +
                     basis_dt_(i, quadrNode) * J_(0,2);
      basis_dy_[i] = basis_dr_(i, quadrNode) * J_(1,0) +
                     basis_ds_(i, quadrNode) * J_(1,1) +
                     basis_dt_(i, quadrNode) * J_(1,2);
      basis_dz_[i] = basis_dr_(i, quadrNode) * J_(2,0) +
                     basis_ds_(i, quadrNode) * J_(2,1) +
                     basis_dt_(i, quadrNode) * J_(2,2);
    }
  }

  inline double getQuadrWeight(const int quadrNode) const
  {
    return(weight_[quadrNode]);
  }

  inline double getDetJacobian(const int quadrNode) const
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
    return(basis_dz_[i]);
  }

  inline int getNumQuadrNodes() const
  {
    return(numQuadrNodes_);
  }

  inline int getNumBasisFunctions() const 
  {
    return(numBasisFunctions_);
  }

  virtual void print(ostream & os) const
  {
    cout << "Number of quadrature nodes = " << numQuadrNodes_ << endl;
    cout << coord_;
  }

protected:

  int numQuadrNodes_;
  int numLocalNodes_;
  int numBasisFunctions_;

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

  mutable Epetra_SerialDenseVector weight_;

  mutable Epetra_SerialDenseVector qr_;
  mutable Epetra_SerialDenseVector qs_;
  mutable Epetra_SerialDenseVector qt_;

  Epetra_SerialDenseMatrix coord_; 

  int ID_;
};

} // namespace quadrature
} // namespace Galeri
#endif
