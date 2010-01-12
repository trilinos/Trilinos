//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_RK_BUTCHER_TABLEAU_BASE_HPP
#define RYTHMOS_RK_BUTCHER_TABLEAU_BASE_HPP

#include "Rythmos_Types.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace Rythmos {

/* \brief . */
template<class Scalar>
class RKButcherTableauBase : 
  virtual public Teuchos::Describable,
  virtual public Teuchos::ParameterListAcceptor,
  virtual public Teuchos::VerboseObject<RKButcherTableauBase<Scalar> >
{
public:
  /** \brief . */
  virtual int numStages() const = 0;
  /** \brief . */
  virtual const Teuchos::SerialDenseMatrix<int,Scalar>& A() const = 0;
  /** \brief . */
  virtual const Teuchos::SerialDenseVector<int,Scalar>& b() const = 0;
  /** \brief . */
  virtual const Teuchos::SerialDenseVector<int,Scalar>& c() const = 0;
  /** \brief . */
  virtual int order() const = 0;
  /** \brief . */
  virtual bool operator== (const RKButcherTableauBase<Scalar>& rkbt) const;
  /** \brief . */
  virtual void setDescription(std::string longDescription) = 0;
};


/* \brief . */
template<class Scalar>
bool RKButcherTableauBase<Scalar>::operator== (const RKButcherTableauBase<Scalar>& rkbt) const
{ 
  if (this->numStages() != rkbt.numStages()) {
    return false;
  }
  if (this->order() != rkbt.order()) {
    return false;
  }
  int N = rkbt.numStages();
  // Check b and c first:
  const Teuchos::SerialDenseVector<int,Scalar> b_ = this->b();
  const Teuchos::SerialDenseVector<int,Scalar> c_ = this->c();
  const Teuchos::SerialDenseVector<int,Scalar> other_b = rkbt.b();
  const Teuchos::SerialDenseVector<int,Scalar> other_c = rkbt.c();
  for (int i=0 ; i<N ; ++i) {
    if (b_(i) != other_b(i)) {
      return false;
    }
    if (c_(i) != other_c(i)) {
      return false;
    }
  }
  // Then check A:
  const Teuchos::SerialDenseMatrix<int,Scalar>& A_ = this->A();
  const Teuchos::SerialDenseMatrix<int,Scalar>& other_A = rkbt.A();
  for (int i=0 ; i<N ; ++i) {
    for (int j=0 ; j<N ; ++j) {
      if (A_(i,j) != other_A(i,j)) {
        return false;
      }
    } 
  }
  return true;
}

} // namespace Rythmos


#endif // RYTHMOS_RK_BUTCHER_TABLEAU_BASE_HPP
