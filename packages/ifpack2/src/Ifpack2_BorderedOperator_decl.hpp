/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

//-----------------------------------------------------
// Ifpack2::BorderedOperator is a translation of the
// LOCA_BorderedSolver_EpetraHouseholder
// implementation written by Eric Phipps.
// DMD.
//------------------------------------------------------


#ifndef IFPACK2_BORDEREDOPERATOR_DECL_HPP
#define IFPACK2_BORDEREDOPERATOR_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

#include "Tpetra_Operator.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <iostream>

namespace Ifpack2 {
//! Ifpack2 bordered operator

/*!
  Ifpack2::BorderedOperator is a concrete class
extending Tpetra::Operator and defining an interface.

  This class extends Tpetra::Operator, providing the additional methods:
  - compute() computes everything required to apply the
    bordered operator, using matrix values  (and assuming that the
    sparsity of the matrix has not been changed);
  - isComputed() should return true if the bordered operator
    has been successfully computed, false otherwise.
  - getRHS() returns a reference to the matrix to be preconditioned.
  - getLHS() returns a reference to the matrix to be preconditioned.

The bordered operator is applied by apply()
(which returns if isComputed() is false). 
Each time compute() is called, the object re-computes the actual values of
the bordered operator.

<b>Title of Method Description</b>
*/

template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType >
class BorderedOperator : virtual public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node > {

  public:

    //! Destructor.
    virtual ~BorderedOperator(){}

    /** \name Methods implementing Tpetra::Operator. */
    //@{

    //! Returns the Map associated with the domain of this operator, which must be compatible with X.getMap().
    virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

    //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
    virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;


    bool hasTransposeApply() const;


    //! Applies the effect of the bordered operator.
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                       Teuchos::ETransp mode = Teuchos::NO_TRANS,
                       Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                       Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

    //@}

  //! constuctor with Tpetra::Operator input.
  BorderedOperator(const Teuchos::RCP< const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node > > &A);

  private:

  //! reference to the operator
  Teuchos::RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;

};//class BorderedOperator

}//namespace Ifpack2

#endif // IFPACK2_BORDEREDOPERATOR_DECL_HPP
