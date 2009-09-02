// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_INVERSEOPERATOR_HPP
#define TPETRA_INVERSEOPERATOR_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_RCP.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Tpetra {

  //! \brief Abstract interface for linear inverse operators accepting Tpetra MultiVector objects.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal, \c GlobalOrdinal and \c Node. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type.  Node is by defult of type Kokkos::DefaultNode::DefaultNodeType.

   A companion class to Tpetra::InverseOperator is Tpetra::Operator.  Both classes support polymorphic behavior for abstract linear operators. 
   However, Tpetra::Operator supports a "forward" operator, typically a matrix-vector multiplication, while Tpetra::InverseOperator supports an inversion
   such as triangular solves or the application of a preconditioner.  Full-featured classes such as Tpetra::CrsMatrix will implement both interfaces.
   Many user-defined classes will only implement one or the other, but there are important cases where both interfaces will be implemented by the same
   derived class. Such is the case for sophisticated multiphysics preconditioners.
   
   */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
	class InverseOperator : public Teuchos::Describable {
	public:

		/** \name Pure virtual functions to be overridden by subclasses. */
    //@{

		//! Returns the Map associated with the domain of this inverse operator, which must be compatible with Y.getMap().
		virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getInverseOperatorDomainMap() const = 0;

		//! Returns the Map associated with the range of this inverse operator, which must be compatible with X.getMap().
		virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getInverseOperatorRangeMap() const = 0;

    //! Computes the inverse operator-multivector operation, given the operator \f$M\f$ and input multivector \f$X\f$, find \f$Y\f$ such that \f$MY = A X\f$.
		virtual void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
						   MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
						   Teuchos::ETransp mode = Teuchos::NO_TRANS) const = 0;

    //@}

	};

} // Tpetra namespace

#endif // TPETRA_INVERSEOPERATOR_HPP
