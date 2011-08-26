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

#ifndef TPETRA_OPERATOR_HPP
#define TPETRA_OPERATOR_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Tpetra {

  //! \brief Abstract interface for linear operators accepting Tpetra MultiVector objects.
  /*!  This class is templated on \c Scalar, \c LocalOrdinal, \c GlobalOrdinal and \c Node. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type  Node is by defult of type Kokkos::DefaultNode::DefaultNodeType.

     A Operator object applies a linear operator to a MultiVector, storing the result in another MultiVector. The scalar type \c Scalar 
     of the Operator specifies the scalar field of the input and output MultiVector objects, not that of the underlying linear operator. Operator is an 
     abstract base class, and interfaces exist for this interface from numerous other classes, including sparse matrices, direct solvers, iterative solvers, 
     and preconditioners.
   */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Operator : virtual public Teuchos::Describable {
  public:

    /** \name Typedefs that give access to the template parameters. */
    //@{

    /// \typedef scalar_type
    /// \brief The type of the entries of the input and output multivectors.
    typedef Scalar scalar_type;

    /// \typedef local_ordinal_type
    /// \brief The local index type.
    typedef LocalOrdinal local_ordinal_type;

    /// \typedef global_ordinal_type
    /// \brief The global index type.
    typedef GlobalOrdinal global_ordinal_type;

    /// \typedef node_type
    /// \brief The Kokkos Node type.
    typedef Node node_type;

    //@}

    /** \name Pure virtual functions to be overridden by subclasses. */
    //@{

    //! Returns the Map associated with the domain of this operator, which must be compatible with X.getMap().
    virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const = 0;

    //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
    virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const = 0;

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
     */
    virtual void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
               Teuchos::ETransp mode = Teuchos::NO_TRANS, 
               Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
               Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const = 0;

    //! Indicates whether this operator supports applying the adjoint operator.
    virtual bool hasTransposeApply() const;

    //@}

  };

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::hasTransposeApply() const {
    return false;
  }

} // Tpetra namespace

#endif // TPETRA_OPERATOR_HPP
