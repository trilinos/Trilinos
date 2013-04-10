// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_OPERATOR_HPP
#define TPETRA_OPERATOR_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace Tpetra {
  //
  // Forward declarations.  The "doxygen" bit simply tells Doxygen
  // (our automatic documentation generation system) to skip forward
  // declarations.
  //
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class Map;

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class MultiVector;
#endif // DOXYGEN_SHOULD_SKIP_THIS

  //! \brief Abstract interface for linear operators accepting Tpetra MultiVector objects.
  /*!  This class is templated on \c Scalar, \c LocalOrdinal, \c GlobalOrdinal and \c Node.
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal
     type, if omitted, defaults to the \c LocalOrdinal type. Node is by default of type Kokkos::DefaultNode::DefaultNodeType.

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
