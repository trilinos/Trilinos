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

#ifndef TPETRA_APPLYOP_HPP
#define TPETRA_APPLYOP_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"


/*! \file Tpetra_ApplyOp_decl.hpp

    The implementations for the class Tpetra::ApplyOp and related non-member constructors.
 */

namespace Tpetra {

  namespace details {

  /// \brief A class for wrapping an Operator apply in a Operator.
  ///
  /// \note Most Tpetra users do not need to use this class.  It will
  ///   be useful to Tpetra users who want to do mixed-precision
  ///   operator apply, where the operator's data
  ///   has a different precision than that of the input and
  ///   output vectors.  If your operator and vectors have the
  ///   same type of entries, then you don't need to use this class.
  ///
  /// This class makes a <tt>Operator<OpScalar, ...></tt> "look
  /// like" an <tt>Operator<Scalar, ...></tt>, where
  /// <tt>OpScalar</tt> and <tt>Scalar</tt> may be different types.
  /// It does so by working around a limitation of C++, namely that
  /// template methods of a class can't be virtual.
  ///

  ///
  /// \tparam Scalar The type of the entries of the input and output
  ///   MultiVector of the apply() method.  Same as the first template
  ///   parameter of Operator.
  ///
  /// \tparam OperatorType The type of the underlying Operator, whose
  ///   first template parameter OperatorType::scalar_type may be different
  ///   than this Operator's Scalar type.
  ///
  template <class Scalar,
	    class OperatorType>
  class ApplyOp :
    public Tpetra::Operator<Scalar,typename OperatorType::local_ordinal_type,typename OperatorType::global_ordinal_type,typename OperatorType::node_type>
  {
  public:
  // \name Typedefs
  //@{

  //! The type of the entries of the input OperatorType.
  typedef typename OperatorType::scalar_type scalar_type;

  //! The type of local indices in the input OperatorType.
  typedef typename OperatorType::local_ordinal_type local_ordinal_type;

  //! The type of global indices in the input OperatorType.
  typedef typename OperatorType::global_ordinal_type global_ordinal_type;

  //! The type of the Kokkos Node used by the input OperatorType.
  typedef typename OperatorType::node_type node_type;

  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //@}

    //! @name Constructor and destructor
    //@{

    /// \brief Constructor
    ///
    /// \param A [in] The Operator to wrap as an
    ///   <tt>Operator<Scalar, ...></tt>.
    ApplyOp(const Teuchos::RCP<const OperatorType> &op)
      : operator_(op) 
    {
    }

    //! Destructor
    virtual ~ApplyOp()
    {
    }

    //@}
    //! @name Methods implementing Operator
    //@{

    /// \brief Compute <tt>Y = beta*Y + alpha*Op(A)*X</tt>, where
    ///   <tt>Op(A)</tt> is either A, \f$A^T\f$, or \f$A^H\f$.
    ///
    /// This method calls the underlying Operator object's
    /// applyTempl<Scalar,Scalar>() method.
    void
    apply (const Tpetra::MultiVector<Scalar,local_ordinal_type,global_ordinal_type,node_type>& X,
	   Tpetra::MultiVector<Scalar,local_ordinal_type,global_ordinal_type,node_type>& Y,
	   Teuchos::ETransp mode = Teuchos::NO_TRANS,
	   Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
	   Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const
    {
      operator_->template applyTempl<Scalar,Scalar>(X, Y, mode, alpha, beta);
    }


    /// \brief Whether this Operator's apply() method can apply the
    ///   transpose or conjugate transpose.
    ///
    /// This is always true, since it is true for the Operator that
    /// this object wraps.
    bool hasTransposeApply() const
    {
      return operator_->hasTransposeApply();
    }

    //! The domain Map of this Operator.
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > & getDomainMap() const
    {
      return operator_->getDomainMap();
    }

    //! The range Map of this Operator.
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > & getRangeMap() const
    {
      return operator_->getRangeMap();
    }

    //@}

  protected:

    //! The underlying Operator object.
    const Teuchos::RCP<const OperatorType> operator_;

  };

} // end of namespace details

} // end of namespace Tpetra

#endif // TPETRA_APPLYOP_HPP
