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

#ifndef TPETRA_MIXEDSCALARMULTIPLYOP_HPP
#define TPETRA_MIXEDSCALARMULTIPLYOP_HPP

/// \file Tpetra_MixedScalarMultiplyOp.hpp
///
/// Declaration and definition of Tpetra::MixedScalarMultiplyOp and its
/// nonmember constructor Tpetra::createMixedScalarMatrixMultiplyOp.

#include "Tpetra_Operator.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_Profiling.hpp"

namespace Tpetra {

  /// \brief A class for wrapping an Operator of one Scalar type into
  ///        an Operator of another Scalar type.
  ///
  /// \tparam Scalar The type of the entries of the input and output
  ///   MultiVector (see apply()).  Same as the first template
  ///   parameter of Operator.
  ///
  /// \tparam OpScalar The type of the entries of the wrapped Operator.
  ///
  /// \tparam LocalOrdinal The second template parameter of Operator.
  ///
  /// \tparam GlobalOrdinal The third template parameter of Operator.
  ///
  /// \tparam Node The fourth template parameter of Operator.
  template <class Scalar,
            class OpScalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class MixedScalarMultiplyOp :
    public Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  public:
    //! The specialization of CrsMatrix which this class wraps.
    using op_type =
      Operator<OpScalar, LocalOrdinal, GlobalOrdinal, Node>;
    //! The specialization of Map which this class uses.
    using map_type = Map<LocalOrdinal, GlobalOrdinal, Node>;

  public:
    //! @name Constructor and destructor
    //@{

    /// \brief Constructor
    ///
    /// \param A [in] The Operator to wrap as an
    ///   <tt>Operator<Scalar, ...></tt>.
    MixedScalarMultiplyOp (const Teuchos::RCP<const op_type>& A) :
      op_ (A)
    {
      Allocate(1);
    }

    //! Destructor (virtual for memory safety of derived classes).
    ~MixedScalarMultiplyOp () override = default;

    //@}
    //! @name Methods implementing Operator
    //@{

    /// \brief Compute <tt>Y = beta*Y + alpha*Op(A)*X</tt>, where
    ///   <tt>Op(A)</tt> is either A, \f$A^T\f$, or \f$A^H\f$.
    void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
           MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one (),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero ()) const override
    {
      TEUCHOS_TEST_FOR_EXCEPTION
        (X.getNumVectors () != Y.getNumVectors (), std::runtime_error,
         Teuchos::typeName (*this) << "::apply(X,Y): X and Y must have the same number of vectors.");

      if (X.getNumVectors() != X_->getNumVectors()) {
        Allocate(X.getNumVectors());
      }

      Tpetra::deep_copy(*X_, X);
      op_->apply (*X_, *Y_, mode, alpha, beta);
      Tpetra::deep_copy(Y, *Y_);
    }

    /// \brief Whether this Operator's apply() method can apply the
    ///   transpose or conjugate transpose.
    ///
    /// This is always true, since it is true for the CrsMatrix that
    /// this object wraps.
    bool hasTransposeApply() const override {
      return true;
    }

    //! The domain Map of this Operator.
    Teuchos::RCP<const map_type> getDomainMap () const override {
      return op_->getDomainMap ();
    }

    //! The range Map of this Operator.
    Teuchos::RCP<const map_type> getRangeMap () const override {
      return op_->getRangeMap ();
    }

    //@}

  protected:
    typedef MultiVector<OpScalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    //! The underlying CrsMatrix object.
    const Teuchos::RCP<const op_type> op_;
    mutable Teuchos::RCP<MV> X_, Y_;

    void
    Allocate(int numVecs) const {
      X_ = rcp(new MV(op_->getDomainMap(),numVecs));
      Y_ = rcp(new MV(op_->getRangeMap(),numVecs));
    }

  };

  /// \brief Non-member function to create a MixedScalarMultiplyOp.
  /// \relatesalso MixedScalarMultiplyOp
  ///
  /// The function has the same template parameters of MixedScalarMultiplyOp.
  ///
  /// \param A [in] The Operator instance to wrap in an MixedScalarMultiplyOp.
  /// \return The MixedScalarMultiplyOp wrapper for the given Operator.
  template <class Scalar,
            class OpScalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  Teuchos::RCP<
    MixedScalarMultiplyOp<Scalar, OpScalar, LocalOrdinal, GlobalOrdinal, Node> >
  createMixedScalarMultiplyOp (const Teuchos::RCP<
    const Operator<OpScalar, LocalOrdinal, GlobalOrdinal, Node> >& A)
  {
    typedef MixedScalarMultiplyOp<Scalar, OpScalar, LocalOrdinal,
      GlobalOrdinal, Node> op_type;
    return Teuchos::rcp (new op_type (A));
  }

} // end of namespace Tpetra

#endif // TPETRA_MIXEDSCALARMULTIPLYOP_HPP
