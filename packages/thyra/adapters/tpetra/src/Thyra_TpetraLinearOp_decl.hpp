// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_TPETRA_LINEAR_OP_DECL_HPP
#define THYRA_TPETRA_LINEAR_OP_DECL_HPP

#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Tpetra_Operator.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

#if defined(HAVE_THYRA_EPETRA) && defined(HAVE_TPETRA_EPETRA)
#  define HAVE_THYRA_TPETRA_EPETRA
#endif

#ifdef HAVE_THYRA_TPETRA_EPETRA
#  include "Thyra_EpetraLinearOpBase.hpp"
#  include "Tpetra_EpetraRowMatrix.hpp"
#endif


namespace Thyra {


/** \brief Concrete Thyra::LinearOpBase subclass for Tpetra::Operator.
 *
 * \todo Finish Documentation
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal=LocalOrdinal,
  class Node=Kokkos::DefaultNode::DefaultNodeType>
class TpetraLinearOp
  : virtual public Thyra::LinearOpDefaultBase<Scalar>
#ifdef HAVE_THYRA_TPETRA_EPETRA
  , virtual public EpetraLinearOpBase
#endif
{
public:

  /** \name Constructors/initializers. */
  //@{

  /** \brief Construct to uninitialized. */
  TpetraLinearOp();

  /** \brief Initialize. */
  void initialize(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator
    );

  /** \brief Initialize. */
  void constInitialize(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator
    );

  /** \brief Get embedded non-const Tpetra::Operator. */
  RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraOperator();

  /** \brief Get embedded const Tpetra::Operator. */
  RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraOperator() const;

  //@}

  /** \name Public Overridden functions from LinearOpBase. */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > range() const;

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > domain() const;

  //@}

#ifdef HAVE_THYRA_TPETRA_EPETRA

  /** \name Overridden from EpetraLinearOpBase */
  //@{

  /** \brief . */
  void getNonconstEpetraOpView(
    const Ptr<RCP<Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    );
  /** \brief . */
  void getEpetraOpView(
    const Ptr<RCP<const Epetra_Operator> > &epetraOp,
    const Ptr<EOpTransp> &epetraOpTransp,
    const Ptr<EApplyEpetraOpAs> &epetraOpApplyAs,
    const Ptr<EAdjointEpetraOp> &epetraOpAdjointSupport
    ) const;

  //@}

#endif // HAVE_THYRA_TPETRA_EPETRA

protected:

  /** \name Protected Overridden functions from LinearOpBase. */
  //@{

  /** \brief . */
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const;

  /** \brief . */
  void applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar> &X_in,
    const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:

  RCP<const VectorSpaceBase<Scalar> >
  rangeSpace_;

  RCP<const VectorSpaceBase<Scalar> >
  domainSpace_;

  Teuchos::ConstNonconstObjectContainer<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  tpetraOperator_;

#ifdef HAVE_THYRA_TPETRA_EPETRA
  mutable RCP<Epetra_Operator> epetraOp_;
#endif

  template<class TpetraOperator_t>
  void initializeImpl(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<TpetraOperator_t> &tpetraOperator
  );

};


/** \brief Nonmmeber constructor for TpetraLinearOp.
 *
 * \relates TpetraLinearOp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
tpetraLinearOp(
  const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
  const RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator
  )
{
  const RCP<TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op =
    Teuchos::rcp(new TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
  op->initialize(rangeSpace, domainSpace, tpetraOperator);
  return op;
}


/** \brief Nonmmeber constructor for TpetraLinearOp.
 *
 * \relates TpetraLinearOp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
constTpetraLinearOp(
  const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
  const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
  const RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraOperator
  )
{
  const RCP<TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op =
    Teuchos::rcp(new TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
  op->constInitialize(rangeSpace, domainSpace, tpetraOperator);
  return op;
}


}  // namespace Thyra


#endif	// THYRA_TPETRA_LINEAR_OP_DECL_HPP
