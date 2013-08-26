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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov)
//
// ***********************************************************************
// @HEADER

#ifndef THYRA_TPETRA_MULTIVECTOR_DECL_HPP
#define THYRA_TPETRA_MULTIVECTOR_DECL_HPP

#include "Thyra_SpmdMultiVectorDefaultBase.hpp"
#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Tpetra_MultiVector.hpp"


namespace Thyra {


/** \brief Concrete implementation of Thyra::MultiVector in terms of
 * Tpetra::MultiVector.
 *
 * \todo Finish documentation!
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraMultiVector
  : virtual public SpmdMultiVectorDefaultBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /// Construct to uninitialized
  TpetraMultiVector();

  /** \brief Initialize.
   */
  void initialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
    const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
    );

  /** \brief Initialize.
   */
  void constInitialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
    const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
    );

  /** \brief Extract the underlying non-const Tpetra::MultiVector object.*/
  RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMultiVector();

  /** \brief Extract the underlying const Tpetra::MultiVector object.*/
  RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVector() const;

  //@}

  /** @name Overridden public functions form MultiVectorAdapterBase */
  //@{
  /** \brief . */
  RCP< const ScalarProdVectorSpaceBase<Scalar> >
  domainScalarProdVecSpc() const;
  //@}

protected:

  /** @name Overridden protected functions from MultiVectorBase */
  //@{
  /** \brief . */
  RCP<const VectorBase<Scalar> > colImpl(Ordinal j) const;
  /** \brief . */
  RCP<VectorBase<Scalar> > nonconstColImpl(Ordinal j);

  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  contigSubViewImpl( const Range1D& colRng ) const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  nonconstContigSubViewImpl( const Range1D& colRng );
//  /** \brief . */
//  RCP<const MultiVectorBase<Scalar> >
//  nonContigSubViewImpl( const ArrayView<const int> &cols ) const;
//  /** \brief . */
//  RCP<MultiVectorBase<Scalar> >
//  nonconstNonContigSubViewImpl( const ArrayView<const int> &cols );
  //@}

  /** @name Overridden protected functions from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpaceImpl() const;
  /** \brief . */
  void getNonconstLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    );
  /** \brief . */
  void getLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<const Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    ) const;
  //@}

private:

  // ///////////////////////////////////////
  // Private data members

  RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetraVectorSpace_;
  RCP<const ScalarProdVectorSpaceBase<Scalar> > domainSpace_;
  Teuchos::ConstNonconstObjectContainer<Tpetra::MultiVector<Scalar, LocalOrdinal,GlobalOrdinal,Node> >
  tpetraMultiVector_;

  // ////////////////////////////////////
  // Private member functions

  template<class TpetraMultiVector_t>
  void initializeImpl(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
    const RCP<TpetraMultiVector_t> &tpetraMultiVector
    );

};


/** \brief Nonmember constructor for TpetraMultiVector.
 *
 * \relates TpetraMultiVector.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
tpetraMultiVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
  )
{
  RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmv =
    Teuchos::rcp(new TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  tmv->initialize(tpetraVectorSpace, domainSpace, tpetraMultiVector);
  return tmv;
}


/** \brief Nonmember constructor for TpetraMultiVector.
 *
 * \relates TpetraMultiVector.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
constTpetraMultiVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<Scalar> > &domainSpace,
  const RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraMultiVector
  )
{
  RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmv =
    Teuchos::rcp(new TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  tmv->constInitialize(tpetraVectorSpace, domainSpace, tpetraMultiVector);
  return tmv;
}


} // end namespace Thyra


#endif // THYRA_TPETRA_MULTIVECTOR_DECL_HPP
