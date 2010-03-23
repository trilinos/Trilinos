// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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

#ifndef THYRA_TPETRA_MULTIVECTOR_DECL_HPP
#define THYRA_TPETRA_MULTIVECTOR_DECL_HPP

#include "Thyra_SpmdMultiVectorBase.hpp"
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
class TpetraMultiVector : virtual public SpmdMultiVectorBase<Scalar> {
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

  /** @name Overridden public functions from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace() const;
  //@}

protected:

  /** @name Overridden protected functions from MultiVectorBase */
  //@{
  /** \brief . */
  RCP<const VectorBase<Scalar> > colImpl(Ordinal j) const;
  /** \brief . */
  RCP<VectorBase<Scalar> > nonconstColImpl(Ordinal j);
//  /** \brief . */
//  RCP<MultiVectorBase<Scalar> >
//  nonconstContigSubViewImpl( const Range1D& colRng );
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
  void getNonconstLocalDataImpl(
    const Ptr<ArrayRCP<Scalar> > &localValues, const Ptr<Ordinal> &leadingDim
    );
  /** \brief . */
  void getLocalDataImpl(
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
