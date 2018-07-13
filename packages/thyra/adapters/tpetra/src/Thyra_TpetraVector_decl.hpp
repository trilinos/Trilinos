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

#ifndef THYRA_TPETRA_VECTOR_DECL_HPP
#define THYRA_TPETRA_VECTOR_DECL_HPP


#include "Thyra_SpmdVectorDefaultBase.hpp"
#include "Thyra_TpetraVectorSpace_decl.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Concrete Thyra::SpmdVectorBase using Tpetra::Vector.
 *
 * \ingroup Tpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraVector
  : virtual public SpmdVectorDefaultBase<Scalar>
{
public:

  /** @name Constructors/initializers */
  //@{

  /** \brief Construct to uninitialized. */
  TpetraVector();

  /** \brief Initialize. */
  void initialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
    );

  /** \brief Initialize. */
  void constInitialize(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
    );

  /** \brief Get the embedded non-const Tpetra::Vector. */
  RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraVector();

  /** \brief Get the embedded non-const Tpetra::Vector. */
  RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraVector() const;

  //@}

  /** @name Overridden from VectorDefaultBase */
  //@{
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > domain() const;
  //@}

  // Should these Impl functions should alsp be protected???
//protected:

  /** @name Overridden from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpaceImpl() const;
  //@}

  /** @name Overridden from SpmdVectorBase */
  //@{
  /** \brief . */
  void getNonconstLocalVectorDataImpl(const Ptr<ArrayRCP<Scalar> > &localValues);
  /** \brief . */
  void getLocalVectorDataImpl(const Ptr<ArrayRCP<const Scalar> > &localValues) const;
  //@}

protected:

  /** @name Overridden protected functions from VectorBase */
  //@{

  /** \brief . */
  virtual void randomizeImpl(Scalar l, Scalar u);

  /** \brief . */
  virtual void absImpl(const VectorBase<Scalar>& x);

  /** \brief . */
  virtual void reciprocalImpl(const VectorBase<Scalar>& x);

  /** \brief . */
  virtual void eleWiseScaleImpl(const VectorBase<Scalar>& x);

  /** \brief . */
  virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  norm2WeightedImpl(const VectorBase<Scalar>& x) const;

  /** \brief . */
  virtual void applyOpImpl(
    const RTOpPack::RTOpT<Scalar> &op,
    const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
    const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const;

  /** \brief . */
  void acquireDetachedVectorViewImpl(
    const Range1D& rng,
    RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;

  /** \brief . */
  void acquireNonconstDetachedVectorViewImpl(
    const Range1D& rng,
    RTOpPack::SubVectorView<Scalar>* sub_vec
    );

  /** \brief . */
  void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<Scalar>* sub_vec
    );

  //@}

  /** @name Overridden protected functions from MultiVectorBase */
  //@{
  /** \brief . */
  virtual void assignImpl(Scalar alpha);

  /** \brief . */
  virtual void assignMultiVecImpl(const MultiVectorBase<Scalar>& mv);

  /** \brief . */
  virtual void scaleImpl(Scalar alpha);

  /** \brief . */
  virtual void updateImpl(
    Scalar alpha,
    const MultiVectorBase<Scalar>& mv
    );

  /** \brief . */
  virtual void linearCombinationImpl(
    const ArrayView<const Scalar>& alpha,
    const ArrayView<const Ptr<const MultiVectorBase<Scalar> > >& mv,
    const Scalar& beta
    );

  /** \brief . */
  virtual void dotsImpl(
    const MultiVectorBase<Scalar>& mv,
    const ArrayView<Scalar>& prods
    ) const;

  /** \brief . */
  virtual void norms1Impl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const;

  /** \brief . */
  virtual void norms2Impl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const;

  /** \brief . */
  virtual void normsInfImpl(
    const ArrayView<typename ScalarTraits<Scalar>::magnitudeType>& norms
    ) const;

  //@}

  /** @name Overridden protected functions from LinearOpBase */
  //@{

  /** \brief . */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:

  // ///////////////////////////////////////
  // Private data members
  
  RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  tpetraVectorSpace_;

  mutable RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  domainSpace_;
  
  Teuchos::ConstNonconstObjectContainer<Tpetra::Vector<Scalar, LocalOrdinal,GlobalOrdinal,Node> >
  tpetraVector_;

  typedef Tpetra::MultiVector<Scalar, LocalOrdinal,GlobalOrdinal,Node> TpetraMultiVector_t;

  // ////////////////////////////////////
  // Private member functions

  template<class TpetraVector_t>
  void initializeImpl(
    const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
    const RCP<TpetraVector_t> &tpetraVector
    );

  // Non-throwing Tpetra Vector/MultiVector extraction methods.
  // Return null if casting failed.
  RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraMultiVector(const RCP<MultiVectorBase<Scalar> >& mv) const;

  RCP<const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraMultiVector(const RCP<const MultiVectorBase<Scalar> >& mv) const;

  RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getTpetraVector(const RCP<VectorBase<Scalar> >& v) const;

  RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getConstTpetraVector(const RCP<const VectorBase<Scalar> >& v) const;

};


/** \brief Nonmember constructor for TpetraVector.
 *
 * \relates TpetraVector
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
inline
RCP<TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
tpetraVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
  )
{
  RCP<TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > v =
    Teuchos::rcp(new TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  v->initialize(tpetraVectorSpace, tpetraVector);
  return v;
}


/** \brief Nonmember constructor for TpetraVector.
 *
 * \relates TpetraVector
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
inline
RCP<const TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
constTpetraVector(
  const RCP<const TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVectorSpace,
  const RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &tpetraVector
  )
{
  RCP<TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > v =
    Teuchos::rcp(new TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>);
  v->constInitialize(tpetraVectorSpace, tpetraVector);
  return v;
}


} // end namespace Thyra


#endif // THYRA_TPETRA_VECTOR_DECL_HPP
