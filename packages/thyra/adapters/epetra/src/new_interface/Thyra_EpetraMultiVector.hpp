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

#ifndef THYRA_EPETRA_MULTIVECTOR_DECL_HPP
#define THYRA_EPETRA_MULTIVECTOR_DECL_HPP

// Not directly needed in this file, but this way we made the
// macro HAVE_THYRA_EPETRA_REFACTOR available to files that include
// this header. This way, they do not need to include the config.h
// header manually. That's nice, because in the future we may deprecate
// and then remove the old interface, making the config.h file pointless.
// If that happens, we may remove it, and at that point all files including
// it would have to be updated. This was, only the adapters headers need to
// be updated.
#include "ThyraEpetraAdapters_config.h"

#include "Thyra_SpmdMultiVectorDefaultBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

class Epetra_MultiVector;

namespace Thyra {

class EpetraVectorSpace;

/** \brief Concrete implementation of Thyra::MultiVector in terms of
 * Epetra_MultiVector.
 *
 * \todo Finish documentation!
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraMultiVector : virtual public SpmdMultiVectorDefaultBase<double>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /// Construct to uninitialized
  EpetraMultiVector();

  /** \brief Initialize.
   */
  void initialize(
    const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
    const RCP<Epetra_MultiVector> &epetraMultiVector
    );

  /** \brief Initialize.
   */
  void constInitialize(
    const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
    const RCP<const Epetra_MultiVector> &epetraMultiVector
    );

  /** \brief Extract the underlying non-const Epetra_MultiVector object.*/
  RCP<Epetra_MultiVector> getEpetraMultiVector();

  /** \brief Extract the underlying const Epetra_MultiVector object.*/
  RCP<const Epetra_MultiVector> getConstEpetraMultiVector() const;

  //@}

  /** @name Overridden public functions form MultiVectorAdapterBase */
  //@{
  /** \brief . */
  RCP< const ScalarProdVectorSpaceBase<double> >
  domainScalarProdVecSpc() const;
  //@}

protected:

  /** @name Overridden protected functions from MultiVectorBase */
  //@{
  /** \brief . */
  virtual void assignImpl(double alpha);

  /** \brief . */
  virtual void assignMultiVecImpl(const MultiVectorBase<double>& mv);

  /** \brief . */
  virtual void scaleImpl(double alpha);

  /** \brief . */
  virtual void updateImpl(
    double alpha,
    const MultiVectorBase<double>& mv
    );

  /** \brief . */
  virtual void linearCombinationImpl(
    const ArrayView<const double>& alpha,
    const ArrayView<const Ptr<const MultiVectorBase<double> > >& mv,
    const double& beta
    );

  /** \brief . */
  virtual void dotsImpl(
    const MultiVectorBase<double>& mv,
    const ArrayView<double>& prods
    ) const;

  /** \brief . */
  virtual void norms1Impl(
    const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
    ) const;

  /** \brief . */
  virtual void norms2Impl(
    const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
    ) const;

  /** \brief . */
  virtual void normsInfImpl(
    const ArrayView<typename ScalarTraits<double>::magnitudeType>& norms
    ) const;

  /** \brief . */
  RCP<const VectorBase<double> > colImpl(Ordinal j) const;
  /** \brief . */
  RCP<VectorBase<double> > nonconstColImpl(Ordinal j);

  /** \brief . */
  RCP<const MultiVectorBase<double> >
  contigSubViewImpl(const Range1D& colRng) const;
  /** \brief . */
  RCP<MultiVectorBase<double> >
  nonconstContigSubViewImpl(const Range1D& colRng);
  /** \brief . */
  RCP<const MultiVectorBase<double> >
  nonContigSubViewImpl(const ArrayView<const int>& cols_in) const;
  /** \brief . */
  RCP<MultiVectorBase<double> >
  nonconstNonContigSubViewImpl(const ArrayView<const int>& cols_in);

  /** \brief . */
  virtual void mvMultiReductApplyOpImpl(
    const RTOpPack::RTOpT<double> &primary_op,
    const ArrayView<const Ptr<const MultiVectorBase<double> > > &multi_vecs,
    const ArrayView<const Ptr<MultiVectorBase<double> > > &targ_multi_vecs,
    const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
    const Ordinal primary_global_offset
    ) const;

  /** \brief . */
  void acquireDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::ConstSubMultiVectorView<double>* sub_mv
    ) const;

  /** \brief . */
  void acquireNonconstDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<double>* sub_mv
    );

  /** \brief . */
  void commitNonconstDetachedMultiVectorViewImpl(
    RTOpPack::SubMultiVectorView<double>* sub_mv
    );

//  /** \brief . */
//  RCP<const MultiVectorBase<double> >
//  nonContigSubViewImpl( const ArrayView<const int> &cols ) const;
//  /** \brief . */
//  RCP<MultiVectorBase<double> >
//  nonconstNonContigSubViewImpl( const ArrayView<const int> &cols );
  //@}

  /** @name Overridden protected functions from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RCP<const SpmdVectorSpaceBase<double> > spmdSpaceImpl() const;
  /** \brief . */
  void getNonconstLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<double> > &localValues, const Ptr<Ordinal> &leadingDim
    );
  /** \brief . */
  void getLocalMultiVectorDataImpl(
    const Ptr<ArrayRCP<const double> > &localValues, const Ptr<Ordinal> &leadingDim
    ) const;

  //@}

  /** @name Overridden protected functions from MultiVectorAdapterBase */
  //@{
  /** \brief . */
  virtual void euclideanApply(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &X,
    const Ptr<MultiVectorBase<double> > &Y,
    const double alpha,
    const double beta
    ) const;

  //@}

private:

  // ///////////////////////////////////////
  // Private data members

  RCP<const EpetraVectorSpace> epetraVectorSpace_;
  RCP<const ScalarProdVectorSpaceBase<double> > domainSpace_;
  Teuchos::ConstNonconstObjectContainer<Epetra_MultiVector> epetraMultiVector_;

  // ////////////////////////////////////
  // Private member functions

  template<class EpetraMultiVector_t>
  void initializeImpl(
    const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
    const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
    const RCP<EpetraMultiVector_t> &epetraMultiVector
    );

  // Non-throwing Epetra MultiVector extraction methods.
  // Return null if casting failed.
  RCP<Epetra_MultiVector>
  getEpetraMultiVector(const RCP<MultiVectorBase<double> >& mv) const;

  RCP<const Epetra_MultiVector>
  getConstEpetraMultiVector(const RCP<const MultiVectorBase<double> >& mv) const;

};


/** \brief Nonmember constructor for EpetraMultiVector.
 *
 * \relates EpetraMultiVector.
 */
inline
RCP<EpetraMultiVector>
epetraMultiVector(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
  const RCP<Epetra_MultiVector> &epetraMultiVector
  )
{
  RCP<EpetraMultiVector> emv = Teuchos::rcp(new EpetraMultiVector());
  emv->initialize(epetraVectorSpace, domainSpace, epetraMultiVector);
  return emv;
}

/** \brief Nonmember constructor for EpetraMultiVector.
 *
 * \relates EpetraMultiVector.
 */
inline
RCP<const EpetraMultiVector>
constEpetraMultiVector(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<const ScalarProdVectorSpaceBase<double> > &domainSpace,
  const RCP<const Epetra_MultiVector> &epetraMultiVector
  )
{
  RCP<EpetraMultiVector> emv = Teuchos::rcp(new EpetraMultiVector());
  emv->constInitialize(epetraVectorSpace, domainSpace, epetraMultiVector);
  return emv;
}


} // end namespace Thyra


#endif // THYRA_EPETRA_MULTIVECTOR_DECL_HPP
