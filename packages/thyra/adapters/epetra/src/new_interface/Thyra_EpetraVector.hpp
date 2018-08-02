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

#ifndef THYRA_EPETRA_VECTOR_HPP
#define THYRA_EPETRA_VECTOR_HPP

// Not directly needed in this file, but this way we made the
// macro HAVE_THYRA_EPETRA_REFACTOR available to files that include
// this header. This way, they do not need to include the config.h
// header manually. That's nice, because in the future we may deprecate
// and then remove the old interface, making the config.h file pointless.
// If that happens, we may remove it, and at that point all files including
// it would have to be updated. This was, only the adapters headers need to
// be updated.
#include "ThyraEpetraAdapters_config.h"

#include "Thyra_SpmdVectorDefaultBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

class Epetra_Vector;
class Epetra_MultiVector;

namespace Thyra {

class EpetraVectorSpace;

/** \brief Concrete Thyra::SpmdVectorBase using Epetra_Vector.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraVector : virtual public SpmdVectorDefaultBase<double>
{
public:

  /** @name Constructors/initializers */
  //@{

  /** \brief Construct to uninitialized. */
  EpetraVector();

  /** \brief Initialize. */
  void initialize(
    const RCP<const VectorSpaceBase<double>>& epetraVectorSpace,
    const RCP<Epetra_Vector>& epetraVector
    );


  /** \brief Initialize. */
  void constInitialize(
    const RCP<const VectorSpaceBase<double>>& epetraVectorSpace,
    const RCP<const Epetra_Vector>& epetraVector
    );

  /** \brief Get the embedded non-const Epetra_Vector. */
  RCP<Epetra_Vector>
  getEpetraVector();

  /** \brief Get the embedded non-const Epetra_Vector. */
  RCP<const Epetra_Vector>
  getConstEpetraVector() const;

  //@}

  /** @name Overridden from VectorDefaultBase */
  //@{
  /** \brief . */
  RCP<const VectorSpaceBase<double>> domain() const;
  //@}

  // Should these Impl functions should alsp be protected???
//protected:

  /** @name Overridden from SpmdMultiVectorBase */
  //@{
  /** \brief . */
  RCP<const SpmdVectorSpaceBase<double>> spmdSpaceImpl() const;
  //@}

  /** @name Overridden from SpmdVectorBase */
  //@{
  /** \brief . */
  void getNonconstLocalVectorDataImpl(const Ptr<ArrayRCP<double>> &localValues);
  /** \brief . */
  void getLocalVectorDataImpl(const Ptr<ArrayRCP<const double>> &localValues) const;
  //@}

protected:

  /** @name Overridden protected functions from VectorBase */
  //@{

  /** \brief . */
  virtual void randomizeImpl(double l, double u);

  /** \brief . */
  virtual void absImpl(const VectorBase<double>& x);

  /** \brief . */
  virtual void reciprocalImpl(const VectorBase<double>& x);

  /** \brief . */
  virtual void eleWiseScaleImpl(const VectorBase<double>& x);

  /** \brief . */
  virtual typename Teuchos::ScalarTraits<double>::magnitudeType
  norm2WeightedImpl(const VectorBase<double>& x) const;

  /** \brief . */
  virtual void applyOpImpl(
    const RTOpPack::RTOpT<double> &op,
    const ArrayView<const Ptr<const VectorBase<double>>> &vecs,
    const ArrayView<const Ptr<VectorBase<double>>> &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const;

  /** \brief . */
  void acquireDetachedVectorViewImpl(
    const Range1D& rng,
    RTOpPack::ConstSubVectorView<double>* sub_vec
    ) const;

  /** \brief . */
  void acquireNonconstDetachedVectorViewImpl(
    const Range1D& rng,
    RTOpPack::SubVectorView<double>* sub_vec
    );

  /** \brief . */
  void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<double>* sub_vec
    );

  //@}

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
    const ArrayView<const Ptr<const MultiVectorBase<double>>>& mv,
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

  //@}

  /** @name Overridden protected functions from LinearOpBase */
  //@{

  /** \brief . */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &X,
    const Ptr<MultiVectorBase<double>> &Y,
    const double alpha,
    const double beta
    ) const;

  //@}

private:

  // ///////////////////////////////////////
  // Private data members
  
  RCP<const EpetraVectorSpace>
  epetraVectorSpace_;

  mutable RCP<const EpetraVectorSpace>
  domainSpace_;
  
  Teuchos::ConstNonconstObjectContainer<Epetra_Vector>
  epetraVector_;

  // ////////////////////////////////////
  // Private member functions

  template<class EpetraVector_t>
  void initializeImpl(
    const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
    const RCP<EpetraVector_t> &epetraVector
    );

  // Non-throwing Epetra Vector/MultiVector extraction methods.
  // Return null if casting failed.
  RCP<Epetra_MultiVector>
  getEpetraMultiVector(const RCP<MultiVectorBase<double>>& mv) const;

  RCP<const Epetra_MultiVector>
  getConstEpetraMultiVector(const RCP<const MultiVectorBase<double>>& mv) const;

  RCP<Epetra_Vector>
  getEpetraVector(const RCP<VectorBase<double>>& v) const;

  RCP<const Epetra_Vector>
  getConstEpetraVector(const RCP<const VectorBase<double>>& v) const;

};


/** \brief Nonmember constructor for EpetraVector.
 *
 * \relates EpetraVector
 */

inline
RCP<EpetraVector>
epetraVector(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<Epetra_Vector> &epetraVector
  )
{
  RCP<EpetraVector> v = Teuchos::rcp(new EpetraVector());
  v->initialize(epetraVectorSpace, epetraVector);
  return v;
}


/** \brief Nonmember constructor for EpetraVector.
 *
 * \relates EpetraVector
 */

inline
RCP<const EpetraVector>
constEpetraVector(
  const RCP<const VectorSpaceBase<double>> &epetraVectorSpace,
  const RCP<const Epetra_Vector> &epetraVector
  )
{
  RCP<EpetraVector> v = Teuchos::rcp(new EpetraVector());
  v->constInitialize(epetraVectorSpace, epetraVector);
  return v;
}


} // end namespace Thyra


#endif // THYRA_EPETRA_VECTOR_HPP
