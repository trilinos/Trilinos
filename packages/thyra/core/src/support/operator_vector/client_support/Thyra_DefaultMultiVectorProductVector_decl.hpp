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

#ifndef THYRA_MULTI_VECTOR_PRODUCT_VECTOR_DECL_HPP
#define THYRA_MULTI_VECTOR_PRODUCT_VECTOR_DECL_HPP


#include "Thyra_ProductVectorBase.hpp" // Interface
#include "Thyra_VectorDefaultBase.hpp" // Implementation
#include "Thyra_DefaultProductVector.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


template<class Scalar> class DefaultMultiVectorProductVectorSpace;


/** \brief Concrete implementation of a product vector which is really
 * composed out of the columns of a multi-vector.
 *
 * Note that clients should almost never be creating objects of this
 * type explicitly and should instead use <tt>DefaultMultiVectorProductVectorSpace</tt>
 * as a factory.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultMultiVectorProductVector
  : virtual public ProductVectorBase<Scalar>,
    virtual protected VectorDefaultBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  DefaultMultiVectorProductVector();

  /** \brief Initialize with a non-const multi-vector. */
  void initialize(
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace,
    const RCP<MultiVectorBase<Scalar> > &multiVec
    );

  /** \brief Initialize with a const multi-vector. */
  void initialize(
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace,
    const RCP<const MultiVectorBase<Scalar> > &multiVec
    );

  // ToDo: Add const version of above function also when needed!

  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  getNonconstMultiVector();

  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  getMultiVector() const;

  /** \brief . */
  void uninitialize();

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

  /** @name Overridden from ProductVectorBase */
  //@{

  /** \brief . */
  RCP<VectorBase<Scalar> >
  getNonconstVectorBlock(const int k); 
  /** \brief . */
  RCP<const VectorBase<Scalar> >
  getVectorBlock(const int k) const;

  //@}

  /** @name Overridden from ProductMultiVectorBase */
  //@{

  /** \brief . */
  RCP<const ProductVectorSpaceBase<Scalar> >
  productSpace() const;
  /** \brief . */
  bool blockIsConst(const int k) const; 
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  getNonconstMultiVectorBlock(const int k);
  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  getMultiVectorBlock(const int k) const;

  //@}

  /** @name Overridden public functions from VectorBase */
  //@{

  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > space() const;

  //@}

protected:

  /** @name Overridden protected functions from VectorBase */
  //@{

  /** \brief . */
  void applyOpImpl(
    const RTOpPack::RTOpT<Scalar> &op,
    const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
    const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
    const Ptr<RTOpPack::ReductTarget> &reduct_obj,
    const Ordinal global_offset
    ) const;
  /** \brief . */
  void acquireDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief . */
  void releaseDetachedVectorViewImpl(
    RTOpPack::ConstSubVectorView<Scalar>* sub_vec
    ) const;
  /** \brief . */
  void acquireNonconstDetachedVectorViewImpl(
    const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec
    );
  /** \brief . */
  void commitNonconstDetachedVectorViewImpl(
    RTOpPack::SubVectorView<Scalar>* sub_vec
    );
  /** \brief . */
  void setSubVectorImpl(
    const RTOpPack::SparseSubVectorT<Scalar>& sub_vec
    );

  //@}

private:

  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<MultiVectorBase<Scalar> > CNMVC;

  // //////////////////////////////
  // Private data members

  int numBlocks_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > productSpace_;
  CNMVC multiVec_;

  // //////////////////////////////
  // Private member functions

  RCP<const DefaultProductVector<Scalar> >
  getDefaultProductVector() const;

};


/** \brief Nonmember constructor that just wraps an existing non-const
 * MultiVector as a non-const product vector.
 *
 * \relates DefaultMultiVectorProductVector
 */
template<class Scalar>
inline
RCP<DefaultMultiVectorProductVector<Scalar> >
multiVectorProductVector(
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace,
  const RCP<MultiVectorBase<Scalar> > &multiVec
  )
{
  RCP<DefaultMultiVectorProductVector<Scalar> > multiVecProdVec
    = Teuchos::rcp(new DefaultMultiVectorProductVector<Scalar>());
  multiVecProdVec->initialize(productSpace,multiVec);
  return multiVecProdVec;
}


/** \brief Nonmember constructor that just wraps an existing const MultiVector
 * as a const product vector.
 *
 * \relates DefaultMultiVectorProductVector
 */
template<class Scalar>
inline
RCP<const DefaultMultiVectorProductVector<Scalar> >
multiVectorProductVector(
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace,
  const RCP<const MultiVectorBase<Scalar> > &multiVec
  )
{
  RCP<DefaultMultiVectorProductVector<Scalar> > multiVecProdVec
    = Teuchos::rcp(new DefaultMultiVectorProductVector<Scalar>());
  multiVecProdVec->initialize(productSpace,multiVec);
  return multiVecProdVec;
}


// ToDo: Add non-const and const versions of the nonmember constructor
// functions to wrap already created multi-vectors once needed!


/** \brief Nonmember constructor that creates a new uninitialized product
 * vector represented underneath as a multi-vector.
 *
 * \relates DefaultMultiVectorProductVector
 */
template<class Scalar>
inline
RCP<DefaultMultiVectorProductVector<Scalar> >
multiVectorProductVector(
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &productSpace
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(productSpace));
#endif
  return multiVectorProductVector(
    productSpace,
    createMembers(productSpace->getBlock(0),productSpace->numBlocks())
    );
}


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_PRODUCT_VECTOR_DECL_HPP
