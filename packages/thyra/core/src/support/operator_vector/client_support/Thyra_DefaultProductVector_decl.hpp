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

#ifndef THYRA_DEFAULT_PRODUCT_VECTOR_DECL_HPP
#define THYRA_DEFAULT_PRODUCT_VECTOR_DECL_HPP

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_VectorDefaultBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_as.hpp"


namespace Thyra {


/** \brief . */
template <class Scalar> class DefaultProductVectorSpace;


/** \brief Concrete implementation of a product vector.
 *
 * Note that clients should almost never be creating objects of this
 * type explicitly and should instead use <tt>DefaultProductVectorSpace</tt>
 * as a factory.
 *
 * ToDo: Finish documentation!
 *
 * The default constructor is made private to avoid accidental default
 * construction.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultProductVector
  : virtual public ProductVectorBase<Scalar>,
    virtual protected VectorDefaultBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  DefaultProductVector();

  /** \brief Constructs to initialized (calls <tt>initialize()</tt>). */
  DefaultProductVector(
    const RCP<const DefaultProductVectorSpace<Scalar> >  &productSpace
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const ArrayView<const RCP<VectorBase<Scalar> > > &vecs
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const ArrayView<const RCP<const VectorBase<Scalar> > > &vecs
    );

  /** \brief Uninitialize.
   *
   * ToDo: Finish documentation.
   */
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

  /** @name Extensions to ProductVectorBase suitable for physically-blocked vectors */
  //@{

  /** \brief . */
  void setBlock(int i, const RCP<const VectorBase<Scalar> >& b); 
  /** \brief . */
  void setNonconstBlock(int i, const RCP<VectorBase<Scalar> >& b); 
  //@}

  /** @name Overridden from ProductVectorBase */
  //@{

  /** \brief . */
  RCP<VectorBase<Scalar> > getNonconstVectorBlock(const int k); 
  /** \brief . */
  RCP<const VectorBase<Scalar> > getVectorBlock(const int k) const;

  //@}

  /** @name Overridden public functions from ProductMultiVectorBase */
  //@{

  /** \brief . */
  RCP<const ProductVectorSpaceBase<Scalar> > productSpace() const;
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

public:

#ifndef THYRA_HIDE_DEPRECATED_CODE
  /** \name Deprecated */
  //@{

  /** \brief Deprecated. */
  THYRA_DEPRECATED DefaultProductVector(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
    const RCP<VectorBase<Scalar> > vecs[]
    )
    :numBlocks_(0)
    { initialize(productSpace_in, Teuchos::arrayView(vecs, productSpace_in->numBlocks())); }

  /** \brief Deprecated. */
  THYRA_DEPRECATED void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
    const RCP<VectorBase<Scalar> > vecs[]
    )
    { initialize(productSpace_in, Teuchos::arrayView(vecs, productSpace_in->numBlocks())); }

  /** \brief Deprecated. */
  THYRA_DEPRECATED void initialize(
    const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace_in,
    const RCP<const VectorBase<Scalar> > vecs[]
    )
    { initialize(productSpace_in, Teuchos::arrayView(vecs, productSpace_in->numBlocks())); }

  //@}

#endif // THYRA_HIDE_DEPRECATED_CODE
private:

  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<VectorBase<Scalar> > CNVC;

  // //////////////////////////////
  // Private data members

  RCP<const DefaultProductVectorSpace<Scalar> > productSpace_;
  Array<CNVC> vecs_;
  // cache
  int numBlocks_;

};


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
inline
RCP<DefaultProductVector<Scalar> >
defaultProductVector()
{
  return Teuchos::rcp(new DefaultProductVector<Scalar>);
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
inline
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace
  )
{
  return Teuchos::rcp(
    new DefaultProductVector<Scalar>(productSpace)
    );
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const ArrayView<const RCP<VectorBase<Scalar> > > &vecs
  )
{
  RCP<DefaultProductVector<Scalar> > pv = defaultProductVector<Scalar>();
  pv->initialize(productSpace, vecs);
  return pv;
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const ArrayView<const RCP<const VectorBase<Scalar> > > &vecs
  )
{
  RCP<DefaultProductVector<Scalar> > pv = defaultProductVector<Scalar>();
  pv->initialize(productSpace, vecs);
  return pv;
}


/** \brief Return a casted non-const ProductVectorBase object or create a new
 * DefaultProductVector object with one component.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<ProductVectorBase<Scalar> >
castOrCreateNonconstProductVectorBase(const RCP<VectorBase<Scalar> > v);


/** \brief Return a casted const ProductVectorBase object or create a new
 * DefaultProductVector object with one component.
 *
 * \relates DefaultProductVector
 */
template<class Scalar>
RCP<const ProductVectorBase<Scalar> >
castOrCreateProductVectorBase(const RCP<const VectorBase<Scalar> > v);

#ifndef THYRA_HIDE_DEPRECATED_CODE
/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const RCP<VectorBase<Scalar> > vecs[]
  )
{
  return defaultProductVector<Scalar>(
    productSpace, Teuchos::arrayView(vecs, productSpace->numBlocks()));
}


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED
RCP<DefaultProductVector<Scalar> >
defaultProductVector(
  const RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const RCP<const VectorBase<Scalar> > vecs[]
  )
{
  return defaultProductVector<Scalar>(productSpace,
    Teuchos::arrayView(vecs, productSpace->numBlocks()));
}
#endif // THYRA_HIDE_DEPRECATED_CODE

} // namespace Thyra


#endif // THYRA_DEFAULT_PRODUCT_VECTOR_DECL_HPP
