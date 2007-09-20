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

#ifndef THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DECL_HPP
#define THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DECL_HPP

#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_MultiVectorDefaultBase.hpp"
#include "Thyra_SingleScalarLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief . */
template <class Scalar> class DefaultProductVectorSpace;


/** \brief Concrete implementation of a product multi-vector.
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
class DefaultProductMultiVector
  : virtual public ProductMultiVectorBase<Scalar>,
    virtual protected MultiVectorDefaultBase<Scalar>,
    virtual protected SingleScalarLinearOpBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief. Constructs to initialized (calls <tt>initialize()</tt>). */
  DefaultProductMultiVector(
    const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const int numMembers
    );

  /** \brief. Constructs to initialized (calls <tt>initialize()</tt>). */
  DefaultProductMultiVector(
    const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const Teuchos::RCP<MultiVectorBase<Scalar> > multiVecs[]
    );

  /** \brief. Constructs to initialized (calls <tt>initialize()</tt>). */
  DefaultProductMultiVector(
    const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const Teuchos::RCP<const MultiVectorBase<Scalar> > multiVecs[]
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const int numMembers
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const Teuchos::RCP<MultiVectorBase<Scalar> > multiVecs[]
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const Teuchos::RCP<const MultiVectorBase<Scalar> > multiVecs[]
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

  /** \name Overridden from ProductMultiVectorBase */
  //@{

  /** \brief . */
  Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
  productSpace() const;
  /** \brief . */
  bool blockIsConst(const int k) const;
  /** \brief . */
  Teuchos::RCP<MultiVectorBase<Scalar> >
  getNonconstMultiVectorBlock(const int k);
  /** \brief . */
  Teuchos::RCP<const MultiVectorBase<Scalar> >
  getMultiVectorBlock(const int k) const;

  //@}

  /** \name Overriden from MultiVectorBase */
  //@{

  /** \brief . */
  Teuchos::RCP<const VectorBase<Scalar> > col(Index j) const;
  /** \brief . */
  Teuchos::RCP<VectorBase<Scalar> > col(Index j);
  /** \brief . */
  Teuchos::RCP<const MultiVectorBase<Scalar> >
  subView( const Range1D& colRng ) const;
  /** \brief . */
  Teuchos::RCP<MultiVectorBase<Scalar> >
  subView( const Range1D& colRng );
  /** \brief . */
  Teuchos::RCP<const MultiVectorBase<Scalar> >
  subView( const int numCols, const int cols[] ) const;
  /** \brief . */
  Teuchos::RCP<MultiVectorBase<Scalar> >
  subView( const int numCols, const int cols[] );
  /** \brief . */
  void mvMultiReductApplyOpImpl(
    const RTOpPack::RTOpT<Scalar> &primary_op,
    const int num_multi_vecs,
    const MultiVectorBase<Scalar>*const multi_vecs[],
    const int num_targ_multi_vecs,
    MultiVectorBase<Scalar>*const targ_multi_vecs[],
    RTOpPack::ReductTarget*const reduct_objs[],
    const Index primary_first_ele_offset,
    const Index primary_sub_dim,
    const Index primary_global_offset,
    const Index secondary_first_ele_offset,
    const Index secondary_sub_dim
    ) const;
  /** \brief . */
  void acquireDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::ConstSubMultiVectorView<Scalar> *sub_mv
    ) const;
  /** \brief . */
  void releaseDetachedMultiVectorViewImpl(
    RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv
    ) const;
  /** \brief . */
  void acquireNonconstDetachedMultiVectorViewImpl(
    const Range1D &rowRng,
    const Range1D &colRng,
    RTOpPack::SubMultiVectorView<Scalar> *sub_mv
    );
  /** \brief . */
  void commitNonconstDetachedMultiVectorViewImpl(
    RTOpPack::SubMultiVectorView<Scalar>* sub_mv
    );
  /** \brief . */
  Teuchos::RCP<MultiVectorBase<Scalar> > clone_mv() const;

  //@}

  /** \name Overriden from LinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RCP< const VectorSpaceBase<Scalar> >
  range() const;
  /** \brief . */
  Teuchos::RCP< const VectorSpaceBase<Scalar> >
  domain() const;

  //@}

protected:

  /** \name Overridden from SingleScalarLinearOpBase */
  //@{

  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  /** \brief . */
  void apply(
    const ETransp M_trans,
    const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:

  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<MultiVectorBase<Scalar> > CNMVC;

  // //////////////////////////////
  // Private data members

  Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > productSpace_;
  Teuchos::Array<CNMVC> multiVecs_;
  // cache
  int numBlocks_;

  // //////////////////////////////
  // Private member functions

  template<class MultiVectorType>
  void initializeImpl(
    const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
    const Teuchos::RCP<MultiVectorType> multiVecs[]
    );
  
  void assertInitialized() const;

  void validateColIndex(const int j) const;

  //  Not defined and not to be called!
  DefaultProductMultiVector();

};


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductMultiVector
 */
template<class Scalar>
inline
Teuchos::RCP<DefaultProductMultiVector<Scalar> >
defaultProductMultiVector(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const int numMembers
  )
{
  return Teuchos::rcp(
    new DefaultProductMultiVector<Scalar>(productSpace,numMembers)
    );
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductMultiVector
 */
template<class Scalar>
Teuchos::RCP<DefaultProductMultiVector<Scalar> >
defaultProductMultiVector(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const Teuchos::RCP<MultiVectorBase<Scalar> > multiVecs[]
  )
{
  return Teuchos::rcp(
    new DefaultProductMultiVector<Scalar>(productSpace,multiVecs)
    );
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultProductMultiVector
 */
template<class Scalar>
Teuchos::RCP<DefaultProductMultiVector<Scalar> >
defaultProductMultiVector(
  const Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > &productSpace,
  const Teuchos::RCP<const MultiVectorBase<Scalar> > multiVecs[]
  )
{
  return Teuchos::rcp(
    new DefaultProductMultiVector<Scalar>(productSpace,multiVecs)
    );
}


// /////////////////////////
// Inline members

#ifndef TEUCHOS_DEBUG


template<class Scalar>
inline
void DefaultProductMultiVector<Scalar>::assertInitialized() const
{}


template<class Scalar>
inline
void DefaultProductMultiVector<Scalar>::validateColIndex(const int j) const
{}



#endif // TEUCHOS_DEBUG



} // namespace Thyra


#endif // THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_DECL_HPP
