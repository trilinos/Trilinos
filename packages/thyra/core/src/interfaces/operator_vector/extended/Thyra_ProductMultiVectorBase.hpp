// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_PRODUCT_MULTI_VECTOR_BASE_HPP
#define THYRA_PRODUCT_MULTI_VECTOR_BASE_HPP

#include "Thyra_MultiVectorBase.hpp"

namespace Thyra {


template<class Scalar> class ProductVectorSpaceBase;


/** \brief Base interface for product multi-vectors.
 *
 * This class defines an abstract interface for a multi-vector that is built
 * out of the one or more other multi-vectors to form a product multi-vector.
 * This class is only an interface.  A standard implementation of this
 * interface that should be sufficient for 99% or so of use cases is provided
 * in the concrete subclass <tt>DefaultProductMultiVector</tt>.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class ProductMultiVectorBase : virtual public MultiVectorBase<Scalar> {
public:

  /** \brief Returns the associated product vector space that represents the
   * range.
   *
   * If <tt>*this</tt> is uninitialized then <tt>return.get()==NULL</tt>.
   */
  virtual Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
  productSpace() const = 0;

  /** \brief Return if the <tt>kth</tt> multi-vector block is const-only.
   *
   * \param k [in] The (zero-based) <tt>kth</tt> block index specifying which
   * multi-vector block to access.
   *
   * Preconditions:<ul>
   * <li> <tt>productSpace().get()!=NULL</tt>
   * <li> <tt>0 <= k && k < productSpace()->numBlocks()</tt>
   * </ul>
   */
  virtual bool blockIsConst(const int k) const = 0;

  /** \brief Returns a non-persisting non-<tt>const</tt> view of the
   * zero-based <tt>kth</tt> block multi-vector.
   *
   * \param k [in] The (zero-based) <tt>kth</tt> block index specifying which
   * multi-vector block to access.
   *
   * Preconditions:<ul>
   * <li> <tt>productSpace().get()!=NULL</tt>
   * <li> <tt>0 <= k && k < productSpace()->numBlocks()</tt>
   * </ul>
   *
   * Note that <tt>*this</tt> is not guaranteed to be modified until the smart
   * pointer returned from this function, as well as any other smart pointers
   * created from this smart pointer, are destroyed.  This requirement allows
   * more flexibility in how this function is implemented.
   *
   * Also note that no further interactions with <tt>*this</tt> should be
   * performed until the view returned from this function is released as
   * described above.
   */
  virtual Teuchos::RCP<MultiVectorBase<Scalar> >
  getNonconstMultiVectorBlock(const int k) = 0;

  /** \brief Returns a non-persisting <tt>const</tt> view of the (zero-based)
   * <tt>kth</tt> block multi-vector.
   *
   * \param k [in] The (zero-based) <tt>kth</tt> block index specifying which
   * multi-vector block to access.
   *
   * Preconditions:<ul>
   * <li> <tt>productSpace().get()!=NULL</tt>
   * <li> <tt>0 <= k && k < productSpace()->numBlocks()</tt>
   * </ul>
   */
  virtual Teuchos::RCP<const MultiVectorBase<Scalar> >
  getMultiVectorBlock(const int k) const = 0;

private:
  
  // Not defined and not to be called
  ProductMultiVectorBase<Scalar>&
  operator=(const ProductMultiVectorBase<Scalar>&);
  
};


} // namespace Thyra


#endif // THYRA_PRODUCT_MULTI_VECTOR_BASE_HPP
