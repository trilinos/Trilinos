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

#ifndef THYRA_PRODUCT_VECTOR_SPACE_BASE_HPP
#define THYRA_PRODUCT_VECTOR_SPACE_BASE_HPP


#include "Thyra_VectorSpaceBase.hpp"
#include "Teuchos_ExpandScalarTypeMacros.hpp"


namespace Thyra {


/** \brief Base interface for product vector spaces.
 *
 * This class defines an abstract interface for a vector space that is
 * built out of the one or more other vector spaces to form what
 * mathematicians like to call a "product space".
 * 
 * For example, one can think of a product space as the concatenation
 * of one or more vector spaces <tt>V[k]</tt> where
 * <tt>k=0,...,numBlocks-1</tt>.  A product space <tt>Z</tt> would
 * then be represented as:
 
 \verbatim

     [ V[0]           ]
 Z = [ V[1]           ]
     [ .              ]
     [ V[numBlocks-1] ]

 \endverbatim
 *
 * The total number of constituent vector spaces is returned by the
 * <tt>numBlocks()</tt> function.  Smart pointers to the constituent
 * vector space blocks themselves are returned using the
 * <tt>getBlock()</tt> function.
 *
 * The vectors created by <tt>this->createMember()</tt> (which is
 * inherited from the <tt>VectorSpaceBase</tt> interface) must support the
 * <tt>ProductVectorBase</tt> interface
 * (i.e. <tt>dynamic_cast<ProductVectorBase<Scalar>*>(&*this->createMember())
 * != NULL</tt>).  Likewise, the multi-vectors created by
 * <tt>this->createMembers()</tt> must support the
 * <tt>ProductMultiVectorBase</tt> interface
 * (i.e. <tt>dynamic_cast<ProductMultiVectorBase<Scalar>*>(&*this->createMember())
 * != NULL</tt>)
 *
 * This class is only an interface.  A standard implementation of this
 * interface that should be sufficient for 99% or so of use cases is
 * provided in the concrete subclass <tt>DefaultProductVectorSpace</tt>.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class ProductVectorSpaceBase : virtual public VectorSpaceBase<Scalar> {
public:

  /** \brief Returns the number of blocks that make up this product space.
   *
   * Preconditions:<ul>
   * <li> <tt>this->dim() > 0</tt>.
   * </ul>
   */
  virtual int numBlocks() const = 0;

  /** \brief Returns a vector space for the <tt>kth</tt> (zero-based) block.
   *
   * Preconditions:<ul>
   * <li> <tt>0 <= k <= numBlocks()-1</tt>
   * </ul>
   */
  virtual Teuchos::RCP<const VectorSpaceBase<Scalar> > getBlock(const int k) const = 0; 

#ifdef DOXYGEN_COMPILE
private:
  const VectorSpaceBase<Scalar> *spaces;
#endif

private:
  
  // Not defined and not to be called
  ProductVectorSpaceBase<Scalar>&
  operator=(const ProductVectorSpaceBase<Scalar>&);

};


/** \brief Dynamic cast from a <tt>VectorSpaceBase</tt> to a
 * <tt>ProductVectorSpaceBase</tt> object and thow exception if this fails.
 *
 * \relates ProductVectorSpaceBase
 */
template<class Scalar>
inline
RCP<ProductVectorSpaceBase<Scalar> >
nonconstProductVectorSpaceBase(
  const RCP<VectorSpaceBase<Scalar> > &v,
  const bool forceSuccess = true
  )
{
  return Teuchos::rcp_dynamic_cast<ProductVectorSpaceBase<Scalar> >(
    v, forceSuccess);
}


/** \brief Dynamic cast from a <tt>const VectorSpaceBase</tt> to a <tt>const
 * ProductVectorSpaceBase</tt> object and thow exception if this fails.
 *
 * \relates ProductVectorSpaceBase
 */
template<class Scalar>
inline
RCP<const ProductVectorSpaceBase<Scalar> >
productVectorSpaceBase(
  const RCP<const VectorSpaceBase<Scalar> > &v,
  const bool forceSuccess = true
  )
{
  return Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<Scalar> >(
    v, forceSuccess);
}


/** \brief Inline overload of nonconstProductVectorSpaceBase<Scalar>(..) for
 * double.
 *
 * \relates ProductVectorSpaceBase
 */
inline
RCP<ProductVectorSpaceBase<double> >
nonconstProductVectorSpaceBase(
  const RCP<VectorSpaceBase<double> > &vs,
  const bool forceSuccess = true
  )
{
  return nonconstProductVectorSpaceBase<double>(vs, forceSuccess);
}


/** \brief Inline overload of productVectorSpaceBase<Scalar>(..) for double.
 *
 * \relates ProductVectorSpaceBase
 */
inline
RCP<const ProductVectorSpaceBase<double> >
productVectorSpaceBase(
  const RCP<const VectorSpaceBase<double> > &vs,
  const bool forceSuccess = true
  )
{
  return productVectorSpaceBase<double>(vs, forceSuccess);
}


} // namespace Thyra

#endif // THYRA_PRODUCT_VECTOR_SPACE_BASE_HPP
