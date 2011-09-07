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

#ifndef THYRA_PRODUCT_VECTOR_BASE_HPP
#define THYRA_PRODUCT_VECTOR_BASE_HPP


#include "Thyra_ProductMultiVectorBase.hpp"


namespace Thyra {


/** \brief Base interface for product vectors.
 *
 * This class defines an abstract interface for a vector that is built
 * out of the one or more other vectors to form what mathematicians
 * like to call a "product vector".
 *
 * A product vector is simply the concatenation of two or more vectors
 * to form a larger "composite" vector.  Specifically, a product vector
 * with <tt>numBlock</tt> constituent block vectors represents the blocked
 * vector

 \verbatim

         [ v[0]           ]
         [ v[1]           ]
 this =  [ .              ]
         [ v[numBlocks-1] ]
 \endverbatim
 
 * The constituent vectors <tt>v[k]</tt> can be accessed through the
 * <tt>const</tt> and non-<tt>const</tt> access functions
 * <tt>getBlock()</tt>.
 *
 * A product vector knows its product space which is returned by the
 * </tt>productSpace()</tt> function.  A <tt>%ProductVectorBase</tt>
 * object is created by a <tt>ProductVectorSpaceBase</tt> object and
 * never directly created by clients.
 *
 * This class is only an interface.  A standard implementation of this
 * interface that should be sufficient for 99% or so of use cases is
 * provided in the concrete subclass <tt>DefaultProductVector</tt>.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class ProductVectorBase
  : virtual public VectorBase<Scalar>
  , virtual public ProductMultiVectorBase<Scalar>
{
public:

  /** \brief Returns a non-persisting non-<tt>const</tt> view of the
   * (zero-based) <tt>k</tt>th block vector.
   *
   * \param k [in] The (zero-based) <tt>k</tt>th block index specifying which
   * vector block to access.
   *
   * Preconditions:<ul>
   * <li> <tt>productSpace().get()!=NULL</tt>
   * <li> <tt>0 <= k && k < productSpace()->numBlocks()</tt>
   * <li> <tt>blockIsConst(k)==false</tt>
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
  virtual RCP<VectorBase<Scalar> >
  getNonconstVectorBlock(const int k) = 0; 

  /** \brief Returns a non-persisting <tt>const</tt> view of the (zero-based)
   * <tt>k</tt>th block vector.
   *
   * \param k [in] The (zero-based) <tt>k</tt>th block index specifying which
   * vectorblock to access.
   *
   * Preconditions:<ul>
   * <li> <tt>productSpace().get()!=NULL</tt>
   * <li> <tt>0 <= k <= productSpace()->numBlocks()-1tt>
   * </ul>
   */
  virtual RCP<const VectorBase<Scalar> >
  getVectorBlock(const int k) const = 0; 

private:
  
  // Not defined and not to be called
  ProductVectorBase<Scalar>&
  operator=(const ProductVectorBase<Scalar>&);

};


/** \brief Dynamic cast from a <tt>VectorBase</tt> to a
 * <tt>ProductVectorBase</tt> object and thow exception if this fails.
 *
 * \relates ProductVectorBase
 */
template<class Scalar>
inline
RCP<Thyra::ProductVectorBase<Scalar> >
nonconstProductVectorBase(
  const RCP<Thyra::VectorBase<Scalar> > &v
  )
{
  return Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(v, true);
}


/** \brief Dynamic cast from a <tt>const VectorBase</tt> to a <tt>const
 * ProductVectorBase</tt> object and thow exception if this fails.
 *
 * \relates ProductVectorBase
 */
template<class Scalar>
inline
RCP<const Thyra::ProductVectorBase<Scalar> >
productVectorBase(
  const RCP<const Thyra::VectorBase<Scalar> > &v
  )
{
  return Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar> >(v, true);
}


} // namespace Thyra

#endif // THYRA_PRODUCT_VECTOR_BASE_HPP
