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

#ifndef THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
#define THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP

#include "Thyra_MultiVectorBaseDecl.hpp"
#include "RTOpPack_RTOpT.hpp"

namespace Thyra {

/** \brief Node subclass that uses a default <tt>MultiVectorBase</tt>
 * implementation</tt> to provide default implementations for as many other
 * functions in <tt>MultiVectorBase</tt> interface the as is reasonable.
 *
 * <b>Notes to subclass developers</b>
 *
 * Only three function overrides are required in order to create a concrete
 * <tt>MultiVectorBase</tt> subclass: <tt>range()</tt>, <tt>domain()</tt> and
 * the non-const version of <tt>col()</tt>.  All of the other functions have
 * default implementations.  However, a good implementation will provide
 * optimized overrides of at least the functions <tt>apply()</tt> and
 * <tt>applyTranspose()</tt>.  The non-const versions of <tt>subView()</tt>
 * should be overridden if subviews are important.  The default implementation
 * will not achieve near-optimal performance in many cases.
 *
 * \ingroup Thyra_Op_Vec_general_adapter_support_code_grp
 */
template<class Scalar>
class MultiVectorDefaultBase : virtual public MultiVectorBase<Scalar>
{
public:

  /** \brief . */
  using MultiVectorBase<Scalar>::describe;
  /** \brief . */
  using MultiVectorBase<Scalar>::applyOp;
  /** \brief . */
  using MultiVectorBase<Scalar>::col; // Inject *all* functions!
  /** \brief . */
  using MultiVectorBase<Scalar>::subView; // Inject *all* functions!

  /** \name Overridden public member functions from MultiVectorBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const Range1D& colRng ) const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const Range1D& colRng );
  /** \brief . */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] ) const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > subView( const int numCols, const int cols[] );
  //@}

};

} // namespace Thyra

#endif // THYRA_MULTI_VECTOR_DEFAULT_BASE_DECL_HPP
