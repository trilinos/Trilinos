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

#ifndef THYRA_VECTOR_SERIAL_BASE_DECL_HPP
#define THYRA_VECTOR_SERIAL_BASE_DECL_HPP

#include "Thyra_VectorDefaultBaseDecl.hpp"

namespace Thyra {

/** \brief Efficient base subclass of serial vectors.
 *
 * This base subclass contains the an implementation of
 * <tt>applyOp()</tt> that relies on implementations of the methods
 * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt> and
 * <tt>commitDetachedView()</tt>.  This class also contains default
 * implementations of <tt>acquireDetachedView()</tt>,
 * <tt>releaseDetachedView()</tt> and <tt>commitDetachedView()</tt>.
 *
 * <b>Notes to subclass developers</b>
 *
 * All that is needed to develop a concrete subclass is to override the pure
 * virtual functions <tt>getData()</tt>, <tt>commitData()</tt> and
 * <tt>space()</tt> (for which a subclass of <tt>SerialVectorSpaceBase</tt>
 * should be used).
 *
 * \ingroup Thyra_Op_Vec_adapters_serial_support_grp
 */
template<class Scalar>
class SerialVectorBase : virtual public VectorDefaultBase<Scalar> {
public:

  /** \brief . */
  using VectorDefaultBase<Scalar>::applyOp;
  /** \brief . */
  using VectorBase<Scalar>::acquireDetachedView;
  /** \brief . */
  using VectorBase<Scalar>::releaseDetachedView;
  /** \brief . */
  using VectorBase<Scalar>::commitDetachedView;

  /** @name Constructors */
  //@{

  /** \brief . */
  SerialVectorBase();

  ///@}

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Sets a non-<tt>const</tt> pointer to the beginning of the
   * vector data (and its stride).
   *
   * @param  values  [out] On output <tt>*values</tt> will point to an array of the values.
   * @param  stride  [out] On output <tt>*stride</tt> will be the stride between elements in <tt>(*values)[]</tt>
   *
   * Preconditions:<ul>
   * <li> <tt>values!=NULL</tt>
   * <li> <tt>stride!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*values!=NULL</tt>
   * <li> <tt>*stride!=0</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be committed
   * back by a call to <tt>this->commitData()</tt> in case dynamic
   * memory allocation had to be used and therefore the pointer
   * returned does not point to internal storage.
   */
  virtual void getData( Scalar** values, Index* stride ) = 0;

  /** \brief Commits updated vector data that was accessed using <tt>this->getData()</tt>.
   *
   * @param  values  [in/out] On input <tt>*values</tt> must have been set by
   *                 a previous call to <tt>this->getData()</tt>.  On output
   *                 <tt>*values==NULL</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>values!=NULL</tt>
   * <li> <tt>*values!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*this</tt> will be updated to the entries in <tt>*values</tt>.
   * <li> <tt>*values==NULL</tt>
   * </ul>
   */
  virtual void commitData( Scalar** values ) = 0;

  //@}

  /** @name Virtual methods with default implementations. */
  //@{

  /** \brief Returns a <tt>const</tt> pointer to the beginning of the vector data.
   *
   * @param  values  [out] On output <tt>*values</tt> will point to an array of the values.
   * @param  stride  [out] On output <tt>*stride</tt> will be the stride between elements in <tt>(*values)[]</tt>
   *
   * Preconditions:<ul>
   * <li> <tt>values!=NULL</tt>
   * <li> <tt>stride!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*values!=NULL</tt>
   * <li> <tt>*stride!=0</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be freed by
   * a call to <tt>this->freeData()</tt> in case dynamic memory
   * allocation had to be used and therefore the pointer returned does
   * not point to internal storage.
   *
   * The default implementation performs a <tt>const_cast</tt> of
   * <tt>this</tt> and then calls the non-<tt>const</tt> version of
   * this function.  An override of this function should only be
   * provided if dynamic memory allocation is used and data copies
   * have to be performed.  If this function is overridden then the
   * function <tt>freeData()</tt> must be overridden as well!
   */
  virtual void getData( const Scalar** values, Index* stride ) const;

  /** \brief Frees a <tt>const</tt> view of vector data that was
   * accessed using <tt>this->getData()</tt>.
   *
   * @param  values  [in/out] On input <tt>*values</tt> must have been set by
   *                 a previous call to <tt>this->getData()</tt>.  On output
   *                 <tt>*values==NULL</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>values!=NULL</tt>
   * <li> <tt>*values!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*values==NULL</tt>
   * </ul>
   *
   * The default implementation performs a <tt>const_cast</tt> of
   * <tt>this</tt> and then calls the non-<tt>const</tt> function
   * <tt>this->commitData()</tt>.  If the <tt>const</tt> version of the
   * <tt>getData()</tt> function is overridden then this function must
   * be overridden also.
   */
  virtual void freeData( const Scalar** values ) const;

  //@}

  /** @name Overridden from VectorBase */
  //@{

  /** \brief Implements this method through the methods
   * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt> and
   * <tt>commitDetachedView()</tt>.
   *
   * Note that if this method is entered again before a call has
   * been completed, then this is an indication that the methods
   * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt> and/or
   * <tt>commitDetachedView()</tt> have not been overridden properly.
   */
  void applyOp(
    const RTOpPack::RTOpT<Scalar>   &op
    ,const int                      num_vecs
    ,const VectorBase<Scalar>*const vecs[]
    ,const int                      num_targ_vecs
    ,VectorBase<Scalar>*const       targ_vecs[]
    ,RTOpPack::ReductTarget         *reduct_obj
    ,const Index                    first_ele
    ,const Index                    sub_dim
    ,const Index                    global_offset
    ) const;
  /// 
  void acquireDetachedView( const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const;
  /** \brief . */
  void releaseDetachedView( RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const;
  /** \brief . */
  void acquireDetachedView( const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec );
  /** \brief . */
  void commitDetachedView( RTOpPack::SubVectorView<Scalar>* sub_vec );
  /** \brief . */
  void setSubVector( const RTOpPack::SparseSubVectorT<Scalar>& sub_vec );

  //@}

private:

  // ///////////////////////////////////////
  // Private data members
  
  mutable bool   in_applyOp_;

}; // end class SerialVectorBase

} // end namespace Thyra

#endif // THYRA_VECTOR_SERIAL_BASE_DECL_HPP
