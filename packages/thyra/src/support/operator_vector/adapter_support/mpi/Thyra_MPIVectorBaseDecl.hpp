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

#ifndef THYRA_MPI_VECTOR_BASE_DECL_HPP
#define THYRA_MPI_VECTOR_BASE_DECL_HPP

#include "Thyra_VectorDefaultBaseDecl.hpp"
#include "Thyra_MPIVectorSpaceDefaultBaseDecl.hpp"

#define THYRA_MPI_VECTOR_BASE_DUMP

namespace Thyra {

/** \brief Base class for MPI-based SPMD vectors.
 *
 * By inheriting from this base class, vector implementations allow their
 * vector objects to be seamlessly combined with other MPI-based vector
 * objects (of different concrete types) in <tt>applyOp()</tt>.  A big part of
 * this protocol is that every vector object can expose an
 * <tt>MPIVectorSpaceBase</tt> object through the virtual function
 * <tt>mpiSpace()</tt>.
 *
 * This base class contains an implementation of <tt>applyOp()</tt> that
 * relies on implementations of the methods (<tt>const</tt>)
 * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt>, (non-<tt>const</tt>)
 * <tt>acquireDetachedView()</tt> and <tt>commitDetachedView()</tt> (which all have
 * default implementations in this subclass).  In essence, this implementation
 * will only call the <tt>acquireDetachedView()</tt> methods using a range of
 * (global) indexes for elements that exist on the local processor.  As long
 * as the number of local elements on each processor is fairly large, the
 * virtual function call overhead will be minimal and this will result in a
 * near optimal implementation.
 *
 * <b>Notes to subclass developers</b>
 *
 * Concrete subclasses must override only two functions: <tt>mpiSpace()</tt>
 * and <tt>getLocalData(Scalar**,Index*)</tt>.  The default implementation of
 * <tt>getLocalData(cons Scalar**,Index*)</tt> should rarely need to be
 * overridden as it just calls the pure-virtual non-<tt>const</tt> version.
 * Note that providing an implementation for <tt>mpiSpace()</tt> of course
 * means having to implement or use a pre-implemented
 * <tt>MPIVectorSpaceBase</tt> subclass.
 *
 * Vector subclasses must also call the protected function
 * <tt>updateMpiState()</tt> whenever the state of <tt>*this->mpiSpace()</tt>
 * vector space changes.  This function gathers some cached data that makes
 * the rest of the class more efficient.  This function must be called in a
 * constructor or any other function that changes the state of the vector
 * space.
 *
 * If the <tt>acquireDetachedView()</tt> methods are ever called with index ranges
 * outside of those of the local processor, then the default implementations
 * in <tt>VectorBase</tt> of all of the methods (<tt>const</tt>)
 * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt>, (non-<tt>const</tt>)
 * <tt>acquireDetachedView()</tt> and <tt>commitDetachedView()</tt> are called instead.
 * Alternatively, a subclass could provide more specialized implementations of
 * these methods (for more efficient gather/scatter operations) if desired but
 * this should not be needed for most use cases.
 *
 * It is interesting to note that in the above use case when the explicit
 * subvector access methods call on its default implementation defined in
 * <tt>VectorBase</tt> (which calls on <tt>applyOp()</tt>), the operator will
 * be properly applied since the version of <tt>applyOp()</tt> implemented in
 * this class will only request local vector data and hence there will only be
 * two levels of recursion for any call to an explicit subvector access
 * method.  This is a truly elegant result.
 *
 * As described in the documentation for <tt>MPIVectorSpaceBase</tt>, it is
 * possible that at runtime it may be discovered that the mapping of vector
 * data to processors does not fall under this design in which case the method
 * <tt>applyOp()</tt> should be overridden to handle this which will of course
 * remove the possibility of interoperability with other MPI-based vector
 * objects.  This, however, should never be the case.
 *
 * Note that vector subclass derived from this node interface class must only
 * be directly used in SPMD mode to work properly.
 *
 * \ingroup Thyra_Op_Vec_adapters_MPI_support_grp
 */
template<class Scalar>
class MPIVectorBase : virtual public VectorDefaultBase<Scalar> {
public:

  /** \brief . */
  using VectorDefaultBase<Scalar>::applyOp;

  /** \brief . */
  MPIVectorBase();

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Returns the MPI-based vector space object for <tt>*this</tt> vector.
   */
  virtual Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;

  /** \brief Returns a non-<tt>const</tt> pointer to the beginning of the
   * local vector data (and its stride).
   *
   * @param  localValues  [out] On output <tt>*localValues</tt> will point to an array of the local values.
   * @param  stride  [out] On output <tt>*stride</tt> will be the stride between elements in <tt>(*localValues)[]</tt>
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>stride!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*stride!=0</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be committed
   * back by a call to <tt>this->commitLocalData()</tt> in case dynamic
   * memory allocation had to be used and therefore the pointer return
   * does not point to internal storage.
   */
  virtual void getLocalData( Scalar** localValues, Index* stride ) = 0;

  /** \brief Commits updated local vector data that was accessed using <tt>this->getLocalData()</tt>.
   *
   * @param  localValues  [in/out] On input <tt>localValues</tt> must have been set by
   *                      a previous call to <tt>this->getData()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*this</tt> will be updated to the entries in <tt>*localValues</tt>.
   * </ul>
   */
  virtual void commitLocalData( Scalar* localValues ) = 0;

  //@}

  /** @name Virtual methods with default implementations. */
  //@{

  /** \brief Returns a <tt>const</tt> pointer to the beginning of the local
   * vector data.
   *
   * @param  localValues  [out] On output <tt>*localValues</tt> will point to an array of the local values.
   * @param  stride  [out] On output <tt>*stride</tt> will be the stride between elements in <tt>(*localValues)[]</tt>
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>stride!=NULL</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*stride!=0</tt>
   * </ul>
   *
   * Note, the data view returned from this function must be freed by
   * a call to <tt>this->freeLocalData()</tt> in case dynamic memory
   * allocation had to be used and therefore the pointer returned does
   * not point to internal storage.
   *
   * The default implementation performs a <tt>const_cast</tt> of
   * <tt>this</tt> and then calls the non-<tt>const</tt> version of
   * this function.  An override of this function should only be
   * provided if dynamic memory allocation is used and data copies
   * have to be performed.  If this function is overridden then the
   * function <tt>freeLocalData()</tt> must be overridden as well!
   */
  virtual void getLocalData( const Scalar** localValues, Index* stride ) const;

  /** \brief Frees a <tt>const</tt> view of local vector data that was
   * accessed using <tt>this->getLocalData()</tt>.
   *
   * @param  localValues  [in/out] On input <tt>localValues</tt> must have been set by
   *                 a previous call to <tt>this->getLocalData()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * </ul>
   *
   * The default implementation performs a <tt>const_cast</tt> of
   * <tt>this</tt> and then calls the non-<tt>const</tt> function
   * <tt>this->commitLocalData()</tt>.  If the <tt>const</tt> version of the
   * <tt>getData()</tt> function is overridden then this function must
   * be overridden also.
   */
  virtual void freeLocalData( const Scalar* localValues ) const;

  /** \brief Implements the <tt>%applyOp()</tt> method through the
   * methods <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt> and
   * <tt>commitDetachedView()</tt> as described above.
   *
   * \param  comm
   *           [in] The MPI communicator to use in the global reduction.
   *           If <tt>comm==MPI_COMM_NULL</tt>, then the local communicator
   *           will be used instead.
   *
   * Note that if this method is entered again before a call has
   * been completed, then this is an indication that the methods
   * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt> and/or
   * <tt>commitDetachedView()</tt> have not been overridden properly
   * and this method will then throw an exception.
   */
  void applyOp(
    MPI_Comm                        comm
    ,const RTOpPack::RTOpT<Scalar>  &op
    ,const int                      num_vecs
    ,const VectorBase<Scalar>*      vecs[]
    ,const int                      num_targ_vecs
    ,VectorBase<Scalar>*            targ_vecs[]
    ,RTOpPack::ReductTarget         *reduct_obj
    ,const Index                    first_ele_offset
    ,const Index                    sub_dim
    ,const Index                    global_offset
    ) const;

  //@}

  /** @name Overridden from MultiVectorBase */
  //@{

  /// Returns <tt>this->mpiSpace()</tt>
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > space() const;
  /** \brief Calls <tt>this->applyOp(MPI_COMM_NULL,op,...)</tt>. */
  void applyOp(
    const RTOpPack::RTOpT<Scalar>   &op
    ,const int                      num_vecs
    ,const VectorBase<Scalar>*      vecs[]
    ,const int                      num_targ_vecs
    ,VectorBase<Scalar>*            targ_vecs[]
    ,RTOpPack::ReductTarget         *reduct_obj
    ,const Index                    first_ele_offset
    ,const Index                    sub_dim
    ,const Index                    global_offset
    ) const;
  /// Implemented through <tt>this->getLocalData()</tt>
  void acquireDetachedView( const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const;
  /// Implemented through <tt>this->freeLocalData()</tt>
  void releaseDetachedView( RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const;
  /// Implemented through <tt>this->getLocalData()</tt>
  void acquireDetachedView( const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec );
  /// Implemented through <tt>this->commitLocalData()</tt>
  void commitDetachedView( RTOpPack::SubVectorView<Scalar>* sub_vec );

  //@}

protected:

  /** \brief Subclasses must call this function whenever the structure of the
   * VectorSpaceBase changes.
   */
  virtual void updateMpiSpace();

private:

  // ///////////////////////////////////////
  // Private data members
  
  mutable bool in_applyOp_;

  // Cached (only on vector space!)
  mutable Index  globalDim_;
  mutable Index  localOffset_;
  mutable Index  localSubDim_;

  // /////////////////////////////////////
  // Private member functions

  Range1D validateRange( const Range1D& rng_in ) const;

#ifdef THYRA_MPI_VECTOR_BASE_DUMP
public:
  static bool show_dump;
#endif // THYRA_MPI_VECTOR_BASE_DUMP

}; // end class MPIVectorBase

} // end namespace Thyra

#endif // THYRA_MPI_VECTOR_BASE_DECL_HPP
