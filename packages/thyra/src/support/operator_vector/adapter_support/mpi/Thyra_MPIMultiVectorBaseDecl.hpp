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

#ifndef THYRA_MPI_MULTI_VECTOR_BASE_DECL_HPP
#define THYRA_MPI_MULTI_VECTOR_BASE_DECL_HPP

#include "Thyra_MultiVectorDefaultBaseDecl.hpp"
#include "Thyra_SingleScalarEuclideanLinearOpBaseDecl.hpp"
#include "Teuchos_BLAS.hpp"

namespace Thyra {

/** \brief . */
template<class Scalar> class MPIVectorSpaceBase;

/** \brief Base class for MPI-based SPMD multi-vectors.
 *
 * By inheriting from this base class, multi-vector implementations
 * allow their multi-vector objects to be seamlessly combined with
 * other MPI-based multi-vector objects (of different concrete types)
 * in <tt>applyOp()</tt> and <tt>apply()</tt>.  A big part of this
 * protocol is that every multi-vector object can expose an
 * <tt>MPIVectorSpaceBase</tt> object through the virtual function
 * <tt>mpiSpace()</tt>.
 *
 * This base class contains an implementation of <tt>applyOp()</tt>
 * that relies on implementations of the methods (<tt>const</tt>)
 * <tt>acquireDetachedView()</tt>, <tt>releaseDetachedView()</tt>,
 * (non-<tt>const</tt>) <tt>acquireDetachedView()</tt> and
 * <tt>commitDetachedView()</tt> (which all have default
 * implementations in this subclass).  In essence, this implementation
 * will only call the <tt>acquireDetachedView()</tt> methods using a
 * range of (global) indexes for elements that exist on the local
 * processor.  As long as the number of local elements on each
 * processor is fairly large, the virtual function call overhead will
 * be minimal and this will result in a near optimal implementation.
 *
 * <b>Notes to subclass developers</b>
 *
 * Concrete subclasses must override only five functions:
 * <tt>mpiSpace()</tt>, <tt>getLocalData(const Scalar**,Index*)</tt>,
 * <tt>freeLocalData(const Scalar**,Index*)</tt>,
 * <tt>getLocalData(Scalar**,Index*)</tt>,
 * <tt>commitLocalData(Scalar**,Index*)</tt>.  Note that overriding
 * the <tt>mpiSpace()</tt> function requires implementing or using a
 * pre-implemented concrete <tt>MPIVectorSpace</tt> object.
 *
 * If the <tt>acquireDetachedView()</tt> methods are ever called with
 * index ranges outside of those of the local processor, then the
 * default implementations in <tt>MultiVectorBase</tt> of all of the
 * methods (<tt>const</tt>) <tt>MultiVectorBase::acquireDetachedView()</tt>,
 * <tt>MultiVectorBase::releaseDetachedView()</tt>, (non-<tt>const</tt>)
 * <tt>MultiVectorBase::acquireDetachedView()</tt> and
 * <tt>MultiVectorBase::commitDetachedView()</tt> are called in instead.
 * Alternatively, a subclass could provide more specialized
 * implementations of these methods (for more efficient gather/scatter
 * operations) if desired but this should not be needed for most use
 * cases.
 *
 * It is interesting to note that in the above use case that the
 * explicit subvector access methods call on its default
 * implementation defined in <tt>MultiVectorBase</tt> (which calls on
 * <tt>applyOp()</tt>) and the operator will be properly applied since
 * the version of <tt>applyOp()</tt> implemented in this class will
 * only request local vector data and hence there will only be two
 * levels of recursion for any call to an explicit subvector access
 * method.  This is a truly elegant result.
 *
 * As described in the documentation for <tt>MPIVectorSpaceBase</tt>,
 * it is possible that at runtime it may be discovered that the
 * mapping of vector data to processors does not fall under this
 * design in which case the method <tt>applyOp()</tt> should be
 * overridden to handle this which will of course remove the
 * possibility of interoperability with other MPI-based vector
 * objects.  As long as ghost data is not included this should never
 * be an issue.
 *
 * Note that multi-vector subclass derived from this base class must
 * only be directly used in SPMD mode for this to work properly.
 *
 * \ingroup Thyra_Op_Vec_adapters_MPI_support_grp
 */
template<class Scalar>
class MPIMultiVectorBase
  : virtual public MultiVectorDefaultBase<Scalar>
  , virtual public SingleScalarEuclideanLinearOpBase<Scalar>
{
public:

  /** \brief . */
  using SingleScalarEuclideanLinearOpBase<Scalar>::euclideanApply;
  /** \brief . */
  using SingleScalarEuclideanLinearOpBase<Scalar>::apply;
  /** \brief . */
  using MultiVectorDefaultBase<Scalar>::applyOp;

  /** @name  Constructors / initializers / accessors */
  //@{

  /** \brief . */
  MPIMultiVectorBase();

  //@}

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Returns the MPI-based vector space object for the range of <tt>*this</tt> multi-vector.
   */
  virtual Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const = 0;

  /** \brief Returns a non-<tt>const</tt> pointer to a Fortran-style view of the local multi-vector data.
   *
   * @param  localValues [out] On output <tt>*localValues</tt> will point to 
   *                     the first element in the first column of the local multi-vector
   *                     stored as a column-major dense Fortran-style matrix.
   * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
   *                     of the Fortran-style local multi-vector.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>leadingDim!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*leadingDim!=0</tt>
   * </ul>
   *
   * The function <tT>commitLocalData()</tt> must be called to
   * commit changes to the data.
   */
  virtual void getLocalData( Scalar **localValues, Index *leadingDim ) = 0;

  /** \brief Commit view of local data that was gotten from <tt>getLocalData()</tt>.
   *
   * @param  localValues [in/out] On input <tt>localValues</tt> must be the pointer set
   *                     by <tt>getLocalData()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*this</tt> will be updated to the entries in <tt>*localValues</tt>.
   * </ul>
   */
  virtual void commitLocalData( Scalar *localValues ) = 0;

  /** \brief Returns a <tt>const</tt>  pointer to a Fortran-style view of the local multi-vector data.
   *
   * @param  localValues [out] On output <tt>*localValues</tt> will point to 
   *                     the first element in the first column of the local multi-vector
   *                     stored as a column-major dense Fortran-style matrix.
   * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
   *                     of the Fortran-style local multi-vector.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * <li> <tt>leadingDim!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*localValues!=NULL</tt>
   * <li> <tt>*leadingDim!=0</tt>
   * </ul>
   */
  virtual void getLocalData( const Scalar **localValues, Index *leadingDim ) const = 0;

  /** \brief Free view of local data that was gotten from <tt>getLocalData()</tt>.
   *
   * @param  localValues [in/out] On input <tt>localValues</tt> must be the pointer set
   *                     by <tt>getLocalData()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>localValues!=NULL</tt>
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>*this</tt> will be updated to the entries in <tt>*localValues</tt>.
   * </ul>
   */
  virtual void freeLocalData( const Scalar *localValues ) const = 0;

  //@}

  /** @name Overridden from EuclideanLinearOpBase */
  //@{

  /// Returns <tt>mpiSpace</tt>.
  Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Calls <tt>EuclideanLinearOpBase::apply()</tt> to disambiguate <tt>apply()</tt>
   */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}

  /** @name Overridden from MultiVectorBase */
  //@{
  /** \brief . */
  void applyOp(
    const RTOpPack::RTOpT<Scalar>         &primary_op
    ,const int                            num_multi_vecs
    ,const MultiVectorBase<Scalar>*const  multi_vecs[]
    ,const int                            num_targ_multi_vecs
    ,MultiVectorBase<Scalar>*const        targ_multi_vecs[]
    ,RTOpPack::ReductTarget*const         reduct_objs[]
    ,const Index                          primary_first_ele
    ,const Index                          primary_sub_dim
    ,const Index                          primary_global_offset
    ,const Index                          secondary_first_ele
    ,const Index                          secondary_sub_dim
    ) const;
  /** \brief . */
  void acquireDetachedView(
    const Range1D                       &rowRng
    ,const Range1D                      &colRng
    ,RTOpPack::ConstSubMultiVectorView<Scalar>  *sub_mv
    ) const;
  /** \brief . */
  void releaseDetachedView( RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv ) const;
  /** \brief . */
  void acquireDetachedView(
    const Range1D                                &rowRng
    ,const Range1D                               &colRng
    ,RTOpPack::SubMultiVectorView<Scalar>    *sub_mv
    );
  /** \brief . */
  void commitDetachedView( RTOpPack::SubMultiVectorView<Scalar>* sub_mv );
  //@}

protected:

  /** @name Overridden from SingleScalarEuclideanLinearOpBase */
  //@{

  /** \brief For complex <tt>Scalar</tt> types returns <tt>true</tt> for
   * <tt>NOTRANS</tt>, <tt>TRANS</tt>, and <tt>CONJTRANS</tt> and for real
   * types returns <tt>true</tt> for all values of <tt>M_trans</tt>.
   */
  bool opSupported(ETransp M_trans) const;

  /** \brief Uses GEMM(...) and MPI_Allreduce(...) to implement.
   *
   * ToDo: Finish documentation!
   */
  void euclideanApply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;

  //@}

  /** @name Protected functions for subclasses to call. */
  //@{

  /** \brief Subclasses should call whenever the structure of the VectorSpaceBase changes.
   *
   * This function can be overridden by subclasses but this
   * particular function implementation must be called from within
   * any override.
   */
  virtual void updateMpiSpace();

  /** \brief Validate and resize the row range.
   *
   * This function throws an exception if the input range is invalid
   */
  Range1D validateRowRange( const Range1D& rowRng ) const;

  /** \brief Validate and resize the column range.
   *
   * This function throws an exception if the input range is invalid
   */
  Range1D validateColRange( const Range1D& rowCol ) const;

  //@}
  
private:
  
  // ///////////////////////////////////////
  // Private data members
  
  mutable bool in_applyOp_;

  mutable Teuchos::BLAS<int,Scalar> blas_;

  // cached
  Index  globalDim_;
  Index  localOffset_;
  Index  localSubDim_;
  Index  numCols_;

  
}; // end class MPIMultiVectorBase

} // end namespace Thyra

#endif // THYRA_MPI_MULTI_VECTOR_BASE_DECL_HPP
