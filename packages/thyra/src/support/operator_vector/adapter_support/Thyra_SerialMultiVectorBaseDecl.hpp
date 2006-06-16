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

#ifndef THYRA_SERIAL_MULTI_VECTOR_BASE_DECL_HPP
#define THYRA_SERIAL_MULTI_VECTOR_BASE_DECL_HPP

#include "Thyra_MultiVectorDefaultBaseDecl.hpp"
#include "Thyra_SingleScalarEuclideanLinearOpBaseDecl.hpp"
#include "Teuchos_BLAS.hpp"

namespace Thyra {

/** \brief Base class for serial shared-memory multi-vectors.
 *
 * By inheriting from this base class, multi-vector implementations
 * allow their multi-vector objects to be seamlessly combined with
 * other serial multi-vector objects (of different concrete types) in
 * <tt>applyOp()</tt> and <tt>apply()</tt>.
 *
 * This base class contains implementations of <tt>applyOp()</tt> and
 * <tt>apply()</tt> that rely on implementations of the methods
 * (<tt>const</tt>) <tt>acquireDetachedView()</tt>,
 * <tt>releaseDetachedView()</tt>, (non-<tt>const</tt>)
 * <tt>acquireDetachedView()</tt> and <tt>commitDetachedView()</tt>
 * (which all have default implementations in this subclass).  As long
 * as the number of elements is fairly large, the virtual function
 * call overhead will be minimal and this will result in a near
 * optimal implementation.
 *
 * <b>Notes to subclass developers</b>
 *
 * Concrete subclasses must override only four functions that give
 * direct access to MultiVectorBase data: <<tt>getData(const
 * Scalar**,Index*)</tt>, <tt>freeData(const
 * Scalar**,Index*)</tt>, <tt>getData(Scalar**,Index*)</tt>,
 * <tt>commitData(Scalar**,Index*)</tt>.
 *
 * \ingroup Thyra_Op_Vec_adapters_serial_support_grp
 */
template<class Scalar>
class SerialMultiVectorBase
  : virtual public MultiVectorDefaultBase<Scalar>
  , virtual public SingleScalarEuclideanLinearOpBase<Scalar>
{
public:

  /** \brief . */
  using SingleScalarEuclideanLinearOpBase<Scalar>::euclideanApply;
  /** \brief . */
  using MultiVectorDefaultBase<Scalar>::applyOp;

  /** @name  Constructors / initializers / accessors */
  //@{

  /** \brief . */
  SerialMultiVectorBase();

  //@}

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Returns a <tt>const</tt>  pointer to a Fortran-style view of the local multi-vector data.
   *
   * @param  values      [out] On output <tt>*values</tt> will point to 
   *                     the first element in the first column of the local multi-vector
   *                     stored as a column-major dense Fortran-style matrix.
   * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
   *                     of the Fortran-style local multi-vector.
   *
   */
  virtual void getData( const Scalar **values, Index *leadingDim ) const = 0;

  /** \brief Free view of local data that was gotten from <tt>getData()</tt>.
   *
   * @param  values      [in/out] On input <tt>values</tt> must be the pointer set
   *                     by <tt>getData()</tt>.
   */
  virtual void freeData( const Scalar *values ) const = 0;

  /** \brief Returns a non-<tt>const</tt> pointer to a Fortran-style view of the local multi-vector data.
   *
   * @param  values      [out] On output <tt>*values</tt> will point to 
   *                     the first element in the first column of the local multi-vector
   *                     stored as a column-major dense Fortran-style matrix.
   * @param  leadingDim  [out] On output <tt>*leadingDim</tt> gives the leading dimension
   *                     of the Fortran-style local multi-vector.
   *
   * The function <tT>commitData()</tt> must be called to
   * commit changes to the data.
   */
  virtual void getData( Scalar **values, Index *leadingDim ) = 0;

  /** \brief Commit view of local data that was gotten from <tt>getData()</tt>.
   *
   * @param  values      [in/out] On input <tt>*values</tt> must be the pointer set
   *                     by <tt>getData()</tt>.
   */
  virtual void commitData( Scalar *values ) = 0;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Calls <tt>EuclideanLinearOpBase::apply()</tt> to disambiguate <tt>apply()</tt>
   */
  /*
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  */

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
    const Range1D                               &rowRng
    ,const Range1D                              &colRng
    ,RTOpPack::ConstSubMultiVectorView<Scalar>  *sub_mv
    ) const;
  /** \brief . */
  void releaseDetachedView( RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv ) const;
  /** \brief . */
  void acquireDetachedView(
    const Range1D                           &rowRng
    ,const Range1D                          &colRng
    ,RTOpPack::SubMultiVectorView<Scalar>   *sub_mv
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

  /** \brief Uses GEMM(...) to implement.
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

  /** @name Miscellaneous functions for subclasses to call */
  //@{

  /** \brief Subclasses should call whenever the structure of any
   * <tt>VectorSpaceBase</tt> changes.
   *
   * This function can be overridden by subclasses but this
   * particular function implementation must be called from within
   * any override.
   */
  virtual void updateSpace();

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
  Index  numRows_;
  Index  numCols_;
  
}; // end class SerialMultiVectorBase

} // end namespace Thyra

#endif // THYRA_SERIAL_MULTI_VECTOR_BASE_DECL_HPP
