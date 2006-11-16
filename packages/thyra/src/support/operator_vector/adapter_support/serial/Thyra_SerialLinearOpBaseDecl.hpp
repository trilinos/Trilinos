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

#ifndef THYRA_SERIAL_LINEAR_OP_BASE_DECL_HPP
#define THYRA_SERIAL_LINEAR_OP_BASE_DECL_HPP

#include "Thyra_SingleScalarEuclideanLinearOpBaseDecl.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace Thyra {

/** \brief Base subclass for simplistic in-core serial linear operators.
 *
 * This subclass defines machinery for developing concrete
 * <tt>LinearOpBase</tt> subclasses for serial vectors where it is assumed
 * that all of the elements in associated vectors and multi-vectors
 * are immediately available.
 *
 * This base subclass derives from <tt>EuclideanLinearOpBase</tt> and
 * therefore any application-specific scalar products can easily be
 * incorporated.
 *
 * <b>Notes to subclass developers</b>
 *
 * The only function that a subclass must override in order to provide
 * a concrete implementation is the explicit single-vector version
 * \ref apply_expl_vec "euclideanApply()".
 *
 * This function is called on the subclass implementation passing in
 * views of explicit data.  The raw pointers to the input and
 * input/output arrays are passed in simple templated classes
 * <tt>RTOpPack::ConstSubVectorView</tt> and
 * <tt>RTOpPack::SubVectorView</tt>.  Getting raw pointers out of
 * these objects is easy.
 *
 * It is very easy to create concrete subclasses of
 * <tt>SerialLinearOpBase</tt> (see <tt>SerialTridiagLinearOp</tt> for
 * a concrete example).  All one has to do is to create a derived base class
 * (<tt>MySerialLinearOp</tt> for example) of the form:
 *
 \code
template<class Scalar>
class MySerialLinearOp : public SerialLinearOpBase<Scalar> {
private:
  // Declare your classes private data
  ...
public:
  // Declare you classes constructors, destructor, and other initialization functions
  ...
protected:
  // Override of the version of euclideanApply() that takes explicit vector data
  void euclideanApply(
    const ETransp                                M_trans
    ,const RTOpPack::ConstSubVectorView<Scalar>          &x_in
    ,const RTOpPack::SubVectorView<Scalar>   *y_out
    ,const Scalar                                alpha
    ,const Scalar                                beta
    ) const
    {
      // Get raw pointers to vector data to make me feel better!
      const Scalar *x     = x_in.values();
      const Index  x_dim  = x_in.subDim();
      Scalar       *y     = y_out->values();
      const Index  y_dim  = y_out->subDim();
      // Perform operation the operation
      if( real_trans(M_trans) == ::Thyra::NOTRANS ) {
        // Perform the non-transposed operation: y = alpha*M*x + beta*y
        ...
      }
      else {
        // Perform the transposed operation:     y = alpha*M'*x + beta*y
        ...
      }
  };
 \endcode
 *
 * If you do not need to handle arbitrary scalar data types then you
 * do not have to support them.  For example, to define a subclass
 * that only supports <b><tt>double</tt></b> you would declare a
 * non-templated version of the form:
 *
 \code
class MySerialLinearOp : public SerialLinearOpBase<double> {
private:
  // Declare your classes private data
  ...
public:
  // Declare you classes constructors, destructor and other initialization functions
  ...
protected:
  // Override of the version of euclideanApply() that takes explicit vector data
  void euclideanApply(
    const ETransp                                M_trans
    ,const RTOpPack::ConstSubVectorView<double>          &x_in
    ,const RTOpPack::SubVectorView<double>   *y_out
    ,const double                                alpha
    ,const double                                beta
    ) const
    {
      // Get raw pointers to vector data to make me feel better!
      const double *x     = x_in.values();
      const Index  x_dim  = x_in.subDim();
      double       *y     = y_out->values();
      const Index  y_dim  = y_out->subDim();
      // Perform operation the operation
      if( real_trans(M_trans) == ::Thyra::NOTRANS ) {
        // Perform the non-transposed operation: y = alpha*M*x + beta*y
        ...
      }
      else {
        // Perform the transposed operation:     y = alpha*M'*x + beta*y
        ...
      }
  };
 \endcode
 *
 * By default, pointers to explicit data returned from
 * <tt>x.values()</tt> and <tt>y->values()</tt> above are forced to
 * have unit stride to simplify things.  However, if your subclass can
 * efficiently handle non-unit stride vector data (as the BLAS can for
 * example) then you can allow this by calling the function
 * <tt>this->forceUnitStride()</tt> and passing in <tt>false</tt>.
 * The function <tt>this->forceUnitStride()</tt> can only be called by
 * your subclasses as it is declared <tt>protected</tt> so do not
 * worry about silly users messing with this, it is none of their
 * business.
 *
 * The explicit multi-vector version of \ref apply_expl_multi_vec "euclideanApply()"
 * has a default implementation that calls the explicit
 * single-vector version (that a subclass must supply) one column at a
 * time.  A subclass should only override this default multi-vector version
 * if it can do something more efficient.
 *
 * \ingroup Thyra_Op_Vec_adapters_serial_support_grp
 */
template<class Scalar>
class SerialLinearOpBase : virtual public SingleScalarEuclideanLinearOpBase<Scalar>
{
public:

  /** \brief . */
  using SingleScalarEuclideanLinearOpBase<Scalar>::euclideanApply;

  /** @name Overridden from EuclideanLinearOpBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;
  /** \brief . */
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;
  /** \brief Calls protected <tt>euclideanApply()</tt> function.
   *
   * \anchor apply_multi_vec
   */
  void euclideanApply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  //@}

protected:

  /** @name Protected constructors/initializers/accessors */
  //@{

  /** \brief Set if unit stride is forced for vector data views or not
   *
   * @param forceUnitStride  [in]
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->forceUnitStride() == forceUnitStride</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, forceUnitStride );

  /** Construct to uninitialized
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->domain().get() == NULL</tt>
   * <li><tt>this->range().get() == NULL</tt>
   * </ul>
   */
  SerialLinearOpBase();

  /** \brief Initialize vector spaces using pre-formed <tt>ScalarProdVectorSpaceBase</tt> objects.
   *
   * @param  range    [in] Smart pointer to range space
   * @param  domain   [in] Smart pointer to domain space
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>range.get()  != NULL</tt>
   * <li><tt>domain.get() != NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range().get()  == range.get()</tt>
   * <li><tt>this->domain().get() == domain.get()</tt>
   * </ul>
   */
  virtual void setSpaces(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >      &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >     &domain
    );

  /** \brief Initialize vector spaces given dimensions (uses <tt>DefaultSerialVectorSpace</tt>).
   *
   * @param  dimRange   [in] The dimension of the serial range space
   * @param  dimDomain  [in] The dimension of the serial domain space
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>dimRange  > 0</tt>
   * <li><tt>dimDomain > 0</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range()->dim()  == dimRange</tt>
   * <li><tt>this->domain()->dim() == dimDomain</tt>
   * <li><tt>dynamic_cast<const DefaultSerialVectorSpace<Scalar>*>(this->range().get())  != NULL</tt>
   * <li><tt>dynamic_cast<const DefaultSerialVectorSpace<Scalar>*>(this->domain().get()) != NULL</tt>
   * </ul>
   */
  virtual void setDimensions(
    const Index                                                 dimRange
    ,const Index                                                dimDomain
    );

  //@}

  /** @name Protected virtual functions to be overridden by subclasses */
  //@{

  /** \brief Apply the operator to explicit vector data.
   *
   * \anchor apply_expl_vec
   *
   * See <tt>LinearOpBase::euclideanApply()</tt> for a discussion of the
   * arguments to this function.  What differentiates this function is
   * that <tt>x</tt> and <tt>y</tt> are passed as objects with
   * explicit pointers to vector data.
   *
   * Since this function is protected and does not get directly called by a client.
   * Instead, this function is called by the vector version of \ref apply_vec "euclideanApply()".
   */
  virtual void euclideanApply(
    const ETransp                                M_trans
    ,const RTOpPack::ConstSubVectorView<Scalar>          &x
    ,const RTOpPack::SubVectorView<Scalar>   *y
    ,const Scalar                                alpha
    ,const Scalar                                beta
    ) const = 0;

  /** \brief Apply the operator to explicit multi-vector data.
   *
   * \anchor apply_expl_multi_vec
   *
   * See <tt>LinearOpBase::euclideanApply()</tt> for a discussion of the
   * arguments to this function.  What differentiates this function is
   * that <tt>X</tt> and <tt>Y</tt> are passed as objects with
   * explicit pointers to multi-vector data.
   *
   * Since this function is protected and does not get directly called by a client.
   * Instead, this function is called by the multi-vector version of \ref apply_multi_vec "euclideanApply()".
   *
   * The default implementation just calls the above vector version
   * one column at a time.  A subclass should only override this
   * function if it can provide a cache-smart version.  At any rate,
   * one can get up and going very quickly by just providing an
   * override for the simpler single-vector version.  Then latter, if
   * profiling data justifies it, one can provide a specialized
   * override for this function in an attempt to improve performance.
   */
  virtual void euclideanApply(
    const ETransp                                     M_trans
    ,const RTOpPack::ConstSubMultiVectorView<Scalar>          &X
    ,const RTOpPack::SubMultiVectorView<Scalar>   *Y
    ,const Scalar                                     alpha
    ,const Scalar                                     beta
    ) const;

  //@}

private:

  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >    range_;
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >    domain_;

};	// end class LinearOpBase

}	// end namespace Thyra

#endif	// THYRA_SERIAL_LINEAR_OP_BASE_DECL_HPP
