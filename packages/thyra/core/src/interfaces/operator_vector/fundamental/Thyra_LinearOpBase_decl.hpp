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

#ifndef THYRA_LINEAR_OP_DECL_HPP
#define THYRA_LINEAR_OP_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ExpandScalarTypeMacros.hpp"
#include "Teuchos_PromotionTraits.hpp"


namespace Thyra {


/** \brief Base class for all linear operators.
 *
 * \section Thyra_LO_outline_sec Outline
 *
 * <ul>
 * <li>\ref Thyra_LO_intro_sec
 * <li>\ref Thyra_LO_spaces_sec
 * <li>\ref Thyra_LO_adjoint_relation_sec
 * <li>\ref Thyra_LO_aliasing_sec
 * <li>\ref Thyra_LO_optional_adjoints_sec
 * <li>\ref Thyra_LO_initialization_sec
 * <li>\ref Thyra_LO_testing_sec
 * <li>\ref Thyra_LO_dev_notes_sec
 * </ul>
 *
 * \section Thyra_LO_intro_sec Introduction
 *
 * A linear operator can optionally perform following operations
 *
 * <ul>
 * <li><b>Forward non-conjugate apply</b> \verbatim Y = alpha*M*X + beta*Y \endverbatim
 * <li><b>Forward conjugate apply</b> \verbatim Y = alpha*conjugate(M)*X + beta*Y \endverbatim
 * <li><b>Transpose non-conjugate apply</b> \verbatim Y = alpha*transpose(M)*X + beta*Y \endverbatim
 * <li><b>Transpose conjugate (i.e. adjoint) apply</b> \verbatim Y = alpha*adjoint(M)*X + beta*Y \endverbatim
 * </ul>
 *
 * through the <tt>apply()</tt> function where <tt>Y</tt> and <tt>X</tt> are
 * <tt>MultiVectorBase</tt> objects.  The reason for the exact form of the
 * above operations is that there are direct BLAS and equivalent versions of
 * these operations and performing a sum-into multiplication is more efficient
 * in general.
 *
 * \section Thyra_LO_spaces_sec Range and domain spaces
 *
 * A linear operator has vector spaces associated with it for the vectors
 * <tt>x</tt> and <tt>y</tt> that lie in the domain and the range spaces of
 * the non-transposed linear operator <tt>y = M*x</tt>.  These spaces are
 * returned by <tt>domain()</tt> and <tt>range()</tt>.
 *
 * \section Thyra_LO_adjoint_relation_sec Scalar products and the adjoint relation
 *
 * Note that the vector spaces returned from <tt>domain()</tt> and
 * <tt>range()</tt> may have specialized implementations of the scalar product
 * \f$<u,w>\f$ (i.e. \f$<u,w> \neq u^H w\f$ in general).  As a result, the
 * operator and adjoint operator must obey the defined scalar products.
 * Specifically, for any two vectors \f$w\in\mathcal{D}\f$ (in the domain
 * space of <tt>A</tt>) and \f$u\in\mathcal{R}\f$ (in the range space of
 * <tt>A</tt>), the adjoint operation must obey the adjoint property
 *
 \f[
  <u,A v>_{\mathcal{R}} =\!= <A^H u, v>_{\mathcal{D}}
 \f]
 *
 * where \f$<.,.>_{\mathcal{R}}\f$ is the scalar product defined by
 * <tt>this->range()->scalarProd()</tt> and \f$<.,.>_{\mathcal{D}}\f$ is the
 * scalar product defined by <tt>this->domain()->scalarProd()</tt>.  This
 * property of the adjoint can be checked numerically, if adjoints are
 * supported, using the testing utility class <tt>LinearOpTester</tt>.
 *
 * \section Thyra_LO_aliasing_sec Aliasing policy
 *
 * It is strictly forbidden to alias the input/output object <tt>Y</tt> with
 * the input object <tt>X</tt> in <tt>apply()</tt>.  Allowing aliasing would
 * greatly complicate the development of concrete subclasses.
 *
 * \section Thyra_LO_optional_adjoints_sec Optional support for specific types of operator applications
 *
 * This interface does not require that a linear operator implementation
 * support all of the different types of operator applications defined in the
 * \ref Thyra_LO_intro_sec "introduction" above.  If a <tt>%LinearOpBase</tt>
 * object can not support a particular type of operator application, then this
 * is determined by the functions <tt>opSupported()</tt>.
 *
 * \section Thyra_LO_testing_sec Testing LinearOpBase objects
 *
 * The concrete class <tt>LinearOpTester</tt> provides a full featured set of
 * tests for any <tt>%LinearOpBase</tt> object.  This testing class can check
 * if the operator is truly "linear", and/or if the adjoint relationship
 * holds, and/or if an operator is symmetric.  All of the tests are controlled
 * by the client, can be turned on and off, and pass/failure is determined by
 * tolerances that the client can specify.  In addition, this testing class
 * can also check if two linear operators are approximately equal.
 *
 * \section Thyra_LO_initialization_sec Initialization states
 *
 * A <tt>%LinearOpBase</tt> object has three different states of
 * initialization.  These three initailziation states, a description of their
 * definition, and non-member helper functions that return these states are
 * given below:
 *
 * <ul>
 *
 * <li><b>Fully Uninitialized</b>:
 *     State: <tt>(is_null(this->range()) && is_null(this->domain()))</tt>,
 *     Nonmember function: <tt>isFullyUninitialized()</tt>
 *
 * <li><b>Partially Initialized</b>:
 *     State: <tt>(!is_null(this->range()) && !is_null(this->domain()))
 *            && (!this->opSupported(M_trans))</tt>
 *            for all values of <tt>M_trans</tt>,
 *     Nonmember function: <tt>isPartiallyInitialized()</tt>
 *
 * <li><b>Fully Initialized</b>:
 *     State: <tt>(!is_null(this->range()) && !is_null(this->domain()))
 *            && (this->opSupported(M_trans)</tt>
 *            for at least one valid value for <tt>M_trans</tt>,
 *     Nonmember function: <tt>isFullyInitialized()</tt>
 *
 * </ul>
 *
 * These three different states of initialization allow for the simplification
 * of the implementation of many different types of use cases.
 *
 * \section Thyra_LO_dev_notes_sec Notes for subclass developers
 *
 * There are only foure functions that a concrete subclass is required to
 * override: <tt>domain()</tt>, <tt>range()</tt> <tt>opSupportedImpl()</tt>,
 * and <tt>applyImpl()</tt>.  Note that the functions <tt>domain()</tt> and
 * <tt>range()</tt> should simply return <tt>VectorSpaceBase</tt> objects for
 * subclasses that are already defined for the vectors that the linear
 * operator interacts with through the function <tt>apply()</tt>.  The
 * function <tt>opSupportedImpl()</tt> just returns what operations are
 * supported and is therefore trivial to implement.  Therefore, given that
 * appropriate <tt>VectorSpaceBase</tt> and <tt>MultiVectorBase</tt> (and/or
 * <tt>VectorBase</tt>) subclasses exist, the only real work involved in
 * implementing a <tt>LinearOpBase</tt> subclass is in defining a single
 * function <tt>applyImpl()</tt>.
 *
 * If possible, the subclass should also override the <tt>clone()</tt>
 * function which allows clients to create copies of a <tt>LinearOpBase</tt>
 * object.  This functionality is useful in some circumstances.  However, this
 * functionality is not required and the default <tt>clone()</tt>
 * implementation returns a null smart pointer object.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
class LinearOpBase : virtual public Teuchos::Describable {
public:

  /** @name Public interface functions */
  //@{

  /** \brief Return a smart pointer for the range space for <tt>this</tt> operator.
   *
   * Note that a return value of <tt>is_null(returnVal)</tt> is a flag that
   * <tt>*this</tt> is not fully initialized.
   *
   * If <tt>nonnull(returnVal)</tt>, it is required that the object referenced
   * by <tt>*returnVal</tt> must have lifetime that extends past the
   * lifetime of the returned smart pointer object.  However, the object
   * referenced by <tt>*returnVal</tt> may change if <tt>*this</tt>
   * modified so this reference should not be maintained for too long.
   *
   * <b>New Behavior!</b> It is required that the <tt>%VectorSpaceBase</tt>
   * object embedded in <tt>return</tt> must be valid past the lifetime of
   * <tt>*this</tt> linear operator object.
   */
  virtual RCP< const VectorSpaceBase<Scalar> > range() const = 0;

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   *
   * Note that a return value of <tt>is_null(returnVal)</tt> is a flag
   * that <tt>*this</tt> is not fully initialized.
   *
   * If <tt>nonnull(returnVal)</tt>, it is required that the object referenced
   * by <tt>*returnVal</tt> must have lifetime that extends past the lifetime
   * of the returned smart pointer object.  However, the object referenced by
   * <tt>*returnVal</tt> may change if <tt>*this</tt> modified so this
   * reference should not be maintained for too long.
   *
   * <b>New Behavior!</b> It is required that the <tt>%VectorSpaceBase</tt>
   * object embedded in <tt>return</tt> must be valid past the lifetime of
   * <tt>*this</tt> linear operator object.
   */
  virtual RCP< const VectorSpaceBase<Scalar> > domain() const = 0;

  /** \brief Return if the <tt>M_trans</tt> operation of <tt>apply()</tt> is
   * supported or not.
   *
   * Preconditions:<ul>
   * <li> <tt>isPartiallyInitialized(*this)</tt>
   * </ul>
   *
   * Note that an operator must support at least one of the values of
   * <tt>ETrans</tt> (i.e. the transposed or the non-transposed operations
   * must be supported, both can not be unsupported)
   */
  bool opSupported(EOpTransp M_trans) const
    {
      return opSupportedImpl(M_trans);
    }

  /** \brief Apply the linear operator to a multi-vector : <tt>Y =
   * alpha*op(M)*X + beta*Y</tt>.
   *
   * \param M_trans [in] Determines whether the operator is applied or the
   * adjoint for <tt>op(M)</tt>.
   *
   * \param X [in] The right hand side multi-vector.
   *
   * \param Y [in/out] The target multi-vector being transformed.  When
   * <tt>beta==0.0</tt>, this multi-vector can have uninitialized elements.
   *
   * \param alpha [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
   * The default value of <tt>alpha</tt> is </tt>1.0</tt>
   *
   * \param beta [in] The multiplier for the target multi-vector <tt>Y</tt>.
   * The default value of <tt>beta</tt> is <tt>0.0</tt>.
   * 
   * <b>Preconditions:</b><ul>
   *
   * <li> <tt>nonnull(this->domain()) && nonnull(this->range())</tt>
   *
   * <li> <tt>this->opSupported(M_trans)==true</tt> (throw
   * <tt>Exceptions::OpNotSupported</tt>)
   *
   * <li> <tt>X.range()->isCompatible(*op(this)->domain()) == true</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li> <tt>Y->range()->isCompatible(*op(this)->range()) == true</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li> <tt>Y->domain()->isCompatible(*X.domain()) == true</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li> <tt>Y</tt> can not alias <tt>X</tt>.  It is up to the client to
   * ensure that <tt>Y</tt> and <tt>X</tt> are distinct since in general this
   * can not be verified by the implementation until, perhaps, it is too late.
   * If possible, an exception will be thrown if aliasing is detected.
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> Is it not obvious?  After the function returns the multi-vector <tt>Y</tt>
   *      is transformed as indicated above.
   * </ul>
   */
  void apply(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const
    {
      applyImpl(M_trans, X, Y, alpha, beta);
    }

  /** \brief Clone the linear operator object (if supported).
   *
   * The primary purpose for this function is to allow a client to capture the
   * current state of a linear operator object and be guaranteed that some
   * other client will not alter its behavior.  A smart implementation will
   * use reference counting and lazy evaluation internally and will not
   * actually copy any large amount of data unless it has to.
   *
   * The default implementation returns <tt>is_null(returnVal)</tt> which is
   * allowable.  A linear operator object is not required to return a non-NULL
   * value but many good matrix-based linear operator implementations will.
   */
  virtual RCP<const LinearOpBase<Scalar> > clone() const;

  //@}

  /** \name Deprecated. */
  //@{

  /** \brief Deprecated. */
  THYRA_DEPRECATED bool applySupports( const EConj conj ) const;

  /** \brief Deprecated. */
  THYRA_DEPRECATED void apply(
    const EConj conj,
    const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
    const Scalar alpha = static_cast<Scalar>(1.0),
    const Scalar beta = static_cast<Scalar>(0.0)
    ) const;

  /** \brief Deprecated. */
  THYRA_DEPRECATED bool applyTransposeSupports( const EConj conj ) const;

  /** \brief Deprecated. */
  THYRA_DEPRECATED void applyTranspose(
    const EConj conj,
    const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
    const Scalar alpha = static_cast<Scalar>(1.0),
    const Scalar beta = static_cast<Scalar>(0.0)
    ) const;

  //@}

protected:

  /** \name Protected virtual functions to be overridden by subclasses. */
  //@{

  /** \brief Override in subclass. */
  virtual bool opSupportedImpl(EOpTransp M_trans) const = 0;

  /** \brief Override in subclass. */
  virtual void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const = 0;

  //@}

private:
  
  // Not defined and not to be called
  LinearOpBase<Scalar>&
  operator=(const LinearOpBase<Scalar>&);

};


/** \brief Determines if a linear operator is in the "Fully Uninitialized"
 * state or not.
 *
 * \relates LinearOpBase
 */
template<class Scalar>
bool isFullyUninitialized( const LinearOpBase<Scalar> &M );


/** \brief Determines if a linear operator is in the "Partially Initialized"
 * state or not.
 *
 * \relates LinearOpBase
 */
template<class Scalar>
bool isPartiallyInitialized( const LinearOpBase<Scalar> &M );


/** \brief Determines if a linear operator is in the "Fully Initialized"
 * state or not.
 *
 * \relates LinearOpBase
 */
template<class Scalar>
bool isFullyInitialized( const LinearOpBase<Scalar> &M );


/** \brief Determines if an operation is supported for a single scalar type.
 *
 * \relates LinearOpBase
 */
template<class Scalar>
inline
bool opSupported( const LinearOpBase<Scalar> &M, EOpTransp M_trans );


/** \brief Non-member function call for <tt>M.apply(...)</tt>.
 *
 * \relates LinearOpBase
 */
template<class Scalar>
void apply(
  const LinearOpBase<Scalar> &M,
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha = static_cast<Scalar>(1.0),
  const Scalar beta = static_cast<Scalar>(0.0)
  );


/** \brief Calls <tt>apply<double>(...)</tt>.
 *
 * Non-tempalted double inlined non-member helper function.
 *
 * \relates LinearOpBase
 */
inline
void apply(
  const LinearOpBase<double> &M,
  const EOpTransp M_trans,
  const MultiVectorBase<double> &X,
  const Ptr<MultiVectorBase<double> > &Y,
  const double alpha = 1.0,
  const double beta = 0.0
  );


// Deprecated


/** \brief Deprecated. */
template<class Scalar>
inline
THYRA_DEPRECATED void apply(
  const LinearOpBase<Scalar> &M,
  const EConj conj,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha = static_cast<Scalar>(1.0),
  const Scalar beta = static_cast<Scalar>(0.0)
  );


/** \brief Deprecated. */
template<class Scalar>
inline
THYRA_DEPRECATED void applyTranspose(
  const LinearOpBase<Scalar> &M,
  const EConj conj,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha = static_cast<Scalar>(1.0),
  const Scalar beta = static_cast<Scalar>(0.0)
  );


/** \brief Deprecated. */
template<class Scalar>
THYRA_DEPRECATED void apply(
  const LinearOpBase<Scalar> &M,
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha = static_cast<Scalar>(1.0),
  const Scalar beta = static_cast<Scalar>(0.0)
  );


}	// end namespace Thyra


//
// Inline and other Template Implementations
//


template<class Scalar>
inline
bool Thyra::isFullyUninitialized( const LinearOpBase<Scalar> &M )
{
  return ( is_null(M.range()) || is_null(M.domain()) );
}


template<class Scalar>
bool Thyra::isPartiallyInitialized( const LinearOpBase<Scalar> &M )
{
  return
    (
      ( !is_null(M.range()) && !is_null(M.domain()) )
      && 
      (
        !opSupported(M,NOTRANS) && !opSupported(M,CONJ)
        && !opSupported(M,TRANS) && !opSupported(M,CONJTRANS)
        )
      );
}


template<class Scalar>
bool Thyra::isFullyInitialized( const LinearOpBase<Scalar> &M )
{
  return
    (
      ( !is_null(M.range()) && !is_null(M.domain()) )
      && 
      (
        opSupported(M,NOTRANS) || opSupported(M,CONJ)
        || opSupported(M,TRANS) || opSupported(M,CONJTRANS)
        )
      );
}


template<class Scalar>
inline
bool Thyra::opSupported( const LinearOpBase<Scalar> &M, EOpTransp M_trans )
{
  return M.opSupported(M_trans);
}


inline
void Thyra::apply(
  const LinearOpBase<double> &M,
  const EOpTransp M_trans,
  const MultiVectorBase<double> &X,
  const Ptr<MultiVectorBase<double> > &Y,
  const double alpha,
  const double beta
  )
{
  apply<double>(M, M_trans, X, Y, alpha, beta);
}


// Deprecated


template<class Scalar>
inline
void Thyra::apply(
  const LinearOpBase<Scalar> &M,
  const EConj conj,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha,
  const Scalar beta
  )
{
  M.apply(conj, X, Y, alpha, beta);
}


template<class Scalar>
inline
void Thyra::applyTranspose(
  const LinearOpBase<Scalar> &M,
  const EConj conj,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha,
  const Scalar beta
  )
{
  M.applyTranspose(conj, X, Y, alpha, beta);
}


template<class Scalar>
inline
void Thyra::apply(
  const LinearOpBase<Scalar> &M,
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  MultiVectorBase<Scalar> *Y,
  const Scalar alpha,
  const Scalar beta
  )
{
  apply(M, M_trans, X, Teuchos::ptr(Y), alpha, beta);
}


#endif	// THYRA_LINEAR_OP_DECL_HPP
