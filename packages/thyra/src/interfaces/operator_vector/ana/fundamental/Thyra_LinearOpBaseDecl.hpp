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
#include "Teuchos_PromotionTraits.hpp"

namespace Thyra {

/** \brief Base class for all linear operators.
 *
 * \section Thyra_LO_outline_sec Outline
 *
 * <ul>
 * <li>\ref Thyra_LO_intro_sec
 * <li>\ref Thyra_LO_spaces_sec
 * <li>\ref Thyra_LO_scalar_types_sec
 * <li>\ref Thyra_LO_single_scalar_sec
 * <li>\ref Thyra_LO_adjoint_relation_sec
 * <li>\ref Thyra_LO_aliasing_sec
 * <li>\ref Thyra_LO_optional_adjoints_sec
 * <li>\ref Thyra_LO_testing_sec
 * <li>\ref Thyra_LO_dev_notes_sec
 * </ul>
 *
 * \section Thyra_LO_intro_sec Introduction
 *
 * A linear operator can optionally perform the forward operations
 *
 * <ul>
 * <li><b>Forward non-conjugate apply</b> \verbatim Y = alpha*M*X + beta*Y \endverbatim
 * <li><b>Forward conjugate apply</b> \verbatim Y = alpha*conjugate(M)*X + beta*Y \endverbatim
 * </ul>
 *
 * through the <tt>apply()</tt> function and the operations
 *
 * <ul>
 * <li><b>Transpose non-conjugate apply</b> \verbatim Y = alpha*transpose(M)*X + beta*Y \endverbatim
 * <li><b>Transpose conjugate (i.e. adjoint) apply</b> \verbatim Y = alpha*adjoint(M)*X + beta*Y \endverbatim
 * </ul>
 *
 * through the <tt>applyTranspose()</tt> function where <tt>Y</tt> and
 * <tt>X</tt> are <tt>MultiVectorBase</tt> objects.  The reason for the exact
 * form of the above operations is that there are direct BLAS and equivalent
 * versions of these operations and performing a sum-into multiplication is
 * more efficient in general.
 *
 * \section Thyra_LO_spaces_sec Range and domain spaces
 *
 * A linear operator has vector spaces associated with it for the vectors
 * <tt>x</tt> and <tt>y</tt> that lie in the domain and the range spaces of
 * the non-transposed linear operator <tt>y = M*x</tt> and these spaces are
 * returned by <tt>domain()</tt> and <tt>range()</tt>.
 *
 * \section Thyra_LO_scalar_types_sec Support for different range and domain scalar types
 *
 * This interface allows for different scalar types for the range and domain
 * spaces, <tt>RangeScalar</tt> and <tt>DomainScalar</tt> respectively.  This
 * is needed to support such things as real-to-complex FFTs (see the example
 * <tt>RealToComplex1DFFTLinearOp</tt> for instance) and real to
 * extended-precision linear operators.
 *
 * \section Thyra_LO_single_scalar_sec Single scalar type linear operators
 *
 * While this interface supports the notion of different range and domain
 * scalar types, many different <tt>%LinearOpBase</tt> implementations and
 * many different ANAs will only support a single scalar type.  There is a lot
 * of support for single-scalar-type linear operators.
 *
 * First, note that there is a forward declaration for this class of the form
 
 \code
  template<class RangeScalar, class DomainScalar = RangeScalar> class LinearOpBase;
 \endcode

 * that allows the class to be refereed to using just one Scalar type as
 * <tt>LinearOpBase<Scalar></tt>.  This is useful for both clients and
 * subclass implementations.
 *
 * When a client ANA can only support a single-scalar type of linear operator,
 * it may be more convenient to use some of the wrapper functions such as
 * <tt>Thyra::apply()</tt> and <tt>Thyra::opSupported()</tt> that are
 * described \ref Thyra_Op_Vec_LinearOpBase_support_grp "here".
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
 * the input object <tt>X</tt> in <tt>apply()</tt> or
 * <tt>applyTranspose()</tt>.  Allowing aliasing would greatly complicate the
 * development of concrete subclasses.
 *
 * \section Thyra_LO_optional_adjoints_sec Optional support for specific types of operator applications
 *
 * This interface does not require that a linear operator implementation
 * support all of the different types of operator applications defined in the
 * \ref Thyra_LO_intro_sec "introduction" above.  If a
 * <tt>%LinearOpBase</tt> object can not support a particular type of operator
 * application, then this is determined by the functions
 * <tt>applySupports()</tt> and <tt>applyTransposeSupports()</tt>.
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
 * \section Thyra_LO_dev_notes_sec Notes for subclass developers
 *
 * There are only three functions that a concrete subclass is required to
 * override: <tt>domain()</tt>, <tt>range()</tt> and <tt>apply()</tt>.  Note
 * that the functions <tt>domain()</tt> and <tt>range()</tt> should simply
 * return <tt>VectorSpaceBase</tt> objects for subclasses that are already
 * defined for the vectors that the linear operator interacts with through the
 * function <tt>apply()</tt>.  Therefore, given that appropriate
 * <tt>VectorSpaceBase</tt> and <tt>MultiVectorBase</tt> (and/or
 * <tt>VectorBase</tt>) subclasses exist, the only real work involved in
 * implementing a <tt>LinearOpBase</tt> subclass is in defining a single
 * function <tt>apply()</tt>.
 *
 * This interface provides default implementations for the functions
 * <tt>applyTranspose()</tt> and <tt>applyTransposeSupports()</tt> where it is
 * assumed that the operator does not support transpose (or adjoint) operator
 * applications.  If transpose (and/or adjoint) operator applications can be
 * supported, then the functions <tt>applyTranspose()</tt> and
 * <tt>applyTransposeSupports()</tt> should be overridden as well.
 *
 * If possible, the subclass should also override the <tt>clone()</tt>
 * function which allows clients to create copies of a <tt>LinearOpBase</tt>
 * object.  This functionality is useful in some circumstances.  However, this
 * functionality is not required and the default <tt>clone()</tt>
 * implementation returns a null smart pointer object.
 *
 * If a concrete subclass can only support a single scalar type, then the
 * concrete subclass should perhaps inherit form the
 * <tt>SingleScalarLinearOpBase</tt> node subclass.  This node subclass
 * provides just a single function <tt>SingleScalarLinearOpBase::apply()</tt>
 * that will support both forward and transpose (adjoint) operator
 * applications.
 *
 * If, in addition to only supporting a single scalar type, a concrete
 * subclass can only support single RHS operator applications, then perhaps
 * the node subclass <tt>SingleRhsLinearOpBase</tt> should be inherited from.
 * This node subclass provides a version of
 * <tt>SingleRhsLinearOpBase::apply()</tt> that takes <tt>VectorBase</tt>
 * objects instead of <tt>MultiVectorBase</tt> objects.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class RangeScalar, class DomainScalar>
class LinearOpBase : virtual public Teuchos::Describable {
public:

  /** @name Public pure virtual functions (must be overridden by subclass) */
  //@{

  /** \brief Return a smart pointer for the range space for <tt>this</tt> operator.
   *
   * Note that a return value of <tt>return.get()==NULL</tt> is a flag
   * that <tt>*this</tt> is not fully initialized.
   *
   * If <tt>return.get()!=NULL</tt>, it is required that the object
   * referenced by <tt>*return.get()</tt> must have lifetime that
   * extends past the lifetime of the returned smart pointer object.
   * However, the object referenced by <tt>*return.get()</tt> may
   * change if <tt>*this</tt> modified so this reference should not
   * be maintained for too long.
   *
   * <b>New Behavior!</b> It is required that the <tt>%VectorSpaceBase</tt>
   * object embedded in <tt>return</tt> must be valid past the lifetime of
   * <tt>*this</tt> linear operator object.
   */
  virtual Teuchos::RefCountPtr< const VectorSpaceBase<RangeScalar> > range() const = 0;

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   *
   * Note that a return value of <tt>return.get()==NULL</tt> is a flag
   * that <tt>*this</tt> is not fully initialized.
   *
   * If <tt>return.get()!=NULL</tt>, it is required that the object
   * referenced by <tt>*return.get()</tt> must have lifetime that
   * extends past the lifetime of the returned smart pointer object.
   * However, the object referenced by <tt>*return.get()</tt> may
   * change if <tt>*this</tt> modified so this reference should not
   * be maintained for too long.
   *
   * <b>New Behavior!</b> It is required that the <tt>%VectorSpaceBase</tt>
   * object embedded in <tt>return</tt> must be valid past the lifetime of
   * <tt>*this</tt> linear operator object.
   */
  virtual Teuchos::RefCountPtr< const VectorSpaceBase<DomainScalar> > domain() const = 0;

  /** \brief Apply the forward non-conjugate or conjugate linear operator to a
   * multi-vector : <tt>Y = alpha*M*X + beta*Y</tt>.
   *
   * @param  conj
   *                [in] Determines whether the elements are non-conjugate or conjugate.
   *                The value <tt>NONCONJ_ELE</tt> gives the standard forward operator
   *                while the value of <tt>CONJ_ELE</tt> gives the forward operator
   *                with the complex conjugate matrix elements.  For a real-valued
   *                operators, this argument is ignored and has no effect.
   * @param  X      [in] The right hand side multi-vector 
   * @param  Y      [in/out] The target multi-vector being transformed
   * @param  alpha  [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
     *                The default value of <tt>alpha</tt> is </tt>1.0</tt>
   * @param  beta   [in] The multiplier for the target multi-vector <tt>Y</tt>.
   *                The default value of <tt>beta</tt> is <tt>0.0</tt>.
   * 
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->applySupports(conj)==true</tt> (throw <tt>Exceptions::OpNotSupported</tt>)
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>X.range()->isCompatible(this->domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>Y->range()->isCompatible(*this->range()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>Y->domain()->isCompatible(*X.domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>Y</tt> can not alias <tt>X</tt>.  It is up to the client to ensure that <tt>Y</tt>
   *      and <tt>X</tt> are distinct since in general this can not be verified by the implementation until,
   *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> Is it not obvious?  After the function returns the multi-vector <tt>Y</tt>
   *      is transformed as indicated above.
   * </ul>
   */
  virtual void apply(
    const EConj                             conj
    ,const MultiVectorBase<DomainScalar>    &X
    ,MultiVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha = Teuchos::ScalarTraits<RangeScalar>::one()
    ,const RangeScalar                      beta  = Teuchos::ScalarTraits<RangeScalar>::zero()
    ) const = 0;

  //@}

  /** @name Public virtual functions with default implementations */
  //@{

  /** \brief Determines if <tt>apply()</tt> supports this <tt>conj</tt> argument.
   *
   * The default implementation returns <tt>true</tt> for real valued scalar types
   * or when <tt>conj==NONCONJ_ELE</tt> for complex valued types.
   */
  virtual bool applySupports( const EConj conj ) const;

  /** \brief Determines if <tt>applyTranspose()</tt> supports this <tt>conj</tt> argument.
   *
   * The default implementation returns <tt>false</tt> which is consistent
   * with the below default implementation for <tt>applyTranspose()</tt>.
   */
  virtual bool applyTransposeSupports( const EConj conj ) const;

  /** \brief Apply the non-conjugate or conjugate transposed linear operator
   * to a multi-vector : <tt>Y = alpha*trans(M)*X + beta*Y</tt>.
   *
   * @param  conj
   *                [in] Determines whether the elements are non-conjugate or conjugate.
   *                The value <tt>NONCONJ_ELE</tt> gives the standard transposed operator
   *                with non-conjugate transposed matrix elements
   *                while the value of <tt>CONJ_ELE</tt> gives the standard adjoint operator
   *                with the complex conjugate transposed matrix elements.  For a real-valued
   *                operators, this argument is ignored and has no effect.
   * @param  X      [in] The right hand side multi-vector 
   * @param  Y      [in/out] The target multi-vector being transformed
   * @param  alpha  [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
     *                The default value of <tt>alpha</tt> is </tt>1.0</tt>
   * @param  beta   [in] The multiplier for the target multi-vector <tt>Y</tt>.
   *                The default value of <tt>beta</tt> is <tt>0.0</tt>.
   * 
   * <b>Preconditions:</b><ul>
   * <li> <tt>this->applyTransposeSupports(conj)==true</tt> (throw <tt>Exceptions::OpNotSupported</tt>)
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>X.range()->isCompatible(this->range()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>Y->range()->isCompatible(*this->domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>Y->domain()->isCompatible(*X.domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *      and <tt>X</tt> are distinct since in general this can not be verified by the implementation until,
   *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
   * <li> <tt>Y</tt> can not alias <tt>X</tt>.  It is up to the client to ensure that <tt>Y</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> Is it not obvious?  After the function returns the multi-vector <tt>Y</tt>
   *      is transformed as indicated above.
   * </ul>
   *
   * The default implementation throws an exception but gives a very good
   * error message.  The assumption here is that most linear operators will
   * not be able to support an transpose apply and that is why this default
   * implementation is provided.
   */
  virtual void applyTranspose(
    const EConj                            conj
    ,const MultiVectorBase<RangeScalar>    &X
    ,MultiVectorBase<DomainScalar>         *Y
    ,const DomainScalar                     alpha = Teuchos::ScalarTraits<DomainScalar>::one()
    ,const DomainScalar                     beta  = Teuchos::ScalarTraits<DomainScalar>::zero()
    ) const;

  /** \brief Clone the linear operator object (if supported).
   *
   * The primary purpose for this function is to allow a client to capture the
   * current state of a linear operator object and be guaranteed that some
   * other client will not alter its behavior.  A smart implementation will
   * use reference counting and lazy evaluation internally and will not
   * actually copy any large amount of data unless it has to.
   *
   * The default implementation returns <tt>return.get()==NULL</tt> which is
   * allowable.  A linear operator object is not required to return a non-NULL
   * value but many good matrix-based linear operator implementations will.
   */
  virtual Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> > clone() const;

  //@}

};	// end class LinearOpBase

/** \defgroup Thyra_Op_Vec_LinearOpBase_support_grp Support functions for LinearOpBase interface

These functions allow a client to use a <tt>LinearOpBase</tt> object more
easily in simpler use cases.

\ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp

*/
//@{

/** \brief Determines if an operation is supported for a single scalar type.
 *
 */
template<class Scalar>
inline bool opSupported( const LinearOpBase<Scalar> &M, ETransp M_trans )
{
  if(real_trans(M_trans)==NOTRANS)
    return M.applySupports(transToConj(M_trans));
  return M.applyTransposeSupports(transToConj(M_trans));
}

/** \brief Call <tt>LinearOpBase::apply()</tt> as a global function call.
 *
 * Calls <tt>M.apply(conj,X,Y,alpha,beta)</tt>.
 *
 */
template<class RangeScalar, class DomainScalar>
inline void apply(
  const LinearOpBase<RangeScalar,DomainScalar>     &M
  ,const EConj                                     conj
  ,const MultiVectorBase<DomainScalar>             &X
  ,MultiVectorBase<RangeScalar>                    *Y
  ,const RangeScalar                               alpha
#ifndef __sun
                                                         = Teuchos::ScalarTraits<RangeScalar>::one()
#endif
  ,const RangeScalar                               beta
#ifndef __sun
                                                         = Teuchos::ScalarTraits<RangeScalar>::zero()
#endif
  )
{
  M.apply(conj,X,Y,alpha,beta);
}

#ifdef __sun

template<class RangeScalar, class DomainScalar>
inline void apply(
  const LinearOpBase<RangeScalar,DomainScalar>     &M
  ,const EConj                                     conj
  ,const MultiVectorBase<DomainScalar>             &X
  ,MultiVectorBase<RangeScalar>                    *Y
  ,const RangeScalar                               alpha
  )
{
  typedef Teuchos::ScalarTraits<RangeScalar> ST;
  apply(M,conj,X,Y,alpha,ST::zero());
}

template<class RangeScalar, class DomainScalar>
inline void apply(
  const LinearOpBase<RangeScalar,DomainScalar>     &M
  ,const EConj                                     conj
  ,const MultiVectorBase<DomainScalar>             &X
  ,MultiVectorBase<RangeScalar>                    *Y
  )
{
  typedef Teuchos::ScalarTraits<RangeScalar> ST;
  apply(M,conj,X,Y,ST::one(),ST::zero());
}

#endif // __sun

/** \brief Call <tt>LinearOpBase::applyTranspose()</tt> as a global function call.
 *
 * Calls <tt>M.applyTranspose(conj,X,Y,alpha,beta)</tt>.
 *
 */
template<class RangeScalar, class DomainScalar>
inline void applyTranspose(
  const LinearOpBase<RangeScalar,DomainScalar>      &M
  ,const EConj                                      conj
  ,const MultiVectorBase<RangeScalar>               &X
  ,MultiVectorBase<DomainScalar>                    *Y
  ,const DomainScalar                               alpha
#ifndef __sun
                                                          = Teuchos::ScalarTraits<DomainScalar>::one()
#endif
  ,const DomainScalar                               beta
#ifndef __sun
                                                          = Teuchos::ScalarTraits<DomainScalar>::zero()
#endif
  )
{
  M.applyTranspose(conj,X,Y,alpha,beta);
}

#ifdef __sun

template<class RangeScalar, class DomainScalar>
inline void applyTranspose(
  const LinearOpBase<RangeScalar,DomainScalar>      &M
  ,const EConj                                      conj
  ,const MultiVectorBase<RangeScalar>               &X
  ,MultiVectorBase<DomainScalar>                    *Y
  ,const DomainScalar                               alpha
  )
{
  typedef Teuchos::ScalarTraits<DomainScalar> ST;
  applyTranspose(M,conj,X,Y,alpha,ST::zero());
}

template<class RangeScalar, class DomainScalar>
inline void applyTranspose(
  const LinearOpBase<RangeScalar,DomainScalar>      &M
  ,const EConj                                      conj
  ,const MultiVectorBase<RangeScalar>               &X
  ,MultiVectorBase<DomainScalar>                    *Y
  )
{
  typedef Teuchos::ScalarTraits<DomainScalar> ST;
  applyTranspose(M,conj,X,Y,ST::one(),ST::zero());
}

#endif

/** \brief Call <tt>LinearOpBase::apply()</tt> or
 *    <tt>LinearOpBase::applyTranspose()</tt> as a global function call (for a
 *    single scalar type).
 *
 * Calls <tt>M.apply(...,X,Y,alpha,beta)</tt> or
 * <tt>M.applyTranspose(...,X,Y,alpha,beta)</tt>.
 *
 */
template<class Scalar>
inline void apply(
  const LinearOpBase<Scalar>        &M
  ,const ETransp                    M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
#ifndef __sun
                                          = Teuchos::ScalarTraits<Scalar>::one()
#endif
  ,const Scalar                     beta
#ifndef __sun
                                          = Teuchos::ScalarTraits<Scalar>::zero()
#endif
  )
{
  if(real_trans(M_trans)==NOTRANS) {
    M.apply(transToConj(M_trans),X,Y,alpha,beta);
  }
  else {
    M.applyTranspose(transToConj(M_trans),X,Y,alpha,beta);
  }
}

#ifdef __sun

template<class Scalar>
inline void apply(
  const LinearOpBase<Scalar>        &M
  ,const ETransp                    M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  )
{
  apply(M,M_trans,X,Y,alpha,Teuchos::ScalarTraits<Scalar>::zero());
}

template<class Scalar>
inline void apply(
  const LinearOpBase<Scalar>        &M
  ,const ETransp                    M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  )
{
  apply(M,M_trans,X,Y,Teuchos::ScalarTraits<Scalar>::one());
}

#endif

//@}

}	// end namespace Thyra

#endif	// THYRA_LINEAR_OP_DECL_HPP
