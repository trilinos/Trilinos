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

#include "Thyra_OpBaseDecl.hpp"

namespace Thyra {

/** \brief Base class for all linear operators.
 *
 * A linear operator can perform the following operations:
 *
 * <ul>
 * <li><tt>y = alpha*op(M)*x + beta*y</tt>  (vector version)
 * <li><tt>Y = alpha*op(M)*X + beta*Y</tt>  (multi-vector version)
 * </ul>
 *
 * through the <tt>apply()</tt> functions where <tt>y</tt> and
 * <tt>x</tt> are <tt>VectorBase</tt> objects while <tt>Y</tt> and
 * <tt>X</tt> are <tt>MultiVectorBase</tt> objects.  The reason for the
 * exact form of the above operations is that there are direct BLAS
 * and equivalent versions of these operations and performing a
 * sum-into multiplication is more efficient in general.
 *
 * A linear operator has vector spaces associated with it for the
 * vectors <tt>x</tt> and <tt>y</tt> that lie in the domain and the
 * range spaces of the non-transposed linear operator <tt>y = M*x</tt>
 * and these spaces are returned by <tt>domain()</tt> and
 * <tt>range()</tt>.
 *
 * Note that the vector spaces returned from <tt>domain()</tt> and
 * <tt>range()</tt> may have specialized implementations of the scalar
 * product \f$<u,w>\f$ (i.e. \f$<u,w> \neq u^T w\f$ in general).  As a
 * result, the operator and adjoint operator must obey the defined
 * scalar product.  Specifically, for any two vectors \f$w\f$ (in the
 * domain space \f$\mathcal{D}\f$) and \f$u\f$ (in the range space
 * \f$\mathcal{R}\f$) the adjoint operation must obey:
 *
 \f[
  <u,A v>_{\mathcal{R}} = <A^T u, v>_{\mathcal{D}}
 \f]
 *
 * where \f$<.,.>_{\mathcal{R}}\f$ is the scalar product defined by
 * <tt>this->range()->scalarProd()</tt> and \f$<.,.>_{\mathcal{D}}\f$
 * is the scalar product defined by
 * <tt>this->domain()->scalarProd()</tt>.  This property of the
 * adjoint can be checked numerically, if adjoints are supported,
 * using the testing utility class <tt>LinearOpTester</tt>.
 *
 * Note that it is strictly forbidden to alias the input/output
 * objects <tt>y</tt> and <tt>Y</tt> with the input objects <tt>x</tt>
 * and <tt>X</tt>.
 *
 * If a <tt>%LinearOpBase</tt> subclass can not support a particular value
 * of <tt>M_tans</tt> in the <tt>apply()</tt> functions, then the function
 * <tt>opSupported()</tt> must return <tt>false</tt> for that
 * particular value of <tt>M_trans</tt>.
 *
 * <b>Notes for subclass developers</b>
 *
 * There are only three functions that a subclass is required to
 * override: <tt>domain()</tt>, <tt>range()</tt> and <tt>apply()</tt>.
 * Note that the functions <tt>domain()</tt> and <tt>range()</tt> should
 * simply return <tt>VectorSpaceBase</tt> objects for subclasses that are
 * already defined for the vectors that the linear operator interacts
 * with through the function <tt>apply()</tt>.  Therefore, given that
 * appropriate <tt>VectorSpaceBase</tt> and <tt>VectorBase</tt> subclasses exist,
 * the only real work involved in implementing a <tt>LinearOpBase</tt> 
 * is defining a single function <tt>apply()</tt>.
 *
 * If a <tt>LinearOpBase</tt> subclass can not support a particular value
 * of the transpose argument <tt>M_trans</tt> in the <tt>apply()</tt>
 * functions, then the function <tt>opSupported(M_trans)</tt> must be
 * overridden to return <tt>false</tt> for this value of
 * <tt>M_trans</tt>.
 *
 * If possible, the subclass should also override the <tt>clone()</tt>
 * function with allows clients to create copies of a
 * <tt>LinearOpBase</tt> object.  This functionality is very important in
 * some circumstances.  However, this functionality is not required
 * and the default <tt>clone()</tt> implementation returns a null
 * smart pointer object.
 *
 * If multi-vectors are supported in general by the application and
 * linear algebra library then, if possible, the subclass should also
 * override the multi-vector version of <tt>apply()</tt>.  In many
 * cases, a specialized multi-vector version will outperform the
 * default implementation (which is based on the single vector
 * version) in this class.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
class LinearOpBase : virtual public OpBase<Scalar> {
public:

  /** \brief . */
  using OpBase<Scalar>::describe;

  /** @name Pure virtual functions (must be overridden by subclass) */
  //@{

  /** \brief Apply the linear operator (or its transpose) to a vector:
   * <tt>y = alpha*op(M)*x + beta*y</tt>.
   *
   * @param  M_trans
   *                [in] Determines whether the transposed or non-transposed
   *                operator is applied as:
   *                <ul>
   *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
   *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
   *                </ul>
   *                where <tt>M == *this</tt>
   * @param  x      [in] The right hand side vector 
   * @param  y      [in/out] The target vector being transformed
   * @param  alpha  [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
     *                The default value of <tt>alpha</tt> is </tt>1.0</tt>
   * @param  beta   [in] The multiplier for the target vector <tt>y</tt>.
   *                The default value of <tt>beta</tt> is <tt>0.0</tt>.
   * 
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>this->opSupported(M_trans)==true</tt> (throw <tt>Exceptions::OpNotSupported</tt>)
   * <li> <tt>y->space()->isCompatible(M_trans==NOTRANS ? *this->range() : *this->domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>x.space()->isCompatible(M_trans==NOTRANS ? *this->domain() : *this->range()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>y</tt> can not alias <tt>x</tt>.  It is up to the client to ensure that <tt>y</tt>
   *      and <tt>x</tt> are distinct since in general this can not be verified by the implementation until,
   *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
   * </ul>
   *
   * Postconditions:<ul>
   * <li> Is it not obvious?  After the function returns the vector <tt>y</tt>
   *      is transformed as indicated above.
   * </ul>
   */
  virtual void apply(
    const ETransp                M_trans
    ,const VectorBase<Scalar>    &x
    ,VectorBase<Scalar>          *y
    ,const Scalar                alpha = Teuchos::ScalarTraits<Scalar>::one()
    ,const Scalar                beta  = Teuchos::ScalarTraits<Scalar>::zero()
    ) const = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief Clone the linear operator object (if supported).
   *
   * The primary purpose for this function is to allow a client to
   * capture the current state of a linear operator object and be
   * guaranteed that some other client will not alter its behavior.
   * A smart implementation will use reference counting and lazy
   * evaluation internally and will not actually copy any large
   * amount of data unless it has to.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * which is allowable by this specification.  A linear operator
   * object is not required to return a non-NULL value but almost
   * every good linear operator implementation should and will.
   */
  virtual Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;

  /** \brief Apply the linear operator (or its transpose) to a multi-vector :
   * <tt>Y = alpha*op(M)*X + beta*Y</tt>.
   *
   * @param  M_trans
   *                [in] Determines whether the transposed or non-transposed
   *                operator is applied as:
   *                <ul>
   *                <li> <tt>op(M) = M</tt>, for <tt>M_trans==NOTRANS</tt>
   *                <li> <tt>op(M) = M'</tt>, for <tt>M_trans==TRANS</tt>
   *                </ul>
   *                where <tt>M == *this</tt>
   * @param  X      [in] The right hand side multi-vector 
   * @param  Y      [in/out] The target multi-vector being transformed
   * @param  alpha  [in] Scalar multiplying <tt>M</tt>, where <tt>M==*this</tt>.
     *                The default value of <tt>alpha</tt> is </tt>1.0</tt>
   * @param  beta   [in] The multiplier for the target multi-vector <tt>y</tt>.
   *                The default value of <tt>beta</tt> is <tt>0.0</tt>.
   * 
   * Preconditions:<ul>
   * <li> <tt>this->domain().get()!=NULL && this->range().get()!=NULL</tt> (throw <tt>std::logic_error</tt>)
   * <li> <tt>this->opSupported(M_trans)==true</tt> (throw <tt>Exceptions::OpNotSupported</tt>)
   * <li> <tt>Y->range()->isCompatible(M_trans==NOTRANS ? *this->range() : *this->domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>Y->domain()->isCompatible(*X.domain()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>X.range()->isCompatible(M_trans==NOTRANS ? *this->domain() : *this->range()) == true</tt>
   *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li> <tt>Y</tt> can not alias <tt>X</tt>.  It is up to the client to ensure that <tt>Y</tt>
   *      and <tt>X</tt> are distinct since in general this can not be verified by the implementation until,
   *      perhaps, it is too late.  If possible, an exception will be thrown if aliasing is detected.
   * </ul>
   *
   * Postconditions:<ul>
   * <li> Is it not obvious?  After the function returns the multi-vector <tt>Y</tt>
   *      is transformed as indicated above.
   * </ul>
   *
   * This function has a default implementation in terms of the
   * <tt>apply()</tt> function for vectors.
   */
  virtual void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha = Teuchos::ScalarTraits<Scalar>::one()
    ,const Scalar                     beta  = Teuchos::ScalarTraits<Scalar>::zero()
    ) const;

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{

  /** \brief Generates a default outputting for all linear operators.
   *
   * Calls on the <tt>this->describe(void)</tt> function for the name
   * of the class (and possibly its instance name) and then if
   * <tt>verbLevel >= VERB_EXTREME</tt>, then the linear operators
   * elements themselves are printed as well.  The format of the
   * output is as follows:
   *
   \verbatim

   type = 'this->describe()', rangeDim = m, domainDim = n
     1:1:a11 1:2:a12 ... 1:n:a1n
     2:1:a21 2:2:a22 ... 1:n:a2n
     .       .           .
     .       .           .
     .       .           .
     m:1:am1 m:2:am2 ... m:n:amn
   \endverbatim
   *
   * Before <tt>type = 'this->describe()'</tt> is printed and after
   * each newline, <tt>leadingIndent</tt> is output.  The
   * <tt>index:value</tt> lines are offset an additional
   * <tt>indentSpacer</tt> amount.  A newline is printed after the
   * last <tt>m:n:amn</tt> entry.
   */
  std::ostream& describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const;

  //@}

};	// end class LinearOpBase

/** \brief Call <tt>LinearOpBase<Scalar>::apply()</tt> as a global function call.
 *
 * Calls <tt>M.apply(M_trans,x,y,alpha,beta);</tt>
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
inline void apply(
  const LinearOpBase<Scalar>   &M
  ,const ETransp               M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
#ifndef __sun
                                          = Teuchos::ScalarTraits<Scalar>::one()
#endif
  ,const Scalar                beta
#ifndef __sun
                                          = Teuchos::ScalarTraits<Scalar>::zero()
#endif
  )
{
  M.apply(M_trans,x,y,alpha,beta);
}

#ifdef __sun

template<class Scalar>
inline void apply(
  const LinearOpBase<Scalar>   &M
  ,const ETransp               M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
  )
{
  apply(M,M_trans,x,y,alpha,Teuchos::ScalarTraits<Scalar>::zero());
}

template<class Scalar>
inline void apply(
  const LinearOpBase<Scalar>   &M
  ,const ETransp               M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  )
{
  apply(M,M_trans,x,y,Teuchos::ScalarTraits<Scalar>::one());
}

#endif

/** \brief Call <tt>LinearOpBase<Scalar>::apply()</tt> as a global function call.
 *
 * Calls <tt>M.apply(M_trans,X,Y,alpha,beta);</tt>
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
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
  M.apply(M_trans,X,Y,alpha,beta);
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

}	// end namespace Thyra

#endif	// THYRA_LINEAR_OP_DECL_HPP
