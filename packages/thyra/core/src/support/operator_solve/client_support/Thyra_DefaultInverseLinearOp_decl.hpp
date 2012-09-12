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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_INVERSE_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_INVERSE_LINEAR_OP_DECL_HPP

#include "Thyra_InverseLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief . */
enum EThrowOnSolveFailure {
  THROW_ON_SOLVE_FAILURE=1 ///< Throw an exception if a solve fails to converge
  ,IGNORE_SOLVE_FAILURE=0  ///< Don't throw an exception if a solve fails to converge
};


/** \brief Concrete <tt>LinearOpBase</tt> subclass that creates an implicit
 * <tt>LinearOpBase</tt> object using the inverse action of a
 * <tt>LinearOpWithSolveBase</tt> object.
 *
 * This class represents an implicit inverse linear operator:

 \verbatim
 
 M = inv(A)
 
 \endverbatim
 
 * where <tt>A</tt> is any <tt>LinearOpWithSolveBase</tt> object.
 * Specifically, the <tt>solve(...)</tt> function <tt>A</tt> is used to
 * implement <tt>this->apply()</tt> and the <tt>solveTranspose(...)</tt>
 * function <tt>A</tt> is used to implement <tt>this->applyTranspose()</tt>.
 *
 * <tt>SolveCriteria</tt> objects can be associated with <tt>A</tt> to define
 * the solve criterion for calling the <tt>A.solve(...,fwdSolveCriteria)</tt>
 * and <tt>A.solveTranspose(...,adjSolveCriteria)</tt>.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultInverseLinearOp : virtual public InverseLinearOpBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized (see postconditions for
   * <tt>uninitialize()</tt>).
   */
  DefaultInverseLinearOp();

  /** Calls <tt>initialize()</tt>.
   */
  DefaultInverseLinearOp(
    const RCP<LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
    const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
    );

  /** Calls <tt>initialize()</tt>.
   *
   * Rather than calling this constructor directly, consider using the non-member helper
   * functions described \ref Thyra_Op_Vec_AddedLinearOp_helpers_grp "here".
   */
  DefaultInverseLinearOp(
    const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
    const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnAdjSolveFailure  = THROW_ON_SOLVE_FAILURE
    );

  /** \brief Initialize given a non-const <tt>LinearOpWithSolveBase</tt>
   * object and an optional <tt>.
   *
   * \param lows [in] The <tt>LinearOpWithSolveBase</tt> object that will
   * <tt>solve(...)</tt> and/or <tt>solveTranspose(...)</tt> will be called
   * on.  Note that <tt>*this</tt> may give up non-const views of
   * <tt>*lows</tt> so that <tt>*lows</tt> may be changed through clients of
   * this object.
   *
   * \param fwdSolveCriteria [in] The criteria used to call
   * <tt>lows->solve(...)</tt>.  If <tt>fwdSolveCriteria==NULL</tt> then the
   * default solve criteria built into <tt>*lows<tt> will be used.  If
   * <tt>fwdSolveCriteria!=NULL</tt> then <tt>*fwdSolveCriteria</tt> will be
   * copied internally.  <b>Warning!</b> If shallow copy is used by any
   * parameters in <tt>fwdSolveCriteria->extraParameter</tt> these these
   * parameters will be "remembered" by <tt>*this</tt>.
   *
   * \param adjSolveCriteria [in] The criteria used to call
   * <tt>lows->solveTranspose(...)</tt>.  If <tt>adjSolveCriteria==NULL</tt>
   * then the default solve criteria built into <tt>*lows<tt> will be used.
   * If <tt>adjSolveCriteria!=NULL</tt> then <tt>*adjSolveCriteria</tt> will
   * be copied internally.  <b>Warning!</b> If shallow copy is used by any
   * parameters in <tt>adjSolveCriteria->extraParameter</tt> these these
   * parameters will be "remembered" by <tt>*this</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>lows.get() != NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->domain().get() == lows->range().get()</tt>
   * <li><tt>this->range().get() == lows->domain().get()</tt>
   * </ul>
   */
  void initialize(
    const RCP<LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
    const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
    );

  /** \brief Initialize given a non-const <tt>LinearOpWithSolveBase</tt>
   * object and an optional <tt>.
   *
   * \param lows [in] The <tt>LinearOpWithSolveBase</tt> object that will
   * <tt>solve(...)</tt> and/or <tt>solveTranspose(...)</tt> will be called
   * on.  Note that <tt>*this</tt> may give up non-const views of
   * <tt>*lows</tt> so that <tt>*lows</tt> may be changed through clients of
   * this object.
   *
   * \param fwdSolveCriteria [in] The criteria used to call
   * <tt>lows->solve(...)</tt>.  If <tt>fwdSolveCriteria==NULL</tt> then the
   * default solve criteria built into <tt>*lows<tt> will be used.  If
   * <tt>fwdSolveCriteria!=NULL</tt> then <tt>*fwdSolveCriteria</tt> will be
   * copied internally.  <b>Warning!</b> If shallow copy is used by any
   * parameters in <tt>fwdSolveCriteria->extraParameter</tt> these these
   * parameters will be "remembered" by <tt>*this</tt>.
   *
   * \param adjSolveCriteria [in] The criteria used to call
   * <tt>lows->solveTranspose(...)</tt>.  If <tt>adjSolveCriteria==NULL</tt>
   * then the default solve criteria built into <tt>*lows<tt> will be used.
   * If <tt>adjSolveCriteria!=NULL</tt> then <tt>*adjSolveCriteria</tt> will
   * be copied internally.  <b>Warning!</b> If shallow copy is used by any
   * parameters in <tt>adjSolveCriteria->extraParameter</tt> these these
   * parameters will be "remembered" by <tt>*this</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>lows.get() != NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->domain().get() == lows->range().get()</tt>
   * <li><tt>this->range().get() == lows->domain().get()</tt>
   * </ul>
   */
  void initialize(
    const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
    const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
    );

  /** \brief Set to uninitialized.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getLows().get()==NULL</tt>
   * <li><tt>this->range().get()==NULL</tt>
   * <li><tt>this->domain().get()==NULL</tt>
   * </ul>
   */
  void uninitialize();

  //@}

  /** @name Overridden from InverseLinearOpBase */
  //@{

  /** \brief . */
  bool isLowsConst() const;
  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> >
  getNonconstLows(); 
  /** \brief . */
  RCP<const LinearOpWithSolveBase<Scalar> >
  getLows() const; 

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Returns <tt>this->getLows()->domain() if
   * <t>this->getLows().get()!=NULL</tt> and returns <tt>Teuchos::null</tt>
   * otherwise.
   */
  RCP< const VectorSpaceBase<Scalar> > range() const;

  /** \brief Returns <tt>this->getLows()->range() if
   * <t>this->getLows().get()!=NULL</tt> and returns <tt>Teuchos::null</tt>
   * otherwise.
   */
  RCP< const VectorSpaceBase<Scalar> > domain() const;

  /** \brief . */
  RCP<const LinearOpBase<Scalar> > clone() const;

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Returns <tt>true</tt> only if all constituent operators support
   * <tt>M_trans</tt>.
   */
  bool opSupportedImpl(EOpTransp M_trans) const;

  /** \brief . */
  void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:

  Teuchos::ConstNonconstObjectContainer<LinearOpWithSolveBase<Scalar> > lows_;
  RCP<SolveCriteria<Scalar> > fwdSolveCriteria_;
  EThrowOnSolveFailure throwOnFwdSolveFailure_;
  RCP<SolveCriteria<Scalar> > adjSolveCriteria_;
  EThrowOnSolveFailure throwOnAdjSolveFailure_;
  
  void assertInitialized() const;

  template<class LOWS>
  void initializeImpl(
    const RCP<LOWS> &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const SolveCriteria<Scalar> *adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure
    );

  // Not defined and not to be called
  DefaultInverseLinearOp(const DefaultInverseLinearOp&);
  DefaultInverseLinearOp& operator=(const DefaultInverseLinearOp&);

};


/** \brief Form a non-const implicit inverse operator <tt>M = inv(A)</tt>.
 *
 * \relates DefaultInverseLinearOp
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
nonconstInverse(
  const RCP<LinearOpWithSolveBase<Scalar> > &A,
  const Ptr<const SolveCriteria<Scalar> > &fwdSolveCriteria = Teuchos::null,
  const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
  const Ptr<const SolveCriteria<Scalar> > &adjSolveCriteria = Teuchos::null,
  const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  );


/** \brief Form a const implicit inverse operator <tt>M = inv(A)</tt>.
 *
 * \relates DefaultInverseLinearOp
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
inverse(
  const RCP<const LinearOpWithSolveBase<Scalar> > &A,
  const Ptr<const SolveCriteria<Scalar> > &fwdSolveCriteria = Teuchos::null,
  const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
  const Ptr<const SolveCriteria<Scalar> > &adjSolveCriteria = Teuchos::null,
  const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  );


// /////////////////////////////////
// Inline members


template<class Scalar>
inline
void DefaultInverseLinearOp<Scalar>::assertInitialized() const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !lows_.getConstObj().get() );
#endif
}


#ifndef THYRA_HIDE_DEPRECATED_CODE
//
// Deprecated
//


/** \brief Deprecated
 *
 * \ingroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED RCP<LinearOpBase<Scalar> >
nonconstInverse(
  const RCP<LinearOpWithSolveBase<Scalar> > &A,
  const SolveCriteria<Scalar> *fwdSolveCriteria,
  const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
  const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
  const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  )
{
  using Teuchos::ptr;
  return nonconstInverse<Scalar>(A, ptr(fwdSolveCriteria), throwOnFwdSolveFailure,
    ptr(adjSolveCriteria), throwOnAdjSolveFailure);
}


/** \brief Deprecated
 *
 * \ingroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED RCP<LinearOpBase<Scalar> >
inverse(
  const RCP<const LinearOpWithSolveBase<Scalar> > &A,
  const SolveCriteria<Scalar> *fwdSolveCriteria,
  const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
  const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
  const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  )
{
  using Teuchos::ptr;
  return inverse<Scalar>(A, ptr(fwdSolveCriteria), throwOnFwdSolveFailure,
    ptr(adjSolveCriteria), throwOnAdjSolveFailure);
}


#endif // THYRA_HIDE_DEPRECATED_CODE
} // end namespace Thyra


#endif	// THYRA_DEFAULT_INVERSE_LINEAR_OP_DECL_HPP
