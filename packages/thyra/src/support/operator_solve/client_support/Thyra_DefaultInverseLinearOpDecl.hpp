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

#ifndef THYRA_DEFAULT_INVERSE_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_INVERSE_LINEAR_OP_DECL_HPP

#include "Thyra_InverseLinearOpBase.hpp"
#include "Thyra_SingleScalarLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_arrayArg.hpp"

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
 * This class represents an implicitly inverse linear operator:

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
class DefaultInverseLinearOp
  : virtual public InverseLinearOpBase<Scalar>,            // Public interface
    virtual protected SingleScalarLinearOpBase<Scalar>     // Implementation detail
{
public:

  /** \brief . */
  using SingleScalarLinearOpBase<Scalar>::apply;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized (see postconditions for
   * <tt>uninitialize()</tt>).
   */
  DefaultInverseLinearOp();

  /** Calls <tt>initialize()</tt>.
   */
  DefaultInverseLinearOp(
    const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &lows,
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
    const Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
    const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnAdjSolveFailure  = THROW_ON_SOLVE_FAILURE
    );

  /** \brief Initialize given a non-const <tt>LinearOpWithSolveBase</tt>
   * object and an optional <tt>.
   *
   * \param  lows
   *           [in] The <tt>LinearOpWithSolveBase</tt> object that will
   *           <tt>solve(...)</tt> and/or <tt>solveTranspose(...)</tt> will be
   *           called on.  Note that <tt>*this</tt> may give up non-const
   *           views of <tt>*lows</tt> so that <tt>*lows</tt> may be changed
   *           through clients of this object.
   * \param  fwdSolveCriteria
   *           [in] The criteria used to call <tt>lows->solve(...)</tt>.  If
   *           <tt>fwdSolveCriteria==NULL</tt> then the default solve criteria
   *           built into <tt>*lows<tt> will be used.  If
   *           <tt>fwdSolveCriteria!=NULL</tt> then <tt>*fwdSolveCriteria</tt>
   *           will be copied internally.  <b>Warning!</b> If shallow copy is
   *           used by any parameters in
   *           <tt>fwdSolveCriteria->extraParameter</tt> these these
   *           parameters will be "remembered" by <tt>*this</tt>.

   * \param  adjSolveCriteria
   *           [in] The criteria used to call
   *           <tt>lows->solveTranspose(...)</tt>.  If
   *           <tt>adjSolveCriteria==NULL</tt> then the default solve criteria
   *           built into <tt>*lows<tt> will be used.  If
   *           <tt>adjSolveCriteria!=NULL</tt> then <tt>*adjSolveCriteria</tt>
   *           will be copied internally.  <b>Warning!</b> If shallow copy is
   *           used by any parameters in
   *           <tt>adjSolveCriteria->extraParameter</tt> these these
   *           parameters will be "remembered" by <tt>*this</tt>.
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
    const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
    const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
    const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
    );

  /** \brief Initialize given a non-const <tt>LinearOpWithSolveBase</tt>
   * object and an optional <tt>.
   *
   * \param  lows
   *           [in] The <tt>LinearOpWithSolveBase</tt> object that will
   *           <tt>solve(...)</tt> and/or <tt>solveTranspose(...)</tt> will be
   *           called on.  Note that <tt>*this</tt> may give up non-const
   *           views of <tt>*lows</tt> so that <tt>*lows</tt> may be changed
   *           through clients of this object.
   * \param  fwdSolveCriteria
   *           [in] The criteria used to call <tt>lows->solve(...)</tt>.  If
   *           <tt>fwdSolveCriteria==NULL</tt> then the default solve criteria
   *           built into <tt>*lows<tt> will be used.  If
   *           <tt>fwdSolveCriteria!=NULL</tt> then <tt>*fwdSolveCriteria</tt>
   *           will be copied internally.  <b>Warning!</b> If shallow copy is
   *           used by any parameters in
   *           <tt>fwdSolveCriteria->extraParameter</tt> these these
   *           parameters will be "remembered" by <tt>*this</tt>.
   * \param  adjSolveCriteria
   *           [in] The criteria used to call
   *           <tt>lows->solveTranspose(...)</tt>.  If
   *           <tt>adjSolveCriteria==NULL</tt> then the default solve criteria
   *           built into <tt>*lows<tt> will be used.  If
   *           <tt>adjSolveCriteria!=NULL</tt> then <tt>*adjSolveCriteria</tt>
   *           will be copied internally.  <b>Warning!</b> If shallow copy is
   *           used by any parameters in
   *           <tt>adjSolveCriteria->extraParameter</tt> these these
   *           parameters will be "remembered" by <tt>*this</tt>.
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
    const Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > &lows,
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
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
  getNonconstLows(); 
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> >
  getLows() const; 

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Returns <tt>this->getLows()->domain() if
   * <t>this->getLows().get()!=NULL</tt> and returns <tt>Teuchos::null</tt>
   * otherwise.
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;

  /** \brief Returns <tt>this->getLows()->range() if
   * <t>this->getLows().get()!=NULL</tt> and returns <tt>Teuchos::null</tt>
   * otherwise.
   */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{

  /** \brief Returns <tt>true</tt> only if all constituent operators support
   * <tt>M_trans</tt>.
   */
  bool opSupported(ETransp M_trans) const;

  /** \brief . */
  void apply(
    const ETransp M_trans,
    const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

private:

  Teuchos::ConstNonconstObjectContainer<LinearOpWithSolveBase<Scalar> > lows_;
  Teuchos::RefCountPtr<SolveCriteria<Scalar> > fwdSolveCriteria_;
  EThrowOnSolveFailure throwOnFwdSolveFailure_;
  Teuchos::RefCountPtr<SolveCriteria<Scalar> > adjSolveCriteria_;
  EThrowOnSolveFailure throwOnAdjSolveFailure_;
  
  void assertInitialized() const;

  template<class LOWS>
  void initializeImpl(
    const Teuchos::RefCountPtr<LOWS> &lows,
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
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
nonconstInverse(
 const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &A,
 const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
 const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
 const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
 const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  );


/** \brief Form a const implicit inverse operator <tt>M = inv(A)</tt>.
 *
 * \relates DefaultInverseLinearOp
 */
template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
inverse(
  const Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > &A,
  const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
  const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
  const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
  const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  );


// /////////////////////////////////
// Inline members


template<class Scalar>
inline
void DefaultInverseLinearOp<Scalar>::assertInitialized() const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !lows_.getConstObj().get() );
#endif
}


} // end namespace Thyra


#endif	// THYRA_DEFAULT_INVERSE_LINEAR_OP_DECL_HPP
