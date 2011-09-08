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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_MULTI_VECTOR_LINEAR_OP_WITH_SOLVE_DECL_HPP


#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Implicit concrete <tt>LinearOpWithSolveBase</tt> subclass that
 * takes a flattended out multi-vector and performs a multi-RHS solve with it.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class DefaultMultiVectorLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  DefaultMultiVectorLinearOpWithSolve();

  /** \brief . */
  void nonconstInitialize(
    const RCP<LinearOpWithSolveBase<Scalar> > &lows,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
    );

  /** \brief . */
  void initialize(
    const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
    );
  
  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> >
  getNonconstLinearOpWithSolve();
  
  /** \brief . */
  RCP<const LinearOpWithSolveBase<Scalar> >
  getLinearOpWithSolve() const;

  // 2007/05/24: rabartl: ToDo: Add a const version of the above function once
  // needed

  /** \brief . */
  void uninitialize();

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > range() const;

  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > domain() const;

  /** \brief . */
  RCP<const LinearOpBase<Scalar> > clone() const;

  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief . */
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

  /** @name Overridden from LinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const;
  /** \brief . */
  SolveStatus<Scalar> solveImpl(
    const EOpTransp transp,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const;
  //@}

private:

  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpWithSolveBase<Scalar> > CNLOWS;

  // //////////////////////////////
  // Private data members

  CNLOWS lows_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecRange_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecDomain_;

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
    const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
    );
  

};


/** \brief Nonmember constructor function.
 *
 * \relates DefaultMultiVectorLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> >
multiVectorLinearOpWithSolve()
{
  return Teuchos::rcp(new DefaultMultiVectorLinearOpWithSolve<Scalar>());
}


/** \brief Nonmember constructor function.
 *
 * \relates DefaultMultiVectorLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> >
nonconstMultiVectorLinearOpWithSolve(
  const RCP<LinearOpWithSolveBase<Scalar> > &lows,
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
  )
{
  RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> >
    mvlows = Teuchos::rcp(new DefaultMultiVectorLinearOpWithSolve<Scalar>());
  mvlows->initialize(lows,multiVecRange,multiVecDomain);
  return mvlows;
}


/** \brief Nonmember constructor function.
 *
 * \relates DefaultMultiVectorLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> >
multiVectorLinearOpWithSolve(
  const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecRange,
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > &multiVecDomain
  )
{
  RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> >
    mvlows = Teuchos::rcp(new DefaultMultiVectorLinearOpWithSolve<Scalar>());
  mvlows->initialize(lows,multiVecRange,multiVecDomain);
  return mvlows;
}


}	// end namespace Thyra


#endif	// THYRA_MULTI_VECTOR_LINEAR_OP_WITH_SOLVE_DECL_HPP
