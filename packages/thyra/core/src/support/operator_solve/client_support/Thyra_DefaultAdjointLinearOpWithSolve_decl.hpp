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


#ifndef THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_DECL_HPP


#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace Thyra {


/** \brief Default concreate decorator subclass for a transpose/adjoint
 * <tt>LinearOpWithSolveBase</tt> object.
 * 
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class DefaultAdjointLinearOpWithSolve : virtual public LinearOpWithSolveBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Constructs to uninitialized. */
  DefaultAdjointLinearOpWithSolve();

  /** \brief Initialize with non-const LOWSB . */
  void initialize( const RCP<LinearOpWithSolveBase<Scalar> > &lows,
    const EOpTransp transp );

  /** \brief Initialize with const LOWSB . */
  void initialize( const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
    const EOpTransp transp );

  /** \brief Get the non-const underlying LOWSB object. */
  const RCP<LinearOpWithSolveBase<Scalar> > getNonconstOp();

  /** \brief Get the const underlying LOWSB object. */
  const RCP<const LinearOpWithSolveBase<Scalar> > getOp() const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > domain() const;

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

  // //////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpWithSolveBase<Scalar> >
  CNCLOWS;

  // /////////////////////////
  // Private data members

  CNCLOWS lows_;
  EOpTransp transp_;

};


/** \brief Nonmember constructor.
 *
 * \brief DefaultAdjointLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultAdjointLinearOpWithSolve<Scalar> >
defaultAdjointLinearOpWithSolve(
  const RCP<const LinearOpWithSolveBase<Scalar> > &lows,
  const EOpTransp transp )
{
  RCP<DefaultAdjointLinearOpWithSolve<Scalar> >
    dalows = Teuchos::rcp(new DefaultAdjointLinearOpWithSolve<Scalar>);
  dalows->initialize(lows, transp);
  return dalows;
}


/** \brief Nonmember constructor.
 *
 * \brief DefaultAdjointLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultAdjointLinearOpWithSolve<Scalar> >
defaultAdjointLinearOpWithSolveNonconst(
  const RCP<LinearOpWithSolveBase<Scalar> > &lows,
  const EOpTransp transp )
{
  RCP<DefaultAdjointLinearOpWithSolve<Scalar> >
    dalows = Teuchos::rcp(new DefaultAdjointLinearOpWithSolve<Scalar>);
  dalows->initialize(lows, transp);
  return dalows;
}


/** \brief Nonmember constructor.
 *
 * \brief DefaultAdjointLinearOpWithSolve
 */
template<class Scalar>
RCP<const LinearOpWithSolveBase<Scalar> >
adjointLows( const RCP<const LinearOpWithSolveBase<Scalar> > &lows )
{
  return defaultAdjointLinearOpWithSolve<Scalar>(lows, CONJTRANS);
}


/** \brief Nonmember constructor.
 *
 * \brief DefaultAdjointLinearOpWithSolve
 */
template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
nonconstAdjointLows( const RCP<LinearOpWithSolveBase<Scalar> > &lows )
{
  return defaultAdjointLinearOpWithSolveNonconst<Scalar>(lows, CONJTRANS);
}



}	// end namespace Thyra


#endif	// THYRA_DEFAULT_ADJOINT_LINEAR_OP_WITH_SOLVE_DECL_HPP
