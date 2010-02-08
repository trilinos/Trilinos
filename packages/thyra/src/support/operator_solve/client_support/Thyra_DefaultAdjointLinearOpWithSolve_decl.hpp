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
