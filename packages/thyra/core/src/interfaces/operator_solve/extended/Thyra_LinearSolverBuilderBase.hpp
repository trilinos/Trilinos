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

#ifndef THYRA_LINEAR_SOLVER_BUILDING_BASE
#define THYRA_LINEAR_SOLVER_BUILDING_BASE

#include "Teuchos_ParameterListAcceptor.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"


namespace Thyra {


/** \brief Abstract interface for an object that can create
 * <tt>LinearOpWithSolveFactoryBase</tt> objects on demand.
 *
 * \ingroup Thyra_Op_Solve_extended_interfaces_code_grp
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class LinearSolverBuilderBase : virtual public Teuchos::ParameterListAcceptor
{
public:
  
  /** \brief Create a new <tt>LinearOpWithSolveFactoryBase</tt> object purely
   * specified by the parameter list.
   *
   * \param  linearSolveStrategyName
   *           [in] The optional name of the linear solve strategy to create.
   *           The most typical use case will pass in an empty string but
   *           there will be cases where I client will want to try to select
   *           a specific linear solver strategy, even for just testing purposes.
   *
   * This function is given no information about the nature of the linear
   * operators to be used.
   */
  virtual Teuchos::RCP<LinearOpWithSolveFactoryBase<Scalar> >
  createLinearSolveStrategy(
    const std::string &linearSolveStrategyName ) const = 0;
  
  /** \brief Create a new <tt>PreconditionerFactoryBase</tt> object purely
   * specified by the parameter list.
   *
   * \param  preconditioningStrategyName
   *           [in] The optional name of the preconditioning strategy to create.
   *           The most typical use case will pass in an empty string but
   *           there will be cases where I client will want to try to select
   *           a specific linear solver strategy, even for just testing purposes.
   *
   * This function is given no information about the nature of the linear
   * operators to be used.
   */
  virtual Teuchos::RCP<PreconditionerFactoryBase<Scalar> >
  createPreconditioningStrategy(
    const std::string &preconditioningStrategyName ) const = 0;

  /* \brief Create a new LinearOpWithSolveFactory object given a typical
   * forward linear operator and a typical solve criteria.
   *
   * \param  typicalFwdOp
   *           [in] A typical forward linear operator that represents the types of
   *           operator that will be used to solve linear system.
   * \param  typicalSolveCriteria
   *           [in] A typical solve criteria that will be used to solve for linear
   *           systems.
   * \param  typicalSolveUse
   *           [in] Determines how the solver will be used.
   * \param  solveStrategy
   *           [out] The LOWSF object that was determined to be the best suited for solving
   *           the typical system given above.
   * \param  initialLOWS
   *           [out] The LOWS object that was created that is consistent with the returned
   *           solve strategy.  If <tt>initialLOWS->get()==NULL</tt> on return then there is no
   *           such object returned.
   * \param  setupTime
   *           [out] The amount of time it took to setup the solver <tt>*initalLOWS</tt> before
   *           a solve was performed.
   * \param  solveTime
   *           [out] The amount of time it took to solve a typical linear system for the
   *           returned <tt>*initalLOWS</tt> object.
   *
   * ToDo: Finish documentation!
   */
  /*
    virtual void createSmartSolveStrategy(
    const Teuchos::RCP<LinearOpBase<Scalar>               &typicalFwdOp
    ,const SolveCritiera<Scalar>                                  &typicalSolveCriteria
    ,const ESupportSolveUse                                       &typicalSolveUse
    ,Teuchos::RCP<LinearOpWithSolveFactoryBase<Scalar> >  *solveStrategy
    ,Teuchos::RCP<Teuchos::ParameterList>                 *solveStrategyParameters
    ,Teuchos::RCP<LinearOpWithSolveBase<Scalar> >         *initialLOWS
    ,double                                                       *setupTime
    ,double                                                       *solveTime
    ) const = 0;
  */

private:
  
  // Not defined and not to be called
  LinearSolverBuilderBase<Scalar>&
  operator=(const LinearSolverBuilderBase<Scalar>&);

};


/** \brief .
 *
 * \relates LinearSolverBuilderBase
 */
template<class Scalar>
Teuchos::RCP<LinearOpWithSolveFactoryBase<Scalar> >
createLinearSolveStrategy(
  const LinearSolverBuilderBase<Scalar> &linearSolverBuilder,
  const std::string &linearSolveStrategyName = ""
  )
{
  return linearSolverBuilder.createLinearSolveStrategy(
    linearSolveStrategyName );
}


/** \brief .
 *
 * \relates LinearSolverBuilderBase
 */
template<class Scalar>
Teuchos::RCP<PreconditionerFactoryBase<Scalar> >
createPreconditioningStrategy(
  const LinearSolverBuilderBase<Scalar> &linearSolverBuilder,
  const std::string &preconditioningStrategyName = ""
  )
{
  return linearSolverBuilder.createPreconditioningStrategy(
    preconditioningStrategyName );
}


} // namespace Thyra

#endif // THYRA_LINEAR_SOLVER_BUILDING_BASE
