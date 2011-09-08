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
