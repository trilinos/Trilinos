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

#ifndef THYRA_NONLINEAR_SOLVER_BUILDER_BASE
#define THYRA_NONLINEAR_SOLVER_BUILDER_BASE

#include "Thyra_NonlinearSolverBase.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"


namespace Thyra {


/** \brief Abstract interface for an object that can create
 * <tt>NonlinearSolverBase</tt> objects on demand.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NonlinearSolverBuilderBase : virtual public Teuchos::ParameterListAcceptor
{
public:
  
  /** \brief Create a new <tt>NonlinearSolverBase</tt> object purely
   * specified by the parameter list.
   *
   * \param nonlinearSolverTypeName [in] The optional name of the nonlinear
   * solver strategy to create.  The most typical use case will pass in an
   * empty string but there will be cases where I client will want to try to
   * select a specific linear solver strategy, even for just testing purposes.
   *
   * This function is given no information about the nature of the nonlinear
   * problem to be solved.
   */
  virtual Teuchos::RCP<NonlinearSolverBase<Scalar> >
  createNonlinearSolver(const std::string &nonlinearSolverTypeName) const = 0;

private:
  
  // Not defined and not to be called
  NonlinearSolverBuilderBase<Scalar>&
  operator=(const NonlinearSolverBuilderBase<Scalar>&);

};


/** \brief .
 *
 * \relates NonlinearSolverBuilderBase
 */
template<class Scalar>
Teuchos::RCP<NonlinearSolverBase<Scalar> >
createNonlinearSolver(
  const NonlinearSolverBuilderBase<Scalar> &nonlinearSolverBuilder,
  const std::string &nonlinearSolverTypeName = ""
  )
{
  return nonlinearSolverBuilder.createNonlinearSolver(nonlinearSolverTypeName);
}


} // namespace Thyra


#endif // THYRA_NONLINEAR_SOLVER_BUILDER_BASE
