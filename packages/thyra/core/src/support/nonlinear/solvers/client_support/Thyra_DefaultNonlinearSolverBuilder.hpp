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

#ifndef THYRA_DEFAULT_NONLINEAR_SOLVER_BUILDER_HPP
#define THYRA_DEFAULT_NONLINEAR_SOLVER_BUILDER_HPP

#include "Thyra_NonlinearSolverBuilderBase.hpp"
#include "Teuchos_AbstractFactory.hpp"


namespace Thyra {


/** \brief Concrete subclass of <tt>Thyra::NonlinearSolverBuilderBase</tt> for
 * creating <tt>NonlinearSolverBase</tt> objects and
 * <tt>PreconditionerFactoryBase</tt> object on demand given configured
 * factory objects.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Nonlin_ME_solvers_grp
 */
class DefaultNonlinearSolverBuilder
  : public Thyra::NonlinearSolverBuilderBase<double>
{
public:

  /** @name Constructors/Initializers/Accessors */
  //@{

  /** \brief . */
  DefaultNonlinearSolverBuilder();

  /** \brief . */
  ~DefaultNonlinearSolverBuilder();

  /** \brief Set a new NonlinearSolverBase factory object. */
  void setNonlinearSolverFactory(
    const RCP<const AbstractFactory<Thyra::NonlinearSolverBase<double> > >
    &nonlinearSolverFactory,
    const std::string &nonlinearSolverTypeName
    );
  
  /** \brief Get the name of the NonlinearSolver type that will be created on
   * the next call to <tt>this->createNonlinearSolver()</tt>.
   */
  std::string getNonlinearSolverName() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}
  
  /** \name Overridden from NonlinearSolverBuilderBase. */
  //@{

  /** \brief . */
  virtual Teuchos::RCP<NonlinearSolverBase<Scalar> >
  createNonlinearSolver(const std::string &nonlinearSolverTypeName) const;

  //@}

private:

  // //////////////////////////////////////
  // Private types

  typedef RCP<const AbstractFactory<Thyra::NonlinearSolverBase<double> > >
  ns_fcty_t;

  // //////////////////////////////////////
  // Private data members
  
  RCP<ParameterList> paramList_;
  mutable RCP<const ParameterList> validParamList_;
  Array<std::string> validNonlinearSolverNames_;
  Array<ns_fcty_t> nonlinearSolverArray_;
  std::string defaultNonlinearSolverName_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults();

};


} // namespace Thyra


/** \brief Inject a new solver type into a DefaultNonlinearSolverBuilder
 * object.
 */
template<class NonlinearSolverType, class Scalar>
void setNonlinearSolverFactory(
  const std::string &nonlinearSolverTypeName,
  const Ptr<DefaultNonlinearSolverBuilder<Scalar> > &defaultNonlinearSolverBuilder
  )
{
  TEST_FOR_EXCEPT(true);
}


#endif // THYRA_DEFAULT_NONLINEAR_SOLVER_BUILDER_HPP
