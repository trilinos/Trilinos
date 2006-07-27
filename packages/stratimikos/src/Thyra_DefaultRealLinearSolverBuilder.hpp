// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDING_BASE
#define THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDING_BASE

#include "Stratimikos_Config.h"
#include "Thyra_LinearSolverBuilderBase.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace Thyra {

/** \brief Concrete subclass of <tt>Thyra::LinearSolverBuilderBase</tt> for
 * creating <tt>LinearOpWithSolveFactoryBase</tt> objects on demand for
 * various Trilinos linear solver packages.
 *
 * The parameters this class accepts are shown below in human readable format
 * and in XML (i.e. machine readable) format.
 *
 * <b>Human readable format for valid parameters (with default values) accepted by this class</b>
 *
 * \verbinclude simple_stratimikos_example.options.readable.out
 *
 * <b>XML format for valid parameters (with default values) accepted by this class</b>
 *
 * \verbinclude simple_stratimikos_example.options.xml.out
 * 
 */
class DefaultRealLinearSolverBuilder : public LinearSolverBuilderBase<double>
{
public:

  /** @name Constructors/Initializers/Accessors */
  //@{

  /** \brief Construct without a parameter list. */
  DefaultRealLinearSolverBuilder();

  /** \brief Construct given a parameter list. */
  DefaultRealLinearSolverBuilder(
    Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
    );

  /** \brief Set a new linear solver strategy factory object. */
  void setLinearSolveStrategyFactory(
    const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolveFactoryBase<double> > >  &solveStrategyFactory
    ,const std::string                                                                                  &solveStrategyName
    );

  /** \brief Set a new preconditioner strategy factory object. */
  void setPreconditioningStrategyFactory(
    const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<PreconditionerFactoryBase<double> > >     &precStrategyFactory
    ,const std::string                                                                                  &precStrategyName
    );

  /** \brief Get the name of the linear solver strategy that will be created. */
  std::string getLinearSolveStrategyName() const;

  /** \brief Get the name of the preconditioner strategy that will be created. */
  std::string getPreconditionerStrategyName() const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}
  
  /** \name Overridden from LinearSolverBuilderBase. */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
  createLinearSolveStrategy(
    const std::string &linearSolveStrategyName
    ) const;
  /** \brief . */
  Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
  createPreconditioningStrategy(
    const std::string &preconditioningStrategyName
    ) const;

  //@}

private:

  // //////////////////////////////////////
  // Private types

  typedef std::map<std::string,Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolveFactoryBase<double> > > >  lowsf_map_t;
  typedef std::map<std::string,Teuchos::RefCountPtr<const Teuchos::AbstractFactory<PreconditionerFactoryBase<double> > > >     pf_map_t;

  // //////////////////////////////////////
  // Private data members
  
  Teuchos::RefCountPtr<Teuchos::ParameterList>                 paramList_;
  mutable Teuchos::RefCountPtr<const Teuchos::ParameterList>   validParamList_;
  lowsf_map_t                                                  lowsf_map_;
  std::vector<std::string>                                     validLowsfNames_;
  std::string                                                  defaultLOWSF_;
  pf_map_t                                                     pf_map_;
  std::vector<std::string>                                     validPfNames_;
  std::string                                                  defaultPF_;

  // //////////////////////////////////////
  // Private member functions

  void initializeDefaults();
  std::string validLinearSolveStrategyNames() const;
  std::string validPreconditioningStrategyNames() const;

};

} // namespace Thyra

#endif // THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDING_BASE
