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

#include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#ifdef HAVE_STRATIMIKOS_AMESOS_THYRA
#  include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_AZTECOO_THYRA
#  include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_BELOS_THYRA
#  include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_IFPACK_THYRA
#  include "Thyra_IfpackPreconditionerFactory.hpp"
#endif
#include "Teuchos_AbstractFactoryStd.hpp"

namespace {

const std::string LinearSolverType_name = "Linear Solver Type";
const std::string PreconditionerType_name = "Preconditioner Type";

} // namespace 

namespace Thyra {

// Constructors/Initializers/Accessors

DefaultRealLinearSolverBuilder::DefaultRealLinearSolverBuilder()
{
  this->initializeDefaults();
}

DefaultRealLinearSolverBuilder::DefaultRealLinearSolverBuilder(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  this->initializeDefaults();
  this->setParameterList(paramList);
}

void DefaultRealLinearSolverBuilder::setLinearSolveStrategyFactory(
  const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolveFactoryBase<double> > >  &solveStrategyFactory
  ,const std::string                                                                                  &solveStrategyName
  )
{
  lowsf_map_[solveStrategyName] = solveStrategyFactory;
  defaultLOWSF_ = solveStrategyName;
  validParamList_ = Teuchos::null;
}

void DefaultRealLinearSolverBuilder::setPreconditioningStrategyFactory(
  const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<PreconditionerFactoryBase<double> > >     &precStrategyFactory
  ,const std::string                                                                                  &precStrategyName
  )
{
  pf_map_[precStrategyName] = precStrategyFactory;
  defaultPF_ = precStrategyName;
  validParamList_ = Teuchos::null;
}


std::string
DefaultRealLinearSolverBuilder::getLinearSolveStrategyName() const
{
  TEST_FOR_EXCEPT(!paramList_.get());
  const std::string &name = paramList_->get(LinearSolverType_name,defaultLOWSF_);
  return name;
}

std::string
DefaultRealLinearSolverBuilder::getPreconditionerStrategyName() const
{
  TEST_FOR_EXCEPT(!paramList_.get());
  const std::string &name = paramList_->get(PreconditionerType_name,defaultPF_);
  return name;
}


// Overridden from ParameterListAcceptor

void DefaultRealLinearSolverBuilder::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(!paramList.get());
  // Only validate this level of parameters and sublists
  paramList->validateParameters(*this->getValidParameters(),0);
  paramList_ = paramList;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
DefaultRealLinearSolverBuilder::getParameterList()
{
  return paramList_;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
DefaultRealLinearSolverBuilder::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
DefaultRealLinearSolverBuilder::getParameterList() const
{
  return paramList_;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
DefaultRealLinearSolverBuilder::getValidParameters() const
{
  if(!validParamList_.get()) {
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      validParamList = Teuchos::rcp(new Teuchos::ParameterList);
    // Linear Solver Type
    validParamList->set(LinearSolverType_name,defaultLOWSF_);
    for(
      lowsf_map_t::const_iterator itr = lowsf_map_.begin();
      itr != lowsf_map_.end();
      ++itr
      )
    {
      const std::string
        &name = itr->first;
      const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
        lowsf = itr->second->create();
      validParamList->sublist(name).setParameters(*lowsf->getValidParameters());
    }
    // Preconditioner Type
    validParamList->set(PreconditionerType_name,defaultPF_);
    for(
      pf_map_t::const_iterator itr = pf_map_.begin();
      itr != pf_map_.end();
      ++itr
      )
    {
      const std::string
        &name = itr->first;
      const Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
        pf = itr->second->create();
      validParamList->sublist(name).setParameters(*pf->getValidParameters());
    }
    validParamList_ = validParamList;
  }
  return validParamList_;
}
  
// Overridden from LinearSolverBuilderBase.

Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
DefaultRealLinearSolverBuilder::createLinearSolveStrategy(
  const std::string &linearSolveStrategyName
  ) const
{
  // Get the name of the linear solve strategy
  const std::string
    &name = ( linearSolveStrategyName.length()
              ? linearSolveStrategyName
              : this->getLinearSolveStrategyName() );
  // Validate that the linear solver strategy with this name exists
  lowsf_map_t::const_iterator itr = lowsf_map_.find(name);
  TEST_FOR_EXCEPTION(
    itr == lowsf_map_.end(), std::invalid_argument
    ,"Error, the value \""<<LinearSolverType_name<<"\"=\""<<name<<"\" is not a valid"
    " linear solver type.  See the list of valid linear solver types from this->getValidParameters()."
    );
  // Create the uninitialized LOWSFB object
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
    lowsf = itr->second->create();
  // First, set the preconditioner factory and its parameters
  if(lowsf->acceptsPreconditionerFactory()) {
    const std::string &pfName = this->getPreconditionerStrategyName();
    Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
      pf = this->createPreconditioningStrategy(pfName);
    if(pf.get())
      lowsf->setPreconditionerFactory(pf,pfName);
  }
  // Now set the parameters for the linear solver (some of which might
  // override some preconditioner factory parameters).
  lowsf->setParameterList(sublist(paramList_,name));
  //
  return lowsf;
}

Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
DefaultRealLinearSolverBuilder::createPreconditioningStrategy(
  const std::string &preconditioningStrategyName
  ) const
{
  // Get the name of the preconditioning strategy
  const std::string
    &name = ( preconditioningStrategyName.length()
              ? preconditioningStrategyName
              : this->getPreconditionerStrategyName() );
  Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
    pf = Teuchos::null;
  if( name != "None" ) {
    pf_map_t::const_iterator itr = pf_map_.find(name);
    TEST_FOR_EXCEPTION(
      itr == pf_map_.end(), std::invalid_argument
      ,"Error, the value \""<<PreconditionerType_name<<"\"=\""<<name<<"\" is not a valid"
      " preconditioner type.  See the list of valid preconditioner types from this->getValidParameters()."
      );
    Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
      pf = itr->second->create();
    pf->setParameterList(sublist(paramList_,name));
  }
  return pf;
}

// private

void DefaultRealLinearSolverBuilder::initializeDefaults()
{
  using Teuchos::rcp;
  using Teuchos::AbstractFactoryStd;
  defaultLOWSF_ = "";
  defaultPF_ = "";
  // Solvers
#ifdef HAVE_STRATIMIKOS_BELOS_THYRA
  setLinearSolveStrategyFactory(
    rcp(new AbstractFactoryStd<LinearOpWithSolveFactoryBase<double>,BelosLinearOpWithSolveFactory<double> >())
    ,"Belos"
    );
#endif
#ifdef HAVE_STRATIMIKOS_AZTECOO_THYRA
  setLinearSolveStrategyFactory(
    rcp(new AbstractFactoryStd<LinearOpWithSolveFactoryBase<double>,AztecOOLinearOpWithSolveFactory>())
    ,"AztecOO"
    );
#endif
#ifdef HAVE_STRATIMIKOS_AMESOS_THYRA
  setLinearSolveStrategyFactory(
    rcp(new AbstractFactoryStd<LinearOpWithSolveFactoryBase<double>,AmesosLinearOpWithSolveFactory>())
    ,"Amesos"
    );
#endif
  // Preconditioners
#ifdef HAVE_STRATIMIKOS_IFPACK_THYRA
  setPreconditioningStrategyFactory(
    rcp(new AbstractFactoryStd<PreconditionerFactoryBase<double>,IfpackPreconditionerFactory>())
    ,"Ifpack"
    );
#endif
}

} // namespace Thyra
