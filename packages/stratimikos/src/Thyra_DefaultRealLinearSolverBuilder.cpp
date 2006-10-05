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

//#define THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDER_DUMP

#include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"

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
#ifdef HAVE_STRATIMIKOS_ML_THYRA
#  include "Thyra_MLPreconditionerFactory.hpp"
#endif

namespace {

const std::string LinearSolverType_name    = "Linear Solver Type";
const std::string LinearSolverTypes_name   = "Linear Solver Types";
const std::string PreconditionerType_name    = "Preconditioner Type";
const std::string PreconditionerTypes_name   = "Preconditioner Types";

} // namespace 

namespace Thyra {

// Constructors/Initializers/Accessors

DefaultRealLinearSolverBuilder::DefaultRealLinearSolverBuilder(
  const std::string    &paramsXmlFileName
  ,const std::string   &extraParamsXmlString
  ,const std::string   &paramsUsedXmlOutFileName
  ,const std::string   &paramsXmlFileNameOption
  ,const std::string   &extraParamsXmlStringOption
  ,const std::string   &paramsUsedXmlOutFileNameOption
  )
  :paramsXmlFileName_(paramsXmlFileName)
   ,extraParamsXmlString_(extraParamsXmlString)
   ,paramsUsedXmlOutFileName_(paramsUsedXmlOutFileName)
   ,paramsXmlFileNameOption_(paramsXmlFileNameOption)
   ,extraParamsXmlStringOption_(extraParamsXmlStringOption)
   ,paramsUsedXmlOutFileNameOption_(paramsUsedXmlOutFileNameOption)
{
  this->initializeDefaults();
}

DefaultRealLinearSolverBuilder::~DefaultRealLinearSolverBuilder()
{
#ifdef TEUCHOS_DEBUG
  // Validate that we read the parameters correctly!
  if(paramList_.get())
    paramList_->validateParameters(*this->getValidParameters(),0);
#endif    
}

void DefaultRealLinearSolverBuilder::setLinearSolveStrategyFactory(
  const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolveFactoryBase<double> > >  &solveStrategyFactory
  ,const std::string                                                                                  &solveStrategyName
  )
{
  lowsf_map_[solveStrategyName] = solveStrategyFactory;
  validLowsfNames_.push_back(solveStrategyName);
  defaultLOWSF_ = solveStrategyName;
  validParamList_ = Teuchos::null;
}

void DefaultRealLinearSolverBuilder::setPreconditioningStrategyFactory(
  const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<PreconditionerFactoryBase<double> > >     &precStrategyFactory
  ,const std::string                                                                                  &precStrategyName
  )
{
  pf_map_[precStrategyName] = precStrategyFactory;
  validPfNames_.push_back(precStrategyName);
  defaultPF_ = precStrategyName;
  validParamList_ = Teuchos::null;
}

void DefaultRealLinearSolverBuilder::setupCLP( Teuchos::CommandLineProcessor *clp )
{
  TEST_FOR_EXCEPT(clp==NULL);
  clp->setOption(
    paramsXmlFileNameOption().c_str(),&paramsXmlFileName_
    ,"Name of an XML file containing parameters for linear solver options to be appended first."
    );
  clp->setOption(
    extraParamsXmlStringOption().c_str(),&extraParamsXmlString_
    ,"An XML string containing linear solver parameters to be appended second."
    );
  clp->setOption(
    paramsUsedXmlOutFileNameOption().c_str(),&paramsUsedXmlOutFileName_
    ,"Name of an XML file that can be written with the parameter list after it has been used on completion of this program."
    );
}

void DefaultRealLinearSolverBuilder::readParameters( std::ostream *out )
{
  if(!paramList_.get())
    paramList_ = Teuchos::rcp(new Teuchos::ParameterList("DefaultRealLinearSolverBuilder"));
  if(paramsXmlFileName().length()) {
    if(out) *out << "\nReading parameters from XML file \""<<paramsXmlFileName()<<"\" ...\n";
    Teuchos::updateParametersFromXmlFile(paramsXmlFileName(),&*paramList_);
  }
  if(extraParamsXmlString().length()) {
    if(out) *out << "\nAppending extra parameters from the XML string \""<<extraParamsXmlString()<<"\" ...\n";
    Teuchos::updateParametersFromXmlString(extraParamsXmlString(),&*paramList_);
  }
}

void DefaultRealLinearSolverBuilder::writeParamsFile(
  const LinearOpWithSolveFactoryBase<double>   &lowsFactory
  ,const std::string                           &outputXmlFileName
  ) const
{
  TEST_FOR_EXCEPT(!paramList_.get());
  std::string xmlOutputFile
    = ( outputXmlFileName.length() ? outputXmlFileName : paramsUsedXmlOutFileName() );
  if(xmlOutputFile.length()) {
    Teuchos::writeParameterListToXmlFile(*paramList_,xmlOutputFile);
  }
}

std::string
DefaultRealLinearSolverBuilder::getLinearSolveStrategyName() const
{
  TEST_FOR_EXCEPT(!paramList_.get());
  return paramList_->get(LinearSolverType_name,defaultLOWSF_);
}

std::string
DefaultRealLinearSolverBuilder::getPreconditionerStrategyName() const
{
  TEST_FOR_EXCEPT(!paramList_.get());
  return paramList_->get(PreconditionerType_name,defaultPF_);
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
    // Linear Solver Types
    validParamList->set(LinearSolverType_name,defaultLOWSF_);
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      linearSolverTypesSL = sublist(validParamList,LinearSolverTypes_name);
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
      linearSolverTypesSL->sublist(name).setParameters(*lowsf->getValidParameters());
    }
    // Preconditioner Type
    validParamList->set(PreconditionerType_name,defaultPF_);
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      precTypesSL = sublist(validParamList,PreconditionerTypes_name);
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
      precTypesSL->sublist(name).setParameters(*pf->getValidParameters());
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
#ifdef THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDER_DUMP
  std::cout << "\nEntering DefaultRealLinearSolverBuilder::createLinearSolveStrategy(...) ...\n";
  std::cout << "\nlinearSolveStrategyName = \"" << linearSolveStrategyName << "\"\n";
  std::cout << "\nlinearSolveStrategyName.length() = " << linearSolveStrategyName.length() << "\n";
  std::cout << "\ndefaultLOWSF_ = \"" << defaultLOWSF_ << "\"\n";
  std::cout << "\nthis->getLinearSolveStrategyName() = \"" << this->getLinearSolveStrategyName() << "\"\n";
#endif
  const std::string
    name = ( linearSolveStrategyName.length()
             ? linearSolveStrategyName
             : this->getLinearSolveStrategyName() );
#ifdef THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDER_DUMP
  std::cout << "\nname = \"" << name << "\"\n";
#endif
  // Validate that the linear solver strategy with this name exists
  lowsf_map_t::const_iterator itr = lowsf_map_.find(name);
  TEST_FOR_EXCEPTION(
    itr == lowsf_map_.end(), std::invalid_argument
    ,"Error, the value \""<<LinearSolverType_name<<"\"=\""<<name<<"\" is not a valid"
    " linear solver type.  Valid linear solve strategy names include "
    <<validLinearSolveStrategyNames()<<"!"
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
  lowsf->setParameterList(sublist(sublist(paramList_,LinearSolverTypes_name),name));
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
    name = ( preconditioningStrategyName.length()
             ? preconditioningStrategyName
             : this->getPreconditionerStrategyName() );
  Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
    pf = Teuchos::null;
  if( name != "None" ) {
    pf_map_t::const_iterator itr = pf_map_.find(name);
    TEST_FOR_EXCEPTION(
      itr == pf_map_.end(), std::invalid_argument
      ,"Error, the value \""<<PreconditionerType_name<<"\"=\""<<name<<"\" is not a valid"
      " preconditioner type.  Valid preconditioning strategy names include "
      <<validPreconditioningStrategyNames()<<"!"
      );
    Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
      pf = itr->second->create();
    pf->setParameterList(sublist(sublist(paramList_,PreconditionerTypes_name),name));
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
#ifdef HAVE_STRATIMIKOS_AMESOS_THYRA
  setLinearSolveStrategyFactory(
    rcp(new AbstractFactoryStd<LinearOpWithSolveFactoryBase<double>,AmesosLinearOpWithSolveFactory>())
    ,"Amesos"
    );
#endif
#ifdef HAVE_STRATIMIKOS_AZTECOO_THYRA
  setLinearSolveStrategyFactory(
    rcp(new AbstractFactoryStd<LinearOpWithSolveFactoryBase<double>,AztecOOLinearOpWithSolveFactory>())
    ,"AztecOO"
    );
#endif
#ifdef HAVE_STRATIMIKOS_AMESOS_THYRA
  if( Teuchos::GlobalMPISession::getNProc() == 1 ) {
    defaultLOWSF_ = "Amesos";
  }
#endif
  // Preconditioners
#ifdef HAVE_STRATIMIKOS_IFPACK_THYRA
  setPreconditioningStrategyFactory(
    rcp(new AbstractFactoryStd<PreconditionerFactoryBase<double>,IfpackPreconditionerFactory>())
    ,"Ifpack"
    );
#endif
#ifdef HAVE_STRATIMIKOS_ML_THYRA
  setPreconditioningStrategyFactory(
    rcp(new AbstractFactoryStd<PreconditionerFactoryBase<double>,MLPreconditionerFactory>())
    ,"ML"
    );
#endif
}

std::string
DefaultRealLinearSolverBuilder::validLinearSolveStrategyNames() const
{
  std::ostringstream oss;
  oss << "{";
  for( int i = 0; i < int(validLowsfNames_.size()); ++i ) {
    oss << "\"" << validLowsfNames_[i] << "\"";
    if( i != int(validLowsfNames_.size()-1) )
      oss << ",";
  }
  oss << "}";
  return oss.str();
}

std::string
DefaultRealLinearSolverBuilder::validPreconditioningStrategyNames() const
{
  std::ostringstream oss;
  oss << "{";
  for( int i = 0; i < int(validPfNames_.size()); ++i ) {
    oss << "\"" << validPfNames_[i] << "\"";
    if( i != int(validPfNames_.size()-1) )
      oss << ",";
  }
  oss << "}";
  return oss.str();
}

} // namespace Thyra
