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
#include "Teuchos_StandardParameterEntryValidators.hpp"

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
const std::string None_name = "None";

Teuchos::RefCountPtr<const Teuchos::StringToIntegralParameterEntryValidator<int> >
lowsfValidator;

Teuchos::RefCountPtr<const Teuchos::StringToIntegralParameterEntryValidator<int> >
pfValidator;

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
    paramList_->validateParameters(*this->getValidParameters(),1);
#endif    
}

void DefaultRealLinearSolverBuilder::setLinearSolveStrategyFactory(
  const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<LinearOpWithSolveFactoryBase<double> > >  &solveStrategyFactory
  ,const std::string                                                                                  &solveStrategyName
  )
{
  validLowsfNames_.push_back(solveStrategyName);
  lowsfArray_.push_back(solveStrategyFactory);
  defaultLOWSF_ = solveStrategyName;
  validParamList_ = Teuchos::null;
}

void DefaultRealLinearSolverBuilder::setPreconditioningStrategyFactory(
  const Teuchos::RefCountPtr<const Teuchos::AbstractFactory<PreconditionerFactoryBase<double> > >     &precStrategyFactory
  ,const std::string                                                                                  &precStrategyName
  )
{
  validPfNames_.push_back(precStrategyName);
  pfArray_.push_back(precStrategyFactory);
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
  if(!lowsfValidator.get())
    this->getValidParameters(); // Make sure lowsfValidator has been initialized!
  return lowsfValidator->getStringValue(*paramList_,LinearSolverType_name,defaultLOWSF_);
}

std::string
DefaultRealLinearSolverBuilder::getPreconditionerStrategyName() const
{
  TEST_FOR_EXCEPT(!paramList_.get());
  if(!pfValidator.get())
    this->getValidParameters(); // Make sure pfValidator has been initialized!
  return pfValidator->getStringValue(*paramList_,PreconditionerType_name,defaultPF_);
}

// Overridden from ParameterListAcceptor

void DefaultRealLinearSolverBuilder::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(!paramList.get());
  // Only validate the zeroth and first level of parameters and sublists as
  // these are all that this class directly controls.  All other parameters
  // and sublusts are handed off to different LOWSFB and PFB objects.
  paramList->validateParameters(*this->getValidParameters(),1);
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
    lowsfValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        validLowsfNames_,LinearSolverType_name
        )
      );
    validParamList->set(
      LinearSolverType_name,defaultLOWSF_
      ,(std::string("Determines the type of linear solver that will be used.\n")
        + "The parameters for each solver type are specified in the sublist \""
        + LinearSolverTypes_name + "\"").c_str()
      ,lowsfValidator
      );
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      linearSolverTypesSL = sublist(validParamList,LinearSolverTypes_name);
    for( int i = 0; i < static_cast<int>(lowsfArray_.size()); ++i ) {
      const std::string
        &lsname = validLowsfNames_[i];
      const Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
        lowsf = lowsfArray_[i]->create();
      linearSolverTypesSL->sublist(lsname).setParameters(*lowsf->getValidParameters());
    }
    // Preconditioner Type
    pfValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        validPfNames_,PreconditionerType_name
        )
      );
    validParamList->set(
      PreconditionerType_name,defaultPF_
      ,(std::string("Determines the type of preconditioner that will be used.\n")
        + "This option is only meaningful for linear solvers that accept preconditioner"
        + " factory objects!\n"
        + "The parameters for each preconditioner are specified in the sublist \""
        + PreconditionerTypes_name + "\"").c_str()
      ,pfValidator
      );
    Teuchos::RefCountPtr<Teuchos::ParameterList>
      precTypesSL = sublist(validParamList,PreconditionerTypes_name);
    for( int i = 0; i < static_cast<int>(pfArray_.size()); ++i ) {
      const std::string
        &pfname = validPfNames_[i+1]; // "None" is the 0th entry!
      const Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
        pf = pfArray_[i]->create();
      precTypesSL->sublist(pfname).setParameters(*pf->getValidParameters());
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
    lsname = ( linearSolveStrategyName.length()
             ? linearSolveStrategyName
             : this->getLinearSolveStrategyName() );
#ifdef THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDER_DUMP
  std::cout << "\nlsname = \"" << lsname << "\"\n";
#endif
  // Get the index of this linear solver strategy (this will validate!)
  const int
    ls_idx = lowsfValidator->getIntegralValue(lsname,LinearSolverType_name);
  // Create the uninitialized LOWSFB object
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
    lowsf = lowsfArray_[ls_idx]->create();
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
  lowsf->setParameterList(sublist(sublist(paramList_,LinearSolverTypes_name),lsname));
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
    pfname = ( preconditioningStrategyName.length()
             ? preconditioningStrategyName
             : this->getPreconditionerStrategyName() );
  Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
    pf = Teuchos::null;
  // Get the index of this preconditioning strategy (this will validate!)
  const int
    pf_idx = pfValidator->getIntegralValue(pfname,PreconditionerType_name);
  if( pf_idx != 0 ) {
    Teuchos::RefCountPtr<PreconditionerFactoryBase<double> >
      pf = pfArray_[pf_idx-1]->create(); // We offset by -1 since "None" is first!
    pf->setParameterList(sublist(sublist(paramList_,PreconditionerTypes_name),pfname));
  }
  return pf;
}

// private

void DefaultRealLinearSolverBuilder::initializeDefaults()
{
  using Teuchos::rcp;
  using Teuchos::AbstractFactoryStd;
  defaultLOWSF_ = "";
  defaultPF_ = None_name;
  validLowsfNames_.resize(0);
  validPfNames_.resize(0);
  validPfNames_.push_back(None_name); // This will offset everything!
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
  // Note: the last LOWSF object set will be the default!
  //
  // Preconditioners
  //
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
  // Note: the last PF object set will be the default!
}

} // namespace Thyra
