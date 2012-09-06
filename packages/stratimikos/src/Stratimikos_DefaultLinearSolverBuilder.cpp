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

#include "Stratimikos_InternalConfig.h"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_DelayedLinearOpWithSolveFactory.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_STRATIMIKOS_AMESOS
#  include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_AZTECOO
#  include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_BELOS
#  include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_IFPACK
#  include "Thyra_IfpackPreconditionerFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_ML
#  include "Thyra_MLPreconditionerFactory.hpp"
#endif


namespace {


const std::string LinearSolverType_name    = "Linear Solver Type";
const std::string LinearSolverTypes_name   = "Linear Solver Types";
const std::string PreconditionerType_name    = "Preconditioner Type";
const std::string PreconditionerTypes_name   = "Preconditioner Types";
const std::string None_name = "None";
const std::string EnableDelayedSolverConstruction_name = "Enable Delayed Solver Construction";
const bool EnableDelayedSolverConstruction_default = false;


} // namespace 


namespace Stratimikos {


// Constructors/Initializers/Accessors


DefaultLinearSolverBuilder::DefaultLinearSolverBuilder(
  const std::string    &paramsXmlFileName_in
  ,const std::string   &extraParamsXmlString_in
  ,const std::string   &paramsUsedXmlOutFileName_in
  ,const std::string   &paramsXmlFileNameOption_in
  ,const std::string   &extraParamsXmlStringOption_in
  ,const std::string   &paramsUsedXmlOutFileNameOption_in
  )
  :paramsXmlFileName_(paramsXmlFileName_in)
  ,extraParamsXmlString_(extraParamsXmlString_in)
  ,paramsUsedXmlOutFileName_(paramsUsedXmlOutFileName_in)
  ,paramsXmlFileNameOption_(paramsXmlFileNameOption_in)
  ,extraParamsXmlStringOption_(extraParamsXmlStringOption_in)
  ,paramsUsedXmlOutFileNameOption_(paramsUsedXmlOutFileNameOption_in)
  ,enableDelayedSolverConstruction_(EnableDelayedSolverConstruction_default)
{
  this->initializeDefaults();
}


DefaultLinearSolverBuilder::~DefaultLinearSolverBuilder()
{
#ifdef TEUCHOS_DEBUG
  // Validate that we read the parameters correctly!
  if (nonnull(paramList_)) {
    paramList_->validateParameters(*this->getValidParameters());
  }
#endif    
}


void DefaultLinearSolverBuilder::setLinearSolveStrategyFactory(
  const RCP<const AbstractFactory<Thyra::LinearOpWithSolveFactoryBase<double> > >
  &solveStrategyFactory,
  const std::string &solveStrategyName,
  const bool makeDefault
  )
{
  validLowsfNames_.push_back(solveStrategyName);
  lowsfArray_.push_back(solveStrategyFactory);
  validParamList_ = Teuchos::null;
  if (makeDefault) {
    setDefaultLinearSolveStrategyFactoryName(solveStrategyName);
  }
}


void DefaultLinearSolverBuilder::setDefaultLinearSolveStrategyFactoryName(
  const std::string &solveStrategyName)
{
  defaultLOWSF_ = solveStrategyName;
}


void DefaultLinearSolverBuilder::setPreconditioningStrategyFactory(
  const RCP<const AbstractFactory<Thyra::PreconditionerFactoryBase<double> > >
  &precStrategyFactory,
  const std::string &precStrategyName,
  const bool makeDefault
  )
{
  validPfNames_.push_back(precStrategyName);
  pfArray_.push_back(precStrategyFactory);
  validParamList_ = Teuchos::null;
  if (makeDefault) {
    setDefaultPreconditioningStrategyFactoryName(precStrategyName);
  }
}


void DefaultLinearSolverBuilder::setDefaultPreconditioningStrategyFactoryName(
  const std::string &precStrategyName)
{
  defaultPF_ = precStrategyName;
}


void DefaultLinearSolverBuilder::setupCLP( Teuchos::CommandLineProcessor *clp )
{
  TEUCHOS_TEST_FOR_EXCEPT(clp==NULL);
  clp->setOption(
    paramsXmlFileNameOption().c_str(),&paramsXmlFileName_
    ,"Name of an XML file containing parameters for linear solver "
    "options to be appended first."
    );
  clp->setOption(
    extraParamsXmlStringOption().c_str(),&extraParamsXmlString_
    ,"An XML string containing linear solver parameters to be appended second."
    );
  clp->setOption(
    paramsUsedXmlOutFileNameOption().c_str(),&paramsUsedXmlOutFileName_
    ,"Name of an XML file that can be written with the parameter list after it "
    "has been used on completion of this program."
    );
}


void DefaultLinearSolverBuilder::readParameters( std::ostream *out )
{
  using Teuchos::parameterList;
  using Teuchos::ptr;
  using Teuchos::updateParametersFromXmlFile;
  using Teuchos::updateParametersFromXmlString;
  using std::endl;
  
  if (!paramList_.get()) {
    paramList_ = parameterList("DefaultLinearSolverBuilder");
  }
  if (paramsXmlFileName().length()) {
    if (out) {
      *out << endl << "Reading parameters from XML file \"" 
           << paramsXmlFileName() << "\" ..." << endl;
    }
    updateParametersFromXmlFile (paramsXmlFileName (), paramList_.ptr());
  }
  if (extraParamsXmlString().length()) {
    if (out) {
      *out << endl << "Appending extra parameters from the XML string \""
           << extraParamsXmlString() << "\" ..." << endl;
    }
    updateParametersFromXmlString (extraParamsXmlString (), paramList_.ptr());
  }
  setParameterList(paramList_);
}


void DefaultLinearSolverBuilder::writeParamsFile(
  const Thyra::LinearOpWithSolveFactoryBase<double> &lowsFactory,
  const std::string &outputXmlFileName
  ) const
{
  justInTimeInitialize();
  const std::string xmlOutputFile =
    ( outputXmlFileName.length() ? outputXmlFileName : paramsUsedXmlOutFileName() );
  if (xmlOutputFile.length()) {
    Teuchos::writeParameterListToXmlFile(*paramList_, xmlOutputFile);
  }
}


std::string
DefaultLinearSolverBuilder::getLinearSolveStrategyName() const
{
  justInTimeInitialize();
  return lowsfValidator_->getStringValue(*paramList_, LinearSolverType_name,
    defaultLOWSF_);
}


std::string
DefaultLinearSolverBuilder::getPreconditionerStrategyName() const
{
  justInTimeInitialize();
  return pfValidator_->getStringValue(*paramList_, PreconditionerType_name,
    defaultPF_);
}


// Overridden from ParameterListAcceptor


void DefaultLinearSolverBuilder::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*this->getValidParameters());
  paramList_ = paramList;
  enableDelayedSolverConstruction_ = paramList_->get(
    EnableDelayedSolverConstruction_name, EnableDelayedSolverConstruction_default );
}


RCP<Teuchos::ParameterList>
DefaultLinearSolverBuilder::getNonconstParameterList()
{
  return paramList_;
}


RCP<Teuchos::ParameterList>
DefaultLinearSolverBuilder::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


RCP<const Teuchos::ParameterList>
DefaultLinearSolverBuilder::getParameterList() const
{
  return paramList_;
}


RCP<const Teuchos::ParameterList>
DefaultLinearSolverBuilder::getValidParameters() const
{
  using Teuchos::rcp_implicit_cast;
  typedef Teuchos::ParameterEntryValidator PEV;
  if (is_null(validParamList_)) {
    RCP<Teuchos::ParameterList>
      validParamList = Teuchos::rcp(new Teuchos::ParameterList);
    // Linear Solver Types
    lowsfValidator_ = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        validLowsfNames_,LinearSolverType_name
        )
      );
    validParamList->set(
      LinearSolverType_name, defaultLOWSF_,
      (std::string("Determines the type of linear solver that will be used.\n")
        + "The parameters for each solver type are specified in the sublist \""
        + LinearSolverTypes_name + "\"").c_str(),
      rcp_implicit_cast<const PEV>(lowsfValidator_)
      );
    Teuchos::ParameterList &linearSolverTypesSL = validParamList->sublist(
      LinearSolverTypes_name,false,
      "Sublists for each of the linear solver types set using the parameter\n"
      "\"" + LinearSolverType_name + "\".  Note that the options for each\n"
      "linear solver type given below will only be used if linear solvers\n"
      "of that type are created.  It is fine to list parameter sublists for\n"
      "linear solver types that are not used."
      );
    for( int i = 0; i < static_cast<int>(lowsfArray_.size()); ++i ) {
      const std::string
        &lsname = validLowsfNames_[i];
      const RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
        lowsf = lowsfArray_[i]->create();
      linearSolverTypesSL.sublist(lsname).setParameters(*lowsf->getValidParameters()
        ).disableRecursiveValidation();
    }
    // Preconditioner Type
    pfValidator_ = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        validPfNames_, PreconditionerType_name ) );
    validParamList->set(
      PreconditionerType_name, defaultPF_,
      (std::string("Determines the type of preconditioner that will be used.\n")
        + "This option is only meaningful for linear solvers that accept preconditioner"
        + " factory objects!\n"
        + "The parameters for each preconditioner are specified in the sublist \""
        + PreconditionerTypes_name + "\"").c_str(),
      rcp_implicit_cast<const PEV>(pfValidator_)
      );
    Teuchos::ParameterList &precTypesSL = validParamList->sublist(
        PreconditionerTypes_name,false,
        "Sublists for each of the preconditioner types set using the parameter\n"
        "\"" + PreconditionerType_name + "\".  Note that the options for each\n"
        "preconditioner type given below will only be used if preconditioners\n"
        "of that type are created.  It is fine to list parameter sublists for\n"
        "preconditioner types that are not used."
        );
    for( int i = 0; i < static_cast<int>(pfArray_.size()); ++i ) {
      const std::string
        &pfname = validPfNames_[i+1]; // "None" is the 0th entry!
      const RCP<Thyra::PreconditionerFactoryBase<double> >
        pf = pfArray_[i]->create();
      precTypesSL.sublist(pfname).setParameters(*pf->getValidParameters()
        ).disableRecursiveValidation();
    }
    // 
    validParamList->set(
      EnableDelayedSolverConstruction_name, EnableDelayedSolverConstruction_default,
      "When this option is set to true, the linear solver factory will be wrapped\n"
      "in a delayed evaluation Decorator factory object.  This results in a delay\n"
      "in the creation of a linear solver (and the associated preconditioner) until\n"
      "the first solve is actually performed.  This helps in cases where it is not\n"
      "known a-priori if a linear solve will be needed on a given linear operator and\n"
      "therefore can significantly improve performance for some types of algorithms\n"
      "such as NOX and LOCA."
      );
    //
    validParamList_ = validParamList;
  }
  return validParamList_;
}

  
// Overridden from LinearSolverBuilderBase.


RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
DefaultLinearSolverBuilder::createLinearSolveStrategy(
  const std::string &linearSolveStrategyName
  ) const
{
  justInTimeInitialize();

  // Get the name of the linear solve strategy
#ifdef THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDER_DUMP
  std::cout << "\nEntering DefaultLinearSolverBuilder"
            << "::createLinearSolveStrategy(...) ...\n";
  std::cout << "\nlinearSolveStrategyName = \""
            << linearSolveStrategyName << "\"\n";
  std::cout << "\nlinearSolveStrategyName.length() = "
            << linearSolveStrategyName.length() << "\n";
  std::cout << "\ndefaultLOWSF_ = \"" << defaultLOWSF_ << "\"\n";
  std::cout << "\nthis->getLinearSolveStrategyName() = \""
            << this->getLinearSolveStrategyName() << "\"\n";
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
    ls_idx = lowsfValidator_->getIntegralValue(lsname, LinearSolverType_name);

  // Create the uninitialized LOWSFB object
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsf = lowsfArray_[ls_idx]->create();

  // First, set the preconditioner factory and its parameters
  if(lowsf->acceptsPreconditionerFactory()) {
    const std::string &pfName = this->getPreconditionerStrategyName();
    RCP<Thyra::PreconditionerFactoryBase<double> >
      pf = this->createPreconditioningStrategy(pfName);
    if(pf.get())
      lowsf->setPreconditionerFactory(pf,pfName);
  }

  // Now set the parameters for the linear solver (some of which might
  // override some preconditioner factory parameters).
  lowsf->setParameterList(
    sublist(sublist(paramList_, LinearSolverTypes_name), lsname));
  //
  if (enableDelayedSolverConstruction_) {
    return Teuchos::rcp(
      new Thyra::DelayedLinearOpWithSolveFactory<double>(lowsf)
      );
  }

  return lowsf;

}


RCP<Thyra::PreconditionerFactoryBase<double> >
DefaultLinearSolverBuilder::createPreconditioningStrategy(
  const std::string &preconditioningStrategyName
  ) const
{
  justInTimeInitialize();

  // Get the name of the preconditioning strategy
  const std::string
    pfname = ( preconditioningStrategyName.length()
             ? preconditioningStrategyName
             : this->getPreconditionerStrategyName() );
  RCP<Thyra::PreconditionerFactoryBase<double> >
    pf = Teuchos::null;

  // Get the index of this preconditioning strategy (this will validate!)
  const int
    pf_idx = pfValidator_->getIntegralValue(pfname, PreconditionerType_name);
  if( pf_idx != 0 ) {
    pf = pfArray_[pf_idx-1]->create(); // We offset by -1 since "None" is first!
    pf->setParameterList(
      sublist(sublist(paramList_, PreconditionerTypes_name), pfname));
  }

  return pf;

}


// private


void DefaultLinearSolverBuilder::initializeDefaults()
{

  using Teuchos::rcp;
  using Teuchos::abstractFactoryStd;

  defaultLOWSF_ = "";
  defaultPF_ = None_name;
  validLowsfNames_.resize(0);
  validPfNames_.resize(0);
  validPfNames_.push_back(None_name); // This will offset everything!

  //
  // Linear Solvers
  //

#ifdef HAVE_STRATIMIKOS_BELOS
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<double>,
    Thyra::BelosLinearOpWithSolveFactory<double> >(),
    "Belos", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_AMESOS
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<double>,
    Thyra::AmesosLinearOpWithSolveFactory>(),
    "Amesos", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_AZTECOO
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<double>,
    Thyra::AztecOOLinearOpWithSolveFactory>(),
    "AztecOO", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_AMESOS
  if (Teuchos::GlobalMPISession::getNProc() == 1) {
    setDefaultLinearSolveStrategyFactoryName("Amesos");
  }
#endif

  // Note: ABove, the last LOWSF object set will be the default unless we are
  // on multiple processors!

  //
  // Preconditioners
  //

#ifdef HAVE_STRATIMIKOS_ML
  setPreconditioningStrategyFactory(
    abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,
    Thyra::MLPreconditionerFactory>(),
    "ML", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_IFPACK
  setPreconditioningStrategyFactory(
    abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,
    Thyra::IfpackPreconditionerFactory>(),
    "Ifpack", true
    );
#endif

  // Note: Above, the last PF object set will be the default!

}


void DefaultLinearSolverBuilder::justInTimeInitialize() const
{
  paramList_.assert_not_null();
  if (is_null(validParamList_)) {
    // Create the validators
    this->getValidParameters();
  }
}



} // namespace Stratimikos
