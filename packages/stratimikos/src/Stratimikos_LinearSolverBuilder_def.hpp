// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STRATIMIKOS_LINEARSOLVERBUILDER_DEF_HPP
#define STRATIMIKOS_LINEARSOLVERBUILDER_DEF_HPP

//#define THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDER_DUMP

#include "Stratimikos_InternalConfig.h"
#include "Thyra_DelayedLinearOpWithSolveFactory.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_STRATIMIKOS_AMESOS
#  include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_AMESOS2
#  include "Thyra_Amesos2LinearOpWithSolveFactory.hpp"
#endif
#if defined(HAVE_STRATIMIKOS_EPETRAEXT) && defined(HAVE_STRATIMIKOS_AZTECOO)
#  include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_BELOS
#  include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_STRATIMIKOS_IFPACK
#  include "Thyra_IfpackPreconditionerFactory.hpp"
#endif
#if defined(HAVE_STRATIMIKOS_IFPACK2) && defined(HAVE_STRATIMIKOS_THYRATPETRAADAPTERS)
#  include "Thyra_Ifpack2PreconditionerFactory.hpp"
#  include "Tpetra_CrsMatrix.hpp"
#  include "Tpetra_RowMatrix.hpp"
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


template<class Scalar>
LinearSolverBuilder<Scalar>::LinearSolverBuilder(
  const std::string &paramsXmlFileName_in,
  const std::string &extraParamsXmlString_in,
  const std::string &paramsUsedXmlOutFileName_in,
  const std::string &paramsXmlFileNameOption_in,
  const std::string &extraParamsXmlStringOption_in,
  const std::string &paramsUsedXmlOutFileNameOption_in,
  const bool &replaceDuplicateFactories_in
  )
  :paramsXmlFileName_(paramsXmlFileName_in)
  ,extraParamsXmlString_(extraParamsXmlString_in)
  ,paramsUsedXmlOutFileName_(paramsUsedXmlOutFileName_in)
  ,paramsXmlFileNameOption_(paramsXmlFileNameOption_in)
  ,extraParamsXmlStringOption_(extraParamsXmlStringOption_in)
  ,paramsUsedXmlOutFileNameOption_(paramsUsedXmlOutFileNameOption_in)
  ,replaceDuplicateFactories_(replaceDuplicateFactories_in)
  ,enableDelayedSolverConstruction_(EnableDelayedSolverConstruction_default)
{
  this->initializeDefaults();
}


template<class Scalar>
LinearSolverBuilder<Scalar>::~LinearSolverBuilder()
{
#ifdef TEUCHOS_DEBUG
  // Validate that we read the parameters correctly!
  if (nonnull(paramList_)) {
    paramList_->validateParameters(*this->getValidParameters());
  }
#endif    
}


template<class Scalar>
void LinearSolverBuilder<Scalar>::setLinearSolveStrategyFactory(
  const RCP<const AbstractFactory<Thyra::LinearOpWithSolveFactoryBase<Scalar> > >
  &solveStrategyFactory,
  const std::string &solveStrategyName,
  const bool makeDefault
  )
{
  const int existingNameIdx =
    this->getAndAssertExistingFactoryNameIdx("setLinearSolveStrategyFactory()",
      validLowsfNames_(), solveStrategyName);
  if (existingNameIdx >= 0) {
    validLowsfNames_[existingNameIdx] = solveStrategyName;
    lowsfArray_[existingNameIdx] = solveStrategyFactory;
  }
  else {
    validLowsfNames_.push_back(solveStrategyName);
    lowsfArray_.push_back(solveStrategyFactory);
  }
  validParamList_ = Teuchos::null;
  if (makeDefault) {
    setDefaultLinearSolveStrategyFactoryName(solveStrategyName);
  }
}


template<class Scalar>
void LinearSolverBuilder<Scalar>::setDefaultLinearSolveStrategyFactoryName(
  const std::string &solveStrategyName)
{
  defaultLOWSF_ = solveStrategyName;
}


template<class Scalar>
void LinearSolverBuilder<Scalar>::setPreconditioningStrategyFactory(
  const RCP<const AbstractFactory<Thyra::PreconditionerFactoryBase<Scalar>>>
    &precStrategyFactory,
  const std::string &precStrategyName,
  const bool makeDefault
  )
{
  const int existingNameIdx =
    this->getAndAssertExistingFactoryNameIdx("setPreconditioningStrategyFactory()",
      validPfNames_(), precStrategyName);
  if (existingNameIdx >= 0) {
    validPfNames_[existingNameIdx] = precStrategyName;
    pfArray_[existingNameIdx-1] = precStrategyFactory; // We offset by -1 since "None" is first!
  }
  else {
    validPfNames_.push_back(precStrategyName);
    pfArray_.push_back(precStrategyFactory);
  }
  validParamList_ = Teuchos::null;
  if (makeDefault) {
    setDefaultPreconditioningStrategyFactoryName(precStrategyName);
  }
}


template<class Scalar>
void LinearSolverBuilder<Scalar>::setDefaultPreconditioningStrategyFactoryName(
  const std::string &precStrategyName)
{
  defaultPF_ = precStrategyName;
}


template<class Scalar>
void LinearSolverBuilder<Scalar>::setupCLP( Teuchos::CommandLineProcessor *clp )
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


template<class Scalar>
void LinearSolverBuilder<Scalar>::readParameters( std::ostream *out )
{
  using Teuchos::parameterList;
  using Teuchos::ptr;
  using Teuchos::updateParametersFromXmlFile;
  using Teuchos::updateParametersFromXmlString;
  using std::endl;
  
  if (!paramList_.get()) {
    paramList_ = parameterList("LinearSolverBuilder");
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


template<class Scalar>
void LinearSolverBuilder<Scalar>::writeParamsFile(
  const Thyra::LinearOpWithSolveFactoryBase<Scalar> &/* lowsFactory */,
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


template<class Scalar>
std::string
LinearSolverBuilder<Scalar>::getLinearSolveStrategyName() const
{
  justInTimeInitialize();
  return lowsfValidator_->getStringValue(*paramList_, LinearSolverType_name,
    defaultLOWSF_);
}


template<class Scalar>
std::string
LinearSolverBuilder<Scalar>::getPreconditionerStrategyName() const
{
  justInTimeInitialize();
  return pfValidator_->getStringValue(*paramList_, PreconditionerType_name,
    defaultPF_);
}


// Overridden from ParameterListAcceptor


template<class Scalar>
void LinearSolverBuilder<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*this->getValidParameters());
  paramList_ = paramList;
  enableDelayedSolverConstruction_ = paramList_->get(
    EnableDelayedSolverConstruction_name, EnableDelayedSolverConstruction_default );
}


template<class Scalar>
RCP<Teuchos::ParameterList>
LinearSolverBuilder<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<Teuchos::ParameterList>
LinearSolverBuilder<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
LinearSolverBuilder<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
LinearSolverBuilder<Scalar>::getValidParameters() const
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
      const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
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
      const RCP<Thyra::PreconditionerFactoryBase<Scalar> >
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


template<class Scalar>
RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
LinearSolverBuilder<Scalar>::createLinearSolveStrategy(
  const std::string &linearSolveStrategyName
  ) const
{
  justInTimeInitialize();

  // Get the name of the linear solve strategy
#ifdef THYRA_DEFAULT_REAL_LINEAR_SOLVER_BUILDER_DUMP
  std::cout << "\nEntering LinearSolverBuilder"
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
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    lowsf = lowsfArray_[ls_idx]->create();

  // First, set the preconditioner factory and its parameters
  if(lowsf->acceptsPreconditionerFactory()) {
    const std::string &pfName = this->getPreconditionerStrategyName();
    RCP<Thyra::PreconditionerFactoryBase<Scalar> >
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
      new Thyra::DelayedLinearOpWithSolveFactory<Scalar>(lowsf)
      );
  }

  return lowsf;

}


template<class Scalar>
RCP<Thyra::PreconditionerFactoryBase<Scalar> >
LinearSolverBuilder<Scalar>::createPreconditioningStrategy(
  const std::string &preconditioningStrategyName
  ) const
{
  justInTimeInitialize();

  // Get the name of the preconditioning strategy
  const std::string
    pfname = ( preconditioningStrategyName.length()
             ? preconditioningStrategyName
             : this->getPreconditionerStrategyName() );
  RCP<Thyra::PreconditionerFactoryBase<Scalar> >
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

template<class Scalar>
void LinearSolverBuilder<Scalar>::initializeDefaults()
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

#ifdef HAVE_STRATIMIKOS_AMESOS2
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<Scalar>,
    Thyra::Amesos2LinearOpWithSolveFactory<Scalar>>(),
    "Amesos2", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_BELOS
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<Scalar>,
    Thyra::BelosLinearOpWithSolveFactory<Scalar> >(),
    "Belos", true
    );
#endif

#if defined(HAVE_STRATIMIKOS_IFPACK2) && defined(HAVE_STRATIMIKOS_THYRATPETRAADAPTERS)
  setPreconditioningStrategyFactory(
    abstractFactoryStd<Thyra::PreconditionerFactoryBase<Scalar>,
    Thyra::Ifpack2PreconditionerFactory<Tpetra::RowMatrix<Scalar>>>(),
    "Ifpack2", true
    );
#endif


  // Note: Above, the last PF object set will be the default!

}

template<>
void LinearSolverBuilder<double>::initializeDefaults()
{
  using Scalar = double;
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

#ifdef HAVE_STRATIMIKOS_AMESOS2
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<Scalar>,
    Thyra::Amesos2LinearOpWithSolveFactory<Scalar>>(),
    "Amesos2", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_BELOS
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<Scalar>,
    Thyra::BelosLinearOpWithSolveFactory<Scalar> >(),
    "Belos", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_AMESOS
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<Scalar>,
    Thyra::AmesosLinearOpWithSolveFactory>(),
    "Amesos", true
    );
#endif

#if defined(HAVE_STRATIMIKOS_EPETRAEXT) && defined(HAVE_STRATIMIKOS_AZTECOO)
  setLinearSolveStrategyFactory(
    abstractFactoryStd<Thyra::LinearOpWithSolveFactoryBase<Scalar>,
    Thyra::AztecOOLinearOpWithSolveFactory>(),
    "AztecOO", true
    );
#endif

  // Note: Above, the last LOWSF object set will be the default!
  // (unless we have only one processor, see below:)

#ifdef HAVE_STRATIMIKOS_AMESOS
  if (Teuchos::GlobalMPISession::getNProc() == 1) {
    setDefaultLinearSolveStrategyFactoryName("Amesos");
  }
#endif

  //
  // Preconditioners
  //

#ifdef HAVE_STRATIMIKOS_ML
  setPreconditioningStrategyFactory(
    abstractFactoryStd<Thyra::PreconditionerFactoryBase<Scalar>,
    Thyra::MLPreconditionerFactory>(),
    "ML", true
    );
#endif

#if defined(HAVE_STRATIMIKOS_IFPACK2) && defined(HAVE_STRATIMIKOS_THYRATPETRAADAPTERS)
  setPreconditioningStrategyFactory(
    abstractFactoryStd<Thyra::PreconditionerFactoryBase<Scalar>,
    Thyra::Ifpack2PreconditionerFactory<Tpetra::RowMatrix<Scalar>>>(),
    "Ifpack2", true
    );
#endif

#ifdef HAVE_STRATIMIKOS_IFPACK
  setPreconditioningStrategyFactory(
    abstractFactoryStd<Thyra::PreconditionerFactoryBase<Scalar>,
    Thyra::IfpackPreconditionerFactory>(),
    "Ifpack", true
    );
#endif

  // Note: Above, the last PF object set will be the default!

}


template<class Scalar>
void LinearSolverBuilder<Scalar>::justInTimeInitialize() const
{
  paramList_.assert_not_null();
  if (is_null(validParamList_)) {
    // Create the validators
    this->getValidParameters();
  }
}


template<class Scalar>
int LinearSolverBuilder<Scalar>::getAndAssertExistingFactoryNameIdx(
  const std::string &setFunctionName, const Teuchos::ArrayView<std::string> namesArray,
  const std::string &name) const
{
  const int existingNameIdx =
    LinearSolverBuilderHelpers::existingNameIndex(namesArray, name);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (!replaceDuplicateFactories_ && existingNameIdx >= 0), std::logic_error,
    "ERROR: "<<setFunctionName<<": the name='"<<name<<"' already exists"
    << " at index="<<existingNameIdx<<"!\n"
    << "\n"
    << "TIP: To allow duplicates, change the property replaceDuplicateFactories from"
    << " false to true by calling <thisObject>.replaceDuplicateFactories(true)"
    << " where <thisObject> is of type Stratimikos::LinearSolverBuilder<Scalar>!");
  return existingNameIdx;
}


//
// Explicit instantiation macro
//
// Must be expanded from within the Stratimikos namespace!
//


#define STRATIMIKOS_LINEARSOLVERBUILDER_INSTANT(SCALAR) \
  \
  template class LinearSolverBuilder<SCALAR >;



} // namespace Stratimikos

#endif // STRATIMIKOS_LINEARSOLVERBUILDER_DEF_HPP
