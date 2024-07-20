// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_IfpackPreconditionerFactory.hpp"
#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Ifpack_ValidParameters.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"


namespace {

Teuchos::RCP<Teuchos::Time> overallTimer, creationTimer, factorizationTimer;

const std::string Ifpack_name = "Ifpack";

const std::string IfpackSettings_name = "Ifpack Settings";

const std::string PrecType_name = "Prec Type";
Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<Ifpack::EPrecType> >
precTypeValidator; // Will be setup below!
const Ifpack::EPrecType PrecType_default = Ifpack::ILU;
const std::string PrecTypeName_default = Ifpack::precTypeNames[PrecType_default];

const std::string Overlap_name = "Overlap";
const int Overlap_default = 0;


TEUCHOS_STATIC_SETUP()
{
  TEUCHOS_ADD_STRINGTOINTEGRALVALIDATOR_CONVERTER(Ifpack::EPrecType);
}


} // namespace

namespace Thyra {

// Constructors/initializers/accessors

IfpackPreconditionerFactory::IfpackPreconditionerFactory()
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
  ,precType_(PrecType_default)
  ,overlap_(Overlap_default)
{
  initializeTimers();
  getValidParameters(); // Make sure validators get created!
}

// Overridden from PreconditionerFactoryBase

bool IfpackPreconditionerFactory::isCompatible(
  const LinearOpSourceBase<double> &fwdOpSrc
  ) const
{
  using Teuchos::outArg;
  Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
  EOpTransp epetraFwdOpTransp;
  EApplyEpetraOpAs epetraFwdOpApplyAs;
  EAdjointEpetraOp epetraFwdOpAdjointSupport;
  double epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOpSrc.getOp(), 
    outArg(epetraFwdOp), outArg(epetraFwdOpTransp),
    outArg(epetraFwdOpApplyAs), outArg(epetraFwdOpAdjointSupport),
    outArg(epetraFwdOpScalar)
    );
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}

bool IfpackPreconditionerFactory::applySupportsConj(EConj /* conj */) const
{
  return true;
}

bool IfpackPreconditionerFactory::applyTransposeSupportsConj(EConj /* conj */) const
{
  return false; // See comment below
}

Teuchos::RCP<PreconditionerBase<double> >
IfpackPreconditionerFactory::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<double>());
}

void IfpackPreconditionerFactory::initializePrec(
  const Teuchos::RCP<const LinearOpSourceBase<double> >    &fwdOpSrc
  ,PreconditionerBase<double>                                      *prec
  ,const ESupportSolveUse                                           /* supportSolveUse */
  ) const
{
  using Teuchos::outArg;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::RCP;
  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::set_extra_data;
  using Teuchos::get_optional_extra_data;
  using Teuchos::implicit_cast;
  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);
#ifdef STRATIMIKOS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor overallTimeMonitor(*overallTimer);
#endif
  const Teuchos::RCP<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && implicit_cast<int>(verbLevel) > implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::IfpackPreconditionerFactory::initializePrec(...) ...\n";
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(prec==NULL);
#endif
  Teuchos::RCP<const LinearOpBase<double> >
    fwdOp = fwdOpSrc->getOp();
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(fwdOp.get()==NULL);
#endif
  //
  // Unwrap and get the forward Epetra_Operator object
  //
  Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
  EOpTransp epetraFwdOpTransp;
  EApplyEpetraOpAs epetraFwdOpApplyAs;
  EAdjointEpetraOp epetraFwdOpAdjointSupport;
  double epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp,
    outArg(epetraFwdOp), outArg(epetraFwdOpTransp),
    outArg(epetraFwdOpApplyAs), outArg(epetraFwdOpAdjointSupport),
    outArg(epetraFwdOpScalar)
    );
  // Validate what we get is what we need
  RCP<const Epetra_RowMatrix>
    epetraFwdRowMat = rcp_dynamic_cast<const Epetra_RowMatrix>(epetraFwdOp,true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    epetraFwdOpApplyAs != EPETRA_OP_APPLY_APPLY, std::logic_error
    ,"Error, incorrect apply mode for an Epetra_RowMatrix"
    );
  //
  // Get the concrete preconditioner object
  //
  DefaultPreconditioner<double>
    *defaultPrec = &Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);
  //
  // Get the EpetraLinearOp object that is used to implement the preconditoner linear op
  //
  RCP<EpetraLinearOp>
    epetra_precOp = rcp_dynamic_cast<EpetraLinearOp>(defaultPrec->getNonconstUnspecifiedPrecOp(),true);
  //
  // Get the embedded Ifpack_Preconditioner object if it exists
  //
  Teuchos::RCP<Ifpack_Preconditioner>
    ifpack_precOp;
  if(epetra_precOp.get())
    ifpack_precOp = rcp_dynamic_cast<Ifpack_Preconditioner>(epetra_precOp->epetra_op(),true);
  //
  // Get the attached forward operator if it exists and make sure that it matches
  //
  if(ifpack_precOp.get()) {
    // ToDo: Get the forward operator and make sure that it matches what is
    // already being used!
  }
  //
  // Permform initialization if needed
  //
  //const bool startingOver = (ifpack_precOp.get() == NULL);
  const bool startingOver = true;
  // ToDo: Comment back in the above original version of startingOver to allow
  // for resuse.  Rob H. just pointed out to me that a memory leak is being
  // created when you just call Ifpack_ILU::Compute() over and over again.
  // Rob H. said that he will check in a fix the the development branch when
  // he can.
  if(startingOver) {
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      *out << "\nCreating the initial Ifpack_Preconditioner object of type \'"<<Ifpack::toString(precType_)<<"\' ...\n";
    timer.start(true);
#ifdef STRATIMIKOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor creationTimeMonitor(*creationTimer);
#endif
    // Create the initial preconditioner
    ifpack_precOp = rcp(
      ::Ifpack::Create(
        precType_
        ,const_cast<Epetra_RowMatrix*>(&*epetraFwdRowMat)
        ,overlap_
        )
      );
    timer.stop();
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      OSTab(out).o() <<"=> Creation time = "<<timer.totalElapsedTime()<<" sec\n";
    // Set parameters if the list exists
    if(paramList_.get()) {
      Teuchos::ParameterList
        &ifpackSettingsPL = paramList_->sublist(IfpackSettings_name);
      // Above will create new sublist if it does not exist!
      TEUCHOS_TEST_FOR_EXCEPT(0!=ifpack_precOp->SetParameters(ifpackSettingsPL));
      // Above, I have not idea how any error messages for a mistake will be
      // reported back to the user!
    }
    // Initialize the structure for the preconditioner
    TEUCHOS_TEST_FOR_EXCEPT(0!=ifpack_precOp->Initialize());
  }
  //
  // Attach the epetraFwdOp to the ifpack_precOp to guarantee that it will not go away
  //
  set_extra_data(epetraFwdOp, "IFPF::epetraFwdOp", Teuchos::inOutArg(ifpack_precOp),
    Teuchos::POST_DESTROY, false);
  //
  // Update the factorization
  //
  {
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      *out << "\nComputing the preconditioner ...\n";
#ifdef STRATIMIKOS_TEUCHOS_TIME_MONITOR
    Teuchos::TimeMonitor factorizationTimeMonitor(*factorizationTimer);
#endif
    timer.start(true);
    TEUCHOS_TEST_FOR_EXCEPT(0!=ifpack_precOp->Compute());
    timer.stop();
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      OSTab(out).o() <<"=> Setup time = "<<timer.totalElapsedTime()<<" sec\n";
  }
  //
  // Compute the conditioner number estimate if asked
  //

  // ToDo: Implement

  //
  // Attach fwdOp to the ifpack_precOp
  //
  set_extra_data(fwdOpSrc, "IFPF::fwdOpSrc", Teuchos::inOutArg(ifpack_precOp),
    Teuchos::POST_DESTROY, false);
  //
  // Initialize the output EpetraLinearOp
  //
  if(startingOver) {
    epetra_precOp = rcp(new EpetraLinearOp);
  }
  epetra_precOp->initialize(
    ifpack_precOp
    ,epetraFwdOpTransp
    ,EPETRA_OP_APPLY_APPLY_INVERSE
    ,EPETRA_OP_ADJOINT_SUPPORTED
    );
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_MEDIUM)) {
    *out << "\nDescription of created preconditioner:\n";
    OSTab tab2(out);
    ifpack_precOp->Print(*out);
  }

  //
  // Initialize the preconditioner
  //
  defaultPrec->initializeUnspecified(
    Teuchos::rcp_implicit_cast<LinearOpBase<double> >(epetra_precOp)
    );
  totalTimer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nTotal time in IfpackPreconditionerFactory = "<<totalTimer.totalElapsedTime()<<" sec\n";
  if(out.get() && implicit_cast<int>(verbLevel) > implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nLeaving Thyra::IfpackPreconditionerFactory::initializePrec(...) ...\n";
}

void IfpackPreconditionerFactory::uninitializePrec(
  PreconditionerBase<double>                                * /* prec */
  ,Teuchos::RCP<const LinearOpSourceBase<double> >  * /* fwdOpSrc */
  ,ESupportSolveUse                                         * /* supportSolveUse */
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement when needed!
}

// Overridden from ParameterListAcceptor

void IfpackPreconditionerFactory::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
{
  TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),2);
  // Note: The above validation will validate right down into the the sublist
  // named IfpackSettings_name!
  paramList_ = paramList;
  overlap_ = paramList_->get(Overlap_name,Overlap_default);
  std::ostringstream oss;
  oss << "(sub)list \""<<paramList->name()<<"\"parameter \"Prec Type\"";
  precType_ =
    ( paramList_.get()
      ? precTypeValidator->getIntegralValue(*paramList_,PrecType_name,PrecTypeName_default)
      : PrecType_default
      );
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
#ifdef TEUCHOS_DEBUG
  // Validate my use of the parameters!
  paramList->validateParameters(*this->getValidParameters(),1);
#endif
}

Teuchos::RCP<Teuchos::ParameterList>
IfpackPreconditionerFactory::getNonconstParameterList()
{
  return paramList_;
}

Teuchos::RCP<Teuchos::ParameterList>
IfpackPreconditionerFactory::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RCP<const Teuchos::ParameterList>
IfpackPreconditionerFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RCP<const Teuchos::ParameterList>
IfpackPreconditionerFactory::getValidParameters() const
{
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  typedef Teuchos::ParameterEntryValidator PEV;
  static Teuchos::RCP<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList(Ifpack_name));
    {
      // Create the validator for the preconditioner type!
      Teuchos::Array<std::string>
        precTypeNames;
      precTypeNames.insert(
        precTypeNames.begin(),
        &Ifpack::precTypeNames[0],
        &Ifpack::precTypeNames[0] + Ifpack::numPrecTypes
        );
      Teuchos::Array<Ifpack::EPrecType>
        precTypeValues;
      precTypeValues.insert(
        precTypeValues.begin(),
        &Ifpack::precTypeValues[0],
        &Ifpack::precTypeValues[0] + Ifpack::numPrecTypes
        );
      precTypeValidator = rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<Ifpack::EPrecType>(
          precTypeNames,precTypeValues,PrecType_name
          )
        );
    }
    validParamList->set(
      PrecType_name, PrecTypeName_default,
      "Type of Ifpack preconditioner to use.",
      rcp_implicit_cast<const PEV>(precTypeValidator)
      );
    validParamList->set(
      Overlap_name, Overlap_default,
      "Number of rows/columns overlapped between subdomains in different"
      "\nprocesses in the additive Schwarz-type domain-decomposition preconditioners."
      );
    validParamList->set(
      IfpackSettings_name, Ifpack_GetValidParameters(),
      "Preconditioner settings that are passed onto the Ifpack preconditioners themselves."
      );
    // Note that in the above setParameterList(...) function that we actually
    // validate down into the first level of this sublist.  Really the
    // Ifpack_Preconditioner objects themselves should do validation but we do
    // it ourselves taking the return from the Ifpack_GetValidParameters()
    // function as gospel!
    Teuchos::setupVerboseObjectSublist(&*validParamList);
  }
  return validParamList;
}

// Public functions overridden from Teuchos::Describable

std::string IfpackPreconditionerFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::IfpackPreconditionerFactory{";
  oss << "precType=\"" << ::Ifpack::toString(precType_) << "\"";
  oss << ",overlap=" << overlap_;
  oss << "}";
  return oss.str();
}

// private

void IfpackPreconditionerFactory::initializeTimers()
{
  if(!overallTimer.get()) {
#ifdef STRATIMIKOS_TEUCHOS_TIME_MONITOR
    overallTimer       = Teuchos::TimeMonitor::getNewTimer("Stratimikos: IfpackPF");
    creationTimer      = Teuchos::TimeMonitor::getNewTimer("Stratimikos: IfpackPF:Creation");
    factorizationTimer = Teuchos::TimeMonitor::getNewTimer("Stratimikos: IfpackPF:Factorization");
#endif
  }
}

} // namespace Thyra
