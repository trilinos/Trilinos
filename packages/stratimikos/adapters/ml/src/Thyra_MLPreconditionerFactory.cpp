// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Thyra_MLPreconditionerFactory.hpp"

#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_ValidateParameters.h"
#include "ml_RefMaxwell.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"
#include "Teuchos_iostream_helpers.hpp"


namespace {


enum EMLProblemType {
  ML_PROBTYPE_NONE,
  ML_PROBTYPE_SMOOTHED_AGGREGATION, 
  ML_PROBTYPE_NONSYMMETRIC_SMOOTHED_AGGREGATION,
  ML_PROBTYPE_DOMAIN_DECOMPOSITION,
  ML_PROBTYPE_DOMAIN_DECOMPOSITION_ML,
  ML_PROBTYPE_MAXWELL,
  ML_PROBTYPE_REFMAXWELL
};
const std::string BaseMethodDefaults_valueNames_none = "none";
const Teuchos::Array<std::string> BaseMethodDefaults_valueNames
= Teuchos::tuple<std::string>(
  BaseMethodDefaults_valueNames_none,
  "SA", 
  "NSSA",
  "DD",
  "DD-ML",
  "maxwell",
  "refmaxwell"
  );


TEUCHOS_ENUM_INPUT_STREAM_OPERATOR(EMLProblemType)


TEUCHOS_STATIC_SETUP()
{
  TEUCHOS_ADD_STRINGTOINTEGRALVALIDATOR_CONVERTER(EMLProblemType);
}

const std::string BaseMethodDefaults_name = "Base Method Defaults";
const std::string BaseMethodDefaults_default = "SA";
Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<EMLProblemType>
  >
BaseMethodDefaults_validator;

const std::string ReuseFineLevelSmoother_name = "Reuse Fine Level Smoother";
const bool ReuseFineLevelSmoother_default = false;
  
const std::string MLSettings_name = "ML Settings";


} // namespace


namespace Thyra {


using Teuchos::RCP;
using Teuchos::ParameterList;


// Constructors/initializers/accessors

  
MLPreconditionerFactory::MLPreconditionerFactory()
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
{}


// Overridden from PreconditionerFactoryBase


bool MLPreconditionerFactory::isCompatible(
  const LinearOpSourceBase<double> &fwdOpSrc
  ) const
{
  using Teuchos::outArg;
  Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
  EOpTransp epetraFwdOpTransp;
  EApplyEpetraOpAs epetraFwdOpApplyAs;
  EAdjointEpetraOp epetraFwdOpAdjointSupport;
  double epetraFwdOpScalar;
  Teuchos::RCP<const LinearOpBase<double> >
    fwdOp = fwdOpSrc.getOp();
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp,
    outArg(epetraFwdOp),outArg(epetraFwdOpTransp),
    outArg(epetraFwdOpApplyAs),
    outArg(epetraFwdOpAdjointSupport),
    outArg(epetraFwdOpScalar)
    );
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}


bool MLPreconditionerFactory::applySupportsConj(EConj conj) const
{
  return true;
}


bool MLPreconditionerFactory::applyTransposeSupportsConj(EConj conj) const
{
  return false; // See comment below
}


RCP<PreconditionerBase<double> >
MLPreconditionerFactory::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<double>());
}


void MLPreconditionerFactory::initializePrec(
  const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
  PreconditionerBase<double> *prec,
  const ESupportSolveUse supportSolveUse
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
  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && implicit_cast<int>(verbLevel) > implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::MLPreconditionerFactory::initializePrec(...) ...\n";

  // Get the problem type
  const EMLProblemType problemType = BaseMethodDefaults_validator->getIntegralValue(*paramList_,BaseMethodDefaults_name,BaseMethodDefaults_default);

  Teuchos::RCP<const LinearOpBase<double> > fwdOp = fwdOpSrc->getOp();
#ifdef _DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(fwdOp.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(prec==NULL);
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
    fwdOp,outArg(epetraFwdOp),outArg(epetraFwdOpTransp),outArg(epetraFwdOpApplyAs),
    outArg(epetraFwdOpAdjointSupport),outArg(epetraFwdOpScalar)
                                             );
  // Validate what we get is what we need
  RCP<const Epetra_RowMatrix>
    epetraFwdRowMat = rcp_dynamic_cast<const Epetra_RowMatrix>(epetraFwdOp,true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    epetraFwdOpApplyAs != EPETRA_OP_APPLY_APPLY, std::logic_error
    ,"Error, incorrect apply mode for an Epetra_RowMatrix"
    );
  RCP<const Epetra_CrsMatrix> epetraFwdCrsMat = rcp_dynamic_cast<const Epetra_CrsMatrix>(epetraFwdRowMat);

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
  // Get the embedded ML_Epetra Preconditioner object if it exists
  //
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_precOp;
  Teuchos::RCP<ML_Epetra::RefMaxwellPreconditioner> rm_precOp;
  if(epetra_precOp.get()) {
    if(problemType == ML_PROBTYPE_REFMAXWELL)
      rm_precOp = rcp_dynamic_cast<ML_Epetra::RefMaxwellPreconditioner>(epetra_precOp->epetra_op(),true);
    else
      ml_precOp = rcp_dynamic_cast<ML_Epetra::MultiLevelPreconditioner>(epetra_precOp->epetra_op(),true);
  }
  //
  // Get the attached forward operator if it exists and make sure that it matches
  //
  if(ml_precOp!=Teuchos::null) {
    // Get the forward operator and make sure that it matches what is
    // already being used!
    const Epetra_RowMatrix & rm = ml_precOp->RowMatrix();
   
    TEUCHOS_TEST_FOR_EXCEPTION(
       &rm!=&*epetraFwdRowMat, std::logic_error
       ,"ML requires Epetra_RowMatrix to be the same for each initialization of the preconditioner"
       );
  }
  // NOTE: No such check exists for RefMaxwell

  //
  // Perform initialization if needed
  //
  const bool startingOver = (ml_precOp.get() == NULL && rm_precOp.get() == NULL);
  if(startingOver) 
  {
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      *out << "\nCreating the initial ML_Epetra::MultiLevelPreconditioner object...\n";
    timer.start(true);
    // Create the initial preconditioner: DO NOT compute it yet

    if(problemType==ML_PROBTYPE_REFMAXWELL)
      rm_precOp = rcp(new ML_Epetra::RefMaxwellPreconditioner(*epetraFwdCrsMat, paramList_->sublist(MLSettings_name),false));
    else
      ml_precOp = rcp(new ML_Epetra::MultiLevelPreconditioner(*epetraFwdRowMat, paramList_->sublist(MLSettings_name),false));

    
    timer.stop();
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      OSTab(out).o() <<"> Creation time = "<<timer.totalElapsedTime()<<" sec\n";
    // RAB: Above, I am just passing a string to ML::Create(...) in order
    // get this code written.  However, in the future, it would be good to
    // copy the contents of what is in ML::Create(...) into a local
    // function and then use switch(...) to create the initial
    // ML_Epetra::MultiLevelPreconditioner object.  This would result in better validation
    // and faster code.
    // Set parameters if the list exists
    if(paramList_.get()) {
      if (problemType==ML_PROBTYPE_REFMAXWELL) {
	TEUCHOS_TEST_FOR_EXCEPT(0!=rm_precOp->SetParameterList(paramList_->sublist(MLSettings_name)));
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPT(0!=ml_precOp->SetParameterList(paramList_->sublist(MLSettings_name)));
      }
    }
  }
  //
  // Attach the epetraFwdOp to the ml_precOp to guarantee that it will not go away
  //
  if (problemType==ML_PROBTYPE_REFMAXWELL)
    set_extra_data(epetraFwdOp, "IFPF::epetraFwdOp", Teuchos::inOutArg(rm_precOp),
		   Teuchos::POST_DESTROY, false);
  else
    set_extra_data(epetraFwdOp, "IFPF::epetraFwdOp", Teuchos::inOutArg(ml_precOp),
		   Teuchos::POST_DESTROY, false);
  //
  // Update the factorization
  //
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the preconditioner ...\n";
  timer.start(true);
  if (problemType==ML_PROBTYPE_REFMAXWELL) {
    if (startingOver) {
      TEUCHOS_TEST_FOR_EXCEPT(0!=rm_precOp->ComputePreconditioner());
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(0!=rm_precOp->ReComputePreconditioner());
    }
  }
  else {
    if (startingOver) {
      TEUCHOS_TEST_FOR_EXCEPT(0!=ml_precOp->ComputePreconditioner());
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(0!=ml_precOp->ReComputePreconditioner(paramList_->get<bool>(ReuseFineLevelSmoother_name)));
    }
  }
  timer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() <<"=> Setup time = "<<timer.totalElapsedTime()<<" sec\n";
  //
  // Compute the conditioner number estimate if asked
  //

  // ToDo: Implement

  //
  // Attach fwdOp to the ml_precOp
  //
  if (problemType==ML_PROBTYPE_REFMAXWELL)
    set_extra_data(fwdOp, "IFPF::fwdOp", Teuchos::inOutArg(rm_precOp),
		   Teuchos::POST_DESTROY, false);
  else
    set_extra_data(fwdOp, "IFPF::fwdOp", Teuchos::inOutArg(ml_precOp),
		   Teuchos::POST_DESTROY, false);
  //
  // Initialize the output EpetraLinearOp
  //
  if(startingOver) {
    epetra_precOp = rcp(new EpetraLinearOp);
  }
  // ToDo: Look into adjoints again.
  if (problemType==ML_PROBTYPE_REFMAXWELL)
    epetra_precOp->initialize(rm_precOp,epetraFwdOpTransp,EPETRA_OP_APPLY_APPLY_INVERSE,EPETRA_OP_ADJOINT_UNSUPPORTED);  
  else
    epetra_precOp->initialize(ml_precOp,epetraFwdOpTransp,EPETRA_OP_APPLY_APPLY_INVERSE,EPETRA_OP_ADJOINT_UNSUPPORTED);
  //
  // Initialize the preconditioner
  //
  defaultPrec->initializeUnspecified(
    Teuchos::rcp_implicit_cast<LinearOpBase<double> >(epetra_precOp)
    );
  totalTimer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nTotal time in MLPreconditionerFactory = "<<totalTimer.totalElapsedTime()<<" sec\n";
  if(out.get() && implicit_cast<int>(verbLevel) > implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nLeaving Thyra::MLPreconditionerFactory::initializePrec(...) ...\n";
}


void MLPreconditionerFactory::uninitializePrec(
  PreconditionerBase<double> *prec,
  Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOp,
  ESupportSolveUse *supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}


// Overridden from ParameterListAcceptor


void MLPreconditionerFactory::setParameterList(
  Teuchos::RCP<ParameterList> const& paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);

  // Do not recurse
  paramList->validateParameters(*this->getValidParameters(),0);
  paramList_ = paramList;

  // set default for reuse of fine level smoother
  if(!paramList_->isType<bool>(ReuseFineLevelSmoother_name))
    paramList_->set<bool>(ReuseFineLevelSmoother_name,ReuseFineLevelSmoother_default);

  const EMLProblemType
    defaultType = BaseMethodDefaults_validator->getIntegralValue(
      *paramList_,BaseMethodDefaults_name,BaseMethodDefaults_default
      );
  if( ML_PROBTYPE_NONE != defaultType ) {
    const std::string defaultTypeStr = BaseMethodDefaults_valueNames[defaultType];
    
    // ML will do validation on its own.  We don't need to duplicate that here.
    Teuchos::ParameterList defaultParams;
    if(defaultType == ML_PROBTYPE_REFMAXWELL) {	
      ML_Epetra::SetDefaultsRefMaxwell(defaultParams);
    }
    else {	
      TEUCHOS_TEST_FOR_EXCEPTION(0!=ML_Epetra::SetDefaults(defaultTypeStr,defaultParams)				
				 ,Teuchos::Exceptions::InvalidParameterValue
				 ,"Error, the ML problem type \"" << defaultTypeStr << "\' is not recognized by ML!"
				 );
    }
    
    // Note, the only way the above exception message could be generated is if
    // a default problem type was removed from ML_Epetra::SetDefaults(...).
    // When a new problem type is added to this function, it must be added to
    // our enum EMLProblemType along with associated objects ...  In other
    // words, this adapter must be maintained as ML is maintained.  An
    // alternative design would be to just pass in whatever string to this
    // function.  This would improve maintainability but it would not generate
    // very good error messages when a bad string was passed in.  Currently,
    // the error message attached to the exception will contain the list of
    // valid problem types.
    paramList_->sublist(MLSettings_name).setParametersNotAlreadySet(
      defaultParams);
  }

}


RCP<ParameterList>
MLPreconditionerFactory::getNonconstParameterList()
{
  return paramList_;
}


RCP<ParameterList>
MLPreconditionerFactory::unsetParameterList()
{
  Teuchos::RCP<ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


RCP<const ParameterList>
MLPreconditionerFactory::getParameterList() const
{
  return paramList_;
}


RCP<const ParameterList>
MLPreconditionerFactory::getValidParameters() const
{
  // NOTE: We're only going to use this function to genrate valid *Stratimikos* parameters.
  // Since ML's parameters can be validated internally, we'll handle those separarely. 


  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::implicit_cast;
  using Teuchos::rcp_implicit_cast;
  typedef Teuchos::ParameterEntryValidator PEV;

  static RCP<const ParameterList> validPL;

  if(is_null(validPL)) {

    RCP<ParameterList>
      pl = rcp(new ParameterList());

    BaseMethodDefaults_validator = rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<EMLProblemType>(
        BaseMethodDefaults_valueNames,
        tuple<std::string>(
          "Do not set any default parameters",
          "Set default parameters for a smoothed aggregation method",
	  "Set default parameters for a nonsymmetric smoothed aggregation method",
          "Set default parameters for a domain decomposition method",
          "Set default parameters for a domain decomposition method special to ML",
          "Set default parameters for a Maxwell-type of preconditioner",
          "Set default parameters for a RefMaxwell-type preconditioner"
          ),
        tuple<EMLProblemType>(
          ML_PROBTYPE_NONE,
          ML_PROBTYPE_SMOOTHED_AGGREGATION,
	  ML_PROBTYPE_NONSYMMETRIC_SMOOTHED_AGGREGATION,
          ML_PROBTYPE_DOMAIN_DECOMPOSITION,
          ML_PROBTYPE_DOMAIN_DECOMPOSITION_ML,
          ML_PROBTYPE_MAXWELL,
	  ML_PROBTYPE_REFMAXWELL
          ),
        BaseMethodDefaults_name
        )
      );

    pl->set(BaseMethodDefaults_name,BaseMethodDefaults_default,
      "Select the default method type which also sets parameter defaults\n"
      "in the sublist \"" + MLSettings_name + "\"!",
      rcp_implicit_cast<const PEV>(BaseMethodDefaults_validator)
      );

    pl->set(ReuseFineLevelSmoother_name,ReuseFineLevelSmoother_default,
      "Enables/disables the reuse of the fine level smoother.");

    ParameterList mlpl;    
    pl->set(MLSettings_name,mlpl,
        "Sampling of the parameters directly accepted by ML\n"
        "This list of parameters is generated by combining all of\n"
        "the parameters set for all of the default problem types supported\n"
        "by ML.  Therefore, do not assume these parameters are at values that\n"
        "are reasonable to ML.  This list is just to give a sense of some of\n"
        "the parameters that ML accepts.  Consult ML documentation on how to\n"
        "set these parameters.  Also, you can print the parameter list after\n"
        "it is used and see what defaults where set for each default problem\n"
        "type.  Warning! the parameters in this sublist are currently *not*\n"
        "being validated by ML!"
        );     

    validPL = pl;

  }
  
  return validPL;

}


// Public functions overridden from Teuchos::Describable


std::string MLPreconditionerFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::MLPreconditionerFactory";
  return oss.str();
}


} // namespace Thyra
