// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Thyra_MueLuPreconditionerFactory.hpp"

#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StaticSetupMacro.hpp"
#include "Teuchos_iostream_helpers.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#include "MueLu_EpetraOperator.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_Level.hpp"

namespace {

const std::string MueLuSettings_name = "MueLu Settings";

} // namespace


namespace Thyra {


using Teuchos::RCP;
using Teuchos::ParameterList;


// Constructors/initializers/accessors

  
MueLuPreconditionerFactory::MueLuPreconditionerFactory()
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
{}


// Overridden from PreconditionerFactoryBase


bool MueLuPreconditionerFactory::isCompatible(
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
  if( !dynamic_cast<const Epetra_CrsMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}


bool MueLuPreconditionerFactory::applySupportsConj(EConj conj) const
{
  return false;
}


bool MueLuPreconditionerFactory::applyTransposeSupportsConj(EConj conj) const
{
  return false;
}


RCP<PreconditionerBase<double> >
MueLuPreconditionerFactory::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<double>());
}


void MueLuPreconditionerFactory::initializePrec(
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

  typedef Kokkos::DefaultNode::DefaultNodeType NO;
  typedef Kokkos::DefaultKernels<double,int,NO>::SparseOps LMO;

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && implicit_cast<int>(verbLevel) > implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::MueLuPreconditionerFactory::initializePrec(...) ...\n";

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

  RCP<const Epetra_CrsMatrix>
    epetraFwdCrsMat = rcp_dynamic_cast<const Epetra_CrsMatrix>(epetraFwdOp,true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    epetraFwdOpApplyAs != EPETRA_OP_APPLY_APPLY, std::logic_error
    ,"Error, incorrect apply mode for an Epetra_CrsMatrix"
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
  // Get the embedded MueLu::EpetraOperator object if it exists
  //

  Teuchos::RCP<MueLu::EpetraOperator> muelu_precOp;
  if(epetra_precOp.get())
    muelu_precOp = rcp_dynamic_cast<MueLu::EpetraOperator>(epetra_precOp->epetra_op(),true);
  //
  // Get the attached forward operator if it exists and make sure that it matches
  //
  if(muelu_precOp!=Teuchos::null) {
    // TODO
//     // Get the forward operator and make sure that it matches what is
//     // already being used!
//     const Epetra_CrsMatrix & rm = muelu_precOp->CrsMatrix();
   
//     TEUCHOS_TEST_FOR_EXCEPTION(
//        &rm!=&*epetraFwdRowMat, std::logic_error
//        ,"MueLu requires Epetra_RowMatrix to be the same for each initialization of the preconditioner"
//        );
  }

  MueLu::ParameterListInterpreter<double, int, int, NO, LMO> mueluFactory(*paramList_);

  //
  // Perform initialization if needed
  //
  const bool startingOver = (muelu_precOp.get() == NULL);
  if(startingOver) 
  {
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      *out << "\nCreating the initial MueLu::EpetraOperator object...\n";
    timer.start(true);
    // Create the initial preconditioner: DO NOT compute it yet

    // Turns a Epetra_CrsMatrix into a Xpetra::Matrix
    RCP<Epetra_CrsMatrix> epetraFwdCrsMatNonConst = rcp_const_cast<Epetra_CrsMatrix>(epetraFwdCrsMat); // !! TODO: MueLu interface should accept const matrix as input.

    RCP<Xpetra::CrsMatrix<double, int, int, NO, LMO> > mueluAcrs = rcp(new Xpetra::EpetraCrsMatrix(epetraFwdCrsMatNonConst)); //TODO: should not be needed
    RCP<Xpetra::Matrix <double, int, int, NO, LMO> >   mueluA    = rcp(new Xpetra::CrsMatrixWrap<double, int, int, NO, LMO>(mueluAcrs));

    const RCP<MueLu::Hierarchy<double,int, int, NO, LMO > > muelu_hierarchy = mueluFactory.CreateHierarchy();
    muelu_hierarchy->GetLevel(0)->Set("A", mueluA);
    muelu_precOp = rcp(new MueLu::EpetraOperator(muelu_hierarchy));
    
    timer.stop();
    if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
      OSTab(out).o() <<"> Creation time = "<<timer.totalElapsedTime()<<" sec\n";

    //     if(paramList_.get())
    //       TEUCHOS_TEST_FOR_EXCEPT(
    //         0!=muelu_precOp->SetParameterList(paramList_->sublist(MueLuSettings_name))
    //         );
  }

  //
  // Attach the epetraFwdOp to the muelu_precOp to guarantee that it will not go away
  //
  set_extra_data(epetraFwdOp, "IFPF::epetraFwdOp", Teuchos::inOutArg(muelu_precOp),
    Teuchos::POST_DESTROY, false);

  //
  // Update the factorization
  //
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the preconditioner ...\n";
  timer.start(true);

  mueluFactory.SetupHierarchy(*muelu_precOp->GetHierarchy());
  
  timer.stop();
  if(out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() <<"=> Setup time = "<<timer.totalElapsedTime()<<" sec\n";
  //
  // Compute the conditioner number estimate if asked
  //

  // ToDo: Implement

  //
  // Attach fwdOp to the muelu_precOp
  //
  set_extra_data(fwdOp, "IFPF::fwdOp", Teuchos::inOutArg(muelu_precOp),
    Teuchos::POST_DESTROY, false);
  //
  // Initialize the output EpetraLinearOp
  //
  if(startingOver) {
    epetra_precOp = rcp(new EpetraLinearOp);
  }
  epetra_precOp->initialize(
    muelu_precOp
    ,epetraFwdOpTransp
    ,EPETRA_OP_APPLY_APPLY_INVERSE
    ,EPETRA_OP_ADJOINT_UNSUPPORTED  // ToDo: Look into adjoints again.
    );
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


void MueLuPreconditionerFactory::uninitializePrec(
  PreconditionerBase<double> *prec,
  Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOp,
  ESupportSolveUse *supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}


// Overridden from ParameterListAcceptor


void MueLuPreconditionerFactory::setParameterList(
  Teuchos::RCP<ParameterList> const& paramList
  )
{
     TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);
//   paramList->validateParameters(*this->getValidParameters(),0);
     paramList_ = paramList;
//   const EMueLuProblemType
//     defaultType = BaseMethodDefaults_validator->getIntegralValue(
//       *paramList_,BaseMethodDefaults_name,BaseMethodDefaults_default
//       );
//   if( MueLu_PROBTYPE_NONE != defaultType ) {
//     const std::string
//       defaultTypeStr = BaseMethodDefaults_valueNames[defaultType];
//     Teuchos::ParameterList defaultParams;
//     TEUCHOS_TEST_FOR_EXCEPTION(
//       0!=ML_Epetra::SetDefaults(defaultTypeStr,defaultParams)
//       ,Teuchos::Exceptions::InvalidParameterValue
//       ,"Error, the ML problem type \"" << defaultTypeStr << "\' is not recongnised by ML!"
//       );
//     // Note, the only way the above exception message could be generated is if
//     // a default problem type was removed from ML_Epetra::SetDefaults(...).
//     // When a new problem type is added to this function, it must be added to
//     // our enum EMLProblemType along with associated objects ...  In other
//     // words, this adapter must be maintained as ML is maintained.  An
//     // alternative design would be to just pass in whatever string to this
//     // function.  This would improve maintainability but it would not generate
//     // very good error messages when a bad string was passed in.  Currenly,
//     // the error message attached to the exception will contain the list of
//     // valid problem types.
//     paramList_->sublist(MueLuSettings_name).setParametersNotAlreadySet(
//       defaultParams);
//   }
// #ifdef TEUCHOS_DEBUG
//   paramList->validateParameters(*this->getValidParameters(),0);
// #endif
}


RCP<ParameterList>
MueLuPreconditionerFactory::getNonconstParameterList()
{
  return paramList_;
}


RCP<ParameterList>
MueLuPreconditionerFactory::unsetParameterList()
{
  Teuchos::RCP<ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


RCP<const ParameterList>
MueLuPreconditionerFactory::getParameterList() const
{
  return paramList_;
}


RCP<const ParameterList>
MueLuPreconditionerFactory::getValidParameters() const
{

//   using Teuchos::rcp;
//   using Teuchos::tuple;
//   using Teuchos::implicit_cast;
//   using Teuchos::rcp_implicit_cast;
//   typedef Teuchos::ParameterEntryValidator PEV;

  static RCP<const ParameterList> validPL;

  if(is_null(validPL)) {

    RCP<ParameterList>
      pl = rcp(new ParameterList());

//     BaseMethodDefaults_validator = rcp(
//       new Teuchos::StringToIntegralParameterEntryValidator<EMueLuProblemType>(
//         BaseMethodDefaults_valueNames,
//         tuple<std::string>(
//           "Do not set any default parameters",
//           "Set default parameters for a smoothed aggregation method",
//           "Set default parameters for a domain decomposition method",
//           "Set default parameters for a domain decomposition method special to ML",
//           "Set default parameters for a Maxwell-type of linear operator"
//           ),
//         tuple<EMueLuProblemType>(
//           MueLu_PROBTYPE_NONE,
//           MueLu_PROBTYPE_SMOOTHED_AGGREGATION,
//           MueLu_PROBTYPE_DOMAIN_DECOMPOSITION,
//           MueLu_PROBTYPE_DOMAIN_DECOMPOSITION_MueLu,
//           MueLu_PROBTYPE_MAXWELL
//           ),
//         BaseMethodDefaults_name
//         )
//       );

//     pl->set(BaseMethodDefaults_name,BaseMethodDefaults_default,
//       "Select the default method type which also sets parameter defaults\n"
//       "in the sublist \"" + MueLuSettings_name + "\"!",
//       rcp_implicit_cast<const PEV>(BaseMethodDefaults_validator)
//       );

// /* 2007/07/02: rabartl:  The statement below should be the correct way to
//  * get the list of valid parameters but it seems to be causing problems so
//  * I am commenting it out for now.
//  */
// /*
//     pl->sublist(
//       MueLuSettings_name, false,
//       "Parameters directly accepted by ML_Epetra interface."
//       ).setParameters(*rcp(ML_Epetra::GetValidMueLuPParameters()));
// */
    
//     {
//       ParameterList &mlSettingsPL = pl->sublist(
//         MueLuSettings_name, false,
//         "Sampling of the parameters directly accepted by ML\n"
//         "This list of parameters is generated by combining all of\n"
//         "the parameters set for all of the default problem types supported\n"
//         "by MueLu.  Therefore, do not assume these parameters are at values that\n"
//         "are reasonable to MueLu.  This list is just to give a sense of some of\n"
//         "the parameters that MueLu accepts.  Consult MueLu documentation on how to\n"
//         "set these parameters.  Also, you can print the parameter list after\n"
//         "it is used and see what defaults where set for each default problem\n"
//         "type.  Warning! the parameters in this sublist are currently *not*\n"
//         "being validated by MueLu!"
//         );
//       //std::cout << "\nMueLuSettings doc before = " << pl->getEntryPtr(MueLuSettings_name)->docString() << "\n";
//       { // Set of valid parameters (but perhaps not accetable values)
//         for (
//           int i = 0;
//           i < implicit_cast<int>(BaseMethodDefaults_valueNames.size());
//           ++i
//           )
//         {
//           ParameterList defaultParams;
//           const std::string defaultTypeStr = BaseMethodDefaults_valueNames[i];
//           if (defaultTypeStr != BaseMethodDefaults_valueNames_none) {
//             TEUCHOS_TEST_FOR_EXCEPTION(
//               0!=ML_Epetra::SetDefaults(defaultTypeStr,defaultParams)
//               ,Teuchos::Exceptions::InvalidParameterValue
//               ,"Error, the ML problem type \"" << defaultTypeStr
//               << "\' is not recongnised by ML!"
//               );
//             mlSettingsPL.setParameters(defaultParams);
//           }
//         }
//       }
//     }

    validPL = pl;

  }

  return validPL;

}


// Public functions overridden from Teuchos::Describable

std::string MueLuPreconditionerFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::MueLuPreconditionerFactory";
  return oss.str();
}

//
//
//

void addMueLuToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder & builder,
                                  const std::string & stratName)
{
  TEUCHOS_TEST_FOR_EXCEPTION(builder.getValidParameters()->sublist("Preconditioner Types").isParameter(stratName),std::logic_error,
                             "MueLu::addMueLuToStratimikosBuilder cannot add \"" + stratName +"\" because it is already included in builder!");
  
  // use default constructor to add Teko::StratimikosFactory
  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,Thyra::MueLuPreconditionerFactory>(),
                                            stratName);
}

} // namespace Thyra
