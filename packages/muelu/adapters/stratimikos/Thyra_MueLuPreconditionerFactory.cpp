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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
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
      : epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
  {}


  // Overridden from PreconditionerFactoryBase

  bool MueLuPreconditionerFactory::isCompatible(const LinearOpSourceBase<double> &fwdOpSrc) const {
    using Teuchos::outArg;
    Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
    EOpTransp epetraFwdOpTransp;
    EApplyEpetraOpAs epetraFwdOpApplyAs;
    EAdjointEpetraOp epetraFwdOpAdjointSupport;
    double epetraFwdOpScalar;
    Teuchos::RCP<const LinearOpBase<double> > fwdOp = fwdOpSrc.getOp();
    epetraFwdOpViewExtractor_->getEpetraOpView(fwdOp,
                                               outArg(epetraFwdOp),
                                               outArg(epetraFwdOpTransp),
                                               outArg(epetraFwdOpApplyAs),
                                               outArg(epetraFwdOpAdjointSupport),
                                               outArg(epetraFwdOpScalar)
                                              );
    return (!dynamic_cast<const Epetra_CrsMatrix*>(&*epetraFwdOp) ? false : true);
  }

  RCP<PreconditionerBase<double> > MueLuPreconditionerFactory::createPrec() const {
    return Teuchos::rcp(new DefaultPreconditioner<double>());
  }


  void MueLuPreconditionerFactory::initializePrec(const Teuchos::RCP<const LinearOpSourceBase<double> >& fwdOpSrc,
                                                  PreconditionerBase<double>* prec,
                                                  const ESupportSolveUse supportSolveUse
                                                 ) const
  {
    using Teuchos::outArg;
    using Teuchos::OSTab;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::set_extra_data;
    using Teuchos::implicit_cast;

    typedef KokkosClassic::DefaultNode::DefaultNodeType NO;

    RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));

    const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
    if (out.get() && implicit_cast<int>(verbLevel) > implicit_cast<int>(Teuchos::VERB_LOW))
      out = this->getOStream();
    Teuchos::OSTab tab(out);

    *out << "\nEntering Thyra::MueLuPreconditionerFactory::initializePrec(...) ..." << std::endl;

    Teuchos::Time totalTimer(""), timer("");
    totalTimer.start(true);

    Teuchos::RCP<const LinearOpBase<double> > fwdOp = fwdOpSrc->getOp();
#ifdef _DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(fwdOp.get() == NULL);
    TEUCHOS_TEST_FOR_EXCEPT(prec        == NULL);
#endif

    // Unwrap and get the forward Epetra_Operator object
    Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
    EOpTransp epetraFwdOpTransp;
    EApplyEpetraOpAs epetraFwdOpApplyAs;
    EAdjointEpetraOp epetraFwdOpAdjointSupport;
    double epetraFwdOpScalar;
    epetraFwdOpViewExtractor_->getEpetraOpView(fwdOp, outArg(epetraFwdOp), outArg(epetraFwdOpTransp), outArg(epetraFwdOpApplyAs),
                                               outArg(epetraFwdOpAdjointSupport), outArg(epetraFwdOpScalar));
    // Validate what we get is what we need
    RCP<const Epetra_CrsMatrix> epetraFwdCrsMat = rcp_dynamic_cast<const Epetra_CrsMatrix>(epetraFwdOp,true);
    TEUCHOS_TEST_FOR_EXCEPTION(epetraFwdOpApplyAs != EPETRA_OP_APPLY_APPLY, std::logic_error,
                               "Error, incorrect apply mode for an Epetra_CrsMatrix");

    // Get the concrete preconditioner object
    DefaultPreconditioner<double> *defaultPrec = &Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);

    // Get the EpetraLinearOp object that is used to implement the preconditoner linear op
    RCP<EpetraLinearOp> epetra_precOp = rcp_dynamic_cast<EpetraLinearOp>(defaultPrec->getNonconstUnspecifiedPrecOp(), true);

    // Get the embedded MueLu::EpetraOperator object if it exists
    Teuchos::RCP<MueLu::EpetraOperator> muelu_precOp;
    if (epetra_precOp.get())
      muelu_precOp = rcp_dynamic_cast<MueLu::EpetraOperator>(epetra_precOp->epetra_op(), true);

    // Get the attached forward operator if it exists and make sure that it matches
    if (muelu_precOp != Teuchos::null) {
      // TODO
      //     // Get the forward operator and make sure that it matches what is
      //     // already being used!
      //     const Epetra_CrsMatrix & rm = muelu_precOp->CrsMatrix();

      //     TEUCHOS_TEST_FOR_EXCEPTION(
      //        &rm!=&*epetraFwdRowMat, std::logic_error
      //        ,"MueLu requires Epetra_RowMatrix to be the same for each initialization of the preconditioner"
      //        );
    }

    MueLu::ParameterListInterpreter<double, int, int, NO> mueluFactory(*paramList_);

    // Perform initialization if needed
    const bool startingOver = (muelu_precOp.get() == NULL);
    if (startingOver) {
      *out << "\nCreating the initial MueLu::EpetraOperator object..." << std::endl;
      timer.start(true);
      // Create the initial preconditioner: DO NOT compute it yet

      // Turns a Epetra_CrsMatrix into a Xpetra::Matrix
      RCP<Epetra_CrsMatrix> epetraFwdCrsMatNonConst = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(epetraFwdCrsMat); // !! TODO: MueLu interface should accept const matrix as input.

      RCP<Xpetra::CrsMatrix<double, int, int, NO> > mueluAcrs = rcp(new Xpetra::EpetraCrsMatrix(epetraFwdCrsMatNonConst)); //TODO: should not be needed
      RCP<Xpetra::Matrix <double, int, int, NO> >   mueluA    = rcp(new Xpetra::CrsMatrixWrap<double, int, int, NO>(mueluAcrs));

      const RCP<MueLu::Hierarchy<double,int, int, NO > > muelu_hierarchy = mueluFactory.CreateHierarchy();
      muelu_hierarchy->GetLevel(0)->Set("A", mueluA);
      muelu_precOp = rcp(new MueLu::EpetraOperator(muelu_hierarchy));

      timer.stop();
      OSTab(out).o() << "> Creation time = " << timer.totalElapsedTime() << " sec\n";
    }

    // Attach the epetraFwdOp to the muelu_precOp to guarantee that it will not go away
    set_extra_data(epetraFwdOp, "IFPF::epetraFwdOp", Teuchos::inOutArg(muelu_precOp), Teuchos::POST_DESTROY, false);

    // Update the factorization
    *out << "\nComputing the preconditioner ..." << std::endl;
    timer.start(true);

    mueluFactory.SetupHierarchy(*muelu_precOp->GetHierarchy());

    timer.stop();
    OSTab(out).o() << "=> Setup time = " << timer.totalElapsedTime() << " sec\n";

    // Compute the conditioner number estimate if asked
    // TODO: Implement

    // Attach fwdOp to the muelu_precOp
    set_extra_data(fwdOp, "IFPF::fwdOp", Teuchos::inOutArg(muelu_precOp), Teuchos::POST_DESTROY, false);

    // Initialize the output EpetraLinearOp
    if (startingOver)
      epetra_precOp = rcp(new EpetraLinearOp);
    epetra_precOp->initialize(muelu_precOp,
                              epetraFwdOpTransp,
                              EPETRA_OP_APPLY_APPLY_INVERSE,
                              EPETRA_OP_ADJOINT_UNSUPPORTED  // ToDo: Look into adjoints again.
                             );

    // Initialize the preconditioner
    defaultPrec->initializeUnspecified(Teuchos::rcp_implicit_cast<LinearOpBase<double> >(epetra_precOp));
    totalTimer.stop();

    *out << "\nTotal time in MLPreconditionerFactory = " << totalTimer.totalElapsedTime() << " sec" << std::endl;
    *out << "\nLeaving Thyra::MLPreconditionerFactory::initializePrec(...) ..." << std::endl;
  }

  void MueLuPreconditionerFactory::uninitializePrec(PreconditionerBase<double> *prec,
                                                    Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOp,
                                                    ESupportSolveUse *supportSolveUse) const
  {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  // Overridden from ParameterListAcceptor
  void MueLuPreconditionerFactory::setParameterList(const Teuchos::RCP<ParameterList>& paramList) {
    TEUCHOS_TEST_FOR_EXCEPT(paramList.get() == NULL);
    paramList_ = paramList;
  }

  RCP<ParameterList> MueLuPreconditionerFactory::unsetParameterList() {
    Teuchos::RCP<ParameterList> paramList = paramList_;
    paramList_ = Teuchos::null;
    return paramList;
  }

  RCP<ParameterList> MueLuPreconditionerFactory::getNonconstParameterList() {
    return paramList_;
  }

  RCP<const ParameterList> MueLuPreconditionerFactory::getParameterList() const {
    return paramList_;
  }

  RCP<const ParameterList> MueLuPreconditionerFactory::getValidParameters() const {
    static RCP<const ParameterList> validPL;

    if (is_null(validPL))
      validPL = rcp(new ParameterList());

    return validPL;
  }

  // Public functions overridden from Teuchos::Describable
  std::string MueLuPreconditionerFactory::description() const {
    std::ostringstream oss;
    oss << "Thyra::MueLuPreconditionerFactory";
    return oss.str();
  }


  void addMueLuToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder & builder, const std::string & stratName) {
    TEUCHOS_TEST_FOR_EXCEPTION(builder.getValidParameters()->sublist("Preconditioner Types").isParameter(stratName), std::logic_error,
                               "MueLu::addMueLuToStratimikosBuilder cannot add \"" + stratName +"\" because it is already included in builder!");

    // use default constructor to add Teko::StratimikosFactory
    builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>, Thyra::MueLuPreconditionerFactory>(), stratName);
  }

} // namespace Thyra
