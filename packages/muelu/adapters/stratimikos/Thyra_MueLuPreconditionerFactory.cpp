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

#include <Epetra_CrsMatrix.h>

#include <Teuchos_AbstractFactoryStd.hpp>
#include <Teuchos_dyn_cast.hpp>
#include <Teuchos_implicit_cast.hpp>
#include <Teuchos_iostream_helpers.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_StaticSetupMacro.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ValidatorXMLConverterDB.hpp>

#include <Thyra_DefaultPreconditioner.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraOperatorViewExtractorStd.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_CreateEpetraPreconditioner.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_Level.hpp"

namespace {

  const std::string MueLuSettings_name = "MueLu Settings";

} // namespace


namespace Thyra {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Constructors/initializers/accessors

  MueLuPreconditionerFactory::MueLuPreconditionerFactory() :
      epetraFwdOpViewExtractor_(rcp(new EpetraOperatorViewExtractorStd())),
      paramList_               (rcp(new ParameterList()))
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
    return rcp(new DefaultPreconditioner<double>());
  }


  void MueLuPreconditionerFactory::initializePrec(const Teuchos::RCP<const LinearOpSourceBase<double> >& fwdOpSrc,
                                                  PreconditionerBase<double>* prec,
                                                  const ESupportSolveUse supportSolveUse
                                                 ) const
  {
    using Teuchos::outArg;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_implicit_cast;
    using Teuchos::rcp_const_cast;
    using Teuchos::set_extra_data;

    typedef KokkosClassic::DefaultNode::DefaultNodeType NO;

    Teuchos::RCP<const LinearOpBase<double> > fwdOp = fwdOpSrc->getOp();
#ifdef _DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(fwdOp.get() == NULL);
    TEUCHOS_TEST_FOR_EXCEPT(prec        == NULL);
#endif

    ParameterList& paramList = *paramList_;

    // Unwrap and get the forward Epetra_Operator object
    Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
    EOpTransp epetraFwdOpTransp;
    EApplyEpetraOpAs epetraFwdOpApplyAs;
    EAdjointEpetraOp epetraFwdOpAdjointSupport;
    double epetraFwdOpScalar;
    epetraFwdOpViewExtractor_->getEpetraOpView(fwdOp, outArg(epetraFwdOp), outArg(epetraFwdOpTransp), outArg(epetraFwdOpApplyAs),
                                               outArg(epetraFwdOpAdjointSupport), outArg(epetraFwdOpScalar));
    // Validate what we get is what we need
    RCP<const Epetra_CrsMatrix> epetraFwdCrsMat = rcp_dynamic_cast<const Epetra_CrsMatrix>(epetraFwdOp, true);
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

    // Do the magic (init/setup/reuse)
    // FIXME: the check for starting over needs more work
    // For instance, what should happen if a user called the first setup with
    // one parameter list, and the second setup with a different one?
    const bool startingOver = (muelu_precOp.is_null() || !paramList.isParameter("reuse: type") || paramList.get<std::string>("reuse: type") == "none");
    if (startingOver)
      muelu_precOp = MueLu::CreateEpetraPreconditioner(rcp_const_cast<Epetra_CrsMatrix>(epetraFwdCrsMat), paramList);
    else
      MueLu::ReuseEpetraPreconditioner(rcp_const_cast<Epetra_CrsMatrix>(epetraFwdCrsMat), *muelu_precOp);

    // Attach epetraFwdOp and fwdOp to the muelu_precOp to guarantee that it will not go away
    set_extra_data(epetraFwdOp, "IFPF::epetraFwdOp", Teuchos::inOutArg(muelu_precOp), Teuchos::POST_DESTROY, false);
    set_extra_data(fwdOp,       "IFPF::fwdOp",       Teuchos::inOutArg(muelu_precOp), Teuchos::POST_DESTROY, false);

    // Compute the conditioner number estimate if asked
    // TODO: Implement

    // Initialize the output EpetraLinearOp
    if (startingOver)
      epetra_precOp = rcp(new EpetraLinearOp);
    epetra_precOp->initialize(muelu_precOp,
                              epetraFwdOpTransp,
                              EPETRA_OP_APPLY_APPLY_INVERSE,
                              EPETRA_OP_ADJOINT_UNSUPPORTED  // TODO: Look into adjoints again.
                             );

    // Initialize the preconditioner
    defaultPrec->initializeUnspecified(rcp_implicit_cast<LinearOpBase<double> >(epetra_precOp));
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
