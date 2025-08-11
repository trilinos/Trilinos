// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_StratimikosFactory.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "Thyra_DefaultPreconditioner.hpp"

#include "Teko_InverseLibrary.hpp"
#include "Teko_Preconditioner.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_ReorderedLinearOp.hpp"

namespace Teko {

using Teuchos::ParameterList;
using Teuchos::RCP;

// hide stuff
namespace {
// Simple preconditioner class that adds a counter
class StratimikosFactoryPreconditioner : public Thyra::DefaultPreconditioner<double> {
 public:
  StratimikosFactoryPreconditioner() : iter_(0) {}

  inline void incrIter() { iter_++; }
  inline std::size_t getIter() { return iter_; }

 private:
  StratimikosFactoryPreconditioner(const StratimikosFactoryPreconditioner &);

  std::size_t iter_;
};

// factory used to initialize the Teko::StratimikosFactory
// user data
class TekoFactoryBuilder
    : public Teuchos::AbstractFactory<Thyra::PreconditionerFactoryBase<double> > {
 public:
  TekoFactoryBuilder(const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> &builder,
                     const Teuchos::RCP<Teko::RequestHandler> &rh)
      : builder_(builder), requestHandler_(rh) {}
  Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > create() const {
    return Teuchos::rcp(new StratimikosFactory(builder_, requestHandler_));
  }

 private:
  Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> builder_;
  Teuchos::RCP<Teko::RequestHandler> requestHandler_;
};
}  // namespace

StratimikosFactory::StratimikosFactory(
    const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> &builder,
    const Teuchos::RCP<Teko::RequestHandler> &rh)
    : builder_(builder) {
  setRequestHandler(rh);
}

// Overridden from PreconditionerFactoryBase
bool StratimikosFactory::isCompatible(
    const Thyra::LinearOpSourceBase<double> & /* fwdOpSrc */) const {
  using Teuchos::outArg;

  TEUCHOS_ASSERT(false);  // what you doing?

  /*
     Teuchos::RCP<const Epetra_Operator> epetraFwdOp;
     Thyra::EOpTransp epetraFwdOpTransp;
     Thyra::EApplyEpetraOpAs epetraFwdOpApplyAs;
     Thyra::EAdjointEpetraOp epetraFwdOpAdjointSupport;
     double epetraFwdOpScalar;

     // check to make sure this is an epetra CrsMatrix
     Teuchos::RCP<const LinearOpBase<double> > fwdOp = fwdOpSrc.getOp();
     epetraFwdOpViewExtractor_->getEpetraOpView(
             fwdOp,
             outArg(epetraFwdOp),outArg(epetraFwdOpTransp),
             outArg(epetraFwdOpApplyAs),
             outArg(epetraFwdOpAdjointSupport),
             outArg(epetraFwdOpScalar));

     if( !dynamic_cast<const Epetra_CrsMatrix*>(&*epetraFwdOp) )
        return false;
  */

  return true;
}

bool StratimikosFactory::applySupportsConj(Thyra::EConj /* conj */) const { return false; }

bool StratimikosFactory::applyTransposeSupportsConj(Thyra::EConj /* conj */) const {
  return false;  // See comment below
}

Teuchos::RCP<Thyra::PreconditionerBase<double> > StratimikosFactory::createPrec() const {
  return Teuchos::rcp(new StratimikosFactoryPreconditioner());
}

void StratimikosFactory::initializePrec(
    const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOpSrc,
    Thyra::PreconditionerBase<double> *prec, const Thyra::ESupportSolveUse supportSolveUse) const {
  Teuchos::RCP<const LinearOpBase<double> > fwdOp = fwdOpSrc->getOp();

  initializePrec_Thyra(fwdOpSrc, prec, supportSolveUse);
}

void StratimikosFactory::initializePrec_Thyra(
    const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > &fwdOpSrc,
    Thyra::PreconditionerBase<double> *prec, const Thyra::ESupportSolveUse /* supportSolveUse */
) const {
  using Teuchos::implicit_cast;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::LinearOpBase;

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const RCP<Teuchos::FancyOStream> out     = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  bool mediumVerbosity =
      (out.get() && implicit_cast<int>(verbLevel) > implicit_cast<int>(Teuchos::VERB_LOW));

  Teuchos::OSTab tab(out);
  if (mediumVerbosity)
    *out << "\nEntering Teko::StratimikosFactory::initializePrec_Thyra(...) ...\n";

  Teuchos::RCP<const LinearOpBase<double> > fwdOp = fwdOpSrc->getOp();

  // Get the concrete preconditioner object
  StratimikosFactoryPreconditioner &defaultPrec =
      Teuchos::dyn_cast<StratimikosFactoryPreconditioner>(*prec);
  Teuchos::RCP<LinearOpBase<double> > prec_Op = defaultPrec.getNonconstUnspecifiedPrecOp();

  // Permform initialization if needed
  const bool startingOver = (prec_Op == Teuchos::null);
  if (startingOver) {
    // start over
    invLib_     = Teuchos::null;
    invFactory_ = Teuchos::null;

    if (mediumVerbosity) *out << "\nCreating the initial Teko Operator object...\n";

    timer.start(true);

    // build library, and set request handler (user defined!)
    invLib_ = Teko::InverseLibrary::buildFromParameterList(
        paramList_->sublist("Inverse Factory Library"), builder_);
    invLib_->setRequestHandler(reqHandler_);

    // build preconditioner factory
    invFactory_ = invLib_->getInverseFactory(paramList_->get<std::string>("Inverse Type"));

    timer.stop();
    if (mediumVerbosity)
      Teuchos::OSTab(out).o() << "> Creation time = " << timer.totalElapsedTime() << " sec\n";
  }

  if (mediumVerbosity) *out << "\nComputing the preconditioner ...\n";

  // setup reordering if required
  std::string reorderType = paramList_->get<std::string>("Reorder Type");
  if (reorderType != "") {
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > blkFwdOp =
        Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<double> >(fwdOp, true);
    RCP<const Teko::BlockReorderManager> brm = Teko::blockedReorderFromString(reorderType);
    Teko::LinearOp blockedFwdOp              = Teko::buildReorderedLinearOp(*brm, blkFwdOp);

    if (prec_Op == Teuchos::null) {
      Teko::ModifiableLinearOp reorderedPrec = Teko::buildInverse(*invFactory_, blockedFwdOp);
      prec_Op = Teuchos::rcp(new ReorderedLinearOp(brm, reorderedPrec));
    } else {
      Teko::ModifiableLinearOp reorderedPrec =
          Teuchos::rcp_dynamic_cast<ReorderedLinearOp>(prec_Op, true)->getBlockedOp();
      Teko::rebuildInverse(*invFactory_, blockedFwdOp, reorderedPrec);
    }
  } else {
    // no reordering required
    if (prec_Op == Teuchos::null)
      prec_Op = Teko::buildInverse(*invFactory_, fwdOp);
    else
      Teko::rebuildInverse(*invFactory_, fwdOp, prec_Op);
  }

  // construct preconditioner
  timer.stop();

  if (mediumVerbosity)
    Teuchos::OSTab(out).o() << "=> Preconditioner construction time = " << timer.totalElapsedTime()
                            << " sec\n";

  // Initialize the preconditioner
  defaultPrec.initializeUnspecified(prec_Op);
  defaultPrec.incrIter();
  totalTimer.stop();

  if (out.get() && implicit_cast<int>(verbLevel) >= implicit_cast<int>(Teuchos::VERB_LOW))
    *out << "\nTotal time in Teko::StratimikosFactory = " << totalTimer.totalElapsedTime()
         << " sec\n";
  if (mediumVerbosity)
    *out << "\nLeaving Teko::StratimikosFactory::initializePrec_Thyra(...) ...\n";
}

void StratimikosFactory::uninitializePrec(
    Thyra::PreconditionerBase<double> * /* prec */,
    Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > * /* fwdOp */,
    Thyra::ESupportSolveUse * /* supportSolveUse */
) const {
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

// Overridden from ParameterListAcceptor

void StratimikosFactory::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const &paramList) {
  TEUCHOS_TEST_FOR_EXCEPT(paramList.get() == NULL);

  paramList->validateParametersAndSetDefaults(*this->getValidParameters(), 0);
  paramList_ = paramList;
}

Teuchos::RCP<Teuchos::ParameterList> StratimikosFactory::getNonconstParameterList() {
  return paramList_;
}

Teuchos::RCP<Teuchos::ParameterList> StratimikosFactory::unsetParameterList() {
  Teuchos::RCP<ParameterList> _paramList = paramList_;
  paramList_                             = Teuchos::null;
  return _paramList;
}

Teuchos::RCP<const Teuchos::ParameterList> StratimikosFactory::getParameterList() const {
  return paramList_;
}

Teuchos::RCP<const Teuchos::ParameterList> StratimikosFactory::getValidParameters() const {
  using Teuchos::implicit_cast;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::tuple;

  static RCP<const ParameterList> validPL;

  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

    pl->set("Test Block Operator", false,
            "If Stratiikos/Teko is used to break an operator into its parts,\n"
            "then setting this parameter to true will compare applications of the\n"
            "segregated operator to the original operator.");
    pl->set("Write Block Operator", false,
            "Write out the segregated operator to disk with the name \"block-?_xx\"");
    pl->set("Strided Blocking", "1",
            "Assuming that the user wants Strided blocking, break the operator into\n"
            "blocks. The syntax can be thought to be associated with the solution\n"
            "vector. For example if your variables are [u v w p T], and we want [u v w]\n"
            "blocked together, and p and T separate then the relevant string is \"3 1 1\".\n"
            "Meaning put the first 3 unknowns per node together and separate the v and w\n"
            "components.");
    pl->set("Reorder Type", "",
            "This specifies how the blocks are reordered for use in the preconditioner.\n"
            "For example, assume the linear system is generated from 3D Navier-Stokes\n"
            "with an energy equation, yielding the unknowns [u v w p T]. If the\n"
            "\"Strided Blocking\" string is \"3 1 1\", then setting this parameter to\n"
            "\"[2 [0 1]]\" will reorder the blocked operator so its nested with the\n"
            "velocity and pressure forming an inner two-by-two block, and then the\n"
            "temperature unknowns forming a two-by-two system with the velocity-pressure\n"
            "block.");
    std::string defaultInverseType = "";
#if defined(Teko_ENABLE_Amesos)
    defaultInverseType = "Amesos";
#elif defined(Teko_ENABLE_Amesos2)
    defaultInverseType = "Amesos2";
#endif
    pl->set("Inverse Type", defaultInverseType,
            "The type of inverse operator the user wants. This can be one of the defaults\n"
            "from Stratimikos, or a Teko preconditioner defined in the\n"
            "\"Inverse Factory Library\".");
    pl->sublist("Inverse Factory Library", false, "Definition of Teko preconditioners.");

    validPL = pl;
  }

  return validPL;
}

std::string StratimikosFactory::description() const {
  std::ostringstream oss;
  oss << "Teko::StratimikosFactory";
  return oss.str();
}

void addTekoToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder &builder,
                                 const std::string &stratName) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      builder.getValidParameters()->sublist("Preconditioner Types").isParameter(stratName),
      std::logic_error,
      "Teko::addTekoToStratimikosBuilder cannot add \"" + stratName +
          "\" because it is already included in builder!");

  Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> builderCopy =
      Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder(builder));

  // use default constructor to add Teko::StratimikosFactory
  Teuchos::RCP<TekoFactoryBuilder> tekoFactoryBuilder =
      Teuchos::rcp(new TekoFactoryBuilder(builderCopy, Teuchos::null));
  builder.setPreconditioningStrategyFactory(tekoFactoryBuilder, stratName);
}

void addTekoToStratimikosBuilder(Stratimikos::DefaultLinearSolverBuilder &builder,
                                 const Teuchos::RCP<Teko::RequestHandler> &rh,
                                 const std::string &stratName) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      builder.getValidParameters()->sublist("Preconditioner Types").isParameter(stratName),
      std::logic_error,
      "Teko::addTekoToStratimikosBuilder cannot add \"" + stratName +
          "\" because it is already included in builder!");

  Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> builderCopy =
      Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder(builder));

  // build an instance of a Teuchos::AbsractFactory<Thyra::PFB> so request handler is passed onto
  // the resulting StratimikosFactory
  Teuchos::RCP<TekoFactoryBuilder> tekoFactoryBuilder =
      Teuchos::rcp(new TekoFactoryBuilder(builderCopy, rh));
  builder.setPreconditioningStrategyFactory(tekoFactoryBuilder, stratName);
}

}  // namespace Teko
