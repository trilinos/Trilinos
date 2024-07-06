#include "MiniEM_Interpolation.hpp"
#include "Panzer_Interpolation.hpp"


void addInterpolationToRequestHandler(const std::string& name,
                                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                                      const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                                      const std::string& lo_basis_name,
                                      const std::string& ho_basis_name,
                                      Intrepid2::EOperator op,
                                      const bool waitForRequest,
                                      const bool dump,
                                      const size_t worksetSize,
                                      const bool matrixFree) {

  // add interpolation callback to request handler
  reqHandler->addRequestCallback(Teuchos::rcp(new InterpolationRequestCallback(name, linObjFactory, lo_basis_name, ho_basis_name, op, waitForRequest, dump, worksetSize, matrixFree)));
}


InterpolationRequestCallback::
InterpolationRequestCallback(const std::string& name,
                             const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                             const std::string& lo_basis_name,
                             const std::string& ho_basis_name,
                             Intrepid2::EOperator op,
                             const bool waitForRequest,
                             const bool dump,
                             const size_t worksetSize,
                             const bool matrixFree)
  : name_(name), linObjFactory_(linObjFactory), lo_basis_name_(lo_basis_name), ho_basis_name_(ho_basis_name), op_(op), dump_(dump), worksetSize_(worksetSize), matrixFree_(matrixFree)
{
  if (!waitForRequest)
    build();
}

void
InterpolationRequestCallback::
build() {
  std::string timerLabel;
  if (!matrixFree_)
    timerLabel = std::string("Mini-EM: assemble ") + name_;
  else
    timerLabel = std::string("Mini-EM: matrix-free setup ") + name_;
  Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer(timerLabel));
  interp_ = panzer::buildInterpolation(linObjFactory_, lo_basis_name_, ho_basis_name_, op_, worksetSize_, matrixFree_);
  if (matrixFree_) {
    auto constTpetraOp = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal> >(interp_, true);
    auto tpetraOp = Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal> >(constTpetraOp)->getTpetraOperator();
    Teuchos::rcp_dynamic_cast<panzer::MatrixFreeInterpolationOp<double,int,panzer::GlobalOrdinal> >(tpetraOp, true)->setName(name_);
  }
  if (dump_ && !matrixFree_) {
    std::string filename = name_ + ".mm";

    mini_em::writeOut(filename, *interp_);
  }
}


bool
InterpolationRequestCallback::
handlesRequest(const Teko::RequestMesg & rm) {
  std::string name = rm.getName();

  return (name==name_);
}


Teko::LinearOp
InterpolationRequestCallback::
request(const Teko::RequestMesg & rm) {
  TEUCHOS_ASSERT(handlesRequest(rm));
  std::string name = rm.getName();

  if(name==name_) {
    if (interp_.is_null()) {
      build();
    }
    return interp_;
  } else
    TEUCHOS_ASSERT(false);
}

void
InterpolationRequestCallback::
preRequest(const Teko::RequestMesg & rm) {
    // checking for its existance is as good as pre requesting
    TEUCHOS_ASSERT(handlesRequest(rm));
  }
