#ifndef _MiniEM_Interpolation_hpp_
#define _MiniEM_Interpolation_hpp_

#include "Teko_Utilities.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestCallback.hpp"

#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif

#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#endif
#include "MiniEM_Utils.hpp"


void addInterpolationToRequestHandler(const std::string& name,
                                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                                      const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                                      const std::string& lo_basis_name,
                                      const std::string& ho_basis_name,
                                      Intrepid2::EOperator op=Intrepid2::OPERATOR_VALUE,
                                      const bool waitForRequest=true,
                                      const bool dump=false,
                                      const size_t worksetSize=1000,
                                      const bool matrixFree=false);


class InterpolationRequestCallback : public Teko::RequestCallback<Teko::LinearOp> {
private:

  std::string name_;
  const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;
  const std::string lo_basis_name_;
  const std::string ho_basis_name_;
  Intrepid2::EOperator op_;
  const bool dump_;
  Teko::LinearOp interp_;
  const size_t worksetSize_;
  const bool matrixFree_;

public:

  InterpolationRequestCallback(const std::string& name,
                               const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                               const std::string& lo_basis_name,
                               const std::string& ho_basis_name,
                               Intrepid2::EOperator op=Intrepid2::OPERATOR_VALUE,
                               const bool waitForRequest=true,
                               const bool dump=false,
                               const size_t worksetSize=1000,
                               const bool matrixFree=false);

  void build();

  bool handlesRequest(const Teko::RequestMesg & rm);

  Teko::LinearOp request(const Teko::RequestMesg & rm);

  void preRequest(const Teko::RequestMesg & rm);
};



#endif
