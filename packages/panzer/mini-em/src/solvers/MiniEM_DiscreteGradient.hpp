#ifndef _MiniEM_DiscreteGradient_hpp_
#define _MiniEM_DiscreteGradient_hpp_

#include "Teko_Utilities.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestCallback.hpp"

#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif

#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidOrientation.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#endif
#include "Tpetra_Import.hpp"

class GradientRequestCallback : public Teko::RequestCallback<Teko::LinearOp> {
private:

  Teko::LinearOp gradient_;

public:

  GradientRequestCallback(const Teko::LinearOp & gradient)
     : gradient_(gradient) {};

  bool handlesRequest(const Teko::RequestMesg & rm)
  {
    std::string name = rm.getName();

    return (name=="Discrete Gradient");
  };

  Teko::LinearOp request(const Teko::RequestMesg & rm)
  {
    TEUCHOS_ASSERT(handlesRequest(rm));
    std::string name = rm.getName();

    if(name=="Discrete Gradient")
      return gradient_;
    else
      TEUCHOS_ASSERT(false);
  };

  void preRequest(const Teko::RequestMesg & rm)
  {
    // checking for its existance is as good as pre requesting
    TEUCHOS_ASSERT(handlesRequest(rm));
  };
};

void addDiscreteGradientToRequestHandler(
       const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
       const Teuchos::RCP<Teko::RequestHandler> & reqHandler);

#endif
