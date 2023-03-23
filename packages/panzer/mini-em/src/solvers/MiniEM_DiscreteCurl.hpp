#ifndef _MiniEM_DiscreteCurl_hpp_
#define _MiniEM_DiscreteCurl_hpp_

#include "Teko_Utilities.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestCallback.hpp"

#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif

#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_IntrepidOrientation.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#endif

class CurlRequestCallback : public Teko::RequestCallback<Teko::LinearOp> {
private:

  Teko::LinearOp curl_;

public:

  CurlRequestCallback(const Teko::LinearOp & curl)
     : curl_(curl) {};

  bool handlesRequest(const Teko::RequestMesg & rm)
  {
    std::string name = rm.getName();

    return (name=="Discrete Curl");
  };

  Teko::LinearOp request(const Teko::RequestMesg & rm)
  {
    TEUCHOS_ASSERT(handlesRequest(rm));
    std::string name = rm.getName();

    if(name=="Discrete Curl")
      return curl_;
    else
      TEUCHOS_ASSERT(false);
  };

  void preRequest(const Teko::RequestMesg & rm)
  {
    // checking for its existance is as good as pre requesting
    TEUCHOS_ASSERT(handlesRequest(rm));
  };
};


void addDiscreteCurlToRequestHandler(
       const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
       const Teuchos::RCP<Teko::RequestHandler> & reqHandler);

#endif
