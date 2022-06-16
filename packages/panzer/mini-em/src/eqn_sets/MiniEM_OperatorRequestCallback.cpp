#include "MiniEM_OperatorRequestCallback.hpp"

#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Panzer_NodeType.hpp"

#include "MiniEM_Utils.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {

OperatorRequestCallback::
OperatorRequestCallback(const Teuchos::RCP<const panzer::GlobalEvaluationDataContainer> & gedc,
                        const bool& dump=false)
  : gedc_(gedc), matrix_output(dump)
{

}

// RequestCallback member functions
///////////////////////////////////////////

Teko::LinearOp OperatorRequestCallback::request(const Teko::RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm));
   std::string name = rm.getName();

   Teuchos::RCP<panzer::LinearObjContainer> loc;
   if(name.substr(0,11)=="Mass Matrix") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("Mass Matrix " + name.substr(12,name.length()-12)+" Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       mini_em::writeOut(name+".mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else if(name.substr(0,9)=="Curl Curl") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("Curl Curl " + name.substr(10,name.length()-10)+" Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       mini_em::writeOut("CurlCurl.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else if(name=="Weak Gradient") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("Weak Gradient Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       mini_em::writeOut("WeakGradient.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else {
     // bad news!
     TEUCHOS_ASSERT(false);
   }

   return Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th();
}


bool OperatorRequestCallback::handlesRequest(const Teko::RequestMesg & rm)
{
   std::string name = rm.getName();

   if(name.substr(0,11)=="Mass Matrix") {
     if(gedc_->containsDataObject("Mass Matrix " + name.substr(12,name.length()-12)+" Scatter Container"))
       return true;
   }
   else if(name.substr(0,9)=="Curl Curl") {
     if(gedc_->containsDataObject("Curl Curl " + name.substr(10,name.length()-10)+" Scatter Container"))
       return true;
   }
   else if(name=="Weak Gradient") {
     if(gedc_->containsDataObject("Weak Gradient Scatter Container"))
       return true;
   }
   return false;
}

void OperatorRequestCallback::preRequest(const Teko::RequestMesg & rm)
{
   // checking for its existence is as good as pre-requesting
   TEUCHOS_ASSERT(handlesRequest(rm));
}

}
