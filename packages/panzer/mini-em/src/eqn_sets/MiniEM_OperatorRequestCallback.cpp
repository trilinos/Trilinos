#include "MiniEM_OperatorRequestCallback.hpp"

#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Panzer_NodeType.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {

OperatorRequestCallback::
OperatorRequestCallback(const Teuchos::RCP<const panzer::GlobalEvaluationDataContainer> & gedc, const bool & matrix_out)
  : gedc_(gedc), matrix_output(matrix_out)
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
       writeOut("MassMatrix.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else if(name=="Weak Gradient") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("Weak Gradient Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       writeOut("WeakGradient.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else {
     // bad news!
     TEUCHOS_ASSERT(false);
   }

   return Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th();
}

void OperatorRequestCallback::writeOut(const std::string & s,const Thyra::LinearOpBase<double> & op) const
{
  using Teuchos::RCP;
  using NT = panzer::TpetraNodeType;
  const RCP<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> >(Teuchos::rcpFromRef(op));
  if(tOp != Teuchos::null) {
    const RCP<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >(tOp->getConstTpetraOperator(),true);
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >::writeSparseFile(s.c_str(),crsOp);
  }
}

bool OperatorRequestCallback::handlesRequest(const Teko::RequestMesg & rm)
{
   std::string name = rm.getName();

   if(name.substr(0,11)=="Mass Matrix") {
     if(gedc_->containsDataObject("Mass Matrix " + name.substr(12,name.length()-12)+" Scatter Container"))
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
