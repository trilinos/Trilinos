// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
   else if(name.substr(0,15)=="SchurComplement") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("SchurComplement " + name.substr(16,name.length()-16)+" Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       mini_em::writeOut("SchurComplement.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else if(name.substr(0,20)=="DarcySchurComplement") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("DarcySchurComplement " + name.substr(21,name.length()-21)+" Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       mini_em::writeOut("DarcySchurComplement.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else if(name.substr(0,24)=="ProjectedSchurComplement") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("ProjectedSchurComplement " + name.substr(25,name.length()-25)+" Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       mini_em::writeOut("ProjectedSchurComplement.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
   }
   else if(name.substr(0,29)=="ProjectedDarcySchurComplement") {
     loc = Teuchos::rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(gedc_->getDataObject("ProjectedDarcySchurComplement " + name.substr(30,name.length()-30)+" Scatter Container"),true)->getGlobalLOC();

     if (matrix_output)
       mini_em::writeOut("ProjectedDarcySchurComplement.mm", *Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(loc,true)->get_A_th());
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
   else if(name.substr(0,15)=="SchurComplement") {
     if(gedc_->containsDataObject("SchurComplement " + name.substr(16,name.length()-16)+" Scatter Container"))
       return true;
   }
   else if(name.substr(0,20)=="DarcySchurComplement") {
     if(gedc_->containsDataObject("DarcySchurComplement " + name.substr(21,name.length()-21)+" Scatter Container"))
       return true;
   }
   else if(name.substr(0,24)=="ProjectedSchurComplement") {
     if(gedc_->containsDataObject("ProjectedSchurComplement " + name.substr(25,name.length()-25)+" Scatter Container"))
       return true;
   }
   else if(name.substr(0,29)=="ProjectedDarcySchurComplement") {
     if(gedc_->containsDataObject("ProjectedDarcySchurComplement " + name.substr(30,name.length()-30)+" Scatter Container"))
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
