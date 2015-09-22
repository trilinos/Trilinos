// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_MODEL_EVALUATOR_FACTORY_IMPL_HPP
#define PANZER_MODEL_EVALUATOR_FACTORY_IMPL_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"

namespace panzer {
  
  template <typename ScalarT, typename LO, typename GO>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > 
  buildModelEvaluator() const
  { 
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > me;
    
    std::string type = getMyParamList.get<std::string>("Model Evaluator Type");
    
    if (type == "Epetra") {
      
      Teuchos::RCP<EpetraExt::ModelEvaluator> epetraModel = 
	Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb, linObjFactory));
    
      Teuchos::RCP<Thyra::EpetraModelEvaluator>
	epetraThyraModel = rcp(new ::Thyra::EpetraModelEvaluator());
      epetraThyraModel->initialize(epetraModel,lowsFactory);
      Teuchos::RCP< ::Thyra::ModelEvaluator<double> > thyraModel = 
	epetraThyraModel;

    }
    else if (type == "Tpetra") {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			 "Tpetra version not supported yet, use Epetra!");
    }
    


  }

}

#endif
