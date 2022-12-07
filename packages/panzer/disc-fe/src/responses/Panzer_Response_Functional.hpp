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

#ifndef __Panzer_Response_Functional_hpp__
#define __Panzer_Response_Functional_hpp__

#include <string>


#include <mpi.h> // need for comm

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "PanzerDiscFE_config.hpp"
#ifdef PANZER_HAVE_EPETRA
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"
#endif

#include "Panzer_ResponseMESupport_Default.hpp"
#include "Panzer_GlobalEvaluationData.hpp"
#include "Panzer_ThyraObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_LinearObjFactory.hpp"


namespace panzer {

/** This class provides a response storage for
  * simple functionals of the solution (i.e. scalar
  * values).
  */
template <typename EvalT>
class Response_Functional : public ResponseMESupport_Default<EvalT>
                          , public GlobalEvaluationData_BCAdjustment {
public:
   typedef typename EvalT::ScalarT ScalarT;

   Response_Functional(const std::string & responseName,MPI_Comm comm,
                       const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & linObjFact=Teuchos::null)
     : ResponseMESupport_Default<EvalT>(responseName,comm), value(0.0), linObjFactory_(linObjFact)
   {
     if(linObjFactory_!=Teuchos::null) {
       // requires thyra object factory
       thyraObjFactory_ = Teuchos::rcp_dynamic_cast<const panzer::ThyraObjFactory<double> >(linObjFactory_,true);
       setSolnVectorSpace(thyraObjFactory_->getThyraDomainSpace());

       // build a ghosted container, with a solution vector
       ghostedContainer_ = linObjFactory_->buildGhostedLinearObjContainer();

       // set ghosted container (work space for assembly)
       linObjFactory_->initializeGhostedContainer(panzer::LinearObjContainer::X,*ghostedContainer_);

       using Teuchos::rcp_dynamic_cast;
     }
   }

   //! provide direct access, this thing is pretty simple
   ScalarT value;

   //! This simply does global summation, then shoves the result into a vector
   virtual void scatterResponse();

   virtual void initializeResponse()
   { value = 0.0; if(ghostedContainer_!=Teuchos::null) ghostedContainer_->initialize(); }

   // from ResponseMESupport_Default

   //! What is the number of values you need locally
   virtual std::size_t localSizeRequired() const { return 1; }

   //! Is the vector distributed (or replicated)
   virtual bool vectorIsDistributed() const { return false; }

   //! Get ghosted responses (this will be filled by the evaluator)
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const
   { return Teuchos::rcp_dynamic_cast<const ThyraObjContainer<double> >(ghostedContainer_)->get_x_th(); }

   void adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows);

private:
   //! Set solution vector space
   void setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs);

   // hide these methods
   Response_Functional();
   Response_Functional(const Response_Functional &);

   Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;
   Teuchos::RCP<const panzer::ThyraObjFactory<double> > thyraObjFactory_;

   Teuchos::RCP<LinearObjContainer> uniqueContainer_;
   Teuchos::RCP<LinearObjContainer> ghostedContainer_;
};

}

#endif
