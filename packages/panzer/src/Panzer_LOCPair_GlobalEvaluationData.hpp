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

#ifndef __Panzer_LOCPair_GlobalEvaluationData_hpp__
#define __Panzer_LOCPair_GlobalEvaluationData_hpp__

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_LinearObjFactory.hpp"

namespace panzer {

/** Class that overides the communication primitives
  * to do nothing. This is used by the <code>LinearObjContainer</code>.
  */
class LOCPair_GlobalEvaluationData : public GlobalEvaluationData {
public:
   LOCPair_GlobalEvaluationData(
                 Teuchos::RCP<const LinearObjFactory<panzer::Traits> > lof,
                 int initParam) : lof_(lof), initParam_(initParam) 
   {
      globalLOC_ = lof_->buildLinearObjContainer();
      ghostedLOC_ = lof_->buildGhostedLinearObjContainer();

      lof_->initializeContainer(initParam,*globalLOC_);
      lof_->initializeGhostedContainer(initParam,*ghostedLOC_);
   }
  
   virtual void ghostToGlobal(int mem) { lof_->ghostToGlobalContainer(*ghostedLOC_,*globalLOC_,mem); }
   virtual void globalToGhost(int mem) { lof_->globalToGhostContainer(*globalLOC_,*ghostedLOC_,mem); }

   virtual void initializeData() { ghostedLOC_->initialize(); } 

   Teuchos::RCP<LinearObjContainer> getGhostedLOC() const { return ghostedLOC_; }
   Teuchos::RCP<LinearObjContainer> getGlobalLOC() const { return globalLOC_; }

private:
   Teuchos::RCP<LinearObjContainer> ghostedLOC_, globalLOC_;

   Teuchos::RCP<const LinearObjFactory<panzer::Traits> > lof_; 
   int initParam_;
};

}

#endif
