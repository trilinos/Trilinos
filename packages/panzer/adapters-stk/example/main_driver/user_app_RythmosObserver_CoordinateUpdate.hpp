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

#ifndef USER_APP_RYTHMOS_OBSERVER_COORDINATEUPDATE_HPP
#define USER_APP_RYTHMOS_OBSERVER_COORDINATEUPDATE_HPP

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeRange.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_WorksetContainer.hpp"

namespace user_app {

  class RythmosObserver_CoordinateUpdate : 
    public Rythmos::IntegrationObserverBase<double> {

  public:
    
    RythmosObserver_CoordinateUpdate(const Teuchos::RCP<panzer::WorksetContainer> & workset_container) :
      m_workset_container(workset_container)
    { }
    
    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    cloneIntegrationObserver() const
    { return Teuchos::rcp(new RythmosObserver_CoordinateUpdate(m_workset_container)); }

    void resetIntegrationObserver(const Rythmos::TimeRange<double>& /* integrationTimeDomain */)
    { }

    void observeCompletedTimeStep(const Rythmos::StepperBase<double>& /* stepper */,
				  const Rythmos::StepControlInfo<double>& /* stepCtrlInfo */,
				  const int /* timeStepIter */)
    { TEUCHOS_ASSERT(m_workset_container!=Teuchos::null); 
      m_workset_container->clear(); }
    
  protected:

    Teuchos::RCP<panzer::WorksetContainer> m_workset_container;
  };

}

#endif
