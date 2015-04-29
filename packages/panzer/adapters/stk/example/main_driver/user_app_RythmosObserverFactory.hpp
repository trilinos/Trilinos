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

#ifndef USER_APP_RYTHMOS_OBSERVER_FACTORY_HPP
#define USER_APP_RYTHMOS_OBSERVER_FACTORY_HPP

#include "Panzer_config.hpp"

#include "Panzer_STK_RythmosObserverFactory.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_Traits.hpp"

#include "Rythmos_CompositeIntegrationObserver.hpp"

// Individual Observers
#include "user_app_RythmosObserver_WriteToExodus.hpp"
#include "user_app_RythmosObserver_CoordinateUpdate.hpp"

namespace user_app {

  class RythmosObserverFactory : public panzer_stk_classic::RythmosObserverFactory {

  public:
    RythmosObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary,
                           const Teuchos::RCP<panzer::WorksetContainer> wkstContainer,
                           bool useCoordinateUpdate)
       : stkIOResponseLibrary_(stkIOResponseLibrary)
       , wkstContainer_(wkstContainer)
       , useCoordinateUpdate_(useCoordinateUpdate) 
    {}

    bool useNOXObserver() const { return false; }
    
    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    buildRythmosObserver(const Teuchos::RCP<panzer_stk_classic::STK_Interface>& mesh,
			 const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & dof_manager,
			 const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      // all done?
      if(!useCoordinateUpdate_) {
        Teuchos::RCP<user_app::RythmosObserver_WriteToExodus> exodus_observer 
            = Teuchos::rcp(new user_app::RythmosObserver_WriteToExodus(mesh,dof_manager,lof,stkIOResponseLibrary_));

        return exodus_observer;
      }

      // note: Composite observer loops in the order added
      Teuchos::RCP<Rythmos::CompositeIntegrationObserver<double> > composite_observer =
        Rythmos::createCompositeIntegrationObserver<double>();

      {
        Teuchos::RCP<user_app::RythmosObserver_WriteToExodus> observer 
            = Teuchos::rcp(new user_app::RythmosObserver_WriteToExodus(mesh,dof_manager,lof,stkIOResponseLibrary_));
        composite_observer->addObserver(observer);
      }

      {
        Teuchos::RCP<user_app::RythmosObserver_CoordinateUpdate> observer 
            = Teuchos::rcp(new user_app::RythmosObserver_CoordinateUpdate(wkstContainer_));
        composite_observer->addObserver(observer);
      }

      return composite_observer;
    }

  private:
    //! Store STK IO response library...be careful, it will be modified externally
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary_;

    Teuchos::RCP<panzer::WorksetContainer> wkstContainer_;

    //! Use the coordinate update observer?
    bool useCoordinateUpdate_;
  };

}

#endif
