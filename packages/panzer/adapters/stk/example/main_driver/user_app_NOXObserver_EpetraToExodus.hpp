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

#ifndef USER_APP_NOX_OBSERVER_EPETRA_TO_EXODUS_HPP
#define USER_APP_NOX_OBSERVER_EPETRA_TO_EXODUS_HPP

#include "NOX_Abstract_PrePostOperator.H"
#include "Teuchos_RCP.hpp"

#include "Panzer_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ResponseAggregator_SolutionWriter.hpp"

#include "NOX_Epetra_Vector.H"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

namespace user_app {
  
  class NOXObserver_EpetraToExodus : public NOX::Abstract::PrePostOperator {
    
  public:
    
    NOXObserver_EpetraToExodus(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
			       const Teuchos::RCP<panzer::UniqueGlobalIndexerBase>& dof_manager,
			       const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof,
                               const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & response_library) :
      m_mesh(mesh),
      m_dof_manager(dof_manager),
      m_lof(lof),
      m_response_library(response_library)
    { 
      TEUCHOS_ASSERT(m_lof!=Teuchos::null);

      // register solution writer response aggregator with this library
      // this is an "Action" only response aggregator
      panzer::ResponseAggregator_Manager<panzer::Traits> & aggMngr = m_response_library->getAggregatorManager();
      panzer_stk::ResponseAggregator_SolutionWriter_Builder builder(mesh);
      builder.setLinearObjFactory(aggMngr.getLinearObjFactory());
      builder.setGlobalIndexer(aggMngr.getGlobalIndexer());
      aggMngr.defineAggregatorTypeFromBuilder("Solution Writer",builder);

      // require a particular "Solution Writer" response
      panzer::ResponseId rid("Main Field Output","Solution Writer");
      std::list<std::string> eTypes;
      eTypes.push_back("Residual");
      #ifdef HAVE_STOKHOS
         eTypes.push_back("SGResidual");
      #endif

      std::list<std::string> eBlocks;
      {
         // get all element blocks and add them to the list
         std::vector<std::string> eBlockNames;
         mesh->getElementBlockNames(eBlockNames);
         for(std::size_t i=0;i<eBlockNames.size();i++)
            eBlocks.push_back(eBlockNames[i]);
      }

      // reserve response guranteeing that we can evaluate it (assuming things are done correctly elsewhere)
      response_library->reserveLabeledBlockAggregatedVolumeResponse("Main Field Output",rid,eBlocks,eTypes);

      // used block LOF, or epetra LOF
      m_isEpetraLOF = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(m_lof)!=Teuchos::null;
    }
      
    void runPreIterate(const NOX::Solver::Generic& solver)
    {

    }
    
    void runPostIterate(const NOX::Solver::Generic& solver)
    {

    }
    
    void runPreSolve(const NOX::Solver::Generic& solver)
    {

    }
    
    void runPostSolve(const NOX::Solver::Generic& solver)
    {
      TEUCHOS_ASSERT(m_lof!=Teuchos::null);

      const NOX::Abstract::Vector& x = solver.getSolutionGroup().getX();
      const NOX::Thyra::Vector* n_th_x = dynamic_cast<const NOX::Thyra::Vector*>(&x);
      TEUCHOS_TEST_FOR_EXCEPTION(n_th_x == NULL, std::runtime_error, "Failed to dynamic_cast to NOX::Thyra::Vector!")
      const ::Thyra::VectorBase<double>& th_x = n_th_x->getThyraVector(); 

      // initialize the assembly container
      panzer::AssemblyEngineInArgs ae_inargs;
      ae_inargs.container_ = m_lof->buildLinearObjContainer();
      ae_inargs.ghostedContainer_ = m_lof->buildGhostedLinearObjContainer();
      ae_inargs.alpha = 0.0;
      ae_inargs.beta = 1.0;
      ae_inargs.evaluate_transient_terms = false;

      // initialize the ghosted container
      m_lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

      // Teuchos::MpiComm<int> comm(Teuchos::opaqueWrapper(dynamic_cast<const Epetra_MpiComm &>(ep_x->Comm()).Comm()));
      Teuchos::MpiComm<int> comm = m_lof->getComm();
      if(m_isEpetraLOF) {
         Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof
            = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(m_lof,true);
         Teuchos::RCP<const Epetra_Vector> ep_x = Thyra::get_Epetra_Vector(*(ep_lof->getMap()), Teuchos::rcp(&th_x, false));

         // initialize the x vector
         const Teuchos::RCP<panzer::EpetraLinearObjContainer> epGlobalContainer
            = Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs.container_,true);
         epGlobalContainer->set_x(Teuchos::rcp_const_cast<Epetra_Vector>(ep_x));
      }
      else {
         // initialize the x vector
         const Teuchos::RCP<panzer::BlockedEpetraLinearObjContainer> blkGlobalContainer
            = Teuchos::rcp_dynamic_cast<panzer::BlockedEpetraLinearObjContainer>(ae_inargs.container_,true);
         blkGlobalContainer->set_x(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(Teuchos::rcpFromRef(th_x)));
      }

      // do import
      m_lof->globalToGhostContainer(*ae_inargs.container_,*ae_inargs.ghostedContainer_,panzer::LinearObjContainer::X);

      // fill STK mesh objects
      m_response_library->evaluateVolumeFieldManagers<panzer::Traits::Residual>(ae_inargs,comm);
      
      // write to disk
      m_mesh->writeToExodus(0.0);
    }
    
  protected:

    void writeToScreen(std::ostream & os,const Thyra::VectorBase<double> & src)
    {
      const Thyra::SpmdVectorBase<double> & spmdSrc =
             Teuchos::dyn_cast<const Thyra::SpmdVectorBase<double> >(src);

      // get access to data
      Teuchos::ArrayRCP<const double> srcData;
      spmdSrc.getLocalData(Teuchos::ptrFromRef(srcData));
      os << "Local Size = " << srcData.size() << std::endl;
      for (int i=0; i < srcData.size(); ++i) {
         os << "   " << srcData[i] << std::endl;
      }
    }

    //! Copy a flat vector into a product vector
    void copyFlatThyraIntoBlockedThyra(const Thyra::VectorBase<double>& src, 
                                       const Teuchos::Ptr<Thyra::VectorBase<double> > & dest) const
    {
      using Teuchos::RCP;
      using Teuchos::ArrayView;
      using Teuchos::rcpFromPtr;
      using Teuchos::rcp_dynamic_cast;
    
      const RCP<Thyra::ProductVectorBase<double> > prodDest =
        Thyra::castOrCreateNonconstProductVectorBase(rcpFromPtr(dest));

      const Thyra::SpmdVectorBase<double> & spmdSrc =
             Teuchos::dyn_cast<const Thyra::SpmdVectorBase<double> >(src);
    
      // get access to flat data
      Teuchos::ArrayRCP<const double> srcData;
      spmdSrc.getLocalData(Teuchos::ptrFromRef(srcData));
    
      std::size_t offset = 0;
      const int numBlocks = prodDest->productSpace()->numBlocks();
      for (int b = 0; b < numBlocks; ++b) {
        const RCP<Thyra::VectorBase<double> > destBlk = prodDest->getNonconstVectorBlock(b);

        // get access to blocked data
        const RCP<Thyra::SpmdVectorBase<double> > spmdBlk =
               rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(destBlk, true);
        Teuchos::ArrayRCP<double> destData;
        spmdBlk->getNonconstLocalData(Teuchos::ptrFromRef(destData));
    
        // perform copy
        for (int i=0; i < destData.size(); ++i) {
          destData[i] = srcData[i+offset];
        }
        offset += destData.size();
      }
    
    }


    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<panzer::UniqueGlobalIndexerBase> m_dof_manager;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > m_lof;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > m_response_library;

    bool m_isEpetraLOF;
  };
}

#endif
