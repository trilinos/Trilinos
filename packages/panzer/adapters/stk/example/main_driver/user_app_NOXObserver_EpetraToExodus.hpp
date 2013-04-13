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
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"

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

      // get all element blocks and add them to the list
      std::vector<std::string> eBlocks;
      mesh->getElementBlockNames(eBlocks);

      panzer_stk::RespFactorySolnWriter_Builder builder;
      builder.mesh = mesh;
      m_response_library->addResponse("Main Field Output",eBlocks,builder);
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
      RCP<const Thyra::VectorBase<double> > th_x = n_th_x->getThyraRCPVector(); 

      // initialize the assembly container
      panzer::AssemblyEngineInArgs ae_inargs;
      ae_inargs.container_ = m_lof->buildLinearObjContainer();
      ae_inargs.ghostedContainer_ = m_lof->buildGhostedLinearObjContainer();
      ae_inargs.alpha = 0.0;
      ae_inargs.beta = 1.0;
      ae_inargs.evaluate_transient_terms = false;

      // initialize the ghosted container
      m_lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

      {
         // initialize the x vector
         const Teuchos::RCP<panzer::ThyraObjContainer<double> > thyraContainer
            = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ae_inargs.container_,true);
         thyraContainer->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(th_x));
      }

      m_response_library->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
      m_response_library->evaluate<panzer::Traits::Residual>(ae_inargs);
      
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
