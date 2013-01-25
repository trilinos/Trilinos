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

#ifndef PANZER_STK_MODEL_EVALUATOR_FACTORY_HPP
#define PANZER_STK_MODEL_EVALUATOR_FACTORY_HPP

#include <iostream>
#include <string>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "Panzer_ConfigDefs.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_ResponseAggregator_Factory.hpp"

#ifdef HAVE_TEKO 
#include "Teko_RequestHandler.hpp"
#endif

namespace Piro {
  template <typename ScalarT> class RythmosSolver;
}

namespace Thyra {
  template<typename ScalarT> class ModelEvaluator;
  template<typename ScalarT> class LinearOpWithSolveFactoryBase;
}

namespace panzer {
  class GlobalData;
  class UniqueGlobalIndexerBase;
  template <typename> class LinearObjFactory;

  template <typename,typename> class BlockedDOFManager;
  template <typename,typename> class DOFManagerFEI;
}

namespace panzer_stk {

  class STKConnManager;
  
  template<typename ScalarT>
  class ModelEvaluatorFactory_Epetra : public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    /** @name Overridden from ParameterListAcceptor */
    //@{
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    //@}

    /** \brief Builds the model evaluators for a panzer assembly
        
	\param comm [in] (Required) Teuchos communicator.  Must be non-null.
	\param global_data [in] (Required) A fully constructed (all members allocated) global data object used to control parameter library and output support. Must be non-null.
	\param eqset_factory [in] (Required) Equation set factory to provide user defined equation sets.
	\param bc_factory [in] (Required) Boundary condition factory to provide user defined boundary conditions.
	\param cm_factory [in] (Required) Closure model factory to provide user defined closure models.
	\param ra_factory [in] (Optional) Response aggregator factory to provide user defined response aggregator types.
    */
    void buildObjects(const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      const panzer::EquationSetFactory & eqset_factory,
                      const panzer::BCStrategyFactory & bc_factory,
		      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
		      const Teuchos::Ptr<const panzer::ResponseAggregatorFactory<panzer::Traits> > ra_factory = Teuchos::null);

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getPhysicsModelEvaluator();
    
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getResponseOnlyModelEvaluator();
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > buildResponseOnlyModelEvaluator(
                                  const Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > & thyra_me,
 		                  const Teuchos::RCP<panzer::GlobalData>& global_data,
                                  const Teuchos::RCP<Piro::RythmosSolver<ScalarT> > & rythmosSolver=Teuchos::null);

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > getResponseLibrary();

    const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & getPhysicsBlocks() const;

    //! Get mesh object used to build model evaluator
    Teuchos::RCP<panzer_stk::STK_Interface> getMesh() const 
    { return m_mesh; }

    //! Get global indexer used to build model evaluator
    Teuchos::RCP<panzer::UniqueGlobalIndexerBase> getGlobalIndexer() const 
    { return m_global_indexer; }

    //! Get linear object factory used to build model evaluator
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > getLinearObjFactory() const
    { return m_lin_obj_factory; }

    #ifdef HAVE_TEKO 
    void setTekoRequestHandler(Teuchos::RCP<Teko::RequestHandler> & reqHandler)
    { m_req_handler = reqHandler; }
    #endif     

    bool isTransient() const  
    { return m_is_transient; }

  protected:
    void addVolumeResponses(panzer::ResponseLibrary<panzer::Traits> & rLibrary,
                            const panzer_stk::STK_Interface & mesh,const Teuchos::ParameterList & pl) const;

    //! build STK mesh factory from a mesh parameter list
    Teuchos::RCP<STK_MeshFactory> buildSTKMeshFactory(const Teuchos::ParameterList & mesh_params) const;

    void finalizeMeshConstruction(const STK_MeshFactory & mesh_factory,
                                  const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                  const Teuchos::MpiComm<int> mpi_comm, 
                                  STK_Interface & mesh) const;

    /** This method determines a coordinate field from the DOF manager.
      * The algorithm is to loop over all the element blocks to find a field
      * that is defined for each. 
      *
      * \returns False if no unique field is found. Otherwise True is returned.
      */
    bool determineCoordinateField(const panzer::DOFManagerFEI<int,int> & globalIndexer,std::string & fieldName) const;

    /** Fill a STL map with the the block ids associated with the pattern for a specific field.
      *
      * \param[in] fieldName Field to fill associative container with
      * \param[out] fieldPatterns A map from element block IDs to field patterns associated with the fieldName
      *                           argument
      */
    void fillFieldPatternMap(const panzer::DOFManagerFEI<int,int> & globalIndexer, const std::string & fieldName, 
                             std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const;

    /** \brief Gets the initial time from either the input parameter list or an exodus file
     *      
     * \param [in] transient_ic_params ParameterList that determines where to get the initial time value.
     * \param [in] mesh STK Mesh database used if the time value should come from the exodus file
    */
    double getInitialTime(Teuchos::ParameterList& transient_ic_params,
			  const panzer_stk::STK_Interface& mesh) const;

    /**
      */
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > initializeSolnWriterResponseLibrary(
                                                                const Teuchos::RCP<panzer::WorksetContainer> & wc,
                                                                const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & ugi,
                                                                const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
                                                                const Teuchos::RCP<panzer_stk::STK_Interface> & mesh) const;

    /**
      */
    void finalizeSolnWriterResponseLibrary(panzer::ResponseLibrary<panzer::Traits> & rl,
                                           const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                           const Teuchos::ParameterList & closure_models,
                                           int workset_size, Teuchos::ParameterList & user_data) const;

    /** Build LOWS factory.
      */
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > buildLOWSFactory(bool blockedAssembly,
                                                                                const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                                                                                const Teuchos::RCP<panzer_stk::STKConnManager> & stkConn_manager,
                                                                                const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                                                                                const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm);

    void writeTopology(const panzer::BlockedDOFManager<int,int> & blkDofs) const;
    void writeTopology(const panzer::DOFManagerFEI<int,int> & dofs,const std::string & block,std::ostream & os) const;

  private:

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_physics_me;
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_rome_me;

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > m_response_library;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > m_physics_blocks;

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<panzer::UniqueGlobalIndexerBase> m_global_indexer;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > m_lin_obj_factory;
    Teuchos::RCP<panzer::GlobalData> m_global_data;
    #ifdef HAVE_TEKO 
    Teuchos::RCP<Teko::RequestHandler> m_req_handler;
    #endif
    bool useDiscreteAdjoint;
    bool m_is_transient;
  };

}

#include "Panzer_STK_ModelEvaluatorFactory_Epetra_impl.hpp"

#endif
