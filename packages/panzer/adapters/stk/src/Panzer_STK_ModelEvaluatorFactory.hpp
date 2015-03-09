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
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_ModelEvaluator.hpp"

#include "Panzer_NodeType.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"

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
  struct GlobalData;
  class UniqueGlobalIndexerBase;
  template <typename> class LinearObjFactory;

  template <typename,typename> class BlockedDOFManager;
#ifdef PANZER_HAVE_FEI
  template <typename,typename> class DOFManagerFEI;
#endif
  template <typename,typename> class DOFManager;
  template <typename> class ConnManagerBase;
}

namespace panzer_stk_classic {

  template <typename GO> class STKConnManager;
  class NOXObserverFactory;
  class RythmosObserverFactory;
  
  template<typename ScalarT>
  class ModelEvaluatorFactory : public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    /** @name Overridden from ParameterListAcceptor */
    //@{
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    //@}

    /** \brief Builds the model evaluators for a panzer assembly
        
        \param[in] comm (Required) Teuchos communicator.  Must be non-null.
        \param[in] global_data (Required) A fully constructed (all members allocated) global data object used to control parameter library and output support. Must be non-null.
        \param[in] eqset_factory (Required) Equation set factory to provide user defined equation sets.
        \param[in] bc_factory (Required) Boundary condition factory to provide user defined boundary conditions.
        \param[in] cm_factory (Required) Closure model factory to provide user defined closure models.
    */
    void buildObjects(const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
                      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
                      const panzer::BCStrategyFactory & bc_factory,
                      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                      bool meConstructionOn=true);

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getPhysicsModelEvaluator();
    
    /** @name Methods for building the solver */
    //@{

    void setNOXObserverFactory(const Teuchos::RCP<const panzer_stk_classic::NOXObserverFactory>& nox_observer_factory);

    void setRythmosObserverFactory(const Teuchos::RCP<const panzer_stk_classic::RythmosObserverFactory>& rythmos_observer_factory);

    template <typename BuilderT>
    int addResponse(const std::string & responseName,const std::vector<panzer::WorksetDescriptor> & wkstDesc,const BuilderT & builder);

    void buildResponses(const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory, 
                        const bool write_graphviz_file=false,
                        const std::string& graphviz_file_prefix="");

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getResponseOnlyModelEvaluator();

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > 
    buildResponseOnlyModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > & thyra_me,
                                    const Teuchos::RCP<panzer::GlobalData>& global_data,
                                    const Teuchos::RCP<Piro::RythmosSolver<ScalarT> > rythmosSolver = Teuchos::null,
                    const Teuchos::Ptr<const panzer_stk_classic::NOXObserverFactory> & in_nox_observer_factory=Teuchos::null,
                    const Teuchos::Ptr<const panzer_stk_classic::RythmosObserverFactory> & in_rythmos_observer_factory=Teuchos::null);

    //@}

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > getResponseLibrary();

    const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & getPhysicsBlocks() const;

    //! Get mesh object used to build model evaluator
    Teuchos::RCP<panzer_stk_classic::STK_Interface> getMesh() const 
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

    /** Clone the internal model evaluator, but use new physics blocks. Note that
      * the physics blocks must be in some sense compatible with the original set.
      */
    Teuchos::RCP<Thyra::ModelEvaluator<double> >
    cloneWithNewPhysicsBlocks(const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<ScalarT> > & solverFactory,
                              const Teuchos::RCP<Teuchos::ParameterList> & physics_block_plist,
                              const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
                              const panzer::BCStrategyFactory & bc_factory,
                              const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & user_cm_factory,
                              bool is_transient,bool is_explicit,
                              const Teuchos::Ptr<const Teuchos::ParameterList> & bc_list=Teuchos::null) const;

    /** \brief Setup the initial conditions in a model evaluator. Note that this
      *        is entirely self contained.
      */
    void setupInitialConditions(Thyra::ModelEvaluator<ScalarT> & model,
                                panzer::WorksetContainer & wkstContainer,
                                const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                const panzer::LinearObjFactory<panzer::Traits> & lof,
                                const Teuchos::ParameterList & initial_cond_pl,
                                const Teuchos::ParameterList & user_data_pl,
                                bool write_dot_files,const std::string & dot_file_prefix) const;

    /** \brief Write the initial conditions to exodus. Note that this
      *        is entirely self contained.
      */
    void writeInitialConditions(const Thyra::ModelEvaluator<ScalarT> & model,
                                const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                const Teuchos::RCP<panzer::WorksetContainer> & wc,
                                const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & ugi,
                                const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
                                const Teuchos::RCP<panzer_stk_classic::STK_Interface> & mesh,
                                const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                const Teuchos::ParameterList & closure_model_pl,
                                const Teuchos::ParameterList & user_data_pl,
                                int workset_size) const;

    /** This method is to assist with construction of the model evaluators.
      */ 
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > 
    buildPhysicsModelEvaluator(bool buildThyraME,
                        const Teuchos::RCP<panzer::FieldManagerBuilder> & fmb,
                        const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & rLibrary,
                                const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
                        const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > > & p_names,
                        const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<ScalarT> > & solverFactory,
                        const Teuchos::RCP<panzer::GlobalData> & global_data,
                        bool is_transient,double t_init) const;


    bool useDynamicCoordinates() const
    { return useDynamicCoordinates_; }

  protected:
 
    Teuchos::RCP<panzer::FieldManagerBuilder> 
    buildFieldManagerBuilder(const Teuchos::RCP<panzer::WorksetContainer> & wc,
                             const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                             const std::vector<panzer::BC> & bcs,
                             const panzer::EquationSetFactory & eqset_factory,
                             const panzer::BCStrategyFactory& bc_factory,
                             const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& volume_cm_factory,
                             const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& bc_cm_factory,
                             const Teuchos::ParameterList& closure_models,
                             const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
                             const Teuchos::ParameterList& user_data,
                             bool writeGraph,const std::string & graphPrefix) const;

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
    bool determineCoordinateField(const panzer::UniqueGlobalIndexerBase & globalIndexer,std::string & fieldName) const;

    /** Fill a STL map with the the block ids associated with the pattern for a specific field.
      *
      * \param[in] fieldName Field to fill associative container with
      * \param[out] fieldPatterns A map from element block IDs to field patterns associated with the fieldName
      *                           argument
      */
    void fillFieldPatternMap(const panzer::UniqueGlobalIndexerBase & globalIndexer, const std::string & fieldName, 
                             std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const;

#ifdef PANZER_HAVE_FEI
    /** Fill a STL map with the the block ids associated with the pattern for a specific field.
      *
      * \param[in] fieldName Field to fill associative container with
      * \param[out] fieldPatterns A map from element block IDs to field patterns associated with the fieldName
      *                           argument
      */
    template <typename GO>
    void fillFieldPatternMap(const panzer::DOFManagerFEI<int,GO> & globalIndexer, const std::string & fieldName, 
                             std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const;
#endif

    /** Fill a STL map with the the block ids associated with the pattern for a specific field.
      *
      * \param[in] fieldName Field to fill associative container with
      * \param[out] fieldPatterns A map from element block IDs to field patterns associated with the fieldName
      *                           argument
      */
    template <typename GO>
    void fillFieldPatternMap(const panzer::DOFManager<int,GO> & globalIndexer, const std::string & fieldName, 
                             std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const;

    /** \brief Gets the initial time from either the input parameter list or an exodus file
     *      
     * \param [in] transient_ic_params ParameterList that determines where to get the initial time value.
     * \param [in] mesh STK Mesh database used if the time value should come from the exodus file
    */
    double getInitialTime(Teuchos::ParameterList& transient_ic_params,
                          const panzer_stk_classic::STK_Interface& mesh) const;

    /**
      */
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > initializeSolnWriterResponseLibrary(
                                                                const Teuchos::RCP<panzer::WorksetContainer> & wc,
                                                                const Teuchos::RCP<panzer::UniqueGlobalIndexerBase> & ugi,
                                                                const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
                                                                const Teuchos::RCP<panzer_stk_classic::STK_Interface> & mesh) const;

    /**
      */
    void finalizeSolnWriterResponseLibrary(panzer::ResponseLibrary<panzer::Traits> & rl,
                                           const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                           const Teuchos::ParameterList & closure_models,
                                           int workset_size, Teuchos::ParameterList & user_data) const;

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
    buildLOWSFactory(bool blockedAssembly,
                     const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                     const Teuchos::RCP<panzer::ConnManagerBase<int> > & conn_manager,
                     const Teuchos::RCP<panzer_stk_classic::STK_Interface> & mesh,
                     const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm);

    /** Build LOWS factory.
      */
    template <typename GO>
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > buildLOWSFactory(bool blockedAssembly,
                                                                                const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer,
                                                                                const Teuchos::RCP<panzer_stk_classic::STKConnManager<GO> > & stkConn_manager,
                                                                                const Teuchos::RCP<panzer_stk_classic::STK_Interface> & mesh,
                                                                                const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm);

    template <typename GO>
    void writeTopology(const panzer::BlockedDOFManager<int,GO> & blkDofs) const;

#ifdef PANZER_HAVE_FEI
    template <typename GO>
    void writeTopology(const panzer::DOFManagerFEI<int,GO> & dofs,const std::string & block,std::ostream & os) const;
#endif

  private:

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_physics_me;
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_rome_me;

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > m_response_library;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > m_physics_blocks;

    Teuchos::RCP<panzer_stk_classic::STK_Interface> m_mesh;
    Teuchos::RCP<panzer::UniqueGlobalIndexerBase> m_global_indexer;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > m_lin_obj_factory;
    Teuchos::RCP<panzer::GlobalData> m_global_data;
    #ifdef HAVE_TEKO 
    Teuchos::RCP<Teko::RequestHandler> m_req_handler;
    #endif
    bool useDiscreteAdjoint;
    bool m_is_transient;
    Teuchos::RCP<const panzer::EquationSetFactory> m_eqset_factory;

    Teuchos::RCP<const panzer_stk_classic::NOXObserverFactory> m_nox_observer_factory;
    Teuchos::RCP<const panzer_stk_classic::RythmosObserverFactory> m_rythmos_observer_factory;
 
    bool useDynamicCoordinates_;
  };

template<typename ScalarT>
template <typename BuilderT>
int ModelEvaluatorFactory<ScalarT>::
addResponse(const std::string & responseName,const std::vector<panzer::WorksetDescriptor> & wkstDesc,const BuilderT & builder)
{
  typedef panzer::ModelEvaluator<double> PanzerME;

  Teuchos::RCP<Thyra::EpetraModelEvaluator> thyra_ep_me = Teuchos::rcp_dynamic_cast<Thyra::EpetraModelEvaluator>(m_physics_me);
  Teuchos::RCP<PanzerME> panzer_me = Teuchos::rcp_dynamic_cast<PanzerME>(m_physics_me);
 
  if(thyra_ep_me!=Teuchos::null && panzer_me==Teuchos::null) {
    // I don't need no const-ness!
    Teuchos::RCP<EpetraExt::ModelEvaluator> ep_me = Teuchos::rcp_const_cast<EpetraExt::ModelEvaluator>(thyra_ep_me->getEpetraModel());
    Teuchos::RCP<panzer::ModelEvaluator_Epetra> ep_panzer_me = Teuchos::rcp_dynamic_cast<panzer::ModelEvaluator_Epetra>(ep_me);

    return ep_panzer_me->addResponse(responseName,wkstDesc,builder);
  }
  else if(panzer_me!=Teuchos::null && thyra_ep_me==Teuchos::null) {
    return panzer_me->addResponse(responseName,wkstDesc,builder);
  }
     
  TEUCHOS_ASSERT(false);
  return -1;
}

}

#endif
