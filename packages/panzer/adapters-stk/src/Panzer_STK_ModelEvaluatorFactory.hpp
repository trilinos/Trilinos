// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_ModelEvaluator.hpp"

#include "Panzer_NodeType.hpp"

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#endif

#ifdef PANZER_HAVE_TEKO
#include "Teko_RequestHandler.hpp"
#endif

namespace Piro {
#ifdef PANZER_HAVE_TEMPUS
  template <typename ScalarT> class TempusSolverForwardOnly;
#endif
}

namespace Thyra {
  template<typename ScalarT> class ModelEvaluator;
  template<typename ScalarT> class LinearOpWithSolveFactoryBase;
}

namespace panzer {
  struct GlobalData;
  class GlobalIndexer;
  template <typename> class LinearObjFactory;

  class BlockedDOFManager;
  class DOFManager;
  class ConnManager;
}

namespace panzer_stk {

  class STKConnManager;
  class NOXObserverFactory;
#ifdef PANZER_HAVE_TEMPUS
  class TempusObserverFactory;
#endif
  class WorksetFactory;

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

    void setNOXObserverFactory(const Teuchos::RCP<const panzer_stk::NOXObserverFactory>& nox_observer_factory);

#ifdef PANZER_HAVE_TEMPUS
    void setTempusObserverFactory(const Teuchos::RCP<const panzer_stk::TempusObserverFactory>& tempus_observer_factory);
#endif

    template <typename BuilderT>
    int addResponse(const std::string & responseName,const std::vector<panzer::WorksetDescriptor> & wkstDesc,const BuilderT & builder);

    void buildResponses(const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                        const bool write_graphviz_file=false,
                        const std::string& graphviz_file_prefix="");

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getResponseOnlyModelEvaluator();

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> >
    buildResponseOnlyModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > & thyra_me,
                                    const Teuchos::RCP<panzer::GlobalData>& global_data,
#ifdef PANZER_HAVE_TEMPUS
                                    const Teuchos::RCP<Piro::TempusSolverForwardOnly<ScalarT> > tempusSolver = Teuchos::null,
#endif
                                    const Teuchos::Ptr<const panzer_stk::NOXObserverFactory> & in_nox_observer_factory=Teuchos::null
#ifdef PANZER_HAVE_TEMPUS
                                    , const Teuchos::Ptr<const panzer_stk::TempusObserverFactory> & in_tempus_observer_factory=Teuchos::null
#endif
                                    );

    //@}

    //! Set user defined workset factory
    void setUserWorksetFactory(Teuchos::RCP<panzer_stk::WorksetFactory>& user_wkst_factory);

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > getResponseLibrary();

    const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & getPhysicsBlocks() const;

    //! Get mesh object used to build model evaluator
    Teuchos::RCP<panzer_stk::STK_Interface> getMesh() const
    { return m_mesh; }

    //! Get global indexer used to build model evaluator
    Teuchos::RCP<panzer::GlobalIndexer> getGlobalIndexer() const
    { return m_global_indexer; }

    //! Get connection manager
    Teuchos::RCP<panzer::ConnManager> getConnManager() const
    { return m_conn_manager; }

    //! Is blocked assembly?
    bool isBlockedAssembly() const
    { return m_blockedAssembly; }

    //! Get linear object factory used to build model evaluator
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > getLinearObjFactory() const
    { return m_lin_obj_factory; }

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
                              const Teuchos::Ptr<const Teuchos::ParameterList> & bc_list=Teuchos::null,
                              const Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > & physics_me=Teuchos::null) const;

    /** \brief Setup the initial conditions in a model evaluator. Note that this
      *        is entirely self contained.
      */
    void setupInitialConditions(Thyra::ModelEvaluator<ScalarT> & model,
                                panzer::WorksetContainer & wkstContainer,
                                const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                const panzer::LinearObjFactory<panzer::Traits> & lof,
                                const Teuchos::ParameterList & closure_pl,
                                const Teuchos::ParameterList & initial_cond_pl,
                                const Teuchos::ParameterList & user_data_pl,
                                bool write_dot_files,const std::string & dot_file_prefix) const;

    /** \brief Write the initial conditions to exodus. Note that this
      *        is entirely self contained.
      */
    void writeInitialConditions(const Thyra::ModelEvaluator<ScalarT> & model,
                                const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                const Teuchos::RCP<panzer::WorksetContainer> & wc,
                                const Teuchos::RCP<const panzer::GlobalIndexer> & ugi,
                                const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
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
                        const std::vector<Teuchos::RCP<Teuchos::Array<double> > > & p_values,
                        const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<ScalarT> > & solverFactory,
                        const Teuchos::RCP<panzer::GlobalData> & global_data,
                        bool is_transient,double t_init) const;


    bool useDynamicCoordinates() const
    { return useDynamicCoordinates_; }

    /** \brief Gets the initial time from either the input parameter list or an exodus file
     *
     * \param [in] transient_ic_params ParameterList that determines where to get the initial time value.
     * \param [in] mesh STK Mesh database used if the time value should come from the exodus file
    */
    double getInitialTime(Teuchos::ParameterList& transient_ic_params,
                          const panzer_stk::STK_Interface& mesh) const;

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
    buildLOWSFactory(bool blockedAssembly,
                     const Teuchos::RCP<const panzer::GlobalIndexer> & globalIndexer,
                     const Teuchos::RCP<panzer::ConnManager> & conn_manager,
                     const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                     const Teuchos::RCP<const Teuchos::MpiComm<int> > & mpi_comm
                     #ifdef PANZER_HAVE_TEKO
                     , const Teuchos::RCP<Teko::RequestHandler> & req_handler=Teuchos::null
                     #endif
                     ) const;

    //! Get the workset container associated with the mesh database.
    Teuchos::RCP<panzer::WorksetContainer> getWorksetContainer() const
    { return m_wkstContainer; }

    //! Add the user fields specified by output_list to the mesh
    void addUserFieldsToMesh(panzer_stk::STK_Interface & mesh,const Teuchos::ParameterList & output_list) const;

    //! build STK mesh factory from a mesh parameter list
    Teuchos::RCP<STK_MeshFactory> buildSTKMeshFactory(const Teuchos::ParameterList & mesh_params) const;

    void finalizeMeshConstruction(const STK_MeshFactory & mesh_factory,
                                  const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                  const Teuchos::MpiComm<int> mpi_comm,
                                  STK_Interface & mesh) const;

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
                             bool writeGraph,const std::string & graphPrefix,
			     bool write_field_managers,const std::string & field_manager_prefix) const;

    /**
      */
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > initializeSolnWriterResponseLibrary(
                                                                const Teuchos::RCP<panzer::WorksetContainer> & wc,
                                                                const Teuchos::RCP<const panzer::GlobalIndexer> & ugi,
                                                                const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                                                const Teuchos::RCP<panzer_stk::STK_Interface> & mesh) const;

    /**
      */
    void finalizeSolnWriterResponseLibrary(panzer::ResponseLibrary<panzer::Traits> & rl,
                                           const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> & cm_factory,
                                           const Teuchos::ParameterList & closure_models,
                                           int workset_size, Teuchos::ParameterList & user_data) const;


  private:

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_physics_me;
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_rome_me;

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > m_response_library;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > m_physics_blocks;

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<panzer::GlobalIndexer> m_global_indexer;
    Teuchos::RCP<panzer::ConnManager> m_conn_manager;
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > m_lin_obj_factory;
    Teuchos::RCP<panzer::GlobalData> m_global_data;
    bool useDiscreteAdjoint;
    bool m_is_transient;
    bool m_blockedAssembly;
    Teuchos::RCP<const panzer::EquationSetFactory> m_eqset_factory;

    Teuchos::RCP<const panzer_stk::NOXObserverFactory> m_nox_observer_factory;
#ifdef PANZER_HAVE_TEMPUS
    Teuchos::RCP<const panzer_stk::TempusObserverFactory> m_tempus_observer_factory;
#endif
    Teuchos::RCP<panzer_stk::WorksetFactory> m_user_wkst_factory;
    Teuchos::RCP<panzer::WorksetContainer> m_wkstContainer;

    bool useDynamicCoordinates_;
  };

template<typename ScalarT>
template <typename BuilderT>
int ModelEvaluatorFactory<ScalarT>::
addResponse(const std::string & responseName,const std::vector<panzer::WorksetDescriptor> & wkstDesc,const BuilderT & builder)
{
  typedef panzer::ModelEvaluator<double> PanzerME;

#ifdef PANZER_HAVE_EPETRA_STACK
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
#else
  Teuchos::RCP<PanzerME> panzer_me = Teuchos::rcp_dynamic_cast<PanzerME>(m_physics_me);
  if(panzer_me!=Teuchos::null) {
    return panzer_me->addResponse(responseName,wkstDesc,builder);
  }
#endif

  TEUCHOS_ASSERT(false);
  return -1;
}

}

#endif
