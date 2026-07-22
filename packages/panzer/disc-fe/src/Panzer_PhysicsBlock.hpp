// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PHYSICS_BLOCK_HPP
#define PANZER_PHYSICS_BLOCK_HPP

#include <string>
#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_EquationSet.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_EvaluatorsRegistrar.hpp"
#include "Panzer_WorksetNeeds.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace shards {
  class CellTopology;
}

namespace panzer {
  class PureBasis;
  class IntegrationRule;
  struct EquationSetFactory;
  struct GlobalData;
  class PhysicsBlock;
}

namespace panzer {


  /** \brief Nonmember function for building the physics blocks from a Teuchos::ParameterList for a given list of element blocks.  A unique physics block object is built for each element block even if multiple element blocks point to the same physics block.
      \relates panzer::PhysicsBlock
  */
  void buildPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
                          const std::map<std::string,Teuchos::RCP<const shards::CellTopology> >& block_ids_to_cell_topo,
                          const Teuchos::RCP<Teuchos::ParameterList>& physics_blocks_plist,
                          const int default_integration_order,
                          const std::size_t workset_size,
                          const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
                          const Teuchos::RCP<panzer::GlobalData>& global_data,
                          const bool build_transient_support,
                          std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                          const std::vector<std::string>& tangent_param_names = std::vector<std::string>());

  /** \brief Nonmember function for reading and constructing physics blocks from a Teuchos::ParameterList for a given list of element blocks.
    *        A unique physics block object is built for each element block even if multiple element blocks point to the same physics block.
    *        The <code>intialize</code> method must be called before the physics blocks are used.
    *
    * \relates panzer::PhysicsBlock
    */
  void readPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
                         const Teuchos::RCP<Teuchos::ParameterList>& physics_blocks_plist,
                         std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks);

  /** \brief Nonmember function for searching and returning a spcific physics block given an element block id. Throws an erro if the physics block is not found.
      \relates panzer::PhysicsBlock

      \param[in] element_block_id The name of the corresponding element block that this function will search for in the physics block vector.
      \param[in] physics_blocks Vector of physics blocks
      ]param[in] throw_on_failure Optional parameter that determines if the function hsould throw on failure.  Default is true.  If set to false and the funtion fails to find the physics block, then a null RCP is returned.
  */
  Teuchos::RCP<panzer::PhysicsBlock> findPhysicsBlock(const std::string element_block_id,
                                                      const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physics_blocks,
                                                      bool throw_on_failure = true);

  //! Object that contains information on the physics and discretization of a block of elements with the SAME topology.
  class PhysicsBlock : public EvaluatorsRegistrar {

  public:
    /** for testing purposes only */
    explicit PhysicsBlock()
      : m_build_transient_support(false), m_global_data(Teuchos::null),
        m_active_evaluation_types(Sacado::mpl::size<panzer::Traits::EvalTypes>::value,true)
    { std::cout << "WARNING: Default constructor for panzer::PhysicsBlock is for testing purposes only!" << std::endl; }

    /** This constructor call initialize.
      */
    PhysicsBlock(const Teuchos::RCP<Teuchos::ParameterList>& physics_block_plist,
                 const std::string & element_block_id,
                 const int default_integration_order,
                 const panzer::CellData & cell_data,
                 const Teuchos::RCP<const panzer::EquationSetFactory>& factory,
                 const Teuchos::RCP<panzer::GlobalData>& global_data,
                 const bool build_transient_support,
                 const std::vector<std::string>& tangent_param_names = std::vector<std::string>());

    /** This constructor allows a bare bones physics block to be initialized
      * that only knows meta data about the mesh (excepting the cell type).
      * This allows a read in from an input deck and then initializtion
      * further into the code.
      */
    PhysicsBlock(const Teuchos::RCP<Teuchos::ParameterList>& physics_block_plist,
                 const std::string & element_block_id);

    PhysicsBlock(const panzer::PhysicsBlock & pb,
                 const panzer::CellData & cell_data);

    /** This constructor builds a bare bones equation set. It will do gather
      * and scatter for a particular field and set of basis functions. It will
      * not have any equation sets associated with it.
      */
    PhysicsBlock(const std::string & element_block_id,
                 const std::string & physics_block_id,
                 const int integration_order,
                 const panzer::CellData & cell_data,
                 const Teuchos::RCP<panzer::GlobalData>& global_data,
                 const Teuchos::RCP<panzer::PureBasis> & fields);

    /** Initialize with cell data, equation set and global data. This is required before
      * the physics blocks can be used.
      */
    void initialize(const int default_integration_order,
                    const bool build_transient_support,
                    const panzer::CellData & cell_data,
                    const Teuchos::RCP<const panzer::EquationSetFactory>& factory,
                    const Teuchos::RCP<panzer::GlobalData>& global_data,
                    const std::vector<std::string>& tangent_param_names = std::vector<std::string>());

    /// Used to save memory by disabling unneeded evaluation types.
    void setActiveEvaluationTypes(const std::vector<bool>& aet);

    /// Used to reactivate all evaluation types if some were temporarily disabled with a call to setActiveEvalautionTypes().
    void activateAllEvaluationTypes();

    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                        const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterDOFProjectionsToIPEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                      const Teuchos::Ptr<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                                      const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                           const panzer::LinearObjFactory<panzer::Traits> & lof,
                                           const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                const Teuchos::ParameterList& models,
                                                const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                    const std::string& model_name,
                                                    const Teuchos::ParameterList& models,
                                                    const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                    const Teuchos::ParameterList& user_data) const;

    void buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                            const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                            const std::string& model_name,
                                            const Teuchos::ParameterList& models,
                                            const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterEquationSetEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                      const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterGatherAndOrientationEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                               const LinearObjFactory<panzer::Traits> & lof,
                                                               const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterDOFProjectionsToIPEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                             const Teuchos::Ptr<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                                             const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterScatterEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                  const LinearObjFactory<panzer::Traits> & lof,
                                                  const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterClosureModelEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                       const Teuchos::ParameterList& models,
                                                       const Teuchos::ParameterList& user_data) const;

    template<typename EvalT>
    void buildAndRegisterInitialConditionEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                           const std::string& model_name,
                                                           const Teuchos::ParameterList& models,
                                                           const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                           const Teuchos::ParameterList& user_data) const;

    const std::vector<std::string>& getDOFNames() const;
    const std::vector<StrPureBasisPair>& getProvidedDOFs() const;

    const std::vector<std::vector<std::string> > & getCoordinateDOFs() const;

    //! Returns list of tangent fields from DOFs and tangent param names
    const std::vector<StrPureBasisPair>& getTangentFields() const;

    /** Build a workset needs object for this physics block.
      */
    WorksetNeeds getWorksetNeeds() const;

    //! Returns the unique set of bases, key is the unique panzer::PureBasis::name() of the basis
    const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& getBases() const;

    //! Returns the unique set of point rules, key is the unique panzer::PointRule::name()
    const std::map<int,Teuchos::RCP<panzer::IntegrationRule> >& getIntegrationRules() const;

    const shards::CellTopology getBaseCellTopology() const;

    std::string physicsBlockID() const;
    std::string elementBlockID() const;

    const panzer::CellData & cellData() const;

    /** Build a copy of this physics block object but use the cell_data
      * passed in by the user. This is useful for creating physics block objects
      * that reside on the boundaries from physics blocks created for the volume.
      */
    Teuchos::RCP<PhysicsBlock> copyWithCellData(const panzer::CellData & cell_data) const;

    Teuchos::RCP<panzer::GlobalData> globalData() const;

    Teuchos::RCP<const FieldLibrary> getFieldLibrary() const
    { return m_field_lib.getConst(); }

    Teuchos::RCP<const FieldLibraryBase> getFieldLibraryBase() const
    { return m_field_lib.getConst(); }

    // return the Physics Block parameter list
    Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const
    { return m_input_parameters; }

  protected:
    void initialize(const Teuchos::RCP<Teuchos::ParameterList>& input_parameters,
                    const int& default_integration_order,
                    const std::string & element_block_id,
                    const panzer::CellData & cell_data,
                    const bool build_transient_support,
                    const std::vector<std::string>& tangent_param_names = std::vector<std::string>());

    std::string m_physics_id;
    std::string m_element_block_id;
    int m_default_integration_order;
    panzer::CellData m_cell_data;
    //! store the input parameter list for copy ctors
    Teuchos::RCP<Teuchos::ParameterList> m_input_parameters;
    bool m_build_transient_support;
    Teuchos::RCP<panzer::GlobalData> m_global_data;

    std::vector<std::string> m_dof_names;
    std::vector<StrPureBasisPair> m_provided_dofs;
    std::vector<StrPureBasisPair> m_tangent_fields;
    std::vector<std::vector<std::string> > m_coordinate_dofs; // coordinate DOFs (defines them)

    //! map of unique bases, key is the panzer::PureBasis::name() corresponding to its value
    std::map<std::string,Teuchos::RCP<panzer::PureBasis> > m_bases;
    //! map of unique integration rules, key is panzer::IntegrationRule::order() corresponding to its value
    std::map<int,Teuchos::RCP<panzer::IntegrationRule> > m_integration_rules;

    std::vector< Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > > m_equation_sets;
    Teuchos::RCP<FieldLibrary> m_field_lib;
    Teuchos::RCP<const panzer::EquationSetFactory> m_eqset_factory;

    /// Returns true for evaluation types that are active.
    std::vector<bool> m_active_evaluation_types;
  };

}

// ************************************************************
// template implementations
// ************************************************************

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterEquationSetEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                                        const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator
    eq_set = m_equation_sets.begin();
  int idx = 0;
  for (;eq_set != m_equation_sets.end(); ++eq_set,++idx) {
    if (m_active_evaluation_types[idx]) {
      EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

      const int di = eqstm.getAsObject<EvalT>()->setDetailsIndex(this->getDetailsIndex());
      eqstm.getAsObject<EvalT>()->buildAndRegisterEquationSetEvaluators(fm, *m_field_lib, user_data);
      eqstm.getAsObject<EvalT>()->setDetailsIndex(di);
    }
  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterGatherAndOrientationEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                                                 const LinearObjFactory<panzer::Traits> & lof,
                                                                                 const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator
    eq_set = m_equation_sets.begin();
  int idx = 0;
  for (;eq_set != m_equation_sets.end(); ++eq_set,++idx) {
    if (m_active_evaluation_types[idx]) {

      EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

      const int di = eqstm.getAsObject<EvalT>()->setDetailsIndex(this->getDetailsIndex());
      eqstm.getAsObject<EvalT>()->buildAndRegisterGatherAndOrientationEvaluators(fm,*m_field_lib,lof,user_data);
      eqstm.getAsObject<EvalT>()->setDetailsIndex(di);
    }
  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterDOFProjectionsToIPEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                                               const Teuchos::Ptr<const panzer::LinearObjFactory<panzer::Traits> > & lof,
                                                                               const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {
    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    // Loop over integration rules
    for (std::map<int,Teuchos::RCP<panzer::IntegrationRule> >::const_iterator ir_iter = m_integration_rules.begin();
         ir_iter != m_integration_rules.end(); ++ ir_iter) {

      Teuchos::RCP<panzer::IntegrationRule> ir = ir_iter->second;

      const int di = eqstm.getAsObject<EvalT>()->setDetailsIndex(this->getDetailsIndex());
      eqstm.getAsObject<EvalT>()->buildAndRegisterDOFProjectionsToIPEvaluators(fm,*m_field_lib->buildFieldLayoutLibrary(*ir),ir,lof,user_data);
      eqstm.getAsObject<EvalT>()->setDetailsIndex(di);
    }

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterScatterEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                                    const LinearObjFactory<panzer::Traits> & lof,
                                                                    const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    const int di = eqstm.getAsObject<EvalT>()->setDetailsIndex(this->getDetailsIndex());
    eqstm.getAsObject<EvalT>()->buildAndRegisterScatterEvaluators(fm,*m_field_lib,lof,user_data);
    eqstm.getAsObject<EvalT>()->setDetailsIndex(di);
  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterClosureModelEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                                         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                                         const Teuchos::ParameterList& models,
                                                                         const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    // Loop over integration rules
    for (std::map<int,Teuchos::RCP<panzer::IntegrationRule> >::const_iterator ir_iter = m_integration_rules.begin();
         ir_iter != m_integration_rules.end(); ++ ir_iter) {

      Teuchos::RCP<panzer::IntegrationRule> ir = ir_iter->second;

      const int di = eqstm.getAsObject<EvalT>()->setDetailsIndex(this->getDetailsIndex());
      eqstm.getAsObject<EvalT>()->buildAndRegisterClosureModelEvaluators(fm,*m_field_lib->buildFieldLayoutLibrary(*ir),ir,factory,models,user_data);
      eqstm.getAsObject<EvalT>()->setDetailsIndex(di);
    }

  }
}

template<typename EvalT>
void panzer::PhysicsBlock::buildAndRegisterInitialConditionEvaluatorsForType(PHX::FieldManager<panzer::Traits>& fm,
                                                                             const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                                                             const std::string& model_name,
                                                                             const Teuchos::ParameterList& models,
                                                                             const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                                             const Teuchos::ParameterList& user_data) const
{
  using std::vector;
  using Teuchos::RCP;
  using panzer::EquationSet_TemplateManager;

  // Loop over equation set template managers
  vector< RCP<EquationSet_TemplateManager<panzer::Traits> > >::const_iterator
    eq_set = m_equation_sets.begin();
  for (;eq_set != m_equation_sets.end(); ++eq_set) {
    std::vector<StrBasisPair> providedDOFs;

    EquationSet_TemplateManager<panzer::Traits> eqstm = *(*eq_set);

    const int di = eqstm.getAsObject<EvalT>()->setDetailsIndex(this->getDetailsIndex());
    eqstm.getAsObject<EvalT>()->buildAndRegisterInitialConditionEvaluators(fm, *m_field_lib, factory, model_name, models, lof, user_data);
    eqstm.getAsObject<EvalT>()->setDetailsIndex(di);
  }
}

#endif
