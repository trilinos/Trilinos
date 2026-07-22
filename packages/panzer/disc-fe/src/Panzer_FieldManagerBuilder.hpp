// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_FIELD_MANAGER_BUILDER_HPP
#define PANZER_FIELD_MANAGER_BUILDER_HPP

#include <iostream>
#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_WorksetContainer.hpp"

// Forward Declarations
namespace panzer {
  struct Traits;
  class Workset;
  struct EquationSetFactory;
  struct BCStrategyFactory;
  class PhysicsBlock;
}

namespace PHX {
  template<typename T> class FieldManager;
}

namespace panzer {

  class GenericEvaluatorFactory {
  public:
    virtual bool registerEvaluators(PHX::FieldManager<panzer::Traits> & fm,const WorksetDescriptor & wd, const PhysicsBlock & pb) const = 0;
  };

  class EmptyEvaluatorFactory : public GenericEvaluatorFactory {
  public:
    bool registerEvaluators(PHX::FieldManager<panzer::Traits> & /* fm */, const WorksetDescriptor & /* wd */, const PhysicsBlock & /* pb */) const
    { return false; }
  };

  class FieldManagerBuilder {

  public:

    typedef std::map<unsigned,panzer::Workset> BCFaceWorksetMap;

    FieldManagerBuilder(bool disablePhysicsBlockScatter=false,bool disablePhysicsBlockGather=false);

    void print(std::ostream& os) const;

    bool physicsBlockScatterDisabled() const
    { return disablePhysicsBlockScatter_; }

    bool physicsBlockGatherDisabled() const
    { return disablePhysicsBlockGather_; }

    void setWorksetContainer(const Teuchos::RCP<WorksetContainer> & wc)
    { worksetContainer_ = wc; }

    Teuchos::RCP<WorksetContainer> getWorksetContainer() const
    { return worksetContainer_; }

    const std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >&
    getVolumeFieldManagers() const {return phx_volume_field_managers_;}

    //! Look up field manager by an element block ID
    Teuchos::RCP< PHX::FieldManager<panzer::Traits> >
    getVolumeFieldManager(const WorksetDescriptor & wd) const
    {
      const std::vector<WorksetDescriptor> & wkstDesc = getVolumeWorksetDescriptors();
      std::vector<WorksetDescriptor>::const_iterator itr = std::find(wkstDesc.begin(),wkstDesc.end(),wd);
      TEUCHOS_ASSERT(itr!=wkstDesc.end());

      // get volume field manager associated with the block ID
      int index = itr - wkstDesc.begin();
      return getVolumeFieldManagers()[index];
    }

    const std::vector<WorksetDescriptor> &
    getVolumeWorksetDescriptors() const { return volume_workset_desc_; }

    const std::map<panzer::BC,
		   std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
		   panzer::LessBC>&
    getBCFieldManagers() const {return bc_field_managers_;}

    // The intention of the next set of functions is to simplify and eventually
    // replace the setup routine above. Its not clear that these functions
    // belong in the field manager builder. Additionally this will add increased
    // flexibility to the field manager build in that the DOFManager will be
    // specified in a more flexable and generic way. Onward.... (ECC - 1/13/11)

    /** Setup the volume field managers. This uses the passed in <code>dofManager</code>
      * and sets it for permenant use.
      */
    void setupVolumeFieldManagers(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
				  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
				  const Teuchos::ParameterList& closure_models,
                                  const LinearObjFactory<panzer::Traits> & lo_factory,
				  const Teuchos::ParameterList& user_data);

    void setupVolumeFieldManagers(const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                  const std::vector<WorksetDescriptor> & wkstDesc,
				  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
				  const Teuchos::ParameterList& closure_models,
                                  const LinearObjFactory<panzer::Traits> & lo_factory,
				  const Teuchos::ParameterList& user_data,
                                  const GenericEvaluatorFactory & gEvalFact,
                                  bool closureModelByEBlock=false);

    /** Build the BC field managers.
      */
    void setupBCFieldManagers(const std::vector<panzer::BC> & bcs,
                              const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
	                      const panzer::EquationSetFactory & eqset_factory,
			      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                              const panzer::BCStrategyFactory& bc_factory,
			      const Teuchos::ParameterList& closure_models,
                              const LinearObjFactory<panzer::Traits> & lo_factory,
			      const Teuchos::ParameterList& user_data)
    { setupBCFieldManagers(bcs,physicsBlocks,Teuchos::ptrFromRef(eqset_factory),cm_factory,bc_factory,closure_models,lo_factory,user_data); }

    void setupBCFieldManagers(const std::vector<panzer::BC> & bcs,
                              const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
			      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                              const panzer::BCStrategyFactory& bc_factory,
			      const Teuchos::ParameterList& closure_models,
                              const LinearObjFactory<panzer::Traits> & lo_factory,
			      const Teuchos::ParameterList& user_data)
    { setupBCFieldManagers(bcs,physicsBlocks,Teuchos::null,cm_factory,bc_factory,closure_models,lo_factory,user_data); }

    void writeVolumeGraphvizDependencyFiles(std::string filename_prefix,
					    const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks) const;

    void writeBCGraphvizDependencyFiles(std::string filename_prefix) const;

    void writeVolumeTextDependencyFiles(std::string filename_prefix,
					const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks) const;

    void writeBCTextDependencyFiles(std::string filename_prefix) const;

    /// Delete all volume field managers, retaining the BC ones.
    void clearVolumeFieldManagers(bool clearVolumeWorksets = true);

    /// Set a vector of active evaluation types to allocate.
    void setActiveEvaluationTypes(const std::vector<bool>& aet);

  private:
    /** Build the BC field managers. This is the real deal, it correclty handles not having an equation set factory.
      */
    void setupBCFieldManagers(const std::vector<panzer::BC> & bcs,
                              const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
	                      const Teuchos::Ptr<const panzer::EquationSetFactory> & eqset_factory,
			      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                              const panzer::BCStrategyFactory& bc_factory,
			      const Teuchos::ParameterList& closure_models,
                              const LinearObjFactory<panzer::Traits> & lo_factory,
			      const Teuchos::ParameterList& user_data);

    void setKokkosExtendedDataTypeDimensions(const std::string & eblock,
                                             const panzer::GlobalIndexer & globalIndexer,
                                             const Teuchos::ParameterList& user_data,
                                             PHX::FieldManager<panzer::Traits> & fm) const;

    //! Phalanx volume field managers for each element block.
    std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >
      phx_volume_field_managers_;

    /** \brief Matches volume field managers so you can determine
      *        the appropriate set of worksets for each field manager.
      */
    std::vector<WorksetDescriptor> volume_workset_desc_;

    /*! \brief Field managers for the boundary conditions

        key is a panzer::BC object.  value is a map of
        field managers where the key is the local side index used by
        intrepid
    */
    std::map<panzer::BC,
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
      panzer::LessBC> bc_field_managers_;

    Teuchos::RCP<WorksetContainer> worksetContainer_;

    /** Set to false by default, enables/disables physics block scattering in
      * newly created field managers.
      */
    bool disablePhysicsBlockScatter_;

    /** Set to false by default, enables/disables physics block scattering in
      * newly created field managers.
      */
    bool disablePhysicsBlockGather_;

    /// Entries correspond to evaluation type mpl vector in traits. A value of true means the evaluation type is active.
    std::vector<bool> active_evaluation_types_;
  };

std::ostream& operator<<(std::ostream& os, const panzer::FieldManagerBuilder & rfd);

} // namespace panzer

#endif
