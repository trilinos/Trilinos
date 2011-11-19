#ifndef PANZER_FIELD_MANAGER_BUILDER_HPP
#define PANZER_FIELD_MANAGER_BUILDER_HPP

#include <iostream>
#include <vector>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_WorksetContainer.hpp"

// Forward Declarations
namespace panzer {
  struct Traits;
  struct Workset;
  template <typename LO, typename GO> class ConnManager;
  template <typename LO, typename GO> class DOFManager;
  class EquationSetFactory;
  class BCStrategyFactory;
  class PhysicsBlock;
}

namespace PHX {
  template<typename T> class FieldManager;
}  

namespace panzer {

  template <typename LO, typename GO>
  class FieldManagerBuilder {

  public:

    typedef std::map<unsigned,panzer::Workset> BCFaceWorksetMap;

    void print(std::ostream& os) const;

    const 
      std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >&
      getVolumeFieldManagers() const {return phx_volume_field_managers_;}

    const std::vector< Teuchos::RCP<std::vector<panzer::Workset> > >& 
      getWorksets() const {return worksets_;}

    const std::map<panzer::BC, 
		   std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
		   panzer::LessBC>& 
      getBCFieldManagers() const {return bc_field_managers_;}

    const std::map<panzer::BC,
		   Teuchos::RCP<std::map<unsigned,panzer::Workset> >,
		   panzer::LessBC>&
      getBCWorksets() const {return bc_worksets_;}

    // The intention of the next set of functions is to simplify and eventually
    // replace the setup routine above. Its not clear that these functions
    // belong in the field manager builder. Additionally this will add increased
    // flexibility to the field manager build in that the DOFManager will be
    // specified in a more flexable and generic way. Onward.... (ECC - 1/13/11)

    /** Setup the volume field managers. This uses the passed in <code>dofManager</code>
      * and sets it for permenant use.
      */
    void setupVolumeFieldManagers(WorksetContainer & wkstContainer, 
                                  const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
				  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
				  const Teuchos::ParameterList& closure_models,
                                  const LinearObjFactory<panzer::Traits> & lo_factory,
				  const Teuchos::ParameterList& user_data);

    /** Build the BC field managers.
      */
    void setupBCFieldManagers(const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>& bc_worksets,
                                                          // boundary condition -> map of (side_id,worksets)
                              const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
	                      const panzer::EquationSetFactory & eqset_factory,
			      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                              const panzer::BCStrategyFactory& bc_factory,
			      const Teuchos::ParameterList& closure_models,
                              const LinearObjFactory<panzer::Traits> & lo_factory,
			      const Teuchos::ParameterList& user_data);

    void writeVolumeGraphvizDependencyFiles(std::string filename_prefix,
					    const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks) const;

  private:

    //! Phalanx volume field managers for each element block.
    std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >
      phx_volume_field_managers_;
    
    //! Volume fill worksets for each element block.
    std::vector< Teuchos::RCP<std::vector<panzer::Workset> > > worksets_;

    /*! \brief Field managers for the boundary conditions

        key is a panzer::BC object.  value is a map of
        field managers where the key is the local side index used by
        intrepid
    */
    std::map<panzer::BC, 
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >,
      panzer::LessBC> bc_field_managers_;

    /*! \brief Worksets for the boundary conditions

        key is a panzer::BC object.  value is a map of
        worksets where the key is the local side index used by
        intrepid.  All elemenst of a boundary are in one workset for
        each local side of the sideset.
    */
    std::map<panzer::BC,
      Teuchos::RCP<std::map<unsigned,panzer::Workset> >,
      panzer::LessBC> bc_worksets_;

  };

template<typename LO, typename GO>
std::ostream& operator<<(std::ostream& os, const panzer::FieldManagerBuilder<LO,GO>& rfd);

} // namespace panzer

#include "Panzer_FieldManagerBuilderT.hpp"

#endif
