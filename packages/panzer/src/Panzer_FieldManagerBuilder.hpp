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

  class GenericEvaluatorFactory {
  public:
    virtual bool registerEvaluators(PHX::FieldManager<panzer::Traits> & fm, const PhysicsBlock & pb) const = 0;
  };

  class EmptyEvaluatorFactory : public GenericEvaluatorFactory {
  public:
    bool registerEvaluators(PHX::FieldManager<panzer::Traits> & fm, const PhysicsBlock & pb) const 
    { return false; }
  };

  class FieldManagerBuilder {

  public:

    typedef std::map<unsigned,panzer::Workset> BCFaceWorksetMap;

    FieldManagerBuilder(bool disablePhysicsBlockScatter=false)
      : disablePhysicsBlockScatter_(disablePhysicsBlockScatter) {}

    void print(std::ostream& os) const;

    bool physicsBlockScatterDisabled() const
    { return disablePhysicsBlockScatter_; }

    void setWorksetContainer(const Teuchos::RCP<WorksetContainer> & wc)
    { worksetContainer_ = wc; }

    Teuchos::RCP<WorksetContainer> getWorksetContainer() const
    { return worksetContainer_; }

    const 
      std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >&
      getVolumeFieldManagers() const {return phx_volume_field_managers_;}

    //! Look up field manager by an element block ID
    Teuchos::RCP< PHX::FieldManager<panzer::Traits> >
    getVolumeFieldManager(const std::string & blockId) const 
    {
       const std::vector<std::string> & blockNames = getElementBlockNames();
       std::vector<std::string>::const_iterator itr = std::find(blockNames.begin(),blockNames.end(),blockId);
       TEUCHOS_ASSERT(itr!=blockNames.end());

       // get volume field manager associated with the block ID
       int index = itr - blockNames.begin();
       return getVolumeFieldManagers()[index];
    }

    const std::vector<std::string> &
      getElementBlockNames() const {return element_block_names_;}

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
				  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
				  const Teuchos::ParameterList& closure_models,
                                  const LinearObjFactory<panzer::Traits> & lo_factory,
				  const Teuchos::ParameterList& user_data,
                                  const GenericEvaluatorFactory & gEvalFact);

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

    //! Phalanx volume field managers for each element block.
    std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >
      phx_volume_field_managers_;

    /** \brief Matches volume field managers so you can determine
      *        the element block name for each field manager.
      */
    std::vector<std::string> element_block_names_;
    
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
  };

std::ostream& operator<<(std::ostream& os, const panzer::FieldManagerBuilder & rfd);

} // namespace panzer

#endif
