#ifndef PANZER_STK_MODEL_EVALUATOR_FACTORY_HPP
#define PANZER_STK_MODEL_EVALUATOR_FACTORY_HPP

#include <string>
#include <map>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "Panzer_ConfigDefs.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

namespace Thyra {
  template<typename ScalarT> class ModelEvaluator;
}

namespace panzer_stk {
  
  template<typename ScalarT>
  class ModelEvaluatorFactory_Epetra : public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    /** @name Overridden from ParameterListAcceptor */
    //@{
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    //@}

    void buildObjects(const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getPhysicsModelEvaluator();
    
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getResponseOnlyModelEvaluator();

  protected:
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
    bool determineCoordinateField(const panzer::DOFManager<int,int> & dofManager,std::string & fieldName) const;

    /** Fill a STL map with the the block ids associated with the pattern for a specific field.
      *
      * \param[in] fieldName Field to fill associative container with
      * \param[out] fieldPatterns A map from element block IDs to field patterns associated with the fieldName
      *                           argument
      */
    void fillFieldPatternMap(const panzer::DOFManager<int,int> & dofManager, const std::string & fieldName, 
                             std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fieldPatterns) const;

  private:

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_physics_me;
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_rome_me;
  };

}

#include "Panzer_STK_ModelEvaluatorFactory_EpetraT.hpp"

#endif
