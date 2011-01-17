#ifndef __Panzer_DOFManagerFactory_hpp__
#define __Panzer_DOFManagerFactory_hpp__

#include "Panzer_UniqueGlobalIndexerFactory.hpp"

namespace panzer {

template <typename LO,typename GO>
class DOFManagerFactory : public virtual UniqueGlobalIndexerFactory<LO,GO,LO,GO> {
public:
   virtual ~DOFManagerFactory() {}


   /** Use the physics block to construct a unique global indexer object.
     * 
     * \param[in] mpiComm MPI communicator to use in the construction
     * \param[in] physicsBlocks A vector of physics block objects that contain
     *                          unknown field information.
     * \param[in] connMngr Connection manager that contains the mesh topology
     *
     * \returns A fully constructed UniqueGlobalIndexer object
     */
   virtual Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO> > 
   buildUniqueGlobalIndexer(MPI_Comm mpiComm,
                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<ConnManager<LO,GO> > & connMngr) const;
};

}

#include "Panzer_DOFManagerFactoryT.hpp"

#endif
