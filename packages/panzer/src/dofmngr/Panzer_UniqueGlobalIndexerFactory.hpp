#ifndef __Panzer_UniqueGlobalIndexerFactory_hpp__
#define __Panzer_UniqueGlobalIndexerFactory_hpp__

// stil includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// panzer includes
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"

namespace panzer {

template <typename LO,GO>
class UniqueGlobalIndexerFactory {

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
                            const Teuchos::RCP<ConnManager<LO,GO> > & connMngr) = 0;
};

}

#endif
