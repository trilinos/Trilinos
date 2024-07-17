// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_GlobalIndexerFactory_hpp__
#define __Panzer_GlobalIndexerFactory_hpp__

// stil includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// panzer includes
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"

namespace panzer {

/** A factory for building GlobalIndexer objects.
  * This is basically a single function that takes an MPI_Comm
  * object, a vector of PhysicsBlocks and a connection manager.
  * The Connection manager can have different local and global
  * index types.  The use case for this is the block assembly
  * functionality.
  */
class GlobalIndexerFactory {
public:

   virtual ~GlobalIndexerFactory() {}

   /** Use the physics block to construct a unique global indexer object.
     * 
     * \param[in] mpiComm MPI communicator to use in the construction
     * \param[in] physicsBlocks A vector of physics block objects that contain
     *                          unknown field information.
     * \param[in] connMngr Connection manager that contains the mesh topology
     * \param[in] fieldOrder Specifies the local ordering of the degrees of
     *            freedom. This is relevant when degrees of freedom are shared
     *            on the same geometric entity. The default is an alphabetical
     *            ordering.
     *
     * \returns A GlobalIndexer object. If buildGlobalUnknowns is true,
     *          the object is fully constructed. If it is false, the caller must
     *          finalize it.
     */
   virtual Teuchos::RCP<panzer::GlobalIndexer> 
   buildGlobalIndexer(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > & mpiComm,
                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<ConnManager> & connMngr,
                            const std::string & fieldOrder="") const = 0;
};

}

#endif
