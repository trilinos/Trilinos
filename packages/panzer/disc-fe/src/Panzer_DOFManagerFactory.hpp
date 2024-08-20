// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_DOFManagerFactory_decl_hpp__
#define __Panzer_DOFManagerFactory_decl_hpp__

#include "PanzerDiscFE_config.hpp"

#include "Panzer_GlobalIndexerFactory.hpp"

namespace panzer {

class DOFManagerFactory : public virtual GlobalIndexerFactory {
public:
   DOFManagerFactory() : useDOFManagerFEI_(false), useTieBreak_(false), useNeighbors_(false) {}

   virtual ~DOFManagerFactory() {}


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
                            const std::string & fieldOrder="") const;

   void setUseDOFManagerFEI(bool flag)
   { 
     useDOFManagerFEI_ = flag; 
   }

   bool getUseDOFManagerFEI() const
   { 
     return false;
   }

   void setUseTieBreak(bool flag) 
   { useTieBreak_ = flag; }

   bool getUseTieBreak() const
   { return useTieBreak_; }

   void setUseNeighbors(bool flag)
   { useNeighbors_ = flag; }

   bool getUseNeighbors() const
   { return useNeighbors_; }

   static void buildFieldOrder(const std::string & fieldOrderStr,std::vector<std::string> & fieldOrder);

protected:

   bool useDOFManagerFEI_;
   bool useTieBreak_;
   bool useNeighbors_;
};

}

#endif
