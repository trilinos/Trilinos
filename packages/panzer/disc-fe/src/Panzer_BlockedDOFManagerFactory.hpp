// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_BlockedDOFManagerFactory_decl_hpp__
#define __Panzer_BlockedDOFManagerFactory_decl_hpp__

#include "Panzer_GlobalIndexerFactory.hpp"

namespace panzer {

class BlockedDOFManagerFactory : public virtual GlobalIndexerFactory {
public:
   BlockedDOFManagerFactory() : useDOFManagerFEI_(false), useTieBreak_(false) {}
   virtual ~BlockedDOFManagerFactory() {}

   /** Does a fieldOrder string require blocking? 
     * A field order is basically stetup like this
     *    blocked: <field 0> <field 1> 
     * where two blocks will be created. To merge fields
     * between blocks use a hyphen, i.e.
     *    blocked: <field 0> <field 1> - <field 2> - <field 3>
     * This will create 2 blocks, the first contains only <field 0>
     * and the second combines <field 1>, <field 2> and <field 3>. Note
     * the spaces before and after the hyphen, these are important!
     */
   static bool requiresBlocking(const std::string & fieldorder);

   /** Does a fieldOrder string require blocking? 
     * A field order is basically stetup like this
     *    blocked: <field 0> <field 1> 
     * where two blocks will be created. To merge fields
     * between blocks use a hyphen, i.e.
     *    blocked: <field 0> <field 1> - <field 2> - <field 3>
     * This will create 2 blocks, the first contains only <field 0>
     * and the second combines <field 1>, <field 2> and <field 3>. Note
     * the spaces before and after the hyphen, these are important!
     */
   static void buildBlocking(const std::string & fieldorder,std::vector<std::vector<std::string> > & blocks);

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
     * \param[in] buildGlobalUnknowns Build the global unknowns before
     *            returning. The default value gives backwards-compatible
     *            behavior. Set this to false if the caller will initialize the
     *            DOF manager in additional ways before issuing the call to
     *            build the global unknowns itself.
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
   { useDOFManagerFEI_ = flag; }

   bool getUseDOFManagerFEI() const
   { return useDOFManagerFEI_; }

   void setUseTieBreak(bool flag) 
   { useTieBreak_ = flag; }

   bool getUseTieBreak()
   { return useTieBreak_; }

private:
   bool useDOFManagerFEI_;
   bool useTieBreak_;
};

}

#endif
