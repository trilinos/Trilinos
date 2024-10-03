// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BlockedDOF_MANAGER_FACTORY_IMPL_HPP
#define PANZER_BlockedDOF_MANAGER_FACTORY_IMPL_HPP

#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_HashUtils.hpp"

namespace panzer {

bool BlockedDOFManagerFactory::
requiresBlocking(const std::string & fieldOrder)
{
   std::vector<std::string> tokens;

   // break it up on spaces
   StringTokenizer(tokens,fieldOrder," ",true);
 
   if(tokens.size()<2) // there has to be at least 2 tokens to block
      return false;

   // check the prefix - must signal "blocked"
   if(tokens[0]!="blocked:")
      return false;

   // loop over tokens
   bool acceptsHyphen = false;
   for(std::size_t i=1;i<tokens.size();i++) {

      // acceptsHyphen can't be false, and then a hyphen accepted
      TEUCHOS_TEST_FOR_EXCEPTION(tokens[i]=="-" && !acceptsHyphen,std::logic_error,
                                 "Blocked assembly: Error \"Field Order\" hyphen error at "
                                 "token " << i); 

      if(acceptsHyphen && tokens[i]=="-")
         acceptsHyphen = false; 
      else { // token must be a field
         acceptsHyphen = true; 
      }
   }

   return true;
}

void BlockedDOFManagerFactory::
buildBlocking(const std::string & fieldOrder,std::vector<std::vector<std::string> > & blocks)
{
   // now we don't have to check
   TEUCHOS_ASSERT(requiresBlocking(fieldOrder));

   std::vector<std::string> tokens;

   // break it up on spaces
   StringTokenizer(tokens,fieldOrder," ",true);

   Teuchos::RCP<std::vector<std::string> > current;
   for(std::size_t i=1;i<tokens.size();i++) {
 
      if(tokens[i]!="-" && tokens[i-1]!="-") {
         // if there is something to add, add it to the blocks
         if(current!=Teuchos::null) 
            blocks.push_back(*current);

         current = Teuchos::rcp(new std::vector<std::string>);
      }

      if(tokens[i]!="-")
         current->push_back(tokens[i]);
   }

   if(current!=Teuchos::null) 
      blocks.push_back(*current);
}

Teuchos::RCP<panzer::GlobalIndexer> 
BlockedDOFManagerFactory::buildGlobalIndexer(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > & mpiComm,
                                                   const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                                   const Teuchos::RCP<ConnManager> & connMngr,
                                                   const std::string & fieldOrder) const
{
   TEUCHOS_ASSERT(requiresBlocking(fieldOrder));

   Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
   pout->setShowProcRank(true);
   pout->setOutputToRootOnly(0);

   // build the DOF manager for the problem
   Teuchos::RCP<panzer::BlockedDOFManager> dofManager 
         = Teuchos::rcp(new panzer::BlockedDOFManager(connMngr,*mpiComm));
   dofManager->enableTieBreak(useTieBreak_);
   dofManager->setUseDOFManagerFEI(useDOFManagerFEI_);

   // by default assume orientations are not required
   bool orientationsRequired = false;

   std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
   for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
      Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
       
      const std::vector<StrPureBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<StrPureBasisPair,StrPureBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end()); 

      // add basis to DOF manager: block specific
      std::set<StrPureBasisPair,StrPureBasisComp>::const_iterator fieldItr; 
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {
         // determine if orientations are required
         orientationsRequired |= fieldItr->second->requiresOrientations();

         Teuchos::RCP< Intrepid2::Basis<PHX::Device::execution_space,double,double> > intrepidBasis 
               = fieldItr->second->getIntrepid2Basis();
         Teuchos::RCP<Intrepid2FieldPattern> fp = Teuchos::rcp(new Intrepid2FieldPattern(intrepidBasis));
         dofManager->addField(pb->elementBlockID(),fieldItr->first,fp);
      }
   }

   // set orientations required flag
   dofManager->setOrientationsRequired(orientationsRequired);

   std::vector<std::vector<std::string> > blocks;
   buildBlocking(fieldOrder,blocks);
   dofManager->setFieldOrder(blocks);

   dofManager->buildGlobalUnknowns();
   // dofManager->printFieldInformation(*pout);

   return dofManager;
}

}

#endif
