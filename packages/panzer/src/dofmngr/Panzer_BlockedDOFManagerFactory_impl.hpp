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

#ifndef PANZER_BlockedDOF_MANAGER_FACTORY_IMPL_HPP
#define PANZER_BlockedDOF_MANAGER_FACTORY_IMPL_HPP

#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_String_Utilities.hpp"

namespace panzer {

template <typename LO,typename GO>
bool BlockedDOFManagerFactory<LO,GO>::
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

template <typename LO,typename GO>
void BlockedDOFManagerFactory<LO,GO>::
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

template <typename LO,typename GO>
Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,std::pair<int,GO> > > 
BlockedDOFManagerFactory<LO,GO>::buildUniqueGlobalIndexer(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > & mpiComm,
                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<ConnManager<LO,GO> > & connMngr,
                            const std::string & fieldOrder) const
{
   TEUCHOS_ASSERT(requiresBlocking(fieldOrder));

   Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
   pout->setShowProcRank(true);
   pout->setOutputToRootOnly(0);

   // build the DOF manager for the problem
   Teuchos::RCP<panzer::BlockedDOFManager<LO,GO> > dofManager 
         = Teuchos::rcp(new panzer::BlockedDOFManager<LO,GO>(connMngr,*mpiComm));
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

         Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepidBasis 
               = fieldItr->second->getIntrepidBasis();
         Teuchos::RCP<IntrepidFieldPattern> fp = Teuchos::rcp(new IntrepidFieldPattern(intrepidBasis));
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
