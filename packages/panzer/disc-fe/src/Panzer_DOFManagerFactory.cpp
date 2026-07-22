// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DOF_MANAGER_FACTORY_IMPL_HPP
#define PANZER_DOF_MANAGER_FACTORY_IMPL_HPP

#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

#include "Panzer_GlobalIndexer_Utilities.hpp"

namespace panzer {

Teuchos::RCP<panzer::GlobalIndexer> 
DOFManagerFactory::buildGlobalIndexer(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > & mpiComm,
                                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                                            const Teuchos::RCP<ConnManager> & connMngr,
                                            const std::string & fieldOrder) const
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManagerFactory::buildUnqueGlobalIndexer",BUGI);

   Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
   pout->setShowProcRank(true);
   pout->setOutputToRootOnly(0);

   // build the DOF manager for the problem
   Teuchos::RCP<panzer::DOFManager> dofManager 
     = Teuchos::rcp(new panzer::DOFManager(connMngr,*mpiComm));
 
   dofManager->enableTieBreak(useTieBreak_);
   dofManager->useNeighbors(useNeighbors_);

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
         // PureBasis::EElementSpace space = fieldItr->second->getElementSpace();
         // orientationsRequired |= ((space==PureBasis::HDIV) || (space==PureBasis::HCURL)); 
         orientationsRequired |= fieldItr->second->requiresOrientations();

         Teuchos::RCP< Intrepid2::Basis<PHX::Device::execution_space,double,double> > intrepidBasis 
               = fieldItr->second->getIntrepid2Basis();
         Teuchos::RCP<Intrepid2FieldPattern> fp = Teuchos::rcp(new Intrepid2FieldPattern(intrepidBasis));
         dofManager->addField(pb->elementBlockID(),fieldItr->first,fp);

         // *pout << "\"" << fieldItr->first << "\" Field Pattern = \n";
         // fp->print(*pout);
      }
   } 

   // set orientations required flag
   dofManager->setOrientationsRequired(orientationsRequired);

   if(fieldOrder!="") {
      std::vector<std::string> fieldOrderV;
      buildFieldOrder(fieldOrder,fieldOrderV);
      dofManager->setFieldOrder(fieldOrderV);
   }

   {
     PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManagerFactory::buildUnqueGlobalIndexer:buildGlobalUnknowns",BGU);
     dofManager->buildGlobalUnknowns();
   }

   // dofManager->printFieldInformation(*pout);

   // print out mesh topology information. Uncomment at your own risk, there will
   // be A LOT of information printed to the screen (scaling with the number of elements)
   // {
   //   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   //   out.setShowProcRank(true);
   //   out.setOutputToRootOnly(-1);
   //   printMeshTopology(out,*dofManager);
   // }

   return dofManager;
}

void 
DOFManagerFactory::
buildFieldOrder(const std::string & fieldOrderStr,std::vector<std::string> & fieldOrder)
{
  // this tokenizes "fieldOrderStr" string 
  // and dumps it into "fieldOrder"
  std::stringstream ss;
  ss << fieldOrderStr;

  // until all tokens are eaten
  while(!ss.eof()) {
     std::string token;
     ss >> token;

     // reorder tokens
     if(token!="")
        fieldOrder.push_back(token);
  }
}

}

#endif
