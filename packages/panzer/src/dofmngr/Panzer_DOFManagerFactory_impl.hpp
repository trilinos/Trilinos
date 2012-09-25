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

#ifndef PANZER_DOF_MANAGER_FACTORY_IMPL_HPP
#define PANZER_DOF_MANAGER_FACTORY_IMPL_HPP

#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManager2.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

namespace panzer {

template <typename LO,typename GO>
Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO> > 
DOFManagerFactory<LO,GO>::buildUniqueGlobalIndexer(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > & mpiComm,
                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<ConnManager<LO,GO> > & connMngr,
                            const std::string & fieldOrder) const
{
   if(useDOFManager2_)
      return buildUniqueGlobalIndexer<panzer::DOFManager2<LO,GO> >(mpiComm,physicsBlocks,connMngr,fieldOrder);
   else
      return buildUniqueGlobalIndexer<panzer::DOFManager<LO,GO> >(mpiComm,physicsBlocks,connMngr,fieldOrder);
}

template <typename LO,typename GO>
template <typename DOFManagerT>
Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO> > 
DOFManagerFactory<LO,GO>::buildUniqueGlobalIndexer(const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > & mpiComm,
                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<ConnManager<LO,GO> > & connMngr,
                            const std::string & fieldOrder) const
{
   Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
   pout->setShowProcRank(true);
   pout->setOutputToRootOnly(0);

   // build the DOF manager for the problem
   Teuchos::RCP<DOFManagerT> dofManager 
         = Teuchos::rcp(new DOFManagerT(connMngr,*mpiComm));

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

         Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepidBasis 
               = fieldItr->second->getIntrepidBasis();
         Teuchos::RCP<IntrepidFieldPattern> fp = Teuchos::rcp(new IntrepidFieldPattern(intrepidBasis));
         dofManager->addField(pb->elementBlockID(),fieldItr->first,fp);

         // *pout << "\"" << fieldItr->first << "\" Field Pattern = \n";
         // fp->print(*pout);
      }
   } 

   // set orientations required flag
   dofManager->setOrientationsRequired(orientationsRequired);

   if(fieldOrder!="") {
      std::vector<std::string> fieldOrderV;

      // this basiclly tokenzies "fieldOrder" string 
      // and dumps it into "fieldOrderV"
      std::stringstream ss;
      ss << fieldOrder;

      // until all tokens are eaten
      while(!ss.eof()) {
         std::string token;
         ss >> token;
 
         // reorder tokens
         if(token!="")
            fieldOrderV.push_back(token);
      }

      // do some stuff columinating in 
      dofManager->setFieldOrder(fieldOrderV);
   }

   dofManager->buildGlobalUnknowns();
   // dofManager->printFieldInformation(*pout);

   return dofManager;
}

}

#endif
