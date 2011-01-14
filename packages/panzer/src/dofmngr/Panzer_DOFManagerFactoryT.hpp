
namespace panzer {

template <typename LO,GO>
Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO> > 
DOFManagerFactory<LO,GO>::buildUniqueGlobalIndexer(MPI_Comm mpiComm,
                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks,
                            const Teuchos::RCP<ConnManager<LO,GO> > & connMngr)
{
   Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
   pout->setShowProcRank(true);
   pout->setOutputToRootOnly(0);

   // build the DOF manager for the problem
   Teuchos::RCP<panzer::DOFManager<LO,GO> > dofManager 
         = Teuchos::rcp(new panzer::DOFManager<LO,GO>(connMngr,mpiComm));

   std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator physIter;
   for(physIter=physicsBlocks.begin();physIter!=physicsBlocks.end();++physIter) {
      Teuchos::RCP<const panzer::PhysicsBlock> pb = *physIter;
       
      const std::vector<StrBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<StrBasisPair,StrBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end()); 

      // add basis to DOF manager: block specific
      std::set<StrBasisPair>::const_iterator fieldItr; 
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {
         Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepidBasis 
               = fieldItr->second->getIntrepidBasis();
         Teuchos::RCP<IntrepidFieldPattern> fp = Teuchos::rcp(new IntrepidFieldPattern(intrepidBasis));
         dofManager->addField(pb->elementBlockID(),fieldItr->first,fp);

         *pout << "\"" << fieldItr->first << "\" Field Pattern = \n";
         fp->print(*pout);
      }
   } 

   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(*pout);

   return dofManager;
}

}
