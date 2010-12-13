#include "write_solution_data.hpp"

#include "Panzer_DOFManager.hpp"

void write_solution_data(const panzer::DOFManager<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x)
{
   typedef Intrepid::FieldContainer<double> FieldContainer;

   // get local IDs
   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds;
   build_local_ids(mesh,localIds);

   // loop over all element blocks
   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > >::const_iterator itr;
   for(itr=localIds.begin();itr!=localIds.end();++itr) {
      std::string blockId = itr->first;
      const std::vector<std::size_t> & localCellIds = *(itr->second);

      std::map<std::string,FieldContainer> data;

      // get all solution data for this block
      gather_in_block(blockId,dofMngr,x,localCellIds,data);

      // write out to stick mesh
      std::map<std::string,FieldContainer>::iterator dataItr;
      for(dataItr=data.begin();dataItr!=data.end();++dataItr) 
         mesh.setSolutionFieldData(dataItr->first,blockId,localCellIds,dataItr->second);
   }
}

void gather_in_block(const std::string & blockId, const panzer::DOFManager<int,int> & dofMngr,
                     const Epetra_Vector & x,const std::vector<std::size_t> & localCellIds,
                     std::map<std::string,Intrepid::FieldContainer<double> > & fc)
{
   panzer::DOFManager<int,int>::const_field_iterator fieldIter;

   for(fieldIter=dofMngr.beginFieldIter();fieldIter!=dofMngr.endFieldIter();++fieldIter) {
      int fieldNum = fieldIter->first;
      std::string fieldStr = fieldIter->second;

      fc[fieldStr].resize(localCellIds.size(),dofMngr.getFieldPattern(blockId,fieldStr)->numberIds());

      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::vector<int> GIDs, LIDs;
         std::size_t cellLocalId = localCellIds[worksetCellIndex];
      
         dofMngr.getElementGIDs(cellLocalId,GIDs);
      
         // caculate the local IDs for this element
         LIDs.resize(GIDs.size());
         for(std::size_t i=0;i<GIDs.size();i++)
            LIDs[i] = x.Map().LID(GIDs[i]);
   
         const std::vector<int> & elmtOffset = dofMngr.getGIDFieldOffsets(blockId,fieldNum);
   
         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            fc[fieldStr](worksetCellIndex,basis) = x[lid];
         }
      }
   }
}

void build_local_ids(const panzer_stk::STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds)
{
   // defines ordering of blocks
   std::vector<std::string> blockIds;
   mesh.getElementBlockNames(blockIds);

   std::vector<std::string>::const_iterator idItr;
   for(idItr=blockIds.begin();idItr!=blockIds.end();++idItr) {
      std::string blockId = *idItr;

      localIds[blockId] = Teuchos::rcp(new std::vector<std::size_t>);
      std::vector<std::size_t> & localBlockIds = *localIds[blockId];

      // grab elements on this block
      std::vector<stk::mesh::Entity*> blockElmts;
      mesh.getMyElements(blockId,blockElmts);

      std::vector<stk::mesh::Entity*>::const_iterator itr;
      for(itr=blockElmts.begin();itr!=blockElmts.end();++itr)
         localBlockIds.push_back(mesh.elementLocalId(*itr));

      std::sort(localBlockIds.begin(),localBlockIds.end());
   }
}
