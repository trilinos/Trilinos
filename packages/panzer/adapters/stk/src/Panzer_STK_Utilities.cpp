#include "Panzer_STK_Utilities.hpp"
#include "Intrepid_FieldContainer.hpp"

#include "Panzer_DOFManager.hpp"

namespace panzer_stk {

static void gather_in_block(const std::string & blockId, const panzer::DOFManager<int,int> & dofMngr,
                            const Epetra_Vector & x,const std::vector<std::size_t> & localCellIds,
                            std::map<std::string,Intrepid::FieldContainer<double> > & fc);

void scatter_to_vector(const std::string & blockId, const panzer::DOFManager<int,int> & dofMngr,
                       const std::map<std::string,Intrepid::FieldContainer<double> > & fc,
                       const std::vector<std::size_t> & localCellIds,
                       Epetra_Vector & x);

static void build_local_ids(const panzer_stk::STK_Interface & mesh,
                            std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds);

void write_solution_data(const panzer::DOFManager<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_MultiVector & x)
{
   write_solution_data(dofMngr,mesh,*x(0));
}

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

      // write out to stk mesh
      std::map<std::string,FieldContainer>::iterator dataItr;
      for(dataItr=data.begin();dataItr!=data.end();++dataItr) 
         mesh.setSolutionFieldData(dataItr->first,blockId,localCellIds,dataItr->second);
   }
}

void read_solution_data(const panzer::DOFManager<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_MultiVector & x)
{
   read_solution_data(dofMngr,mesh,*x(0));
}

void read_solution_data(const panzer::DOFManager<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_Vector & x)
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
      const std::set<int> & fieldNums = dofMngr.getFields(blockId);

      // write out to stk mesh
      std::set<int>::const_iterator fieldItr;
      for(fieldItr=fieldNums.begin();fieldItr!=fieldNums.end();++fieldItr) {
         std::string fieldStr = dofMngr.getFieldString(*fieldItr);
         mesh.getSolutionFieldData(fieldStr,blockId,localCellIds,data[fieldStr]);
      }

      // get all solution data for this block
      scatter_to_vector(blockId,dofMngr,data,localCellIds,x);
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

      // grab the field
      Teuchos::RCP<const panzer::FieldPattern> fp = dofMngr.getFieldPattern(blockId,fieldNum);
      if(fp==Teuchos::null) // this field is not in this block
         continue;

      fc[fieldStr].resize(localCellIds.size(),fp->numberIds());

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

void scatter_to_vector(const std::string & blockId, const panzer::DOFManager<int,int> & dofMngr,
                       const std::map<std::string,Intrepid::FieldContainer<double> > & fc,
                       const std::vector<std::size_t> & localCellIds,
                       Epetra_Vector & x)
{
   
   std::map<std::string,Intrepid::FieldContainer<double> >::const_iterator fieldItr;
   for(fieldItr=fc.begin();fieldItr!=fc.end();++fieldItr) {
      std::string fieldStr = fieldItr->first;
      int fieldNum = dofMngr.getFieldNum(fieldStr);
      const Intrepid::FieldContainer<double> & data = fieldItr->second; 

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
         for(int basis=0;basis<data.dimension(1);basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            x[lid] = data(worksetCellIndex,basis);
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

}
