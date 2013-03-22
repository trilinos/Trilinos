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

#include "Panzer_STK_Utilities.hpp"
#include "Intrepid_FieldContainer.hpp"

#include "Panzer_DOFManagerFEI.hpp"

#include <stk_mesh/base/FieldData.hpp>

namespace panzer_stk {

static void gather_in_block(const std::string & blockId, const panzer::UniqueGlobalIndexer<int,int> & dofMngr,
                            const Epetra_Vector & x,const std::vector<std::size_t> & localCellIds,
                            std::map<std::string,Intrepid::FieldContainer<double> > & fc);

void scatter_to_vector(const std::string & blockId, const panzer::DOFManagerFEI<int,int> & dofMngr,
                       const std::map<std::string,Intrepid::FieldContainer<double> > & fc,
                       const std::vector<std::size_t> & localCellIds,
                       Epetra_Vector & x);

static void build_local_ids(const panzer_stk::STK_Interface & mesh,
                            std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds);

void write_cell_data(panzer_stk::STK_Interface & mesh,const std::vector<double> & data,const std::string & fieldName)
{
   std::vector<std::string> blocks;
   mesh.getElementBlockNames(blocks);

   // loop over element blocks
   for(std::size_t eb=0;eb<blocks.size();eb++) {
      const std::string & blockId = blocks[eb];
      panzer_stk::STK_Interface::SolutionFieldType * field = mesh.getCellField(fieldName,blockId);

      std::vector<stk::mesh::Entity*> elements;
      mesh.getMyElements(blockId,elements);

      // loop over elements in this block
      for(std::size_t el=0;el<elements.size();el++) {
         std::size_t localId = mesh.elementLocalId(elements[el]);
         double * solnData = stk::mesh::field_data(*field,*elements[el]);
         TEUCHOS_ASSERT(solnData!=0); // sanity check
         solnData[0] = data[localId];
      }
   }
}

void write_solution_data(const panzer::UniqueGlobalIndexer<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_MultiVector & x,const std::string & prefix,const std::string & postfix)
{
   write_solution_data(dofMngr,mesh,*x(0),prefix,postfix);
}

void write_solution_data(const panzer::UniqueGlobalIndexer<int,int> & dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x,const std::string & prefix,const std::string & postfix)
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
         mesh.setSolutionFieldData(prefix+dataItr->first+postfix,blockId,localCellIds,dataItr->second);
   }
}

void read_solution_data(const panzer::DOFManagerFEI<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_MultiVector & x)
{
   read_solution_data(dofMngr,mesh,*x(0));
}

void read_solution_data(const panzer::DOFManagerFEI<int,int> & dofMngr,const panzer_stk::STK_Interface & mesh,Epetra_Vector & x)
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

void gather_in_block(const std::string & blockId, const panzer::UniqueGlobalIndexer<int,int> & dofMngr,
                     const Epetra_Vector & x,const std::vector<std::size_t> & localCellIds,
                     std::map<std::string,Intrepid::FieldContainer<double> > & fc)
{
   const std::vector<int> & fieldNums = dofMngr.getBlockFieldNumbers(blockId);

   for(std::size_t fieldIndex=0;fieldIndex<fieldNums.size();fieldIndex++) {
      int fieldNum = fieldNums[fieldIndex];
      std::string fieldStr = dofMngr.getFieldString(fieldNum);

      // grab the field
      const std::vector<int> & elmtOffset = dofMngr.getGIDFieldOffsets(blockId,fieldNum);
      fc[fieldStr].resize(localCellIds.size(),elmtOffset.size());

      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::vector<int> GIDs, LIDs;
         std::size_t cellLocalId = localCellIds[worksetCellIndex];
      
         dofMngr.getElementGIDs(cellLocalId,GIDs);
      
         // caculate the local IDs for this element
         LIDs.resize(GIDs.size());
         for(std::size_t i=0;i<GIDs.size();i++)
            LIDs[i] = x.Map().LID(GIDs[i]);
   
         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            fc[fieldStr](worksetCellIndex,basis) = x[lid];
         }
      }
   }
}

void scatter_to_vector(const std::string & blockId, const panzer::DOFManagerFEI<int,int> & dofMngr,
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
