// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerAdaptersSTK_config.hpp"

#ifdef PANZER_HAVE_EPETRA_STACK

#include "Panzer_STK_Utilities.hpp"
#include "Panzer_GlobalIndexer.hpp"

#include "Kokkos_DynRankView.hpp"

#include <stk_mesh/base/FieldBase.hpp>

namespace panzer_stk {

static void gather_in_block(const std::string & blockId, const panzer::GlobalIndexer& dofMngr,
                            const Epetra_Vector & x,const std::vector<std::size_t> & localCellIds,
                            std::map<std::string,Kokkos::DynRankView<double,PHX::Device> > & fc);

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

      std::vector<stk::mesh::Entity> elements;
      mesh.getMyElements(blockId,elements);

      // loop over elements in this block
      for(std::size_t el=0;el<elements.size();el++) {
         std::size_t localId = mesh.elementLocalId(elements[el]);
         double * solnData = stk::mesh::field_data(*field,elements[el]);
         TEUCHOS_ASSERT(solnData!=0); // sanity check
         solnData[0] = data[localId];
      }
   }
}

void write_solution_data(const panzer::GlobalIndexer& dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_MultiVector & x,const std::string & prefix,const std::string & postfix)
{
   write_solution_data(dofMngr,mesh,*x(0),prefix,postfix);
}

void write_solution_data(const panzer::GlobalIndexer& dofMngr,panzer_stk::STK_Interface & mesh,const Epetra_Vector & x,const std::string & prefix,const std::string & postfix)
{
   typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

   // get local IDs
   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds;
   build_local_ids(mesh,localIds);

   // loop over all element blocks
   for(const auto & itr : localIds) {
      const auto blockId = itr.first;
      const auto & localCellIds = *(itr.second);

      std::map<std::string,FieldContainer> data;

      // get all solution data for this block
      gather_in_block(blockId,dofMngr,x,localCellIds,data);

      // write out to stk mesh
      for(const auto & dataItr : data)
         mesh.setSolutionFieldData(prefix+dataItr.first+postfix,blockId,localCellIds,dataItr.second);
   }
}

void gather_in_block(const std::string & blockId, const panzer::GlobalIndexer& dofMngr,
                     const Epetra_Vector & x,const std::vector<std::size_t> & localCellIds,
                     std::map<std::string,Kokkos::DynRankView<double,PHX::Device> > & fc)
{
   const std::vector<int> & fieldNums = dofMngr.getBlockFieldNumbers(blockId);

   for(std::size_t fieldIndex=0;fieldIndex<fieldNums.size();fieldIndex++) {
      int fieldNum = fieldNums[fieldIndex];
      std::string fieldStr = dofMngr.getFieldString(fieldNum);

      // grab the field
      const std::vector<int> & elmtOffset = dofMngr.getGIDFieldOffsets(blockId,fieldNum);
      fc[fieldStr] = Kokkos::DynRankView<double,PHX::Device>("fc",localCellIds.size(),elmtOffset.size());
      auto field = Kokkos::create_mirror_view(fc[fieldStr]);


      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::vector<panzer::GlobalOrdinal> GIDs;
         std::vector<int> LIDs;
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
            field(worksetCellIndex,basis) = x[lid];
         }
      }
      Kokkos::deep_copy(fc[fieldStr], field);
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
      std::vector<stk::mesh::Entity> blockElmts;
      mesh.getMyElements(blockId,blockElmts);

      std::vector<stk::mesh::Entity>::const_iterator itr;
      for(itr=blockElmts.begin();itr!=blockElmts.end();++itr)
         localBlockIds.push_back(mesh.elementLocalId(*itr));

      std::sort(localBlockIds.begin(),localBlockIds.end());
   }
}

}

#endif // PANZER_HAVE_EPETRA_STACK