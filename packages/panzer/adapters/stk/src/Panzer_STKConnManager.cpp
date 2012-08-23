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

#include "Panzer_STKConnManager.hpp"

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include "Intrepid_FieldContainer.hpp"

#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Teuchos_FancyOStream.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {


// Object describing how to sort a vector of elements using
// local ID as the key
class LocalIdCompare {
public:
  LocalIdCompare(const RCP<const STK_Interface> & mesh) : mesh_(mesh) {}

  // Compares two stk mesh entities based on local ID
  bool operator() (stk::mesh::Entity * a,stk::mesh::Entity * b) 
  { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b);}

private:
  RCP<const STK_Interface> mesh_;
};

STKConnManager::STKConnManager(const Teuchos::RCP<STK_Interface> & stkMeshDB)
   : stkMeshDB_(stkMeshDB)
{
}

void STKConnManager::clearLocalElementMapping()
{
   elements_ = Teuchos::null; 

   elementBlocks_.clear();
   elmtLidToConn_.clear();
   connSize_.clear();
}

void STKConnManager::buildLocalElementMapping()
{ 
   clearLocalElementMapping(); // forget the past

   // build element block information
   //////////////////////////////////////////////
   elements_ = Teuchos::rcp(new std::vector<stk::mesh::Entity*>);

   // defines ordering of blocks
   std::vector<std::string> blockIds;
   stkMeshDB_->getElementBlockNames(blockIds);

   std::size_t blockIndex=0;
   std::vector<std::string>::const_iterator idItr;
   for(idItr=blockIds.begin();idItr!=blockIds.end();++idItr,++blockIndex) {
      std::string blockId = *idItr;

      // grab elements on this block
      std::vector<stk::mesh::Entity*> blockElmts;
      stkMeshDB_->getMyElements(blockId,blockElmts); 

      // concatenate them into element LID lookup table
      elements_->insert(elements_->end(),blockElmts.begin(),blockElmts.end());

      // build block to LID map
      elementBlocks_[blockId] = Teuchos::rcp(new std::vector<LocalOrdinal>);
      for(std::size_t i=0;i<blockElmts.size();i++)
         elementBlocks_[blockId]->push_back(stkMeshDB_->elementLocalId(blockElmts[i]));
   }

   // this expensive operation gurantees ordering of local IDs
   std::sort(elements_->begin(),elements_->end(),LocalIdCompare(stkMeshDB_));

   // allocate space for element LID to Connectivty map
   // connectivity size
   elmtLidToConn_.clear();
   elmtLidToConn_.resize(elements_->size(),0);

   connSize_.clear();
   connSize_.resize(elements_->size(),0);
}

void STKConnManager::buildOffsetsAndIdCounts(const panzer::FieldPattern & fp,
                                LocalOrdinal & nodeIdCnt, LocalOrdinal & edgeIdCnt,
                                LocalOrdinal & faceIdCnt, LocalOrdinal & cellIdCnt,
                                GlobalOrdinal & nodeOffset, GlobalOrdinal & edgeOffset,
                                GlobalOrdinal & faceOffset, GlobalOrdinal & cellOffset) const
{
   // get the global counts for all the nodes, faces, edges and cells
   GlobalOrdinal maxNodeId = stkMeshDB_->getMaxEntityId(stkMeshDB_->getNodeRank());
   GlobalOrdinal maxEdgeId = stkMeshDB_->getMaxEntityId(stkMeshDB_->getEdgeRank());
   GlobalOrdinal maxFaceId = stkMeshDB_->getMaxEntityId(stkMeshDB_->getFaceRank());

   // compute ID counts for each sub cell type
   int patternDim = fp.getDimension();
   switch(patternDim) {
   case 3:
     faceIdCnt = fp.getSubcellIndices(2,0).size();
   case 2:
     edgeIdCnt = fp.getSubcellIndices(1,0).size();
   case 1:
     nodeIdCnt = fp.getSubcellIndices(0,0).size();
     cellIdCnt = fp.getSubcellIndices(patternDim,0).size();
     break;
   case 0:
   default:
      TEUCHOS_ASSERT(false);
   };

   // compute offsets for each sub cell type
   nodeOffset = 0;
   edgeOffset = nodeOffset+(maxNodeId+1)*nodeIdCnt;
   faceOffset = edgeOffset+(maxEdgeId+1)*edgeIdCnt;
   cellOffset = faceOffset+(maxFaceId+1)*faceIdCnt;

   // sanity check
   TEUCHOS_ASSERT(nodeOffset <= edgeOffset 
               && edgeOffset <= faceOffset 
               && faceOffset <= cellOffset);
}

STKConnManager::LocalOrdinal STKConnManager::addSubcellConnectivities(
             stk::mesh::Entity * element,unsigned subcellRank,LocalOrdinal idCnt,GlobalOrdinal offset)
{
   if(idCnt<=0) 
      return 0 ;

   // loop over all relations of specified type
   LocalOrdinal numIds = 0;
   stk::mesh::PairIterRelation relations = element->relations(subcellRank);
   for(std::size_t sc=0;sc<relations.size();++sc) {
      stk::mesh::Entity * subcell = relations[sc].entity();

      // add connectivities: adjust for STK indexing craziness
      for(LocalOrdinal i=0;i<idCnt;i++) 
         connectivity_.push_back(offset+idCnt*(subcell->identifier()-1)+i);

      numIds += idCnt;
   }

   return numIds;
}

void STKConnManager::modifySubcellConnectivities(const panzer::FieldPattern & fp, stk::mesh::Entity * element,
                                                 unsigned subcellRank,unsigned subcellId,GlobalOrdinal newId,
                                                 GlobalOrdinal offset)
{
   LocalOrdinal elmtLID = stkMeshDB_->elementLocalId(element);
   GlobalOrdinal * conn = this->getConnectivity(elmtLID);
   const std::vector<int> & subCellIndices = fp.getSubcellIndices(subcellRank,subcellId);

   // add connectivities: adjust for STK indexing craziness
   for(std::size_t i=0;i<subCellIndices.size();i++) {
      conn[subCellIndices[i]] = offset+subCellIndices.size()*(newId-1)+i;
   }
}

void STKConnManager::buildConnectivity(const panzer::FieldPattern & fp)
{
   // get element info from STK_Interface
   // object and build a local element mapping.
   buildLocalElementMapping();

   // Build sub cell ID counts and offsets
   //    ID counts = How many IDs belong on each subcell (number of mesh DOF used)
   //    Offset = What is starting index for subcell ID type?
   //             Global numbering goes like [node ids, edge ids, face ids, cell ids]
   LocalOrdinal nodeIdCnt=0, edgeIdCnt=0, faceIdCnt=0, cellIdCnt=0;
   GlobalOrdinal nodeOffset=0, edgeOffset=0, faceOffset=0, cellOffset=0;
   buildOffsetsAndIdCounts(fp, nodeIdCnt,  edgeIdCnt,  faceIdCnt,  cellIdCnt,
                               nodeOffset, edgeOffset, faceOffset, cellOffset);

   // std::cout << "node: count = " << nodeIdCnt << ", offset = " << nodeOffset << std::endl;
   // std::cout << "edge: count = " << edgeIdCnt << ", offset = " << edgeOffset << std::endl;
   // std::cout << "face: count = " << faceIdCnt << ", offset = " << faceOffset << std::endl;
   // std::cout << "cell: count = " << cellIdCnt << ", offset = " << cellOffset << std::endl;

   // loop over elements and build global connectivity 
   for(std::size_t elmtLid=0;elmtLid!=elements_->size();++elmtLid) {
      GlobalOrdinal numIds = 0;
      stk::mesh::Entity * element = (*elements_)[elmtLid];

      // get index into connectivity array
      elmtLidToConn_[elmtLid] = connectivity_.size();

      // add connecviities for sub cells
      numIds += addSubcellConnectivities(element,stkMeshDB_->getNodeRank(),nodeIdCnt,nodeOffset);
      numIds += addSubcellConnectivities(element,stkMeshDB_->getEdgeRank(),edgeIdCnt,edgeOffset);
      numIds += addSubcellConnectivities(element,stkMeshDB_->getFaceRank(),faceIdCnt,faceOffset);

      // add connectivity for parent cells
      if(cellIdCnt>0) {
         // add connectivities: adjust for STK indexing craziness
         for(LocalOrdinal i=0;i<cellIdCnt;i++) 
            connectivity_.push_back(cellOffset+cellIdCnt*(element->identifier()-1));
      
         numIds += cellIdCnt;
      }

      connSize_[elmtLid] = numIds;
   }

   applyPeriodicBCs( fp, nodeOffset, edgeOffset, faceOffset, cellOffset);
}

std::string STKConnManager::getBlockId(STKConnManager::LocalOrdinal localElmtId) const
{
   // walk through the element blocks and figure out which this ID belongs to
   stk::mesh::Entity * element = (*elements_)[localElmtId];

   return stkMeshDB_->containingBlockId(element);
}

void STKConnManager::applyPeriodicBCs( const panzer::FieldPattern & fp, GlobalOrdinal nodeOffset, GlobalOrdinal edgeOffset, 
                                                                        GlobalOrdinal faceOffset, GlobalOrdinal cellOffset)
{
   using Teuchos::RCP;
   using Teuchos::rcp;


   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > matchedNodes
            = stkMeshDB_->getPeriodicNodePairing();

   // no matchedNodes means nothing to do!
   if(matchedNodes==Teuchos::null) return;
      
   for(std::size_t m=0;m<matchedNodes->size();m++) {
      stk::mesh::EntityId oldNodeId = (*matchedNodes)[m].first;
      std::size_t newNodeId = (*matchedNodes)[m].second;

      std::vector<stk::mesh::Entity*> elements;
      std::vector<int> localIds;

      // get relevent elements and node IDs
      stkMeshDB_->getOwnedElementsSharingNode(oldNodeId,elements,localIds);

      // modify global numbering already built for each element
      for(std::size_t e=0;e<elements.size();e++)
         modifySubcellConnectivities(fp,elements[e],0,localIds[e],newNodeId,nodeOffset);
   }
}

/** Get the coordinates for a specified element block and field pattern.
  */
void STKConnManager::getDofCoords(const std::string & blockId,
                                  const panzer::IntrepidFieldPattern & coordProvider,
                                  std::vector<std::size_t> & localCellIds,
                                  Intrepid::FieldContainer<double> & points) const
{
   int dim = coordProvider.getDimension();
   int numIds = coordProvider.numberIds();

   // grab element vertices
   Intrepid::FieldContainer<double> vertices;
   workset_utils::getIdsAndVertices(*stkMeshDB_,blockId,localCellIds,vertices);

   // setup output array
   points.resize(localCellIds.size(),numIds,dim);
   coordProvider.getInterpolatoryCoordinates(vertices,points);
}

}
