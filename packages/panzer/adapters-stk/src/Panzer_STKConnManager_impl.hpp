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

#include "Kokkos_DynRankView.hpp"

#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Teuchos_FancyOStream.hpp"

namespace panzer_stk {

using Teuchos::RCP;
using Teuchos::rcp;

// Object describing how to sort a vector of elements using
// local ID as the key
class LocalIdCompare {
public:
  LocalIdCompare(const RCP<const STK_Interface> & mesh) : mesh_(mesh) {}

  // Compares two stk mesh entities based on local ID
  bool operator() (stk::mesh::Entity a,stk::mesh::Entity b) 
  { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b);}

private:
  RCP<const STK_Interface> mesh_;
};

template <typename GO>
STKConnManager<GO>::STKConnManager(const Teuchos::RCP<const STK_Interface> & stkMeshDB)
   : stkMeshDB_(stkMeshDB), ownedElementCount_(0)
{
}

template <typename GO>
Teuchos::RCP<panzer::ConnManagerBase<int> > 
STKConnManager<GO>::
noConnectivityClone() const
{
  return Teuchos::rcp(new STKConnManager<GO>(stkMeshDB_));
}

template <typename GO>
void STKConnManager<GO>::clearLocalElementMapping()
{
   elements_ = Teuchos::null; 

   elementBlocks_.clear();
   elmtLidToConn_.clear();
   connSize_.clear();
   elmtToAssociatedElmts_.clear();
}

template <typename GO>
void STKConnManager<GO>::buildLocalElementMapping()
{ 
   clearLocalElementMapping(); // forget the past

   // build element block information
   //////////////////////////////////////////////
   elements_ = Teuchos::rcp(new std::vector<stk::mesh::Entity>);

   // defines ordering of blocks
   std::vector<std::string> blockIds;
   stkMeshDB_->getElementBlockNames(blockIds);

   std::size_t blockIndex=0;
   for(std::vector<std::string>::const_iterator idItr=blockIds.begin();
       idItr!=blockIds.end();++idItr,++blockIndex) {
      std::string blockId = *idItr;

      // grab elements on this block
      std::vector<stk::mesh::Entity> blockElmts;
      stkMeshDB_->getMyElements(blockId,blockElmts); 

      // concatenate them into element LID lookup table
      elements_->insert(elements_->end(),blockElmts.begin(),blockElmts.end());

      // build block to LID map
      elementBlocks_[blockId] = Teuchos::rcp(new std::vector<LocalOrdinal>);
      for(std::size_t i=0;i<blockElmts.size();i++)
         elementBlocks_[blockId]->push_back(stkMeshDB_->elementLocalId(blockElmts[i]));
   }

   ownedElementCount_ = elements_->size(); 

   blockIndex=0;
   for(std::vector<std::string>::const_iterator idItr=blockIds.begin();
       idItr!=blockIds.end();++idItr,++blockIndex) {
      std::string blockId = *idItr;

      // grab elements on this block
      std::vector<stk::mesh::Entity> blockElmts;
      stkMeshDB_->getNeighborElements(blockId,blockElmts); 

      // concatenate them into element LID lookup table
      elements_->insert(elements_->end(),blockElmts.begin(),blockElmts.end());

      // build block to LID map
      neighborElementBlocks_[blockId] = Teuchos::rcp(new std::vector<LocalOrdinal>);
      for(std::size_t i=0;i<blockElmts.size();i++)
         neighborElementBlocks_[blockId]->push_back(stkMeshDB_->elementLocalId(blockElmts[i]));
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

template <typename GO>
void STKConnManager<GO>::buildOffsetsAndIdCounts(const panzer::FieldPattern & fp,
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
     // Intentional fall-through.
   case 2:
     edgeIdCnt = fp.getSubcellIndices(1,0).size();
     // Intentional fall-through.
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

template <typename GO>
typename STKConnManager<GO>::LocalOrdinal STKConnManager<GO>::addSubcellConnectivities(
             stk::mesh::Entity element,unsigned subcellRank,LocalOrdinal idCnt,GlobalOrdinal offset)
{
   if(idCnt<=0) 
      return 0 ;

   // loop over all relations of specified type
   LocalOrdinal numIds = 0;
   stk::mesh::BulkData& bulkData = *stkMeshDB_->getBulkData();
   const stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(subcellRank);
   const size_t num_rels = bulkData.num_connectivity(element, rank);
   stk::mesh::Entity const* relations = bulkData.begin(element, rank);
   for(std::size_t sc=0; sc<num_rels; ++sc) {
     stk::mesh::Entity subcell = relations[sc];

     // add connectivities: adjust for STK indexing craziness
     for(LocalOrdinal i=0;i<idCnt;i++)
       connectivity_.push_back(offset+idCnt*(bulkData.identifier(subcell)-1)+i);

     numIds += idCnt;
   }
   return numIds;
}

template <typename GO>
void STKConnManager<GO>::modifySubcellConnectivities(const panzer::FieldPattern & fp, stk::mesh::Entity element,
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

template <typename GO>
void STKConnManager<GO>::buildConnectivity(const panzer::FieldPattern & fp)
{
#ifdef HAVE_EXTRA_TIMERS
  using Teuchos::TimeMonitor;
  RCP<Teuchos::TimeMonitor> tM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("panzer_stk::STKConnManager::buildConnectivity"))));
#endif

   stk::mesh::BulkData& bulkData = *stkMeshDB_->getBulkData();

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
      stk::mesh::Entity element = (*elements_)[elmtLid];

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
            connectivity_.push_back(cellOffset+cellIdCnt*(bulkData.identifier(element)-1));
      
         numIds += cellIdCnt;
      }

      connSize_[elmtLid] = numIds;
   }

   applyPeriodicBCs( fp, nodeOffset, edgeOffset, faceOffset, cellOffset);

   // This method does not modify connectivity_. But it should be called here
   // because the data it initializes should be available at the same time as
   // connectivity_.
   if (hasAssociatedNeighbors())
     applyInterfaceConditions();
}

template <typename GO>
std::string STKConnManager<GO>::getBlockId(STKConnManager::LocalOrdinal localElmtId) const
{
   // walk through the element blocks and figure out which this ID belongs to
   stk::mesh::Entity element = (*elements_)[localElmtId];

   return stkMeshDB_->containingBlockId(element);
}

template <typename GO>
void STKConnManager<GO>::applyPeriodicBCs( const panzer::FieldPattern & fp, GlobalOrdinal nodeOffset, GlobalOrdinal edgeOffset, 
                                                                        GlobalOrdinal faceOffset, GlobalOrdinal /* cellOffset */)
{
   using Teuchos::RCP;
   using Teuchos::rcp;

#ifdef HAVE_EXTRA_TIMERS
  using Teuchos::TimeMonitor;
  RCP<Teuchos::TimeMonitor> tM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("panzer_stk::STKConnManager::applyPeriodicBCs"))));
#endif

   std::pair<Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >, Teuchos::RCP<std::vector<unsigned int> > > matchedValues
            = stkMeshDB_->getPeriodicNodePairing();

   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > matchedNodes
            = matchedValues.first;
   Teuchos::RCP<std::vector<unsigned int> > matchTypes
            = matchedValues.second;

   // no matchedNodes means nothing to do!
   if(matchedNodes==Teuchos::null) return;

   for(std::size_t m=0;m<matchedNodes->size();m++) {
      stk::mesh::EntityId oldNodeId = (*matchedNodes)[m].first;
      std::size_t newNodeId = (*matchedNodes)[m].second;

      std::vector<stk::mesh::Entity> elements;
      std::vector<int> localIds;

      GlobalOrdinal offset0 = 0; // to make numbering consistent with that in PeriodicBC_Matcher
      GlobalOrdinal offset1 = 0; // offset for dof indexing
      if((*matchTypes)[m] == 0)
        offset1 = nodeOffset-offset0;
      else if((*matchTypes)[m] == 1){
        offset0 = stkMeshDB_->getMaxEntityId(stkMeshDB_->getNodeRank());
        offset1 = edgeOffset-offset0;
      } else if((*matchTypes)[m] == 2){
        offset0 = stkMeshDB_->getMaxEntityId(stkMeshDB_->getNodeRank())+stkMeshDB_->getMaxEntityId(stkMeshDB_->getEdgeRank());
        offset1 = faceOffset-offset0;
      } else
        TEUCHOS_ASSERT(false);

      // get relevent elements and node IDs
      stkMeshDB_->getOwnedElementsSharingNode(oldNodeId-offset0,elements,localIds,(*matchTypes)[m]);

      // modify global numbering already built for each element
      for(std::size_t e=0;e<elements.size();e++){
         modifySubcellConnectivities(fp,elements[e],(*matchTypes)[m],localIds[e],newNodeId,offset1);
      }

   }
}

/** Get the coordinates for a specified element block and field pattern.
  */
template <typename GO>
void STKConnManager<GO>::getDofCoords(const std::string & blockId,
                                  const panzer::Intrepid2FieldPattern & coordProvider,
                                  std::vector<std::size_t> & localCellIds,
                                  Kokkos::DynRankView<double,PHX::Device> & points) const
{
   int dim = coordProvider.getDimension();
   int numIds = coordProvider.numberIds();

   // grab element vertices
   Kokkos::DynRankView<double,PHX::Device> vertices;
   workset_utils::getIdsAndVertices(*stkMeshDB_,blockId,localCellIds,vertices);

   // setup output array
   points = Kokkos::DynRankView<double,PHX::Device>("points",localCellIds.size(),numIds,dim);
   coordProvider.getInterpolatoryCoordinates(vertices,points);
}

template <typename GO>
bool STKConnManager<GO>::hasAssociatedNeighbors() const
{
  return ! sidesetsToAssociate_.empty();
}

template <typename GO>
void STKConnManager<GO>::associateElementsInSideset(const std::string sideset_id)
{
  sidesetsToAssociate_.push_back(sideset_id);
  sidesetYieldedAssociations_.push_back(false);
}

inline std::size_t
getElementIdx(const std::vector<stk::mesh::Entity>& elements,
              stk::mesh::Entity const e)
{
  return static_cast<std::size_t>(
    std::distance(elements.begin(), std::find(elements.begin(), elements.end(), e)));
}

template <typename GO>
void STKConnManager<GO>::applyInterfaceConditions()
{  
  stk::mesh::BulkData& bulkData = *stkMeshDB_->getBulkData();
  elmtToAssociatedElmts_.resize(elements_->size());
  for (std::size_t i = 0; i < sidesetsToAssociate_.size(); ++i) {
    std::vector<stk::mesh::Entity> sides;
    stkMeshDB_->getAllSides(sidesetsToAssociate_[i], sides);
    sidesetYieldedAssociations_[i] = ! sides.empty();
    for (std::vector<stk::mesh::Entity>::const_iterator si = sides.begin();
         si != sides.end(); ++si) {
      stk::mesh::Entity side = *si;
      const size_t num_elements = bulkData.num_elements(side);
      stk::mesh::Entity const* elements = bulkData.begin_elements(side);
      if (num_elements != 2) {
        // If relations.size() != 2 for one side in the sideset, then it's true
        // for all, including the first.
        TEUCHOS_ASSERT(si == sides.begin());
        sidesetYieldedAssociations_[i] = false;
        break;
      }
      const std::size_t ea_id = getElementIdx(*elements_, elements[0]),
        eb_id = getElementIdx(*elements_, elements[1]);
      elmtToAssociatedElmts_[ea_id].push_back(eb_id);
      elmtToAssociatedElmts_[eb_id].push_back(ea_id);
    }
  }
}

template <typename GO>
std::vector<std::string> STKConnManager<GO>::
checkAssociateElementsInSidesets(const Teuchos::Comm<int>& comm) const
{
  std::vector<std::string> sidesets;
  for (std::size_t i = 0; i < sidesetYieldedAssociations_.size(); ++i) {
    int sya, my_sya = sidesetYieldedAssociations_[i] ? 1 : 0;
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &my_sya, &sya);
    if (sya == 0)
      sidesets.push_back(sidesetsToAssociate_[i]);
  }
  return sidesets;
}

template <typename GO>
const std::vector<typename STKConnManager<GO>::LocalOrdinal>&
STKConnManager<GO>::getAssociatedNeighbors(const LocalOrdinal& el) const
{
  return elmtToAssociatedElmts_[el];
}

}
