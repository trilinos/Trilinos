// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerDofMgr_config.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

namespace panzer {


std::vector<Teuchos::RCP<const GlobalIndexer>>
nc2c_vector(const std::vector<Teuchos::RCP<GlobalIndexer > > & ugis)
{
  std::vector<Teuchos::RCP<const GlobalIndexer>> vec;

  for(std::size_t blk=0;blk<ugis.size();blk++) 
    vec.push_back(ugis[blk]);

  return vec;
}

int getFieldBlock(const std::string & fieldName,
                  const std::vector<Teuchos::RCP<const GlobalIndexer>> & ugis)
{
  int fieldNum = -1;
  for(std::size_t blk=0;blk<ugis.size();blk++) {
    fieldNum = ugis[blk]->getFieldNum(fieldName);
    if(fieldNum>=0)
      return blk;
  }

  return fieldNum;
}

int getFieldBlock(const std::string & fieldName,
                  const std::vector<Teuchos::RCP<GlobalIndexer>> & ugis)
{
  int fieldNum = -1;
  for(std::size_t blk=0;blk<ugis.size();blk++) {
    fieldNum = ugis[blk]->getFieldNum(fieldName);
    if(fieldNum>=0)
      return fieldNum;
  }

  return fieldNum;
}

void computeBlockOffsets(const std::string & blockId,
                         const std::vector<Teuchos::RCP<GlobalIndexer>> & ugis,
                         std::vector<int> & blockOffsets)
{
  blockOffsets.resize(ugis.size()+1); // number of fields, plus a sentinnel

  int offset = 0;
  for(std::size_t blk=0;blk<ugis.size();blk++) {
    blockOffsets[blk] = offset;
    offset += ugis[blk]->getElementBlockGIDCount(blockId);
  }
  blockOffsets[ugis.size()] = offset;
}

void computeBlockOffsets(const std::string & blockId,
                         const std::vector<Teuchos::RCP<const GlobalIndexer>> & ugis,
                         std::vector<int> & blockOffsets)
{
  blockOffsets.resize(ugis.size()+1); // number of fields, plus a sentinnel

  int offset = 0;
  for(std::size_t blk=0;blk<ugis.size();blk++) {
    blockOffsets[blk] = offset;
    offset += ugis[blk]->getElementBlockGIDCount(blockId);
  }
  blockOffsets[ugis.size()] = offset;
}

std::string 
printUGILoadBalancingInformation(const GlobalIndexer & ugi)
{
  std::size_t myOwnedCount = static_cast<std::size_t>(ugi.getNumOwned()); 
  std::size_t sum=0,min=0,max=0;

  // get min,max and sum
  Teuchos::reduceAll(*ugi.getComm(),Teuchos::REDUCE_SUM,1,&myOwnedCount,&sum);
  Teuchos::reduceAll(*ugi.getComm(),Teuchos::REDUCE_MIN,1,&myOwnedCount,&min);
  Teuchos::reduceAll(*ugi.getComm(),Teuchos::REDUCE_MAX,1,&myOwnedCount,&max);

  // compute mean and variance
  double dev2 = (double(myOwnedCount)-double(sum)/double(ugi.getComm()->getSize()));
  dev2 *= dev2;

  double variance = 0.0;
  Teuchos::reduceAll(*ugi.getComm(),Teuchos::REDUCE_SUM,1,&dev2,&variance);
 
  double mean = sum / double(ugi.getComm()->getSize());
  variance = variance / double(ugi.getComm()->getSize());

  // now print to a string stream
  std::stringstream ss;
  ss << "Max, Min, Mean, StdDev = " << max << ", " << min << ", " << mean << ", " << std::sqrt(variance);

  return ss.str();
}

void
printMeshTopology(std::ostream & os,const panzer::GlobalIndexer & ugi)
{
  std::vector<std::string> block_ids;

  ugi.getElementBlockIds(block_ids);
  for(std::size_t b=0;b<block_ids.size();b++) {
    // extract the elemnts in each element block
    const std::vector<panzer::LocalOrdinal> & elements = ugi.getElementBlock(block_ids[b]);

    os << "Element Block: \"" << block_ids[b] << "\"" << std::endl;
 
    // loop over element in this element block, write out to 
    for(std::size_t e=0;e<elements.size();e++) {
      // extract LIDs, this is terribly inefficient for certain
      // devices but only used for debugging
      PHX::View<const int*> lids = ugi.getElementLIDs(elements[e]);
      auto lids_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),lids);

      // extract GIDs, this array is filled
      std::vector<panzer::GlobalOrdinal> gids;
      ugi.getElementGIDs(elements[e],gids);

      os << "   local element id = " << elements[e] << ", ";

      os << "  gids =";
      for(std::size_t i=0;i<gids.size();i++)
        os << " " << gids[i];

      os << ",  lids =";
      for(std::size_t i=0;i<gids.size();i++)
        os << " " << lids_host(i);
      os << std::endl;
    }
  }
}

Teuchos::RCP<Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
buildGhostedFieldReducedVector(const GlobalIndexer & ugi)
{
   using Map = Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
   using IntVector = Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>;

   std::vector<panzer::GlobalOrdinal> indices;
   std::vector<std::string> blocks;

   ugi.getOwnedAndGhostedIndices(indices);
   ugi.getElementBlockIds(blocks);

   std::vector<int> fieldNumbers(indices.size(),-1);

   Teuchos::RCP<Map> ghostedMap 
         = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(), Teuchos::arrayViewFromVector(indices),
                                Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::zero(), ugi.getComm()));

   // build a map from local ids to a field number
   for(std::size_t blk=0;blk<blocks.size();blk++) {
      std::string blockId = blocks[blk];

      const std::vector<panzer::LocalOrdinal> & elements = ugi.getElementBlock(blockId);
      const std::vector<int> & fields = ugi.getBlockFieldNumbers(blockId);
 
      // loop over all elements, and set field number in output array
      std::vector<panzer::GlobalOrdinal> gids(fields.size());
      for(std::size_t e=0;e<elements.size();e++) {
         ugi.getElementGIDs(elements[e],gids);

         for(std::size_t f=0;f<fields.size();f++) {
            int fieldNum = fields[f];
            panzer::GlobalOrdinal gid = gids[f];
            std::size_t lid = ghostedMap->getLocalElement(gid); // hash table lookup

            fieldNumbers[lid] = fieldNum; 
         }
      }
   }

   // produce a reduced vector containing only fields known by this processor
   std::vector<panzer::GlobalOrdinal> reducedIndices;
   std::vector<int> reducedFieldNumbers;
   for(std::size_t i=0;i<fieldNumbers.size();i++) {
      if(fieldNumbers[i]>-1) {
         reducedIndices.push_back(indices[i]);
         reducedFieldNumbers.push_back(fieldNumbers[i]);
      }
   }

   Teuchos::RCP<Map> reducedMap 
      = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(), Teuchos::arrayViewFromVector(reducedIndices),
                             Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::zero(), ugi.getComm()));
   return Teuchos::rcp(new IntVector(reducedMap,Teuchos::arrayViewFromVector(reducedFieldNumbers)));
}

void buildGhostedFieldVector(const GlobalIndexer & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > & reducedVec)
{
   using IntVector = Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
   Teuchos::RCP<const IntVector> dest = buildGhostedFieldVector(ugi,reducedVec);

   auto host_values = dest->getLocalViewHost(Tpetra::Access::ReadOnly);

   fieldNumbers.resize(dest->getLocalLength());
   for (size_t i=0; i < fieldNumbers.size(); ++i)
     fieldNumbers[i] = host_values(i,0);
}

Teuchos::RCP<const Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
buildGhostedFieldVector(const GlobalIndexer & ugi,
                        const Teuchos::RCP<const Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > & reducedVec)
{
   using Map = Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
   using IntVector = Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
   using Importer = Tpetra::Import<int,panzer::GlobalOrdinal,panzer::TpetraNodeType>;

   // first step: get a reduced field number vector and build a map to 
   // contain the full field number vector
   ///////////////////////////////////////////////////////////////////////////////

   Teuchos::RCP<Map> destMap;
   {
      std::vector<panzer::GlobalOrdinal> indices;
      ugi.getOwnedAndGhostedIndices(indices);
      destMap = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(), Teuchos::arrayViewFromVector(indices),
                                     Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::zero(), ugi.getComm()));
   }

   Teuchos::RCP<const IntVector> source = reducedVec;
   if(source==Teuchos::null)
     source = buildGhostedFieldReducedVector(ugi);
   Teuchos::RCP<const Map> sourceMap = source->getMap();

   // second step: perform the global communciation required to fix the
   // interface conditions (where this processor doesn't know what field  
   // some indices are)
   ///////////////////////////////////////////////////////////////////////////////
   Teuchos::RCP<IntVector> dest = Teuchos::rcp(new IntVector(destMap));
   Importer importer(sourceMap,destMap);

   dest->doImport(*source,importer,Tpetra::INSERT);
   PHX::Device::execution_space().fence();

   return dest;
}

/** Construct a map that only uses a certain field. */
Teuchos::RCP<const Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
getFieldMap(int fieldNum,const Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> & fieldTVector)
{
   Teuchos::RCP<const Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > origMap = fieldTVector.getMap();
   std::vector<int> fieldVector(fieldTVector.getLocalLength());
   fieldTVector.get1dCopy(Teuchos::arrayViewFromVector(fieldVector));

   std::vector<panzer::GlobalOrdinal> mapVector;
   for(std::size_t i=0;i<fieldVector.size();i++) { 
      if(fieldVector[i]==fieldNum)
         mapVector.push_back(origMap->getGlobalElement(i));
   }

   Teuchos::RCP<Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > finalMap 
      = Teuchos::rcp(new Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(
                                Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(), Teuchos::arrayViewFromVector(mapVector),
                                Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::zero(), origMap->getComm()));

   return finalMap;
}

ArrayToFieldVector::ArrayToFieldVector(const Teuchos::RCP<const GlobalIndexer> & ugi)
      : ugi_(ugi)
{
   gh_reducedFieldVector_ = buildGhostedFieldReducedVector(*ugi_);
   gh_fieldVector_ = buildGhostedFieldVector(*ugi_,gh_reducedFieldVector_);
}


void ArrayToFieldVector::buildFieldVector(const Tpetra::Vector<int,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> & source) const
{
   // build (unghosted) vector and map
   std::vector<panzer::GlobalOrdinal> indices;
   ugi_->getOwnedIndices(indices);

   Teuchos::RCP<const Map> destMap
         = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(), Teuchos::arrayViewFromVector(indices),
                                Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::zero(), ugi_->getComm()));
   Teuchos::RCP<IntVector> localFieldVector = Teuchos::rcp(new IntVector(destMap));

   Tpetra::Import<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> importer(source.getMap(),destMap);
   localFieldVector->doImport(source,importer,Tpetra::INSERT);

   fieldVector_ = localFieldVector;
}

Teuchos::RCP<const Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
ArrayToFieldVector::getFieldMap(const std::string & fieldName) const
{
   return getFieldMap(ugi_->getFieldNum(fieldName));
}

Teuchos::RCP<const Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
ArrayToFieldVector::getFieldMap(int fieldNum) const
{
   if(fieldMaps_[fieldNum]==Teuchos::null) {
      // if neccessary build field vector
      if(fieldVector_==Teuchos::null)
         buildFieldVector(*gh_fieldVector_);

      fieldMaps_[fieldNum] = panzer::getFieldMap(fieldNum,*fieldVector_);
   }

   return fieldMaps_[fieldNum];
}

// *****************************************
// namespace orientation_helpers
// *****************************************

namespace orientation_helpers {


void computeCellEdgeOrientations(const std::vector<std::pair<int,int> > & topEdgeIndices,
                                 const std::vector<panzer::GlobalOrdinal> & topology,
                                 const FieldPattern & fieldPattern, 
                                 std::vector<signed char> & orientation)
{
   // LOCAL element orientations are always set so that they flow in the positive
   // direction along an edge from node 0 to node 1. As a result if the GID of
   // node 0 is larger then node 1 then the GLOBAL orientation is -1 (and positive
   // otherwise). The local definition of the edge direction is defined by 
   // the shards cell topology.

   TEUCHOS_ASSERT(orientation.size()==std::size_t(fieldPattern.numberIds()));

   int edgeDim = 1;

   for(std::size_t e=0;e<topEdgeIndices.size();e++) {
      // grab topological nodes
      const std::pair<int,int> nodes = topEdgeIndices[e]; 

      // extract global values of topological nodes
      panzer::GlobalOrdinal v0 = topology[nodes.first];
      panzer::GlobalOrdinal v1 = topology[nodes.second];

      // using simple rule make a decision about orientation
      signed char edgeOrientation = 1;
      if(v1>v0)
         edgeOrientation = 1; 
      else if(v0>v1)
         edgeOrientation = -1; 
      else
      { TEUCHOS_ASSERT(false); }
      
      // grab edgeIndices to be set to compute orientation
      const std::vector<int> & edgeIndices = fieldPattern.getSubcellIndices(edgeDim,e);
      for(std::size_t s=0;s<edgeIndices.size();s++)
         orientation[edgeIndices[s]] = edgeOrientation;
   }
}

void computeCellFaceOrientations(const std::vector<std::vector<int> > & topFaceIndices,
                                 const std::vector<panzer::GlobalOrdinal> & topology,
                                 const FieldPattern & fieldPattern, 
                                 std::vector<signed char> & orientation)
{
   // LOCAL element orientations are always set so that they flow in the positive
   // direction away from the cell (counter clockwise rotation of the face). To determine
   // the face orientation we use the fact that Shards (and thus the field pattern) always
   // locally orders faces in a counter clockwise direction. A local rule for each element
   // will take advantage of this ensuring both elements agree on the orientation. This rule
   // is to first find the smallest node GID on the face. Then look at the GIDs of the nodes
   // immediately preceding and following that one. If the node following has smaller GID than
   // the preceding node then the face is oriented counter clockwise and thus the orientation
   // is +1. If the node preceding is larger, then the orientation is clockwise and the set to
   // a value of -1.

   // this only works for 3D field patterns
   TEUCHOS_ASSERT(fieldPattern.getDimension()==3);

   TEUCHOS_ASSERT(orientation.size()==std::size_t(fieldPattern.numberIds()));

   int faceDim = 2; 

   for(std::size_t f=0;f<topFaceIndices.size();f++) {
      // grab topological nodes
      const std::vector<int> & nodes = topFaceIndices[f]; 
      std::vector<panzer::GlobalOrdinal> globals(nodes.size());
      for(std::size_t n=0;n<nodes.size();n++)
         globals[n] = topology[nodes[n]]; 

      typename std::vector<panzer::GlobalOrdinal>::const_iterator itr 
          = std::min_element(globals.begin(),globals.end()); 

      TEUCHOS_TEST_FOR_EXCEPTION(itr==globals.end(),std::out_of_range,
                                 "panzer::orientation_helpers::computeCellFaceOrientations: A face index array "
                                 "was empty.");

      // extract global values of topological nodes
      // The face nodes go in counter clockwise order, let v_min be the
      // value with the minimum element then
      //         vbefore => itr => vafter
      // note that the nonsense with the beginning and end has to do with
      // if this iterator was the first or last in the array
      panzer::GlobalOrdinal vbefore = itr==globals.begin() ? *(globals.end()-1) : *(itr-1);
      panzer::GlobalOrdinal vafter = (itr+1)==globals.end() ? *globals.begin() : *(itr+1);

/*
      // sanity check in debug mode (uncomment these lines)
      TEUCHOS_ASSERT(std::find(globals.begin(),globals.end(),vbefore)!=globals.end());
      TEUCHOS_ASSERT(std::find(globals.begin(),globals.end(),vafter)!=globals.end());

      // print out information about the found nodes and also what 
      // order they were in originally
      std::cout << "\nFace Order = ";
      for(std::size_t l=0;l<globals.size();l++)
         std::cout << globals[l] << " ";
      std::cout << std::endl;
      std::cout << "(before,min,after) " << f << ": " << vbefore << " => " << *itr << " => " << vafter << std::endl;
*/

      // note by assumption
      // vbefore < *itr  and *itr < vafter

      // Based on the next lowest global id starting from the minimum
      signed char faceOrientation = 1;
      if(vafter>vbefore) // means smaller in clockwise direction
         faceOrientation = -1; 
      else if(vbefore>vafter) // means smaller in counter clockwise direction
         faceOrientation = 1; 
      else
      { TEUCHOS_ASSERT(false); } // we got an equality somehow!
      
      // grab faceIndices to be set to compute orientation
      const std::vector<int> & faceIndices = fieldPattern.getSubcellIndices(faceDim,f);
      for(std::size_t s=0;s<faceIndices.size();s++)
         orientation[faceIndices[s]] = faceOrientation;
   }
}

void computePatternEdgeIndices(const FieldPattern & pattern,std::vector<std::pair<int,int> > & edgeIndices)
{
   unsigned dim = 1;
   shards::CellTopology cellTopo = pattern.getCellTopology();
   for(unsigned e=0;e<cellTopo.getEdgeCount();e++) {
      // get local vertex ids for a this edge
      unsigned local_v0 = cellTopo.getNodeMap(dim,e,0);
      unsigned local_v1 = cellTopo.getNodeMap(dim,e,1);

      // get sub cell indices for geometric pattern
      const std::vector<int> & v0_indices = pattern.getSubcellIndices(0,local_v0);
      const std::vector<int> & v1_indices = pattern.getSubcellIndices(0,local_v1);

      TEUCHOS_ASSERT(v0_indices.size()>0); // there must be a node
      TEUCHOS_ASSERT(v1_indices.size()>0); // there must be a node

      // take the first index on each vertex and make a edge lookup
      edgeIndices.push_back(std::make_pair(v0_indices[0],v1_indices[0]));
   }
}

void computePatternFaceIndices(const FieldPattern & pattern,std::vector<std::vector<int> > & faceIndices)
{
   // this only works for 3D field patterns
   // TEUCHOS_ASSERT(pattern.getDimension()==3);
   //
   unsigned node_dim = 0; // by assumption
   unsigned subcell_dim = 2;

   if(pattern.getDimension()==3) {
      shards::CellTopology cellTopo = pattern.getCellTopology();
   
      faceIndices.resize(cellTopo.getSubcellCount(subcell_dim));
   
      for(unsigned f=0;f<cellTopo.getSubcellCount(subcell_dim);f++) {
         shards::CellTopology faceTopo(cellTopo.getBaseCellTopologyData(subcell_dim,f));
   
         for(unsigned v=0;v<faceTopo.getNodeCount();v++) {
            // get local vertex ids for a this edge
            unsigned local_v = cellTopo.getNodeMap(subcell_dim,f,v);
   
            // get sub cell indices for geometric pattern
            const std::vector<int> & v_indices = pattern.getSubcellIndices(node_dim,local_v);
      
            TEUCHOS_ASSERT(v_indices.size()>0); // there must be a node
      
            // take the first index on each vertex and make a edge lookup
            faceIndices[f].push_back(v_indices[0]);
         }
      }
   }
   else if(pattern.getDimension()==2) {
      shards::CellTopology cellTopo = pattern.getCellTopology();
   
      faceIndices.resize(1);
   
      for(unsigned v=0;v<cellTopo.getNodeCount();v++)
        faceIndices[0].push_back(v);
   }
}

} // end namespace orientation_helpers
} // end namespace panzer
