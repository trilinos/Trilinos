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

#include <vector>
#include <map>

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<Tpetra::Vector<int,int,GlobalOrdinalT,Node> >
buildGhostedFieldReducedVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi)
{
   typedef Tpetra::Map<int,GlobalOrdinalT,Node> Map;
   typedef Tpetra::Vector<int,int,GlobalOrdinalT,Node> IntVector;

   std::vector<GlobalOrdinalT> indices;
   std::vector<std::string> blocks;

   ugi.getOwnedAndSharedIndices(indices);
   ugi.getElementBlockIds(blocks);

   std::vector<int> fieldNumbers(indices.size(),-1);

   Teuchos::RCP<Map> sharedMap 
         = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(), Teuchos::arrayViewFromVector(indices),
                                Teuchos::OrdinalTraits<GlobalOrdinalT>::zero(), ugi.getComm()));

   // build a map from local ids to a field number
   for(std::size_t blk=0;blk<blocks.size();blk++) {
      std::string blockId = blocks[blk];

      const std::vector<LocalOrdinalT> & elements = ugi.getElementBlock(blockId);
      const std::vector<int> & fields = ugi.getBlockFieldNumbers(blockId);
 
      // loop over all elements, and set field number in output array
      std::vector<GlobalOrdinalT> gids(fields.size());
      for(std::size_t e=0;e<elements.size();e++) {
         ugi.getElementGIDs(elements[e],gids);

         for(std::size_t f=0;f<fields.size();f++) {
            int fieldNum = fields[f];
            GlobalOrdinalT gid = gids[f];
            std::size_t lid = sharedMap->getLocalElement(gid); // hash table lookup

            fieldNumbers[lid] = fieldNum; 
         }
      }
   }

   // produce a reduced vector containing only fields known by this processor
   std::vector<GlobalOrdinalT> reducedIndices;
   std::vector<int> reducedFieldNumbers;
   for(std::size_t i=0;i<fieldNumbers.size();i++) {
      if(fieldNumbers[i]>-1) {
         reducedIndices.push_back(indices[i]);
         reducedFieldNumbers.push_back(fieldNumbers[i]);
      }
   }

   Teuchos::RCP<Map> reducedMap 
      = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(), Teuchos::arrayViewFromVector(reducedIndices),
                             Teuchos::OrdinalTraits<GlobalOrdinalT>::zero(), ugi.getComm()));
   return Teuchos::rcp(new IntVector(reducedMap,Teuchos::arrayViewFromVector(reducedFieldNumbers)));
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,Node> > & reducedVec)
{
   typedef Tpetra::Vector<int,int,GlobalOrdinalT,Node> IntVector;

   Teuchos::RCP<const IntVector> dest = buildGhostedFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>(ugi,reducedVec);

   fieldNumbers.resize(dest->getLocalLength());
   dest->get1dCopy(Teuchos::arrayViewFromVector(fieldNumbers));
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,Node> >
buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                        const Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,Node> > & reducedVec)
{
   typedef Tpetra::Map<int,GlobalOrdinalT,Node> Map;
   typedef Tpetra::Vector<int,int,GlobalOrdinalT,Node> IntVector;
   typedef Tpetra::Import<int,GlobalOrdinalT,Node> Importer;

   // first step: get a reduced field number vector and build a map to 
   // contain the full field number vector
   ///////////////////////////////////////////////////////////////////////////////

   Teuchos::RCP<Map> destMap;
   {
      std::vector<GlobalOrdinalT> indices;
      ugi.getOwnedAndSharedIndices(indices);
      destMap = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(), Teuchos::arrayViewFromVector(indices),
                                     Teuchos::OrdinalTraits<GlobalOrdinalT>::zero(), ugi.getComm()));
   }

   Teuchos::RCP<const IntVector> source = reducedVec;
   if(source==Teuchos::null)
      source = buildGhostedFieldReducedVector<LocalOrdinalT,GlobalOrdinalT,Node>(ugi);
   Teuchos::RCP<const Map> sourceMap = source->getMap();

   // second step: perform the global communciation required to fix the
   // interface conditions (where this processor doesn't know what field  
   // some indices are)
   ///////////////////////////////////////////////////////////////////////////////
   Teuchos::RCP<IntVector> dest = Teuchos::rcp(new IntVector(destMap));
   Importer importer(sourceMap,destMap);

   dest->doImport(*source,importer,Tpetra::INSERT);

   return dest;
}

template <typename ScalarT,typename ArrayT,typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void updateGhostedDataReducedVector(const std::string & fieldName,const std::string blockId,
                                    const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                                    const ArrayT & data,Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> & dataVector)
{
   typedef Tpetra::Map<int,GlobalOrdinalT,Node> Map;
   typedef Tpetra::Vector<int,int,GlobalOrdinalT,Node> IntVector;
   typedef Tpetra::Import<int,GlobalOrdinalT,Node> Importer;

   TEUCHOS_TEST_FOR_EXCEPTION(!ugi.fieldInBlock(fieldName,blockId),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: field name = \""+fieldName+"\" is not in element block = \"" +blockId +"\"!");

   Teuchos::RCP<const Map> dataMap = dataVector.getMap();

   int fieldNum = ugi.getFieldNum(fieldName);
   const std::vector<LocalOrdinalT> & elements = ugi.getElementBlock(blockId);
   const std::vector<int> & fieldOffsets = ugi.getGIDFieldOffsets(blockId,fieldNum);
   
   TEUCHOS_TEST_FOR_EXCEPTION(data.dimension(0)!=(int) elements.size(),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: data cell dimension does not match up with block cell count");

   int rank = data.rank();

   if(rank==2) {
      // loop over elements distributing relevent data to vector
      std::vector<GlobalOrdinalT> gids;
      for(std::size_t e=0;e<elements.size();e++) { 
         ugi.getElementGIDs(elements[e],gids);
   
         for(std::size_t f=0;f<fieldOffsets.size();f++) {
            std::size_t localIndex = dataMap->getLocalElement(gids[fieldOffsets[f]]); // hash table lookup
            dataVector.replaceLocalValue(localIndex,0,data(e,f));
         }
      }
   }
   else if(rank==3) {
      std::size_t entries = data.dimension(2);
 
      TEUCHOS_TEST_FOR_EXCEPTION(dataVector.getNumVectors()!=entries,std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: number of columns in data vector inconsistent with data array");

      // loop over elements distributing relevent data to vector
      std::vector<GlobalOrdinalT> gids;
      for(std::size_t e=0;e<elements.size();e++) { 
         ugi.getElementGIDs(elements[e],gids);
   
         for(std::size_t f=0;f<fieldOffsets.size();f++) {
            std::size_t localIndex = dataMap->getLocalElement(gids[fieldOffsets[f]]); // hash table lookup
            for(std::size_t v=0;v<entries;v++)
               dataVector.replaceLocalValue(localIndex,v,data(e,f,v));
         }
      }
   }
   else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: data array rank must be 2 or 3");
}

/** Construct a map that only uses a certain field.
 */
template <typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Map<int,GlobalOrdinalT,Node> >
getFieldMap(int fieldNum,const Tpetra::Vector<int,int,GlobalOrdinalT,Node> & fieldTVector)
{
   Teuchos::RCP<const Tpetra::Map<int,GlobalOrdinalT,Node> > origMap = fieldTVector.getMap();
   std::vector<int> fieldVector(fieldTVector.getLocalLength());
   fieldTVector.get1dCopy(Teuchos::arrayViewFromVector(fieldVector));

   std::vector<GlobalOrdinalT> mapVector;
   for(std::size_t i=0;i<fieldVector.size();i++) { 
      if(fieldVector[i]==fieldNum)
         mapVector.push_back(origMap->getGlobalElement(i));
   }

   Teuchos::RCP<Tpetra::Map<int,GlobalOrdinalT,Node> > finalMap 
      = Teuchos::rcp(new Tpetra::Map<int,GlobalOrdinalT,Node>(
                                Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(), Teuchos::arrayViewFromVector(mapVector),
                                Teuchos::OrdinalTraits<GlobalOrdinalT>::zero(), origMap->getComm()));

   return finalMap;
}

namespace orientation_helpers {

template <typename GlobalOrdinalT>
void computeCellEdgeOrientations(const std::vector<std::pair<int,int> > & topEdgeIndices,
                                 const std::vector<GlobalOrdinalT> & topology,
                                 const FieldPattern & fieldPattern, 
                                 std::vector<char> & orientation)
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
      GlobalOrdinalT v0 = topology[nodes.first];
      GlobalOrdinalT v1 = topology[nodes.second];

      // using simple rule make a decision about orientation
      char edgeOrientation = 1; 
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

}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
   ArrayToFieldVector(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & ugi)
      : ugi_(ugi)
{
   gh_reducedFieldVector_ = buildGhostedFieldReducedVector<LocalOrdinalT,GlobalOrdinalT,Node>(*ugi_);
   gh_fieldVector_ = buildGhostedFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>(*ugi_,gh_reducedFieldVector_);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> >
ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
   getGhostedDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
{
   int fieldNum = ugi_->getFieldNum(fieldName);
   std::vector<std::string> blockIds;
   ugi_->getElementBlockIds(blockIds);

   // get rank of first data array, determine column count
   int rank = data.begin()->second.rank();
   int numCols = 0;
   if(rank==2)
      numCols = 1;
   else if(rank==3)
      numCols = data.begin()->second.dimension(2);
   else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                          "ArrayToFieldVector::getGhostedDataVector: data array must have rank 2 or 3");


   // first build and fill in final reduced field vector
   /////////////////////////////////////////////////////////////////

   // build field maps as needed
   Teuchos::RCP<const Map> reducedMap = gh_reducedFieldMaps_[fieldNum];
   if(gh_reducedFieldMaps_[fieldNum]==Teuchos::null) {
      reducedMap = panzer::getFieldMap(fieldNum,*gh_reducedFieldVector_);
      gh_reducedFieldMaps_[fieldNum] = reducedMap;
   }

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> > finalReducedVec
      = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node>(reducedMap,numCols));
   for(std::size_t b=0;b<blockIds.size();b++) {
      std::string block = blockIds[b];

      // make sure field is in correct block
      if(!ugi_->fieldInBlock(fieldName,block))
         continue; 

      // extract data vector
      typename std::map<std::string,ArrayT>::const_iterator blockItr = data.find(block);
     TEUCHOS_TEST_FOR_EXCEPTION(blockItr==data.end(),std::runtime_error,
                        "ArrayToFieldVector::getDataVector: can not find block \""+block+"\".");

     const ArrayT & d = blockItr->second;
     updateGhostedDataReducedVector<ScalarT,ArrayT,LocalOrdinalT,GlobalOrdinalT,Node>(fieldName,block,*ugi_,d,*finalReducedVec); 
   }

   // build final (not reduced vector)
   /////////////////////////////////////////////

   Teuchos::RCP<const Map> map = gh_fieldMaps_[fieldNum];
   if(gh_fieldMaps_[fieldNum]==Teuchos::null) {
      map = panzer::getFieldMap(fieldNum,*gh_fieldVector_);
      gh_fieldMaps_[fieldNum] = map;
   }

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> > finalVec
      = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node>(map,numCols));

   // do import from finalReducedVec
   Tpetra::Import<int,GlobalOrdinalT,Node> importer(reducedMap,map);
   finalVec->doImport(*finalReducedVec,importer,Tpetra::INSERT);

   return finalVec;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> >
ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
   getDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
{
   // if neccessary build field vector
   if(fieldVector_==Teuchos::null)
      buildFieldVector(*gh_fieldVector_);

   Teuchos::RCP<const Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> > sourceVec
         = getGhostedDataVector<ScalarT,ArrayT>(fieldName,data);

   // use lazy construction for each field
   int fieldNum = ugi_->getFieldNum(fieldName);
   Teuchos::RCP<const Map> destMap = fieldMaps_[fieldNum];
   if(fieldMaps_[fieldNum]==Teuchos::null) {
      destMap = panzer::getFieldMap(fieldNum,*fieldVector_);
      fieldMaps_[fieldNum] = destMap;
   }

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> > destVec
         = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node>(destMap,sourceVec->getNumVectors()));
   
   // do import
   Tpetra::Import<int,GlobalOrdinalT> importer(sourceVec->getMap(),destMap);
   destVec->doImport(*sourceVec,importer,Tpetra::INSERT); 

   return destVec;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
        buildFieldVector(const Tpetra::Vector<int,int,GlobalOrdinalT,Node> & source) const
{
   // build (unghosted) vector and map
   std::vector<GlobalOrdinalT> indices;
   ugi_->getOwnedIndices(indices);

   Teuchos::RCP<const Map> destMap
         = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(), Teuchos::arrayViewFromVector(indices),
                                Teuchos::OrdinalTraits<GlobalOrdinalT>::zero(), ugi_->getComm()));
   Teuchos::RCP<IntVector> localFieldVector = Teuchos::rcp(new IntVector(destMap));

   Tpetra::Import<int,GlobalOrdinalT> importer(source.getMap(),destMap);
   localFieldVector->doImport(source,importer,Tpetra::INSERT);

   fieldVector_ = localFieldVector;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Map<int,GlobalOrdinalT,Node> >
ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
getFieldMap(const std::string & fieldName) const
{
   int fieldNum = ugi_->getFieldNum(fieldName);

   if(fieldMaps_[fieldNum]==Teuchos::null) {
      // if neccessary build field vector
      if(fieldVector_==Teuchos::null)
         buildFieldVector(*gh_fieldVector_);

      fieldMaps_[fieldNum] = panzer::getFieldMap(fieldNum,*fieldVector_);
   }

   return fieldMaps_[fieldNum];
}
                                   
} // end namspace panzer
