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
#ifndef PANZER_UNIQUEGLOBALINDEXER_EPETRAUTILITIES_IMPL_HPP
#define PANZER_UNIQUEGLOBALINDEXER_EPETRAUTILITIES_IMPL_HPP

#include <vector>
#include <map>

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_MpiComm.h"

#include "Epetra_CombineMode.h"
#include "Epetra_DataAccess.h"

#include <sstream>
#include <cmath>

namespace panzer {


template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<Epetra_IntVector>
buildGhostedFieldReducedVectorEpetra(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi)
{
   typedef Epetra_BlockMap Map;
   typedef Epetra_IntVector IntVector;

   std::vector<GlobalOrdinalT> indices;
   std::vector<std::string> blocks;

   ugi.getOwnedAndGhostedIndices(indices);
   ugi.getElementBlockIds(blocks);

   std::vector<int> fieldNumbers(indices.size(),-1);

   const Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(ugi.getComm());
   Teuchos::RCP<Epetra_MpiComm> comm;
   if (mpiComm != Teuchos::null)
     comm = Teuchos::rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()));

   Teuchos::RCP<Map> ghostedMap
     = Teuchos::rcp(new Map(-1, static_cast<int>(indices.size()), Teuchos::arrayViewFromVector(indices).getRawPtr(),
                            1, Teuchos::OrdinalTraits<int>::zero(), *comm));

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
            std::size_t lid = ghostedMap->LID(gid); // hash table lookup

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
     = Teuchos::rcp(new Map(-1, static_cast<int>(reducedIndices.size()), Teuchos::arrayViewFromVector(reducedIndices).getRawPtr(),
                            1, Teuchos::OrdinalTraits<int>::zero(), *comm));
   return Teuchos::rcp(new IntVector(Copy,*reducedMap,Teuchos::arrayViewFromVector(reducedFieldNumbers).getRawPtr()));
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void buildGhostedFieldVectorEpetra(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Epetra_IntVector> & reducedVec)
{
   typedef Epetra_IntVector IntVector;

   Teuchos::RCP<const IntVector> dest = buildGhostedFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>(ugi,reducedVec);

   fieldNumbers.resize(dest->MyLength());
   Teuchos::ArrayView<int> av = Teuchos::arrayViewFromVector(fieldNumbers);
   dest->ExtractCopy(av.getRawPtr());
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Epetra_IntVector>
buildGhostedFieldVectorEpetra(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                              const Teuchos::RCP<const Epetra_IntVector> & reducedVec)
{
   typedef Epetra_BlockMap Map;
   typedef Epetra_IntVector IntVector;
   typedef Epetra_Import Importer;

   // first step: get a reduced field number vector and build a map to
   // contain the full field number vector
   ///////////////////////////////////////////////////////////////////////////////

   Teuchos::RCP<Map> destMap;
   {
      std::vector<GlobalOrdinalT> indices;
      ugi.getOwnedAndGhostedIndices(indices);

      const Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(ugi.getComm());
      Teuchos::RCP<Epetra_MpiComm> comm;
      if (mpiComm != Teuchos::null)
        comm = Teuchos::rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()));

      destMap = Teuchos::rcp(new Map(-1, static_cast<int>(indices.size()), Teuchos::arrayViewFromVector(indices).getRawPtr(),
                                     1, Teuchos::OrdinalTraits<int>::zero(), *comm));
   }

   Teuchos::RCP<const IntVector> source = reducedVec;
   if(source==Teuchos::null)
      source = buildGhostedFieldReducedVectorEpetra<LocalOrdinalT,GlobalOrdinalT>(ugi);
   Teuchos::RCP<const Map> sourceMap = Teuchos::rcpFromRef(source->Map());

   // second step: perform the global communciation required to fix the
   // interface conditions (where this processor doesn't know what field
   // some indices are)
   ///////////////////////////////////////////////////////////////////////////////
   Teuchos::RCP<IntVector> dest = Teuchos::rcp(new IntVector(*destMap));
   Importer importer(*destMap,*sourceMap);

   dest->Import(*source,importer,Insert);

   return dest;
}

template <typename ScalarT,typename ArrayT,typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void updateGhostedDataReducedVectorEpetra(const std::string & fieldName,const std::string blockId,
                                          const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                                          const ArrayT & data,Epetra_MultiVector & dataVector)
{
   typedef Epetra_BlockMap Map;

   TEUCHOS_TEST_FOR_EXCEPTION(!ugi.fieldInBlock(fieldName,blockId),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: field name = \""+fieldName+"\" is not in element block = \"" +blockId +"\"!");

   Teuchos::RCP<const Map> dataMap = Teuchos::rcpFromRef(dataVector.Map());

   int fieldNum = ugi.getFieldNum(fieldName);
   const std::vector<LocalOrdinalT> & elements = ugi.getElementBlock(blockId);
   const std::vector<int> & fieldOffsets = ugi.getGIDFieldOffsets(blockId,fieldNum);
   
   TEUCHOS_TEST_FOR_EXCEPTION(data.extent(0)!=elements.size(),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: data cell dimension does not match up with block cell count");

   int rank = data.rank();

   if(rank==2) {
      // loop over elements distributing relevent data to vector
      std::vector<GlobalOrdinalT> gids;
      for(std::size_t e=0;e<elements.size();e++) {
         ugi.getElementGIDs(elements[e],gids);
   
         for(std::size_t f=0;f<fieldOffsets.size();f++) {
            std::size_t localIndex = dataMap->LID(gids[fieldOffsets[f]]); // hash table lookup
            dataVector.ReplaceMyValue(localIndex,0,data(e,f));
         }
      }
   }
   else if(rank==3) {
      std::size_t entries = data.extent(2);
 
      TEUCHOS_TEST_FOR_EXCEPTION(dataVector.NumVectors()!=static_cast<int>(entries),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: number of columns in data vector inconsistent with data array");

      // loop over elements distributing relevent data to vector
      std::vector<GlobalOrdinalT> gids;
      for(std::size_t e=0;e<elements.size();e++) {
         ugi.getElementGIDs(elements[e],gids);
   
         for(std::size_t f=0;f<fieldOffsets.size();f++) {
            std::size_t localIndex = dataMap->LID(gids[fieldOffsets[f]]); // hash table lookup
            for(std::size_t v=0;v<entries;v++)
               dataVector.ReplaceMyValue(localIndex,v,data(e,f,v));
         }
      }
   }
   else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: data array rank must be 2 or 3");
}


template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
ArrayToFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>::
   ArrayToFieldVectorEpetra(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & ugi)
      : ugi_(ugi)
{
   gh_reducedFieldVector_ = buildGhostedFieldReducedVectorEpetra<LocalOrdinalT,GlobalOrdinalT>(*ugi_);
   gh_fieldVector_ = buildGhostedFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>(*ugi_,gh_reducedFieldVector_);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Epetra_MultiVector>
ArrayToFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>::
   getGhostedDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
{
   TEUCHOS_ASSERT(data.size()>0); // there must be at least one "data" item

   int fieldNum = ugi_->getFieldNum(fieldName);
   std::vector<std::string> blockIds;
   ugi_->getElementBlockIds(blockIds);

   // get rank of first data array, determine column count
   int rank = data.begin()->second.rank();
   int numCols = 0;
   if(rank==2)
      numCols = 1;
   else if(rank==3)
      numCols = data.begin()->second.extent(2);
   else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                          "ArrayToFieldVectorEpetra::getGhostedDataVector: data array must have rank 2 or 3. This array has rank " << rank << ".");
   }


   // first build and fill in final reduced field vector
   /////////////////////////////////////////////////////////////////

   // build field maps as needed
   Teuchos::RCP<const Map> reducedMap = gh_reducedFieldMaps_[fieldNum];
   if(gh_reducedFieldMaps_[fieldNum]==Teuchos::null) {
      reducedMap = panzer::getFieldMapEpetra(fieldNum,*gh_reducedFieldVector_);
      gh_reducedFieldMaps_[fieldNum] = reducedMap;
   }

   Teuchos::RCP<MultiVector> finalReducedVec
      = Teuchos::rcp(new MultiVector(*reducedMap,numCols));
   for(std::size_t b=0;b<blockIds.size();b++) {
      std::string block = blockIds[b];

      // make sure field is in correct block
      if(!ugi_->fieldInBlock(fieldName,block))
         continue; 

      // extract data vector
      typename std::map<std::string,ArrayT>::const_iterator blockItr = data.find(block);
     TEUCHOS_TEST_FOR_EXCEPTION(blockItr==data.end(),std::runtime_error,
                        "ArrayToFieldVectorEpetra::getDataVector: can not find block \""+block+"\".");

     const ArrayT & d = blockItr->second;
     updateGhostedDataReducedVectorEpetra<ScalarT,ArrayT,LocalOrdinalT,GlobalOrdinalT,Node>(fieldName,block,*ugi_,d,*finalReducedVec);
   }

   // build final (not reduced vector)
   /////////////////////////////////////////////

   Teuchos::RCP<const Map> map = gh_fieldMaps_[fieldNum];
   if(gh_fieldMaps_[fieldNum]==Teuchos::null) {
      map = panzer::getFieldMapEpetra(fieldNum,*gh_fieldVector_);
      gh_fieldMaps_[fieldNum] = map;
   }

   Teuchos::RCP<MultiVector> finalVec
      = Teuchos::rcp(new MultiVector(*map,numCols));

   // do import from finalReducedVec
   Epetra_Import importer(*map,*reducedMap);
   finalVec->Import(*finalReducedVec,importer,Insert);

   return finalVec;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Epetra_MultiVector>
ArrayToFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>::
   getDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
{
   // if neccessary build field vector
   if(fieldVector_==Teuchos::null)
      buildFieldVector(*gh_fieldVector_);

   Teuchos::RCP<const MultiVector> sourceVec
         = getGhostedDataVector<ScalarT,ArrayT>(fieldName,data);

   // use lazy construction for each field
   int fieldNum = ugi_->getFieldNum(fieldName);
   Teuchos::RCP<const Map> destMap = fieldMaps_[fieldNum];
   if(fieldMaps_[fieldNum]==Teuchos::null) {
      destMap = panzer::getFieldMapEpetra(fieldNum,*fieldVector_);
      fieldMaps_[fieldNum] = destMap;
   }

   Teuchos::RCP<MultiVector> destVec
         = Teuchos::rcp(new MultiVector(*destMap,sourceVec->NumVectors()));
   
   // do import
   Epetra_Import importer(*destMap,sourceVec->Map());
   destVec->Import(*sourceVec,importer,Insert);

   return destVec;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ArrayToFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>::
        buildFieldVector(const Epetra_IntVector & source) const
{
   // build (unghosted) vector and map
   std::vector<GlobalOrdinalT> indices;
   ugi_->getOwnedIndices(indices);

   const Teuchos::RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(ugi_->getComm(),true);
   Teuchos::RCP<Epetra_MpiComm> comm;
   if (mpiComm != Teuchos::null)
     comm = Teuchos::rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()));

   Teuchos::RCP<Map> destMap
     = Teuchos::rcp(new Map(-1, static_cast<int>(indices.size()),
                            Teuchos::arrayViewFromVector(indices).getRawPtr(),
                            // &indices.begin(),
                            1, Teuchos::OrdinalTraits<int>::zero(), *comm));

   Teuchos::RCP<IntVector> localFieldVector = Teuchos::rcp(new IntVector(*destMap));

   Epetra_Import importer(*destMap,source.Map());
   localFieldVector->Import(source,importer,Insert);

   fieldVector_ = localFieldVector;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Epetra_Map>
ArrayToFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>::
getFieldMap(const std::string & fieldName) const
{
   return getFieldMap(ugi_->getFieldNum(fieldName));
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Epetra_Map>
ArrayToFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>::
getFieldMap(int fieldNum) const
{
   if(fieldMaps_[fieldNum]==Teuchos::null) {
      // if neccessary build field vector
      if(fieldVector_==Teuchos::null)
         buildFieldVector(*gh_fieldVector_);

      fieldMaps_[fieldNum] = panzer::getFieldMapEpetra(fieldNum,*fieldVector_);
   }

   return fieldMaps_[fieldNum];
}
                                   
} // end namspace panzer

#endif
