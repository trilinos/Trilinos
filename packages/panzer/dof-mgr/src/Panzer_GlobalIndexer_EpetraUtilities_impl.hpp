// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

template <typename ScalarT,typename ArrayT>
void updateGhostedDataReducedVectorEpetra(const std::string & fieldName,const std::string blockId,
                                          const GlobalIndexer & ugi,
                                          const ArrayT & data,Epetra_MultiVector & dataVector)
{
   typedef Epetra_BlockMap Map;

   TEUCHOS_TEST_FOR_EXCEPTION(!ugi.fieldInBlock(fieldName,blockId),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: field name = \""+fieldName+"\" is not in element block = \"" +blockId +"\"!");

   Teuchos::RCP<const Map> dataMap = Teuchos::rcpFromRef(dataVector.Map());

   int fieldNum = ugi.getFieldNum(fieldName);
   const std::vector<panzer::LocalOrdinal> & elements = ugi.getElementBlock(blockId);
   const std::vector<int> & fieldOffsets = ugi.getGIDFieldOffsets(blockId,fieldNum);
   
   TEUCHOS_TEST_FOR_EXCEPTION(data.extent(0)!=elements.size(),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: data cell dimension does not match up with block cell count");

   int rank = data.rank();

   if(rank==2) {
      // loop over elements distributing relevent data to vector
      std::vector<panzer::GlobalOrdinal> gids;
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
      std::vector<panzer::GlobalOrdinal> gids;
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

template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Epetra_MultiVector>
ArrayToFieldVectorEpetra::getGhostedDataVector(const std::string & fieldName,
                                               const std::map<std::string,ArrayT> & data) const
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
     updateGhostedDataReducedVectorEpetra<ScalarT,ArrayT>(fieldName,block,*ugi_,d,*finalReducedVec);
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

template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Epetra_MultiVector>
ArrayToFieldVectorEpetra::getDataVector(const std::string & fieldName,
                                        const std::map<std::string,ArrayT> & data) const
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
                                   
} // end namspace panzer

#endif
