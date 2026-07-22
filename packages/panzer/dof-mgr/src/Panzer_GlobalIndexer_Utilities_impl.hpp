// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <vector>
#include <map>

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"

#include "Kokkos_DynRankView.hpp"

#include <sstream>
#include <cmath>

namespace panzer {

template <typename ScalarT,typename ArrayT>
void updateGhostedDataReducedVector(const std::string & fieldName,const std::string blockId,
                                    const GlobalIndexer & ugi,
                                    const ArrayT & data,Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> & dataVector)
{
   typedef Tpetra::Map<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> Map;

   TEUCHOS_TEST_FOR_EXCEPTION(!ugi.fieldInBlock(fieldName,blockId),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: field name = \""+fieldName+"\" is not in element block = \"" +blockId +"\"!");

   Teuchos::RCP<const Map> dataMap = dataVector.getMap();

   int fieldNum = ugi.getFieldNum(fieldName);
   const std::vector<panzer::LocalOrdinal> & elements = ugi.getElementBlock(blockId);
   const std::vector<int> & fieldOffsets = ugi.getGIDFieldOffsets(blockId,fieldNum);

   TEUCHOS_TEST_FOR_EXCEPTION(data.extent(0)!=elements.size(),std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: data cell dimension does not match up with block cell count");

   int rank = data.rank();
   auto dataVector_host = dataVector.getLocalViewHost(Tpetra::Access::ReadWrite);
   auto data_host = Kokkos::create_mirror_view(data);
   Kokkos::deep_copy(data_host,data);

   if(rank==2) {
      // loop over elements distributing relevent data to vector
      std::vector<panzer::GlobalOrdinal> gids;
      for(std::size_t e=0;e<elements.size();e++) { 
         ugi.getElementGIDs(elements[e],gids);
   
         for(std::size_t f=0;f<fieldOffsets.size();f++) {
            std::size_t localIndex = dataMap->getLocalElement(gids[fieldOffsets[f]]); // hash table lookup
            dataVector_host(localIndex,0) = data_host(e,f);
         }
      }
   }
   else if(rank==3) {
      std::size_t entries = data.extent(2);
 
      TEUCHOS_TEST_FOR_EXCEPTION(dataVector.getNumVectors()!=entries,std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: number of columns in data vector inconsistent with data array");

      // loop over elements distributing relevent data to vector
      std::vector<panzer::GlobalOrdinal> gids;
      for(std::size_t e=0;e<elements.size();e++) { 
         ugi.getElementGIDs(elements[e],gids);
   
         for(std::size_t f=0;f<fieldOffsets.size();f++) {
            std::size_t localIndex = dataMap->getLocalElement(gids[fieldOffsets[f]]); // hash table lookup
            for(std::size_t v=0;v<entries;v++)
              dataVector_host(localIndex,v) = data_host(e,f,v);
         }
      }
   }
   else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                      "panzer::updateGhostedDataReducedVector: data array rank must be 2 or 3");
}

template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
ArrayToFieldVector::getGhostedDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
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
                          "ArrayToFieldVector::getGhostedDataVector: data array must have rank 2 or 3. This array has rank " << rank << ".");
   }


   // first build and fill in final reduced field vector
   /////////////////////////////////////////////////////////////////

   // build field maps as needed
   Teuchos::RCP<const Map> reducedMap = gh_reducedFieldMaps_[fieldNum];
   if(gh_reducedFieldMaps_[fieldNum]==Teuchos::null) {
      reducedMap = panzer::getFieldMap(fieldNum,*gh_reducedFieldVector_);
      gh_reducedFieldMaps_[fieldNum] = reducedMap;
   }

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > finalReducedVec
      = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(reducedMap,numCols));
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
     updateGhostedDataReducedVector<ScalarT,ArrayT>(fieldName,block,*ugi_,d,*finalReducedVec); 
   }

   // build final (not reduced vector)
   /////////////////////////////////////////////

   Teuchos::RCP<const Map> map = gh_fieldMaps_[fieldNum];
   if(gh_fieldMaps_[fieldNum]==Teuchos::null) {
      map = panzer::getFieldMap(fieldNum,*gh_fieldVector_);
      gh_fieldMaps_[fieldNum] = map;
   }

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > finalVec
      = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(map,numCols));

   // do import from finalReducedVec
   Tpetra::Import<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> importer(reducedMap,map);
   finalVec->doImport(*finalReducedVec,importer,Tpetra::INSERT);
   PHX::Device::execution_space().fence();

   return finalVec;
}

template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
ArrayToFieldVector::getDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
{
   // if neccessary build field vector
   if(fieldVector_==Teuchos::null)
      buildFieldVector(*gh_fieldVector_);

   Teuchos::RCP<const Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > sourceVec
         = getGhostedDataVector<ScalarT,ArrayT>(fieldName,data);

   // use lazy construction for each field
   int fieldNum = ugi_->getFieldNum(fieldName);
   Teuchos::RCP<const Map> destMap = fieldMaps_[fieldNum];
   if(fieldMaps_[fieldNum]==Teuchos::null) {
      destMap = panzer::getFieldMap(fieldNum,*fieldVector_);
      fieldMaps_[fieldNum] = destMap;
   }

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > destVec
         = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,int,panzer::GlobalOrdinal,panzer::TpetraNodeType>(destMap,sourceVec->getNumVectors()));
   
   // do import
   Tpetra::Import<int,panzer::GlobalOrdinal,panzer::TpetraNodeType> importer(sourceVec->getMap(),destMap);
   destVec->doImport(*sourceVec,importer,Tpetra::INSERT); 
   PHX::Device::execution_space().fence();

   return destVec;
}
                                   
} // end namspace panzer
