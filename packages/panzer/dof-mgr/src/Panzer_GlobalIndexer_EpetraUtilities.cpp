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

#ifdef PANZER_HAVE_EPETRA

#include "Panzer_GlobalIndexer_EpetraUtilities.hpp"

namespace panzer {

Teuchos::RCP<Epetra_IntVector>
buildGhostedFieldReducedVectorEpetra(const GlobalIndexer & ugi)
{
   typedef Epetra_BlockMap Map;
   typedef Epetra_IntVector IntVector;

   std::vector<panzer::GlobalOrdinal> indices;
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

      const std::vector<panzer::LocalOrdinal> & elements = ugi.getElementBlock(blockId);
      const std::vector<int> & fields = ugi.getBlockFieldNumbers(blockId);
 
      // loop over all elements, and set field number in output array
      std::vector<panzer::GlobalOrdinal> gids(fields.size());
      for(std::size_t e=0;e<elements.size();e++) {
         ugi.getElementGIDs(elements[e],gids);

         for(std::size_t f=0;f<fields.size();f++) {
            int fieldNum = fields[f];
            panzer::GlobalOrdinal gid = gids[f];
            std::size_t lid = ghostedMap->LID(gid); // hash table lookup

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
     = Teuchos::rcp(new Map(-1, static_cast<int>(reducedIndices.size()), Teuchos::arrayViewFromVector(reducedIndices).getRawPtr(),
                            1, Teuchos::OrdinalTraits<int>::zero(), *comm));
   return Teuchos::rcp(new IntVector(Copy,*reducedMap,Teuchos::arrayViewFromVector(reducedFieldNumbers).getRawPtr()));
}

void buildGhostedFieldVectorEpetra(const GlobalIndexer & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Epetra_IntVector> & reducedVec)
{
   typedef Epetra_IntVector IntVector;

   Teuchos::RCP<const IntVector> dest = buildGhostedFieldVectorEpetra(ugi,reducedVec);

   fieldNumbers.resize(dest->MyLength());
   Teuchos::ArrayView<int> av = Teuchos::arrayViewFromVector(fieldNumbers);
   dest->ExtractCopy(av.getRawPtr());
}

Teuchos::RCP<const Epetra_IntVector>
buildGhostedFieldVectorEpetra(const GlobalIndexer & ugi,
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
      std::vector<panzer::GlobalOrdinal> indices;
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
      source = buildGhostedFieldReducedVectorEpetra(ugi);
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

ArrayToFieldVectorEpetra::ArrayToFieldVectorEpetra(const Teuchos::RCP<const GlobalIndexer> & ugi)
      : ugi_(ugi)
{
   gh_reducedFieldVector_ = buildGhostedFieldReducedVectorEpetra(*ugi_);
   gh_fieldVector_ = buildGhostedFieldVectorEpetra(*ugi_,gh_reducedFieldVector_);
}


void ArrayToFieldVectorEpetra::buildFieldVector(const Epetra_IntVector & source) const
{
   // build (unghosted) vector and map
   std::vector<panzer::GlobalOrdinal> indices;
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

Teuchos::RCP<const Epetra_Map>
ArrayToFieldVectorEpetra::getFieldMap(const std::string & fieldName) const
{
   return getFieldMap(ugi_->getFieldNum(fieldName));
}

Teuchos::RCP<const Epetra_Map>
ArrayToFieldVectorEpetra::getFieldMap(int fieldNum) const
{
   if(fieldMaps_[fieldNum]==Teuchos::null) {
      // if neccessary build field vector
      if(fieldVector_==Teuchos::null)
         buildFieldVector(*gh_fieldVector_);

      fieldMaps_[fieldNum] = panzer::getFieldMapEpetra(fieldNum,*fieldVector_);
   }

   return Teuchos::rcp_dynamic_cast<const Epetra_Map>(fieldMaps_[fieldNum]);
}

Teuchos::RCP<const Epetra_BlockMap>
getFieldMapEpetra(int fieldNum,const Epetra_IntVector & fieldTVector)
{
   Teuchos::RCP<const Epetra_BlockMap> origMap = Teuchos::rcpFromRef(fieldTVector.Map());
   std::vector<int> fieldVector(fieldTVector.MyLength());
   Teuchos::ArrayView<int> av = Teuchos::arrayViewFromVector(fieldVector);
   fieldTVector.ExtractCopy(av.getRawPtr());

   std::vector<int> mapVector;
   for(std::size_t i=0;i<fieldVector.size();i++) {
      if(fieldVector[i]==fieldNum)
         mapVector.push_back(origMap->GID(i));
   }

   Teuchos::RCP<Epetra_BlockMap> finalMap
     = Teuchos::rcp(new Epetra_BlockMap(-1, static_cast<int>(mapVector.size()), Teuchos::arrayViewFromVector(mapVector).getRawPtr(), 1, Teuchos::OrdinalTraits<int>::zero(), origMap->Comm()));

   return finalMap;
}

} // end namespace panzer

#endif //end PANZER_HAVE_EPETRA
