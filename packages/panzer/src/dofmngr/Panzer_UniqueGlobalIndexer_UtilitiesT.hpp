#include <vector>
#include <map>

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> >
buildGhostedFieldReducedVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi)
{
   typedef Tpetra::Map<std::size_t,GlobalOrdinalT,Node> Map;
   typedef Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> IntVector;

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
                             const Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> > & reducedVec)
{
   typedef Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> IntVector;

   Teuchos::RCP<const IntVector> dest = buildGhostedFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>(ugi,reducedVec);

   fieldNumbers.resize(dest->getLocalLength());
   dest->get1dCopy(Teuchos::arrayViewFromVector(fieldNumbers));
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> >
buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                        const Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> > & reducedVec)
{
   typedef Tpetra::Map<std::size_t,GlobalOrdinalT,Node> Map;
   typedef Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> IntVector;
   typedef Tpetra::Import<std::size_t,GlobalOrdinalT,Node> Importer;

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
                                    const ArrayT & data,Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node> & dataVector)
{
   typedef Tpetra::Map<std::size_t,GlobalOrdinalT,Node> Map;
   typedef Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> IntVector;
   typedef Tpetra::Import<std::size_t,GlobalOrdinalT,Node> Importer;

   TEST_FOR_EXCEPTION(!ugi.fieldInBlock(fieldName,blockId),std::runtime_error,
                      "panzer::buildGhostedDataReducedVector: field name = \""+fieldName+"\" is not in element block = \"" +blockId +"\"!");

   Teuchos::RCP<const Map> dataMap = dataVector.getMap();

   int fieldNum = ugi.getFieldNum(fieldName);
   const std::vector<LocalOrdinalT> & elements = ugi.getElementBlock(blockId);
   const std::vector<int> & fieldOffsets = ugi.getGIDFieldOffsets(blockId,fieldNum);
   
   TEST_FOR_EXCEPTION(data.dimension(0)!=(int) elements.size(),std::runtime_error,
                      "panzer::buildGhostedDataReducedVector: data cell dimension does not match up with block cell count");

   // loop over elements distributing relevent data to vector
   std::vector<GlobalOrdinalT> gids;
   for(std::size_t e=0;e<elements.size();e++) { 
      ugi.getElementGIDs(elements[e],gids);

      for(std::size_t f=0;f<fieldOffsets.size();f++) {
         std::size_t localIndex = dataMap->getLocalElement(gids[fieldOffsets[f]]); // hash table lookup
         dataVector.replaceLocalValue(localIndex,data(e,f));
      }
   }
}

/** Construct a map that only uses a certain field.
 */
template <typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Map<std::size_t,GlobalOrdinalT,Node> >
getFieldMap(int fieldNum,const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> & fieldTVector)
{
   Teuchos::RCP<const Tpetra::Map<std::size_t,GlobalOrdinalT,Node> > origMap = fieldTVector.getMap();
   std::vector<int> fieldVector(fieldTVector.getLocalLength());
   fieldTVector.get1dCopy(Teuchos::arrayViewFromVector(fieldVector));

   std::vector<GlobalOrdinalT> mapVector;
   for(std::size_t i=0;i<fieldVector.size();i++) { 
      if(fieldVector[i]==fieldNum)
         mapVector.push_back(origMap->getGlobalElement(i));
   }

   Teuchos::RCP<Tpetra::Map<std::size_t,GlobalOrdinalT,Node> > finalMap 
      = Teuchos::rcp(new Tpetra::Map<std::size_t,GlobalOrdinalT,Node>(
                                Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(), Teuchos::arrayViewFromVector(mapVector),
                                Teuchos::OrdinalTraits<GlobalOrdinalT>::zero(), origMap->getComm()));

   return finalMap;
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
Teuchos::RCP<Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node> >
ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
   getGhostedDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
{
   int fieldNum = ugi_->getFieldNum(fieldName);
   std::vector<std::string> blockIds;
   ugi_->getElementBlockIds(blockIds);

   // first build and fill in final reduced field vector
   /////////////////////////////////////////////////////////////////

   // build field maps as needed
   Teuchos::RCP<const Map> reducedMap = gh_reducedFieldMaps_[fieldNum];
   if(gh_reducedFieldMaps_[fieldNum]==Teuchos::null) {
      reducedMap = getFieldMap(fieldNum,*gh_reducedFieldVector_);
      gh_reducedFieldMaps_[fieldNum] = reducedMap;
   }

   Teuchos::RCP<Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node> > finalReducedVec
      = Teuchos::rcp(new Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node>(reducedMap));
   for(std::size_t b=0;b<blockIds.size();b++) {
      std::string block = blockIds[b];

      // make sure field is in correct block
      if(!ugi_->fieldInBlock(fieldName,block))
         continue; 

      // extract data vector
      typename std::map<std::string,ArrayT>::const_iterator blockItr = data.find(block);
     TEST_FOR_EXCEPTION(blockItr==data.end(),std::runtime_error,
                        "ArrayToFieldVector::getDataVector: can not find block \""+block+"\".");

     const ArrayT & d = blockItr->second;
     updateGhostedDataReducedVector<ScalarT,ArrayT,LocalOrdinalT,GlobalOrdinalT,Node>(fieldName,block,*ugi_,d,*finalReducedVec); 
   }

   // build final (not reduced vector)
   /////////////////////////////////////////////

   Teuchos::RCP<const Map> map = gh_fieldMaps_[fieldNum];
   if(gh_fieldMaps_[fieldNum]==Teuchos::null) {
      map = getFieldMap(fieldNum,*gh_fieldVector_);
      gh_fieldMaps_[fieldNum] = map;
   }

   Teuchos::RCP<Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node> > finalVec
      = Teuchos::rcp(new Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node>(map));

   // do import from finalReducedVec
   Tpetra::Import<std::size_t,GlobalOrdinalT,Node> importer(reducedMap,map);
   finalVec->doImport(*finalReducedVec,importer,Tpetra::INSERT);

   return finalVec;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
template <typename ScalarT,typename ArrayT>
Teuchos::RCP<Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node> >
ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
   getDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const
{
   // if neccessary build field vector
   if(fieldVector_==Teuchos::null)
      buildFieldVector(*gh_fieldVector_);

   Teuchos::RCP<const Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node> > sourceVec
         = getGhostedDataVector<ScalarT,ArrayT>(fieldName,data);

   // use lazy construction for each field
   int fieldNum = ugi_->getFieldNum(fieldName);
   Teuchos::RCP<const Map> destMap = fieldMaps_[fieldNum];
   if(fieldMaps_[fieldNum]==Teuchos::null) {
      destMap = getFieldMap(fieldNum,*fieldVector_);
      fieldMaps_[fieldNum] = destMap;
   }

   Teuchos::RCP<Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node> > destVec
         = Teuchos::rcp(new Tpetra::Vector<ScalarT,std::size_t,GlobalOrdinalT,Node>(destMap));
   
   // do import
   Tpetra::Import<std::size_t,GlobalOrdinalT> importer(sourceVec->getMap(),destMap);
   destVec->doImport(*sourceVec,importer,Tpetra::INSERT); 

   return destVec;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>::
        buildFieldVector(const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> & source) const
{
   // build (unghosted) vector and map
   std::vector<GlobalOrdinalT> indices;
   ugi_->getOwnedIndices(indices);

   Teuchos::RCP<const Map> destMap
         = Teuchos::rcp(new Map(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(), Teuchos::arrayViewFromVector(indices),
                                Teuchos::OrdinalTraits<GlobalOrdinalT>::zero(), ugi_->getComm()));
   Teuchos::RCP<IntVector> localFieldVector = Teuchos::rcp(new IntVector(destMap));

   Tpetra::Import<std::size_t,GlobalOrdinalT> importer(source.getMap(),destMap);
   localFieldVector->doImport(source,importer,Tpetra::INSERT);

   fieldVector_ = localFieldVector;
}
                                   
} // end namspace panzer
