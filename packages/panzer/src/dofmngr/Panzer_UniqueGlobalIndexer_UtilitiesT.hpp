#include <vector>

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> >
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

   std::size_t indicesSz = 0;
   Teuchos::RCP<Map> destMap;
   {
      std::vector<GlobalOrdinalT> indices;
      ugi.getOwnedAndSharedIndices(indices);
      indicesSz = indices.size();
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
                                   
} // end namspace panzer
