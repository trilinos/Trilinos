#include "Panzer_UniqueGlobalIndexer.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class EpetraLinearObjFactory
// ************************************************************

template <typename LocalOrdinalT>
EpetraLinearObjFactory<LocalOrdinalT>::EpetraLinearObjFactory(const Teuchos::RCP<Epetra_Comm> & comm,
                                                              const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider)
   : comm_(comm), gidProvider_(gidProvider)
{ }

template <typename LocalOrdinalT>
EpetraLinearObjFactory<LocalOrdinalT>::~EpetraLinearObjFactory()
{ }

// "Get" functions
/////////////////////////////////////////////////////////////////////

// get the map from the matrix
template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<LocalOrdinalT>::getMap() const
{
   if(map_==Teuchos::null) map_ = buildMap();

   return map_;
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<LocalOrdinalT>::getOverlapMap() const
{
   if(overlappedMap_==Teuchos::null) overlappedMap_ = buildOverlapMap();

   return overlappedMap_;
}

// get the graph of the crs matrix
template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<LocalOrdinalT>::getGraph() const
{
   if(graph_==Teuchos::null) graph_ = buildGraph();

   return graph_;
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<LocalOrdinalT>::getOverlapGraph() const
{
   if(overlappedGraph_==Teuchos::null) overlappedGraph_ = buildOverlapGraph();

   return overlappedGraph_;
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import> EpetraLinearObjFactory<LocalOrdinalT>::getOverlapImport() const
{
   return Teuchos::rcp(new Epetra_Import(*getOverlapMap(),*getMap()));
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export> EpetraLinearObjFactory<LocalOrdinalT>::getOverlapExport() const
{
   return Teuchos::rcp(new Epetra_Export(*getOverlapMap(),*getMap()));
}

// "Build" functions
/////////////////////////////////////////////////////////////////////

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<LocalOrdinalT>::buildMap() const
{
   Teuchos::RCP<Epetra_Map> map; // result
   std::vector<int> indices;

   // get the global indices
   gidProvider_->getOwnedIndices(indices);

   map = Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));

   return map;
}

// build the overlapped map
template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<LocalOrdinalT>::buildOverlapMap() const
{
   std::vector<int> indices;

   // get the global indices
   gidProvider_->getOwnedAndSharedIndices(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));
}

// get the graph of the crs matrix
template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<LocalOrdinalT>::buildGraph() const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the overlapped graph
   RCP<Epetra_Map> map = getMap();
   RCP<Epetra_CrsGraph> graph  = rcp(new Epetra_CrsGraph(Copy,*map,0));
   RCP<Epetra_CrsGraph> oGraph = getOverlapGraph();

   // perform the communication to finish building graph
   RCP<Epetra_Export> exporter = getOverlapExport();
   graph->Export( *oGraph, *exporter, Insert );
   graph->FillComplete();
   return graph;
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<LocalOrdinalT>::buildOverlapGraph() const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> map = getOverlapMap();
   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*map,0));

   std::vector<std::string> elementBlockIds;
   
   gidProvider_->getElementBlockIds(elementBlockIds);

   // graph information about the mesh
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;

      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = gidProvider_->getElementBlock(blockId);

      // get information about number of indicies
      std::vector<int> gids;

      // loop over the elemnts
      for(std::size_t i=0;i<elements.size();i++) {

         gidProvider_->getElementGIDs(elements[i],gids);
         for(std::size_t j=0;j<gids.size();j++)
            graph->InsertGlobalIndices(gids[j],gids.size(),&gids[0]);
      }
   }

   // finish filling the graph
   graph->FillComplete();

   return graph;
}

}
