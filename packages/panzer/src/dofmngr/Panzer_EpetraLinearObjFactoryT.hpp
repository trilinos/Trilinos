#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class EpetraLinearObjFactory
// ************************************************************

template <typename Traits,typename LocalOrdinalT>
EpetraLinearObjFactory<Traits,LocalOrdinalT>::EpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                                                              const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider)
   : comm_(comm), gidProvider_(gidProvider)
{ 
   // build and register the gather/scatter evaluators with 
   // the base class.
   buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename LocalOrdinalT>
EpetraLinearObjFactory<Traits,LocalOrdinalT>::~EpetraLinearObjFactory()
{ }

// LinearObjectFactory functions 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildLinearObjContainer() const
{
   Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer);
   container->x    = getEpetraVector();
   container->dxdt = getEpetraVector();
   container->f    = getEpetraVector();
   container->A    = getEpetraMatrix();

   return container;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedLinearObjContainer() const
{
   Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer);
   container->x    = getGhostedEpetraVector();
   container->dxdt = getGhostedEpetraVector();
   container->f    = getGhostedEpetraVector();
   container->A    = getGhostedEpetraMatrix();

   return container;
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::globalToGhostContainer(const LinearObjContainer & in,
                                                                          LinearObjContainer & out) const
{
  using Teuchos::is_null;

  const EpetraLinearObjContainer & e_in = Teuchos::dyn_cast<const EpetraLinearObjContainer>(in); 
  EpetraLinearObjContainer & e_out = Teuchos::dyn_cast<EpetraLinearObjContainer>(out); 
  
  // Operations occur if the GLOBAL container has the correct targets!
  // Users set the GLOBAL continer arguments
  if ( !is_null(e_in.x) )
    globalToGhostEpetraVector(*e_in.x,*e_out.x);
  
  if ( !is_null(e_in.dxdt) )
    globalToGhostEpetraVector(*e_in.dxdt,*e_out.dxdt);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalContainer(const LinearObjContainer & in,
                                                                          LinearObjContainer & out) const
{
   const EpetraLinearObjContainer & e_in = Teuchos::dyn_cast<const EpetraLinearObjContainer>(in); 
   EpetraLinearObjContainer & e_out = Teuchos::dyn_cast<EpetraLinearObjContainer>(out); 

  // Operations occur if the GLOBAL container has the correct targets!
  // Users set the GLOBAL continer arguments
   if ( !is_null(e_out.f) )
     ghostToGlobalEpetraVector(*e_in.f,*e_out.f);

   if ( !is_null(e_out.A) )
     ghostToGlobalEpetraMatrix(*e_in.A,*e_out.A);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalEpetraVector(const Epetra_Vector in,Epetra_Vector & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = getGhostedExport();
   out.PutScalar(0.0);
   out.Export(in,*exporter,Add);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalEpetraMatrix(const Epetra_CrsMatrix in,Epetra_CrsMatrix & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = getGhostedExport();
   out.PutScalar(0.0);
   out.Export(in,*exporter,Add);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::globalToGhostEpetraVector(const Epetra_Vector in,Epetra_Vector & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Import> importer = getGhostedImport();
   out.PutScalar(0.0);
   out.Import(in,*importer,Insert);
}

// "Get" functions
/////////////////////////////////////////////////////////////////////

// get the map from the matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getMap() const
{
   if(map_==Teuchos::null) map_ = buildMap();

   return map_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedMap() const
{
   if(ghostedMap_==Teuchos::null) ghostedMap_ = buildGhostedMap();

   return ghostedMap_;
}

// get the graph of the crs matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGraph() const
{
   if(graph_==Teuchos::null) graph_ = buildGraph();

   return graph_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedGraph() const
{
   if(ghostedGraph_==Teuchos::null) ghostedGraph_ = buildGhostedGraph();

   return ghostedGraph_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedImport() const
{
   return Teuchos::rcp(new Epetra_Import(*getGhostedMap(),*getMap()));
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedExport() const
{
   return Teuchos::rcp(new Epetra_Export(*getGhostedMap(),*getMap()));
}

// "Build" functions
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildMap() const
{
   Teuchos::RCP<Epetra_Map> map; // result
   std::vector<int> indices;

   // get the global indices
   gidProvider_->getOwnedIndices(indices);

   map = Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));

   return map;
}

// build the ghosted map
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedMap() const
{
   std::vector<int> indices;

   // get the global indices
   gidProvider_->getOwnedAndSharedIndices(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));
}

// get the graph of the crs matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGraph() const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the ghosted graph
   RCP<Epetra_Map> map = getMap();
   RCP<Epetra_CrsGraph> graph  = rcp(new Epetra_CrsGraph(Copy,*map,0));
   RCP<Epetra_CrsGraph> oGraph = getGhostedGraph();

   // perform the communication to finish building graph
   RCP<Epetra_Export> exporter = getGhostedExport();
   graph->Export( *oGraph, *exporter, Insert );
   graph->FillComplete();
   return graph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedGraph() const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> map = getGhostedMap();
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

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_Vector> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedEpetraVector() const
{
   Teuchos::RCP<const Epetra_Map> eMap = getGhostedMap(); 
   return Teuchos::rcp(new Epetra_Vector(*eMap));
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_Vector> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getEpetraVector() const
{
   Teuchos::RCP<const Epetra_Map> eMap = getMap(); 
   return Teuchos::rcp(new Epetra_Vector(*eMap));
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_CrsMatrix> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getEpetraMatrix() const
{
   Teuchos::RCP<Epetra_CrsGraph> eGraph = getGraph();
   return Teuchos::rcp(new Epetra_CrsMatrix(Copy, *eGraph));
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_CrsMatrix> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedEpetraMatrix() const
{
   Teuchos::RCP<Epetra_CrsGraph> eGraph = getGhostedGraph(); 
   return Teuchos::rcp(new Epetra_CrsMatrix(Copy, *eGraph));
}

}
