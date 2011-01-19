#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class EpetraLinearObjFactory
// ************************************************************

template <typename LocalOrdinalT>
EpetraLinearObjFactory<LocalOrdinalT>::EpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                                                              const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider)
   : comm_(comm), gidProvider_(gidProvider)
{ }

template <typename LocalOrdinalT>
EpetraLinearObjFactory<LocalOrdinalT>::~EpetraLinearObjFactory()
{ }

// LinearObjectFactory functions 
/////////////////////////////////////////////////////////////////////

// build 
template <typename LocalOrdinalT>
Teuchos::RCP<Thyra::MultiVectorBase<double> > EpetraLinearObjFactory<LocalOrdinalT>::getGhostedVector() const
{
   Teuchos::RCP<const Epetra_Map> eMap = getGhostedMap(); // use ghosted map
   Teuchos::RCP<Epetra_MultiVector> eVec = Teuchos::rcp(new Epetra_Vector(*eMap));
 
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs = Thyra::create_VectorSpace(eMap);
   Teuchos::RCP<Thyra::MultiVectorBase<double> > vec = Thyra::create_MultiVector(eVec,vs,Teuchos::null);

   Teuchos::set_extra_data(eMap,"epetra_map",Teuchos::inOutArg(vec));
   return vec;
}

template <typename LocalOrdinalT>
Teuchos::RCP<Thyra::LinearOpBase<double> > EpetraLinearObjFactory<LocalOrdinalT>::getGhostedMatrix() const
{
   Teuchos::RCP<Epetra_CrsGraph> eGraph = getGhostedGraph(); // use ghosted graph
   Teuchos::RCP<Epetra_CrsMatrix> eMat = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *eGraph));

   return Thyra::nonconstEpetraLinearOp(eMat);
}

template <typename LocalOrdinalT>
Teuchos::RCP<Thyra::MultiVectorBase<double> > EpetraLinearObjFactory<LocalOrdinalT>::getVector() const
{
   Teuchos::RCP<const Epetra_Map> eMap = getMap(); // use fully distributed map
   Teuchos::RCP<Epetra_MultiVector> eVec = Teuchos::rcp(new Epetra_Vector(*eMap));
 
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs = Thyra::create_VectorSpace(eMap);
   Teuchos::RCP<Thyra::MultiVectorBase<double> > vec = Thyra::create_MultiVector(eVec,vs,Teuchos::null);

   Teuchos::set_extra_data(eMap,"epetra_map",Teuchos::inOutArg(vec));
   return vec;
}

template <typename LocalOrdinalT>
Teuchos::RCP<Thyra::LinearOpBase<double> > EpetraLinearObjFactory<LocalOrdinalT>::getMatrix() const
{
   Teuchos::RCP<Epetra_CrsGraph> eGraph = getGraph(); // use fully distributed graph
   Teuchos::RCP<Epetra_CrsMatrix> eMat = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *eGraph));

   return Thyra::nonconstEpetraLinearOp(eMat);
}


template <typename LocalOrdinalT>
void EpetraLinearObjFactory<LocalOrdinalT>::ghostToGlobalMatrix(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & ghostA, 
                                                                const Teuchos::RCP<Thyra::LinearOpBase<double> > & A) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // grab CrsMatrix objects
   RCP<Epetra_Export> exporter = getGhostedExport();
   RCP<const Epetra_CrsMatrix> eGhostA = rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*ghostA));
   RCP<Epetra_CrsMatrix> eA = rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*A));

   // do the global distribution
   eA->PutScalar(0.0);
   eA->Export(*eGhostA,*exporter,Add);
}

template <typename LocalOrdinalT>
void EpetraLinearObjFactory<LocalOrdinalT>::ghostToGlobalVector(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & ghostV, 
                                                                const Teuchos::RCP<Thyra::MultiVectorBase<double> > & V) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // grab maps for construction epetra vectors
   Teuchos::RCP<const Epetra_Map> eGhostMap = getGhostedMap();
   Teuchos::RCP<const Epetra_Map> eMap = getMap();

   // build vectors and exporter
   RCP<Epetra_Export> exporter = getGhostedExport();
   RCP<const Epetra_MultiVector> eGhostV = Thyra::get_Epetra_MultiVector(*eGhostMap,ghostV);
   RCP<Epetra_MultiVector> eV = Thyra::get_Epetra_MultiVector(*eMap,V);

   // do the global distribution
   eV->PutScalar(0.0);
   eV->Export(*eGhostV,*exporter,Add);
}

template <typename LocalOrdinalT>
void EpetraLinearObjFactory<LocalOrdinalT>::globalToGhostVector(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & V, 
                                                                const Teuchos::RCP<Thyra::MultiVectorBase<double> > & ghostV) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // grab maps for construction epetra vectors
   Teuchos::RCP<const Epetra_Map> eGhostMap = getGhostedMap();
   Teuchos::RCP<const Epetra_Map> eMap = getMap();

   // build vectors and importer
   RCP<Epetra_Import> importer = getGhostedImport();
   RCP<const Epetra_MultiVector> eV = Thyra::get_Epetra_MultiVector(*eMap,V);
   RCP<Epetra_MultiVector> eGhostV = Thyra::get_Epetra_MultiVector(*eGhostMap,ghostV);

   // do the global distribution
   eGhostV->PutScalar(0.0);
   eGhostV->Import(*eV,*importer,Insert);
}

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
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<LocalOrdinalT>::getGhostedMap() const
{
   if(ghostedMap_==Teuchos::null) ghostedMap_ = buildGhostedMap();

   return ghostedMap_;
}

// get the graph of the crs matrix
template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<LocalOrdinalT>::getGraph() const
{
   if(graph_==Teuchos::null) graph_ = buildGraph();

   return graph_;
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<LocalOrdinalT>::getGhostedGraph() const
{
   if(ghostedGraph_==Teuchos::null) ghostedGraph_ = buildGhostedGraph();

   return ghostedGraph_;
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import> EpetraLinearObjFactory<LocalOrdinalT>::getGhostedImport() const
{
   return Teuchos::rcp(new Epetra_Import(*getGhostedMap(),*getMap()));
}

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export> EpetraLinearObjFactory<LocalOrdinalT>::getGhostedExport() const
{
   return Teuchos::rcp(new Epetra_Export(*getGhostedMap(),*getMap()));
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

// build the ghosted map
template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<LocalOrdinalT>::buildGhostedMap() const
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

template <typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<LocalOrdinalT>::buildGhostedGraph() const
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

}
