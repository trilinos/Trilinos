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

#ifndef PANZER_EPETRA_LINEAR_OBJ_FACTORY_IMPL_HPP
#define PANZER_EPETRA_LINEAR_OBJ_FACTORY_IMPL_HPP

#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class EpetraLinearObjFactory
// ************************************************************

template <typename Traits,typename LocalOrdinalT>
EpetraLinearObjFactory<Traits,LocalOrdinalT>::EpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                                                                     const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider,
                                                                     bool useDiscreteAdjoint)
   : comm_(comm), gidProvider_(gidProvider), useDiscreteAdjoint_(useDiscreteAdjoint)
{ 
   hasColProvider_ = colGidProvider_!=Teuchos::null;

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename LocalOrdinalT>
EpetraLinearObjFactory<Traits,LocalOrdinalT>::EpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                                                     const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider,
                                                                     bool useDiscreteAdjoint)
   : comm_(Teuchos::null), gidProvider_(gidProvider), rawMpiComm_(comm->getRawMpiComm()), useDiscreteAdjoint_(useDiscreteAdjoint)
{ 
   comm_ = Teuchos::rcp(new Epetra_MpiComm(*rawMpiComm_));
   hasColProvider_ = colGidProvider_!=Teuchos::null;

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename LocalOrdinalT>
EpetraLinearObjFactory<Traits,LocalOrdinalT>::EpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                                                     const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider,
                                                                     const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & colGidProvider,
                                                                     bool useDiscreteAdjoint)
   : comm_(Teuchos::null), gidProvider_(gidProvider), colGidProvider_(colGidProvider), rawMpiComm_(comm->getRawMpiComm()), useDiscreteAdjoint_(useDiscreteAdjoint)
{ 
   comm_ = Teuchos::rcp(new Epetra_MpiComm(*rawMpiComm_));
   hasColProvider_ = colGidProvider_!=Teuchos::null;

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename LocalOrdinalT>
EpetraLinearObjFactory<Traits,LocalOrdinalT>::~EpetraLinearObjFactory()
{ }

// LinearObjectFactory functions 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildLinearObjContainer() const
{
   Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer(getMap(),getMap()));

   return container;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedLinearObjContainer() const
{
   Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer(getGhostedMap(),getGhostedMap()));

   return container;
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::globalToGhostContainer(const LinearObjContainer & in,
                                                                          LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;
   const EpetraLinearObjContainer & e_in = Teuchos::dyn_cast<const EpetraLinearObjContainer>(in); 
   EpetraLinearObjContainer & e_out = Teuchos::dyn_cast<EpetraLinearObjContainer>(out); 
  
   // Operations occur if the GLOBAL container has the correct targets!
   // Users set the GLOBAL continer arguments
   if ( !is_null(e_in.get_x()) && !is_null(e_out.get_x()) && ((mem & LOC::X)==LOC::X))
     globalToGhostEpetraVector(*e_in.get_x(),*e_out.get_x());
  
   if ( !is_null(e_in.get_dxdt()) && !is_null(e_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
     globalToGhostEpetraVector(*e_in.get_dxdt(),*e_out.get_dxdt());

   if ( !is_null(e_in.get_f()) && !is_null(e_out.get_f()) && ((mem & LOC::F)==LOC::F))
      globalToGhostEpetraVector(*e_in.get_f(),*e_out.get_f());
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalContainer(const LinearObjContainer & in,
                                                                          LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;
   const EpetraLinearObjContainer & e_in = Teuchos::dyn_cast<const EpetraLinearObjContainer>(in); 
   EpetraLinearObjContainer & e_out = Teuchos::dyn_cast<EpetraLinearObjContainer>(out); 

  // Operations occur if the GLOBAL container has the correct targets!
  // Users set the GLOBAL continer arguments
   if ( !is_null(e_in.get_x()) && !is_null(e_out.get_x()) && ((mem & LOC::X)==LOC::X))
     ghostToGlobalEpetraVector(*e_in.get_x(),*e_out.get_x());

   if ( !is_null(e_in.get_f()) && !is_null(e_out.get_f()) && ((mem & LOC::F)==LOC::F))
     ghostToGlobalEpetraVector(*e_in.get_f(),*e_out.get_f());

   if ( !is_null(e_in.get_A()) && !is_null(e_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
     ghostToGlobalEpetraMatrix(*e_in.get_A(),*e_out.get_A());
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalEpetraVector(const Epetra_Vector & in,Epetra_Vector & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = getGhostedExport();
   out.PutScalar(0.0);
   out.Export(in,*exporter,Add);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalEpetraMatrix(const Epetra_CrsMatrix & in,Epetra_CrsMatrix & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = getGhostedExport();
   out.PutScalar(0.0);
   out.Export(in,*exporter,Add);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::globalToGhostEpetraVector(const Epetra_Vector & in,Epetra_Vector & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Import> importer = getGhostedImport();
   out.PutScalar(0.0);
   out.Import(in,*importer,Insert);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::
adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                             const LinearObjContainer & globalBCRows,
                             LinearObjContainer & ghostedObjs) const
{
   TEUCHOS_ASSERT(!hasColProvider_); // not implemented

   const EpetraLinearObjContainer & e_localBCRows = Teuchos::dyn_cast<const EpetraLinearObjContainer>(localBCRows); 
   const EpetraLinearObjContainer & e_globalBCRows = Teuchos::dyn_cast<const EpetraLinearObjContainer>(globalBCRows); 
   EpetraLinearObjContainer & e_ghosted = Teuchos::dyn_cast<EpetraLinearObjContainer>(ghostedObjs); 

   TEUCHOS_ASSERT(!Teuchos::is_null(e_localBCRows.get_x()));
   TEUCHOS_ASSERT(!Teuchos::is_null(e_globalBCRows.get_x()));
   
   // pull out jacobian and vector
   Teuchos::RCP<Epetra_CrsMatrix> A = e_ghosted.get_A();
   Teuchos::RCP<Epetra_Vector> f = e_ghosted.get_f();

   const Epetra_Vector & local_bcs  = *(e_localBCRows.get_x());
   const Epetra_Vector & global_bcs = *(e_globalBCRows.get_x());

   TEUCHOS_ASSERT(local_bcs.MyLength()==global_bcs.MyLength());
   for(int i=0;i<local_bcs.MyLength();i++) {
      if(global_bcs[i]==0.0)
         continue;

      int numEntries = 0;
      double * values = 0;
      int * indices = 0;

      if(local_bcs[i]==0.0) { 
         // this boundary condition was NOT set by this processor

         // if they exist put 0.0 in each entry
         if(!Teuchos::is_null(f))
            (*f)[i] = 0.0;
         if(!Teuchos::is_null(A)) {
            A->ExtractMyRowView(i,numEntries,values,indices);
            for(int c=0;c<numEntries;c++) 
               values[c] = 0.0;
         }
      }
      else {
         // this boundary condition was set by this processor

         double scaleFactor = global_bcs[i];

         // if they exist scale linear objects by scale factor
         if(!Teuchos::is_null(f))
            (*f)[i] /= scaleFactor;
         if(!Teuchos::is_null(A)) {
            A->ExtractMyRowView(i,numEntries,values,indices);
            for(int c=0;c<numEntries;c++) 
               values[c] /= scaleFactor;
         }
      }
   }
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::MpiComm<int> EpetraLinearObjFactory<Traits,LocalOrdinalT>::
getComm() const
{
   return Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(dynamic_cast<const Epetra_MpiComm &>(*getEpetraComm()).Comm()));
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > 
EpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraDomainSpace() const
{
   if(domainSpace_==Teuchos::null);
      domainSpace_ = Thyra::create_VectorSpace(getMap());

   return domainSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > 
EpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraRangeSpace() const
{
   if(rangeSpace_==Teuchos::null);
      rangeSpace_ = Thyra::create_VectorSpace(getMap());

   return rangeSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::LinearOpBase<double> > 
EpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraMatrix() const
{
   return Thyra::nonconstEpetraLinearOp(getEpetraMatrix());
}

// Functions for initalizing a container
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::initializeContainer(int mem,LinearObjContainer & loc) const
{
   EpetraLinearObjContainer & eloc = Teuchos::dyn_cast<EpetraLinearObjContainer>(loc);
   initializeContainer(mem,eloc);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::initializeContainer(int mem,EpetraLinearObjContainer & loc) const
{
   typedef EpetraLinearObjContainer ELOC;

   loc.clear();

   if((mem & ELOC::X) == ELOC::X)
      loc.set_x(getEpetraVector());

   if((mem & ELOC::DxDt) == ELOC::DxDt)
      loc.set_dxdt(getEpetraVector());
    
   if((mem & ELOC::F) == ELOC::F)
      loc.set_f(getEpetraVector());

   if((mem & ELOC::Mat) == ELOC::Mat)
      loc.set_A(getEpetraMatrix());
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::initializeGhostedContainer(int mem,LinearObjContainer & loc) const
{
   EpetraLinearObjContainer & eloc = Teuchos::dyn_cast<EpetraLinearObjContainer>(loc);
   initializeGhostedContainer(mem,eloc);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::initializeGhostedContainer(int mem,EpetraLinearObjContainer & loc) const
{
   typedef EpetraLinearObjContainer ELOC;

   loc.clear();

   if((mem & ELOC::X) == ELOC::X)
      loc.set_x(getGhostedEpetraVector());

   if((mem & ELOC::DxDt) == ELOC::DxDt)
      loc.set_dxdt(getGhostedEpetraVector());
    
   if((mem & ELOC::F) == ELOC::F)
      loc.set_f(getGhostedEpetraVector());

   if((mem & ELOC::Mat) == ELOC::Mat)
      loc.set_A(getGhostedEpetraMatrix());
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

// get the map from the matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getColMap() const
{
   if(cMap_==Teuchos::null) cMap_ = buildColMap();

   return cMap_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedMap() const
{
   if(ghostedMap_==Teuchos::null) ghostedMap_ = buildGhostedMap();

   return ghostedMap_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedColMap() const
{
   if(cGhostedMap_==Teuchos::null) cGhostedMap_ = buildGhostedColMap();

   return cGhostedMap_;
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
   if(importer_==Teuchos::null)
      importer_ = Teuchos::rcp(new Epetra_Import(*getGhostedMap(),*getMap()));

   return importer_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedExport() const
{
   if(exporter_==Teuchos::null)
      exporter_ = Teuchos::rcp(new Epetra_Export(*getGhostedMap(),*getMap()));

   return exporter_;
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

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildColMap() const
{
   if(!hasColProvider_)  
     return buildMap();

   std::vector<int> indices;

   // get the global indices
   colGidProvider_->getOwnedIndices(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));
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

// build the ghosted map
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedColMap() const
{
   if(!hasColProvider_)  
     return buildGhostedMap();

   std::vector<int> indices;

   // get the global indices
   colGidProvider_->getOwnedAndSharedIndices(indices);

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
   RCP<Epetra_Map> rMap = getMap();
   RCP<Epetra_Map> cMap = getColMap();
   RCP<Epetra_CrsGraph> graph  = rcp(new Epetra_CrsGraph(Copy,*rMap,0));
   RCP<Epetra_CrsGraph> oGraph = getGhostedGraph();

   // perform the communication to finish building graph
   RCP<Epetra_Export> exporter = getGhostedExport();
   graph->Export( *oGraph, *exporter, Insert );
   graph->FillComplete(*cMap,*rMap);
   graph->OptimizeStorage();
   return graph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedGraph() const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> rMap = getGhostedMap();
   Teuchos::RCP<Epetra_Map> cMap = getGhostedColMap();
   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rMap,0));

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
      std::vector<int> col_gids;

      // loop over the elemnts
      for(std::size_t i=0;i<elements.size();i++) {

         gidProvider_->getElementGIDs(elements[i],gids);
         if(hasColProvider_)
           colGidProvider_->getElementGIDs(elements[i],col_gids);
         else
           col_gids = gids;

         for(std::size_t j=0;j<gids.size();j++)
            graph->InsertGlobalIndices(gids[j],col_gids.size(),&col_gids[0]);
      }
   }

   // finish filling the graph
   graph->FillComplete(*cMap,*rMap);
   graph->OptimizeStorage();

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

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<const Epetra_Comm> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getEpetraComm() const
{
   return comm_;
}

}

#endif
