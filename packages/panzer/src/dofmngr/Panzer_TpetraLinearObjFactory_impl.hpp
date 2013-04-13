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

#ifndef PANZER_TPETRA_LINEAR_OBJ_FACTORY_IMPL_HPP
#define PANZER_TPETRA_LINEAR_OBJ_FACTORY_IMPL_HPP

#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Thyra_TpetraVectorSpace.hpp"
#include "Thyra_TpetraLinearOp.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class TpetraLinearObjFactory
// ************************************************************

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
TpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                       const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & gidProvider)
   : comm_(comm), gidProvider_(gidProvider)
{ 
   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
~TpetraLinearObjFactory()
{ }

// LinearObjectFactory functions 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<LinearObjContainer> 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildLinearObjContainer() const
{
   Teuchos::RCP<ContainerType> container = Teuchos::rcp(new ContainerType(getMap(),getMap()));

   return container;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<LinearObjContainer> 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildGhostedLinearObjContainer() const
{
   Teuchos::RCP<ContainerType> container = Teuchos::rcp(new ContainerType(getGhostedMap(),getGhostedMap()));

   return container;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
globalToGhostContainer(const LinearObjContainer & in,
                       LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;
   typedef LinearObjContainer LOC;

   const ContainerType & t_in = Teuchos::dyn_cast<const ContainerType>(in); 
   ContainerType & t_out = Teuchos::dyn_cast<ContainerType>(out); 
  
   // Operations occur if the GLOBAL container has the correct targets!
   // Users set the GLOBAL continer arguments
   if ( !is_null(t_in.get_x()) && !is_null(t_out.get_x()) && ((mem & LOC::X)==LOC::X))
     globalToGhostTpetraVector(*t_in.get_x(),*t_out.get_x());
  
   if ( !is_null(t_in.get_dxdt()) && !is_null(t_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
     globalToGhostTpetraVector(*t_in.get_dxdt(),*t_out.get_dxdt());

   if ( !is_null(t_in.get_f()) && !is_null(t_out.get_f()) && ((mem & LOC::F)==LOC::F))
      globalToGhostTpetraVector(*t_in.get_f(),*t_out.get_f());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalContainer(const LinearObjContainer & in,
                       LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;

   const ContainerType & t_in = Teuchos::dyn_cast<const ContainerType>(in); 
   ContainerType & t_out = Teuchos::dyn_cast<ContainerType>(out); 

  // Operations occur if the GLOBAL container has the correct targets!
  // Users set the GLOBAL continer arguments
   if ( !is_null(t_in.get_x()) && !is_null(t_out.get_x()) && ((mem & LOC::X)==LOC::X))
     ghostToGlobalTpetraVector(*t_in.get_x(),*t_out.get_x());

   if ( !is_null(t_in.get_f()) && !is_null(t_out.get_f()) && ((mem & LOC::F)==LOC::F))
     ghostToGlobalTpetraVector(*t_in.get_f(),*t_out.get_f());

   if ( !is_null(t_in.get_A()) && !is_null(t_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
     ghostToGlobalTpetraMatrix(*t_in.get_A(),*t_out.get_A());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                          Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<ExportType> exporter = getGhostedExport();
   out.putScalar(0.0);
   out.doExport(in,*exporter,Tpetra::ADD);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalTpetraMatrix(const Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                          Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<ExportType> exporter = getGhostedExport();
   
   out.resumeFill();
   out.setAllToScalar(0.0);
   out.doExport(in,*exporter,Tpetra::ADD);
   out.fillComplete();
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
globalToGhostTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                          Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<ImportType> importer = getGhostedImport();
   out.putScalar(0.0);
   out.doImport(in,*importer,Tpetra::INSERT);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                             const LinearObjContainer & globalBCRows,
                             LinearObjContainer & ghostedObjs) const
{
   typedef Teuchos::ArrayRCP<const double>::Ordinal Ordinal;

   const ContainerType & t_localBCRows = Teuchos::dyn_cast<const ContainerType>(localBCRows); 
   const ContainerType & t_globalBCRows = Teuchos::dyn_cast<const ContainerType>(globalBCRows); 
   ContainerType & t_ghosted = Teuchos::dyn_cast<ContainerType>(ghostedObjs); 

   TEUCHOS_ASSERT(!Teuchos::is_null(t_localBCRows.get_x()));
   TEUCHOS_ASSERT(!Teuchos::is_null(t_globalBCRows.get_x()));
   
   // pull out jacobian and vector
   Teuchos::RCP<CrsMatrixType> A = t_ghosted.get_A();
   Teuchos::RCP<VectorType> f = t_ghosted.get_f();
   Teuchos::ArrayRCP<double> f_array = f!=Teuchos::null ? f->get1dViewNonConst() : Teuchos::null;

   const VectorType & local_bcs  = *(t_localBCRows.get_x());
   const VectorType & global_bcs = *(t_globalBCRows.get_x());
   Teuchos::ArrayRCP<const double> local_bcs_array = local_bcs.get1dView();
   Teuchos::ArrayRCP<const double> global_bcs_array = global_bcs.get1dView();

   TEUCHOS_ASSERT(local_bcs_array.size()==global_bcs_array.size());
   for(Ordinal i=0;i<local_bcs_array.size();i++) {
      if(global_bcs_array[i]==0.0)
         continue;

      std::size_t numEntries = 0;
      std::size_t sz = A->getNumEntriesInLocalRow(i);
      Teuchos::Array<LocalOrdinalT> indices(sz);
      Teuchos::Array<double> values(sz);

      if(local_bcs_array[i]==0.0) { 
         // this boundary condition was NOT set by this processor

         // if they exist put 0.0 in each entry
         if(!Teuchos::is_null(f))
            f_array[i] = 0.0;
         if(!Teuchos::is_null(A)) {
            A->getLocalRowCopy(i,indices,values,numEntries);

            for(std::size_t c=0;c<numEntries;c++) 
               values[c] = 0.0;

            A->replaceLocalValues(i,indices,values);
         }
      }
      else {
         // this boundary condition was set by this processor

         double scaleFactor = global_bcs_array[i];

         // if they exist scale linear objects by scale factor
         if(!Teuchos::is_null(f))
            f_array[i] /= scaleFactor;
         if(!Teuchos::is_null(A)) {
            A->getLocalRowCopy(i,indices,values,numEntries);

            for(std::size_t c=0;c<numEntries;c++) 
               values[c] /= scaleFactor;

            A->replaceLocalValues(i,indices,values);
         }
      }
   }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::MpiComm<int> TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getComm() const
{
   return *Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(getTeuchosComm());
}

//! Get the domain space
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraDomainSpace() const
{
   if(domainSpace_==Teuchos::null);
      domainSpace_ = Thyra::tpetraVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getMap());

   return domainSpace_;
}

//! Get the range space
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraRangeSpace() const
{
   if(rangeSpace_==Teuchos::null);
      rangeSpace_ = Thyra::tpetraVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getMap());

   return rangeSpace_;
}

//! Get a matrix operator
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraMatrix() const
{
   return Thyra::tpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getThyraRangeSpace(),getThyraDomainSpace(),getTpetraMatrix());
}

// Functions for initalizing a container
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeContainer(int mem,LinearObjContainer & loc) const
{
   ContainerType & tloc = Teuchos::dyn_cast<ContainerType>(loc);
   initializeContainer(mem,tloc);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeContainer(int mem,TpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & loc) const
{
   typedef LinearObjContainer LOC;

   loc.clear();

   if((mem & LOC::X) == LOC::X)
      loc.set_x(getTpetraVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt(getTpetraVector());
    
   if((mem & LOC::F) == LOC::F)
      loc.set_f(getTpetraVector());

   if((mem & LOC::Mat) == LOC::Mat)
      loc.set_A(getTpetraMatrix());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeGhostedContainer(int mem,LinearObjContainer & loc) const
{
   ContainerType & tloc = Teuchos::dyn_cast<ContainerType>(loc);
   initializeGhostedContainer(mem,tloc);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeGhostedContainer(int mem,TpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & loc) const
{
   typedef LinearObjContainer LOC;

   loc.clear();

   if((mem & LOC::X) == LOC::X)
      loc.set_x(getGhostedTpetraVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt(getGhostedTpetraVector());
    
   if((mem & LOC::F) == LOC::F)
      loc.set_f(getGhostedTpetraVector());

   if((mem & LOC::Mat) == LOC::Mat)
      loc.set_A(getGhostedTpetraMatrix());
}

// "Get" functions
/////////////////////////////////////////////////////////////////////

// get the map from the matrix
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getMap() const
{
   if(map_==Teuchos::null) map_ = buildMap();

   return map_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedMap() const
{
   if(ghostedMap_==Teuchos::null) ghostedMap_ = buildGhostedMap();

   return ghostedMap_;
}

// get the graph of the crs matrix
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGraph() const
{
   if(graph_==Teuchos::null) graph_ = buildGraph();

   return graph_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedGraph() const
{
   if(ghostedGraph_==Teuchos::null) ghostedGraph_ = buildGhostedGraph();

   return ghostedGraph_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedImport() const
{
   return Teuchos::rcp(new ImportType(getMap(),getGhostedMap()));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedExport() const
{
   return Teuchos::rcp(new ExportType(getGhostedMap(),getMap()));
}

// "Build" functions
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildMap() const
{
   std::vector<GlobalOrdinalT> indices;

   // get the global indices
   gidProvider_->getOwnedIndices(indices);

   return Teuchos::rcp(new MapType(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(),indices,0,comm_));
}

// build the ghosted map
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildGhostedMap() const
{
   std::vector<GlobalOrdinalT> indices;

   // get the global indices
   gidProvider_->getOwnedAndSharedIndices(indices);

   return Teuchos::rcp(new MapType(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(),indices,0,comm_));
}

// get the graph of the crs matrix
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildGraph() const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the ghosted graph
   RCP<MapType> map = getMap();
   RCP<CrsGraphType> graph  = rcp(new CrsGraphType(map,0));
   RCP<CrsGraphType> oGraph = getGhostedGraph();

   // perform the communication to finish building graph
   RCP<ExportType> exporter = getGhostedExport();
   graph->doExport( *oGraph, *exporter, Tpetra::INSERT );
   graph->fillComplete();

   return graph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildGhostedGraph() const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<MapType> map = getGhostedMap();
   Teuchos::RCP<CrsGraphType> graph = Teuchos::rcp(new CrsGraphType(map,0));

   std::vector<std::string> elementBlockIds;
   
   gidProvider_->getElementBlockIds(elementBlockIds);

   // graph information about the mesh
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;

      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = gidProvider_->getElementBlock(blockId);

      // get information about number of indicies
      std::vector<GlobalOrdinalT> gids;

      // loop over the elemnts
      for(std::size_t i=0;i<elements.size();i++) {

         gidProvider_->getElementGIDs(elements[i],gids);
         for(std::size_t j=0;j<gids.size();j++)
            graph->insertGlobalIndices(gids[j],gids);
      }
   }

   // finish filling the graph
   graph->fillComplete();

   return graph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedTpetraVector() const
{
   Teuchos::RCP<const MapType> tMap = getGhostedMap(); 
   return Teuchos::rcp(new VectorType(tMap));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraVector() const
{
   Teuchos::RCP<const MapType> tMap = getMap(); 
   return Teuchos::rcp(new VectorType(tMap));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraMatrix() const
{
   Teuchos::RCP<CrsGraphType> tGraph = getGraph();
   Teuchos::RCP<CrsMatrixType> tMat =  Teuchos::rcp(new CrsMatrixType(tGraph));
   tMat->fillComplete();

   return tMat;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedTpetraMatrix() const
{
   Teuchos::RCP<CrsGraphType> tGraph = getGhostedGraph(); 
   Teuchos::RCP<CrsMatrixType> tMat =  Teuchos::rcp(new CrsMatrixType(tGraph));
   tMat->fillComplete();

   return tMat;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<const Teuchos::Comm<int> > 
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTeuchosComm() const
{
   return comm_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
beginFill(LinearObjContainer & loc) const
{
  ContainerType & tloc = Teuchos::dyn_cast<ContainerType>(loc);
  Teuchos::RCP<CrsMatrixType> A = tloc.get_A();
  if(A!=Teuchos::null) 
    A->resumeFill();
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
endFill(LinearObjContainer & loc) const
{
  ContainerType & tloc = Teuchos::dyn_cast<ContainerType>(loc);
  Teuchos::RCP<CrsMatrixType> A = tloc.get_A();
  if(A!=Teuchos::null) 
    A->fillComplete();
}

}

#endif
