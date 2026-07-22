// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_TpetraLinearObjFactory_impl_hpp__
#define   __Panzer_TpetraLinearObjFactory_impl_hpp__

// Panzer
#include "Panzer_ConnManager.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_EpetraVector_Write_GlobalEvaluationData.hpp"                    // JMG:  Remove this eventually.
#endif
#include "Panzer_TpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GlobalIndexer.hpp"

// Thyra
#include "Thyra_TpetraVectorSpace.hpp"
#include "Thyra_TpetraLinearOp.hpp"

// Tpetra
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace panzer {

using Teuchos::RCP;

// ************************************************************
// class TpetraLinearObjFactory
// ************************************************************

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
TpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                       const Teuchos::RCP<const GlobalIndexer> & gidProvider)
   : comm_(comm), gidProvider_(gidProvider)
{
   hasColProvider_ = colGidProvider_!=Teuchos::null;

   // build and register the gather/scatter evaluators with
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
TpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                       const Teuchos::RCP<const GlobalIndexer> & gidProvider,
                       const Teuchos::RCP<const GlobalIndexer> & colGidProvider)
   : comm_(comm), gidProvider_(gidProvider), colGidProvider_(colGidProvider)
{
   hasColProvider_ = colGidProvider_!=Teuchos::null;

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
   Teuchos::RCP<ContainerType> container = Teuchos::rcp(new ContainerType(getColMap(),getMap()));

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
     globalToGhostTpetraVector(*t_in.get_x(),*t_out.get_x(),true);

   if ( !is_null(t_in.get_dxdt()) && !is_null(t_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
     globalToGhostTpetraVector(*t_in.get_dxdt(),*t_out.get_dxdt(),true);

   if ( !is_null(t_in.get_f()) && !is_null(t_out.get_f()) && ((mem & LOC::F)==LOC::F))
      globalToGhostTpetraVector(*t_in.get_f(),*t_out.get_f(),false);
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
     ghostToGlobalTpetraVector(*t_in.get_x(),*t_out.get_x(),true);

   if ( !is_null(t_in.get_f()) && !is_null(t_out.get_f()) && ((mem & LOC::F)==LOC::F))
     ghostToGlobalTpetraVector(*t_in.get_f(),*t_out.get_f(),false);

   if ( !is_null(t_in.get_A()) && !is_null(t_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
     ghostToGlobalTpetraMatrix(*t_in.get_A(),*t_out.get_A());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                          Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out, bool col) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<ExportType> exporter = col ? getGhostedColExport() : getGhostedExport();
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
   out.fillComplete(out.getDomainMap(),out.getRangeMap());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
globalToGhostTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                          Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out, bool col) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<ImportType> importer = col ? getGhostedColImport() : getGhostedImport();
   out.putScalar(0.0);
   out.doImport(in,*importer,Tpetra::INSERT);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                             const LinearObjContainer & globalBCRows,
                             LinearObjContainer & ghostedObjs,
                             bool zeroVectorRows, bool adjustX) const
{
   typedef Teuchos::ArrayRCP<const double>::Ordinal Ordinal;

   const ContainerType & t_localBCRows = Teuchos::dyn_cast<const ContainerType>(localBCRows);
   const ContainerType & t_globalBCRows = Teuchos::dyn_cast<const ContainerType>(globalBCRows);
   ContainerType & t_ghosted = Teuchos::dyn_cast<ContainerType>(ghostedObjs);

   TEUCHOS_ASSERT(!Teuchos::is_null(t_localBCRows.get_f()));
   TEUCHOS_ASSERT(!Teuchos::is_null(t_globalBCRows.get_f()));

   // pull out jacobian and vector
   Teuchos::RCP<CrsMatrixType> A = t_ghosted.get_A();
   Teuchos::RCP<VectorType> f = t_ghosted.get_f();
   if(adjustX) f = t_ghosted.get_x();
   Teuchos::ArrayRCP<double> f_array = f!=Teuchos::null ? f->get1dViewNonConst() : Teuchos::null;

   const VectorType & local_bcs  = *(t_localBCRows.get_f());
   const VectorType & global_bcs = *(t_globalBCRows.get_f());
   Teuchos::ArrayRCP<const double> local_bcs_array = local_bcs.get1dView();
   Teuchos::ArrayRCP<const double> global_bcs_array = global_bcs.get1dView();

   TEUCHOS_ASSERT(local_bcs_array.size()==global_bcs_array.size());
   for(Ordinal i=0;i<local_bcs_array.size();i++) {
      if(global_bcs_array[i]==0.0)
         continue;

      if(local_bcs_array[i]==0.0 || zeroVectorRows) {
         // this boundary condition was NOT set by this processor

         // if they exist put 0.0 in each entry
         if(!Teuchos::is_null(f))
            f_array[i] = 0.0;
         if(!Teuchos::is_null(A)) {
            std::size_t numEntries = 0;
            std::size_t sz = A->getNumEntriesInLocalRow(i);
	    typename CrsMatrixType::nonconst_local_inds_host_view_type indices("indices", sz);
	    typename CrsMatrixType::nonconst_values_host_view_type values("values", sz);

            A->getLocalRowCopy(i,indices,values,numEntries);

            for(std::size_t c=0;c<numEntries;c++)
	      values(c) = 0.0;

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
            std::size_t numEntries = 0;
            std::size_t sz = A->getNumEntriesInLocalRow(i);
	    typename CrsMatrixType::nonconst_local_inds_host_view_type indices("indices", sz);
	    typename CrsMatrixType::nonconst_values_host_view_type values("values", sz);

            A->getLocalRowCopy(i,indices,values,numEntries);

            for(std::size_t c=0;c<numEntries;c++)
	      values(c) /= scaleFactor;

            A->replaceLocalValues(i,indices,values);
         }
      }
   }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
applyDirichletBCs(const LinearObjContainer & /* counter */,
                  LinearObjContainer & /* result */) const
{
  TEUCHOS_ASSERT(false); // not yet implemented
}

///////////////////////////////////////////////////////////////////////////////
//
//  buildReadOnlyDomainContainer()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename ScalarT, typename LocalOrdinalT,
  typename GlobalOrdinalT, typename NodeT>
Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>
TpetraLinearObjFactory<Traits, ScalarT, LocalOrdinalT, GlobalOrdinalT, NodeT>::
buildReadOnlyDomainContainer() const
{
  using Teuchos::rcp;
  using TVROGED = TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,
    LocalOrdinalT, GlobalOrdinalT, NodeT>;
  auto ged = rcp(new TVROGED);
  ged->initialize(getGhostedImport(), getGhostedColMap(), getColMap());
  return ged;
} // end of buildReadOnlyDomainContainer()

#ifdef PANZER_HAVE_EPETRA_STACK
///////////////////////////////////////////////////////////////////////////////
//
//  buildWriteDomainContainer()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename ScalarT, typename LocalOrdinalT,
  typename GlobalOrdinalT, typename NodeT>
Teuchos::RCP<WriteVector_GlobalEvaluationData>
TpetraLinearObjFactory<Traits, ScalarT, LocalOrdinalT, GlobalOrdinalT, NodeT>::
buildWriteDomainContainer() const
{
  using std::logic_error;
  using Teuchos::rcp;
  using EVWGED = panzer::EpetraVector_Write_GlobalEvaluationData;
  auto ged = rcp(new EVWGED);
  TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "NOT IMPLEMENTED YET")
  return ged;
} // end of buildWriteDomainContainer()
#endif // PANZER_HAVE_EPETRA_STACK

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
   if(domainSpace_==Teuchos::null) {
     if(!hasColProvider_)
       domainSpace_ = Thyra::tpetraVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getMap());
     else
       domainSpace_ = Thyra::tpetraVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getColMap());
   }

   return domainSpace_;
}

//! Get the range space
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraRangeSpace() const
{
   if(rangeSpace_==Teuchos::null)
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
      loc.set_x(getTpetraColVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt(getTpetraColVector());

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
      loc.set_x(getGhostedTpetraColVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt(getGhostedTpetraColVector());

   if((mem & LOC::F) == LOC::F) {
      loc.set_f(getGhostedTpetraVector());
      loc.setRequiresDirichletAdjustment(true);
   }

   if((mem & LOC::Mat) == LOC::Mat) {
      loc.set_A(getGhostedTpetraMatrix());
      loc.setRequiresDirichletAdjustment(true);
   }
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

// get the map from the matrix
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getColMap() const
{
   if(cMap_==Teuchos::null) cMap_ = buildColMap();

   return cMap_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedMap() const
{
   if(ghostedMap_==Teuchos::null) ghostedMap_ = buildGhostedMap();

   return ghostedMap_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedColMap() const
{
   if(cGhostedMap_==Teuchos::null) cGhostedMap_ = buildGhostedColMap();

   return cGhostedMap_;
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
   if(ghostedImporter_==Teuchos::null)
      ghostedImporter_ = Teuchos::rcp(new ImportType(getMap(),getGhostedMap()));

   return ghostedImporter_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedColImport() const
{
   if(!hasColProvider_)
      ghostedColImporter_ = getGhostedImport(); // they are the same in this case

   if(ghostedColImporter_==Teuchos::null)
      ghostedColImporter_ = Teuchos::rcp(new ImportType(getColMap(),getGhostedColMap()));

   return ghostedColImporter_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedExport() const
{
   if(ghostedExporter_==Teuchos::null)
      ghostedExporter_ = Teuchos::rcp(new ExportType(getGhostedMap(),getMap()));

   return ghostedExporter_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedColExport() const
{
   if(!hasColProvider_)
      ghostedColExporter_ = getGhostedExport(); // they are the same in this case

   if(ghostedColExporter_==Teuchos::null)
      ghostedColExporter_ = Teuchos::rcp(new ExportType(getGhostedColMap(),getColMap()));

   return ghostedColExporter_;
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

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildColMap() const
{
   if(!hasColProvider_)
     return buildMap();

   std::vector<GlobalOrdinalT> indices;

   // get the global indices
   colGidProvider_->getOwnedIndices(indices);

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
   gidProvider_->getOwnedAndGhostedIndices(indices);

   return Teuchos::rcp(new MapType(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(),indices,0,comm_));
}

// build the ghosted map
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildGhostedColMap() const
{
   if(!hasColProvider_)
     return buildGhostedMap();

   std::vector<GlobalOrdinalT> indices;

   // get the global indices
   colGidProvider_->getOwnedAndGhostedIndices(indices);

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
   RCP<MapType> rMap = getMap();
   RCP<MapType> cMap = getColMap();
   RCP<CrsGraphType> graph  = rcp(new CrsGraphType(rMap,0));
   RCP<CrsGraphType> oGraph = getGhostedGraph();

   // perform the communication to finish building graph
   RCP<ExportType> exporter = getGhostedExport();
   graph->doExport( *oGraph, *exporter, Tpetra::INSERT );
   graph->fillComplete(cMap,rMap);

   return graph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildGhostedGraph() const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<MapType> rMap = getGhostedMap();
   Teuchos::RCP<MapType> cMap = getGhostedColMap();

   std::vector<std::string> elementBlockIds;
   gidProvider_->getElementBlockIds(elementBlockIds);

   const Teuchos::RCP<const GlobalIndexer>
     colGidProvider = hasColProvider_ ? colGidProvider_ : gidProvider_;
   const Teuchos::RCP<const ConnManager> conn_mgr = colGidProvider->getConnManager();
   const bool han = conn_mgr.is_null() ? false : conn_mgr->hasAssociatedNeighbors();

   // graph information about the mesh
   // Count number of entries per graph row; needed for graph constructor
   std::vector<size_t> nEntriesPerRow(rMap->getLocalNumElements(), 0);

   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;

      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = gidProvider_->getElementBlock(blockId);

      // get information about number of indicies
      std::vector<GlobalOrdinalT> gids;
      std::vector<GlobalOrdinalT> col_gids;

      // loop over the elemnts
      for(std::size_t i=0;i<elements.size();i++) {
         gidProvider_->getElementGIDs(elements[i],gids);

         colGidProvider->getElementGIDs(elements[i],col_gids);
         if (han) {
           const std::vector<LocalOrdinalT>& aes = conn_mgr->getAssociatedNeighbors(elements[i]);
           for (typename std::vector<LocalOrdinalT>::const_iterator eit = aes.begin();
                eit != aes.end(); ++eit) {
             std::vector<GlobalOrdinalT> other_col_gids;
             colGidProvider->getElementGIDs(*eit, other_col_gids);
             col_gids.insert(col_gids.end(), other_col_gids.begin(), other_col_gids.end());
           }
         }

         for(std::size_t j=0;j<gids.size();j++){
            LocalOrdinalT lid = rMap->getLocalElement(gids[j]);
            nEntriesPerRow[lid] += col_gids.size();
         }
      }
   }

   Teuchos::ArrayView<const size_t> nEntriesPerRowView(nEntriesPerRow);
   Teuchos::RCP<CrsGraphType> graph = Teuchos::rcp(new CrsGraphType(rMap,cMap,
                                                          nEntriesPerRowView));

   // Now insert entries into the graph
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;

      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = gidProvider_->getElementBlock(blockId);

      // get information about number of indicies
      std::vector<GlobalOrdinalT> gids;
      std::vector<GlobalOrdinalT> col_gids;

      // loop over the elemnts
      for(std::size_t i=0;i<elements.size();i++) {
         gidProvider_->getElementGIDs(elements[i],gids);

         colGidProvider->getElementGIDs(elements[i],col_gids);
         if (han) {
           const std::vector<LocalOrdinalT>& aes = conn_mgr->getAssociatedNeighbors(elements[i]);
           for (typename std::vector<LocalOrdinalT>::const_iterator eit = aes.begin();
                eit != aes.end(); ++eit) {
             std::vector<GlobalOrdinalT> other_col_gids;
             colGidProvider->getElementGIDs(*eit, other_col_gids);
             col_gids.insert(col_gids.end(), other_col_gids.begin(), other_col_gids.end());
           }
         }

         for(std::size_t j=0;j<gids.size();j++)
            graph->insertGlobalIndices(gids[j],col_gids);
      }
   }

   // finish filling the graph
   graph->fillComplete(cMap,rMap);

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
getGhostedTpetraColVector() const
{
   Teuchos::RCP<const MapType> tMap = getGhostedColMap();
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
Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraColVector() const
{
   Teuchos::RCP<const MapType> tMap = getColMap();
   return Teuchos::rcp(new VectorType(tMap));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraMatrix() const
{
   Teuchos::RCP<CrsGraphType> tGraph = getGraph();
   Teuchos::RCP<CrsMatrixType> tMat =  Teuchos::rcp(new CrsMatrixType(tGraph));
   tMat->fillComplete(tMat->getDomainMap(),tMat->getRangeMap());

   return tMat;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
TpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedTpetraMatrix() const
{
   Teuchos::RCP<CrsGraphType> tGraph = getGhostedGraph();
   Teuchos::RCP<CrsMatrixType> tMat =  Teuchos::rcp(new CrsMatrixType(tGraph));
   tMat->fillComplete(tMat->getDomainMap(),tMat->getRangeMap());

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
    A->fillComplete(A->getDomainMap(),A->getRangeMap());
}

}

#endif // __Panzer_TpetraLinearObjFactory_impl_hpp__
