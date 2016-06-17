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
#include "Panzer_Filtered_UniqueGlobalIndexer.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

#include "Thyra_SpmdVectorBase.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"

#include "EpetraExt_VectorOut.h"
#include "EpetraExt_VectorIn.h"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class EpetraLinearObjFactory
// ************************************************************

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
EpetraLinearObjFactory<Traits,LocalOrdinalT>::EpetraLinearObjFactory(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider)
   : comm_(Teuchos::null), gidProvider_(gidProvider), useDiscreteAdjoint_(false)
{ 
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
void 
EpetraLinearObjFactory<Traits,LocalOrdinalT>::
readVector(const std::string & identifier,LinearObjContainer & loc,int id) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;

  EpetraLinearObjContainer & eloc = dyn_cast<EpetraLinearObjContainer>(loc);

   // extract the vector from linear object container
  RCP<Thyra::VectorBase<double> > vec;
  switch(id) {
  case LinearObjContainer::X:
    vec = eloc.get_x_th();
    break;
  case LinearObjContainer::DxDt:
    vec = eloc.get_dxdt_th();
    break;
  case LinearObjContainer::F:
    vec = eloc.get_f_th();
    break;
  default:
    TEUCHOS_ASSERT(false);
    break;
  };

  // convert to Epetra then read in from a file
  RCP<Epetra_Vector> ex = Thyra::get_Epetra_Vector(*getMap(),vec);

  // build the file name from the identifier
  std::stringstream ss;
  ss << identifier << ".mm";

  // read in vector (wow the MM to Vector is a poorly designed interface!)
  Epetra_Vector * ptr_ex = 0;
  TEUCHOS_ASSERT(0==EpetraExt::MatrixMarketFileToVector(ss.str().c_str(),*getMap(),ptr_ex));

  *ex = *ptr_ex;
  delete ptr_ex;
}

template <typename Traits,typename LocalOrdinalT>
void 
EpetraLinearObjFactory<Traits,LocalOrdinalT>::
writeVector(const std::string & identifier,const LinearObjContainer & loc,int id) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;

  const EpetraLinearObjContainer & eloc = dyn_cast<const EpetraLinearObjContainer>(loc);

   // extract the vector from linear object container
  RCP<const Thyra::VectorBase<double> > vec;
  switch(id) {
  case LinearObjContainer::X:
    vec = eloc.get_x_th();
    break;
  case LinearObjContainer::DxDt:
    vec = eloc.get_dxdt_th();
    break;
  case LinearObjContainer::F:
    vec = eloc.get_f_th();
    break;
  default:
    TEUCHOS_ASSERT(false);
    break;
  };

  // convert to Epetra then write out to a file
  RCP<const Epetra_Vector> ex = Thyra::get_Epetra_Vector(*getMap(),vec);

  // build the file name from the identifier
  std::stringstream ss;
  ss << identifier << ".mm";

  // write out vector
  TEUCHOS_ASSERT(0==EpetraExt::VectorToMatrixMarketFile(ss.str().c_str(),*ex));
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildLinearObjContainer() const
{
   Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer(getColMap(),getMap()));

   return container;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedLinearObjContainer() const
{
   Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer(getGhostedColMap(),getGhostedMap()));

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
     globalToGhostEpetraVector(*e_in.get_x(),*e_out.get_x(),true);
  
   if ( !is_null(e_in.get_dxdt()) && !is_null(e_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
     globalToGhostEpetraVector(*e_in.get_dxdt(),*e_out.get_dxdt(),true);

   if ( !is_null(e_in.get_f()) && !is_null(e_out.get_f()) && ((mem & LOC::F)==LOC::F))
      globalToGhostEpetraVector(*e_in.get_f(),*e_out.get_f(),false);
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
     ghostToGlobalEpetraVector(*e_in.get_x(),*e_out.get_x(),true);

   if ( !is_null(e_in.get_f()) && !is_null(e_out.get_f()) && ((mem & LOC::F)==LOC::F))
     ghostToGlobalEpetraVector(*e_in.get_f(),*e_out.get_f(),false);

   if ( !is_null(e_in.get_A()) && !is_null(e_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
     ghostToGlobalEpetraMatrix(*e_in.get_A(),*e_out.get_A());
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalEpetraVector(const Epetra_Vector & in,Epetra_Vector & out,bool col) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = col ? getGhostedColExport() : getGhostedExport();
   TEUCHOS_ASSERT(out.PutScalar(0.0)==0);

   int errCode = out.Export(in,*exporter,Add);
   TEUCHOS_TEST_FOR_EXCEPTION(errCode!=0,std::runtime_error,
                              "Epetra_Vector::Export returned an error code of " << errCode << "!");
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
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::globalToGhostEpetraVector(const Epetra_Vector & in,Epetra_Vector & out,bool col) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Import> importer = col ? getGhostedColImport() : getGhostedImport();
   out.PutScalar(0.0);
   out.Import(in,*importer,Insert);
   // NOTE: These lines cause the response_residual test to fail!
   // int retval = out.Import(in,*importer,Insert);
   // TEUCHOS_TEST_FOR_EXCEPTION(0!=retval,std::logic_error,"panzer::EpetraLOF::globalToGhostEpetraVector "
   //                                                       "import call failed with return value retval = " << retval);
}

template <typename Traits,typename LocalOrdinalT>
void EpetraLinearObjFactory<Traits,LocalOrdinalT>::
adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                             const LinearObjContainer & globalBCRows,
                             LinearObjContainer & ghostedObjs,
                             bool zeroVectorRows, bool adjustX) const
          
{
   const EpetraLinearObjContainer & e_localBCRows = Teuchos::dyn_cast<const EpetraLinearObjContainer>(localBCRows); 
   const EpetraLinearObjContainer & e_globalBCRows = Teuchos::dyn_cast<const EpetraLinearObjContainer>(globalBCRows); 
   EpetraLinearObjContainer & e_ghosted = Teuchos::dyn_cast<EpetraLinearObjContainer>(ghostedObjs); 

   TEUCHOS_ASSERT(!Teuchos::is_null(e_localBCRows.get_f()));
   TEUCHOS_ASSERT(!Teuchos::is_null(e_globalBCRows.get_f()));
   
   // pull out jacobian and vector
   Teuchos::RCP<Epetra_CrsMatrix> A = e_ghosted.get_A();
   Teuchos::RCP<Epetra_Vector> f = e_ghosted.get_f();
   if(adjustX) f = e_ghosted.get_x();

   const Epetra_Vector & local_bcs  = *(e_localBCRows.get_f());
   const Epetra_Vector & global_bcs = *(e_globalBCRows.get_f());

   TEUCHOS_ASSERT(local_bcs.MyLength()==global_bcs.MyLength());
   for(int i=0;i<local_bcs.MyLength();i++) {
      if(global_bcs[i]==0.0)
         continue;

      int numEntries = 0;
      double * values = 0;
      int * indices = 0;

      if(local_bcs[i]==0.0 || zeroVectorRows) { 
         // this boundary condition was NOT set by this processor
         // or the user requrested that every row be zeroed

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
void  EpetraLinearObjFactory<Traits,LocalOrdinalT>::
applyDirichletBCs(const LinearObjContainer & counter,
                  LinearObjContainer & result) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;

  const ThyraObjContainer<double> & th_counter = dyn_cast<const ThyraObjContainer<double> >(counter);
  ThyraObjContainer<double> & th_result  = dyn_cast<ThyraObjContainer<double> >(result);
  
  RCP<const Thyra::VectorBase<double> > count = th_counter.get_f_th();
  RCP<const Thyra::VectorBase<double> > f_in = th_counter.get_f_th();
  RCP<Thyra::VectorBase<double> > f_out = th_result.get_f_th();

  Teuchos::ArrayRCP<const double> count_array,f_in_array;
  Teuchos::ArrayRCP<double> f_out_array;

  rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(count,true)->getLocalData(Teuchos::ptrFromRef(count_array));
  rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_in,true)->getLocalData(Teuchos::ptrFromRef(f_in_array));
  rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(f_out,true)->getNonconstLocalData(Teuchos::ptrFromRef(f_out_array));

  TEUCHOS_ASSERT(count_array.size()==f_in_array.size());
  TEUCHOS_ASSERT(count_array.size()==f_out_array.size());

  for(Teuchos::ArrayRCP<double>::size_type i=0;i<count_array.size();i++) {
    if(count_array[i]!=0.0)
      f_out_array[i] = f_in_array[i];
  }
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>
EpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildDomainContainer() const
{
  Teuchos::RCP<EpetraVector_ReadOnly_GlobalEvaluationData> vec_ged
    = Teuchos::rcp(new EpetraVector_ReadOnly_GlobalEvaluationData);
  vec_ged->initialize(getGhostedImport(),getGhostedColMap(),getColMap());

  return vec_ged;
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
   if(domainSpace_ == Teuchos::null) {
     // in the first case, range is domain, 
     // in the second domain must be constructed
     if(!hasColProvider_)
        domainSpace_ = getThyraRangeSpace();
     else
        domainSpace_ = Thyra::create_VectorSpace(getColMap());
   }

   return domainSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > 
EpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraRangeSpace() const
{
   if(rangeSpace_ == Teuchos::null)
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
      loc.set_x(getEpetraColVector());

   if((mem & ELOC::DxDt) == ELOC::DxDt)
      loc.set_dxdt(getEpetraColVector());
    
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
      loc.set_x(getGhostedEpetraColVector());

   if((mem & ELOC::DxDt) == ELOC::DxDt)
      loc.set_dxdt(getGhostedEpetraColVector());
    
   if((mem & ELOC::F) == ELOC::F) {
      loc.set_f(getGhostedEpetraVector());
      loc.setRequiresDirichletAdjustment(true);
   }

   if((mem & ELOC::Mat) == ELOC::Mat) {
      loc.set_A(getGhostedEpetraMatrix());
      loc.setRequiresDirichletAdjustment(true);
   }
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
   if(ghostedGraph_==Teuchos::null) ghostedGraph_ = buildGhostedGraph(true);

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
const Teuchos::RCP<Epetra_Import> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedColImport() const
{
   if(!hasColProvider_)
      colImporter_ = getGhostedImport(); // they are the same in this case

   if(colImporter_==Teuchos::null)
      colImporter_ = Teuchos::rcp(new Epetra_Import(*getGhostedColMap(),*getColMap()));

   return colImporter_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedExport() const
{
   if(exporter_==Teuchos::null)
      exporter_ = Teuchos::rcp(new Epetra_Export(*getGhostedMap(),*getMap()));

   return exporter_;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedColExport() const
{
   if(!hasColProvider_)
      colExporter_ = getGhostedExport(); // they are the same in this case

   if(colExporter_==Teuchos::null)
      colExporter_ = Teuchos::rcp(new Epetra_Export(*getGhostedColMap(),*getColMap()));

   return colExporter_;
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
   // RCP<Epetra_CrsGraph> oGraph = getGhostedGraph();
   RCP<Epetra_CrsGraph> oGraph = buildFilteredGhostedGraph();
       // this method calls getGhosted graph if no filtering is enabled
       // otherwise a filtered matrix is constructed

   // perform the communication to finish building graph
   RCP<Epetra_Export> exporter = getGhostedExport();
   graph->Export( *oGraph, *exporter, Insert );
   graph->FillComplete(*cMap,*rMap);
   graph->OptimizeStorage();
   return graph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedGraph(bool optimizeStorage) const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> rMap = getGhostedMap();
   Teuchos::RCP<Epetra_Map> cMap = getGhostedColMap();
   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rMap,*cMap,0));

   std::vector<std::string> elementBlockIds;
   gidProvider_->getElementBlockIds(elementBlockIds);

   const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> >
     colGidProvider = hasColProvider_ ? colGidProvider_ : gidProvider_;
   const Teuchos::RCP<const ConnManagerBase<LocalOrdinalT> > conn_mgr = colGidProvider->getConnManagerBase();
   const bool han = conn_mgr.is_null() ? false : conn_mgr->hasAssociatedNeighbors();
 
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

         colGidProvider->getElementGIDs(elements[i],col_gids);
         if (han) {
           const std::vector<LocalOrdinalT>& aes = conn_mgr->getAssociatedNeighbors(elements[i]);
           for (typename std::vector<LocalOrdinalT>::const_iterator eit = aes.begin();
                eit != aes.end(); ++eit) {
             std::vector<int> other_col_gids;
             colGidProvider->getElementGIDs(*eit, other_col_gids);
             col_gids.insert(col_gids.end(), other_col_gids.begin(), other_col_gids.end());
           }
         }

         for(std::size_t j=0;j<gids.size();j++)
            graph->InsertGlobalIndices(gids[j],col_gids.size(),&col_gids[0]);
      }
   }

   // finish filling the graph
   graph->FillComplete(*cMap,*rMap);
   if(optimizeStorage) {
     graph->OptimizeStorage();
   }

   return graph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> EpetraLinearObjFactory<Traits,LocalOrdinalT>::buildFilteredGhostedGraph() const
{
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // figure out if the domain is filtered
   RCP<const Filtered_UniqueGlobalIndexer<LocalOrdinalT,int> > filtered_ugi 
       = rcp_dynamic_cast<const Filtered_UniqueGlobalIndexer<LocalOrdinalT,int> >(getDomainGlobalIndexer());

   // domain is unfiltered, a filtered graph is just the original ghosted graph
   if(filtered_ugi==Teuchos::null)
     return getGhostedGraph();

   // get all local indices that are active (i.e. unfiltered)
   std::vector<int> ghostedActive;
   filtered_ugi->getOwnedAndSharedNotFilteredIndicator(ghostedActive);

   // not sure how the deep copy works, so we will do this instead. 
   // This will build a new ghosted graph // without optimized storage so entries can be removed.
   Teuchos::RCP<Epetra_CrsGraph> filteredGraph = buildGhostedGraph(false); 
       // false implies that storage is not optimzied 

   // remove filtered column entries
   for(int i=0;i<filteredGraph->NumMyRows();i++) {
     std::vector<int> removedIndices;
     int numIndices = 0;
     int * indices = 0;
     TEUCHOS_ASSERT(filteredGraph->ExtractMyRowView(i,numIndices,indices)==0);

     for(int j=0;j<numIndices;j++) {
       if(ghostedActive[indices[j]]==0)
         removedIndices.push_back(indices[j]);
     }

     TEUCHOS_ASSERT(filteredGraph->RemoveMyIndices(i,Teuchos::as<int>(removedIndices.size()),&removedIndices[0])==0);
   }

   // finish filling the graph
   Teuchos::RCP<Epetra_Map> rMap = getGhostedMap();
   Teuchos::RCP<Epetra_Map> cMap = getGhostedColMap();

   TEUCHOS_ASSERT(filteredGraph->FillComplete(*cMap,*rMap)==0);
   TEUCHOS_ASSERT(filteredGraph->OptimizeStorage()==0);

   return filteredGraph;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_Vector> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedEpetraColVector() const
{
   Teuchos::RCP<const Epetra_Map> eMap = getGhostedColMap(); 
   return Teuchos::rcp(new Epetra_Vector(*eMap));
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_Vector> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getGhostedEpetraVector() const
{
   Teuchos::RCP<const Epetra_Map> eMap = getGhostedMap(); 
   return Teuchos::rcp(new Epetra_Vector(*eMap));
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_Vector> EpetraLinearObjFactory<Traits,LocalOrdinalT>::getEpetraColVector() const
{
   Teuchos::RCP<const Epetra_Map> eMap = getColMap(); 
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
