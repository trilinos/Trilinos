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

#ifndef __Panzer_EpetraLinearObjFactory_decl_hpp__
#define __Panzer_EpetraLinearObjFactory_decl_hpp__

#define PANZER_USE_BLOCKED_EPETRA_LOF

#ifndef PANZER_USE_BLOCKED_EPETRA_LOF
#include <map>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_ScatterResidual_Epetra.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra.hpp"
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_GatherTangent_Epetra.hpp"
#include "Panzer_GatherOrientation.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_ThyraObjFactory.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

namespace panzer {

template <typename Traits,typename LocalOrdinalT>
class EpetraLinearObjFactory : public LinearObjFactory<Traits>
                             , public ThyraObjFactory<double> {
public:

   EpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                          const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider,
                          bool useDiscreteAdjoint=false);

   EpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                          const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & rowProvider,
                          const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & colProvider,
                          bool useDiscreteAdjoint=false);

   /** This has been added for the case when you want an epetra LOF but don't have a comm. It simply can build
     * scatter and gather evaluators. But not necessarily the objects that go into them.
     */
   EpetraLinearObjFactory(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & rowProvider);

   virtual ~EpetraLinearObjFactory();

/*************** Linear object factory methods *******************/

   virtual void readVector(const std::string & identifier,LinearObjContainer & loc,int id) const;

   virtual void writeVector(const std::string & identifier,const LinearObjContainer & loc,int id) const;

   virtual Teuchos::RCP<LinearObjContainer> buildLinearObjContainer() const;

   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveLinearObjContainer() const 
   { return buildLinearObjContainer(); }

   virtual Teuchos::RCP<LinearObjContainer> buildGhostedLinearObjContainer() const;

   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveGhostedLinearObjContainer() const 
   { return buildGhostedLinearObjContainer(); }

   virtual void globalToGhostContainer(const LinearObjContainer & container,
                                       LinearObjContainer & ghostContainer,int) const;
   virtual void ghostToGlobalContainer(const LinearObjContainer & ghostContainer,
                                       LinearObjContainer & container,int) const;

   /** Adjust the residual vector and Jacobian matrix (if they exist) for applied
     * dirichlet conditions. The adjustment considers if a boundary condition was
     * set globally and locally and based on that result adjust the ghosted matrix
     * and residual vector so that when they are summed across processors they resulting
     * Dirichlet condition is correct.
     */
   virtual void adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                                             const LinearObjContainer & globalBCRows,
                                             LinearObjContainer & ghostedObjs,
                                             bool zeroVectorRows=false, bool adjustX=false) const;

   /** Adjust a vector by replacing selected rows with the value of the evaluated
     * dirichlet conditions. This is handled through the standard container mechanism.
     */
   virtual void applyDirichletBCs(const LinearObjContainer & counter,
                                  LinearObjContainer & result) const;

   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> buildDomainContainer() const;

   /** Acess to the MPI Comm used in constructing this LOF.
     */
   virtual Teuchos::MpiComm<int> getComm() const;

   //! Get the domain space
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraDomainSpace() const;

   //! Get the range space
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraRangeSpace() const;

   //! Get a matrix operator
   Teuchos::RCP<Thyra::LinearOpBase<double> > getThyraMatrix() const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   { return Teuchos::rcp(new ScatterResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_,colGidProvider_,useDiscreteAdjoint_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherTangent() const
   { return Teuchos::rcp(new GatherTangent_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherDomain() const
   { if(colGidProvider_!=Teuchos::null)
       return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits,LocalOrdinalT,int>(colGidProvider_));
     return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   { return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return Teuchos::rcp(new ScatterDirichletResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_,colGidProvider_)); }

/*************** Generic helper functions for container setup *******************/
   
   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeContainer(int,LinearObjContainer & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeContainer(int mem,EpetraLinearObjContainer & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeGhostedContainer(int,LinearObjContainer & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeGhostedContainer(int mem,EpetraLinearObjContainer & loc) const;

/*************** Epetra based methods *******************/

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getMap() const;

   //! get the map from the matrix (included for migration to blocked epetra LOF)
   const Teuchos::RCP<Epetra_Map> getMap(int) const { return getMap(); }

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getColMap() const;

   virtual const Teuchos::RCP<Epetra_Map> getColMap(int) const { return getColMap(); }

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedMap() const;
   virtual const Teuchos::RCP<Epetra_Map> getGhostedMap2() const;

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedMap(int) const { return getGhostedMap(); }

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedColMap() const;
   virtual const Teuchos::RCP<Epetra_Map> getGhostedColMap2() const;

   virtual const Teuchos::RCP<Epetra_Map> getGhostedColMap(int) const { return getGhostedColMap(); }

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGraph() const;

   //! get the graph of the crs matrix (included for migration to blocked epetra LOF)
   const Teuchos::RCP<Epetra_CrsGraph> getGraph(int,int) const { return getGraph(); }

   //! get the ghosted graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGhostedGraph() const;

   const Teuchos::RCP<Epetra_CrsGraph> getGhostedGraph(int,int) const { return getGhostedGraph(); }

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getGhostedImport() const;
   virtual const Teuchos::RCP<Epetra_Import> getGhostedImport2() const;

   virtual const Teuchos::RCP<Epetra_Import> getGhostedImport(int ) const { return getGhostedImport(); }

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getGhostedColImport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getGhostedExport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getGhostedColExport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<const Epetra_Comm> getEpetraComm() const;

   //! Get the domain unique global indexer this factory was created with.
   Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> getDomainGlobalIndexer() const
   { return colGidProvider_!=Teuchos::null ? colGidProvider_ : gidProvider_; }

   //! Get the range unique global indexer this factory was created with.
   Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> getRangeGlobalIndexer() const
   { return gidProvider_; }

protected:
   Teuchos::RCP<Epetra_Vector> getGhostedEpetraVector() const;
   Teuchos::RCP<Epetra_Vector> getGhostedEpetraColVector() const;
   Teuchos::RCP<Epetra_Vector> getEpetraVector() const;
   Teuchos::RCP<Epetra_Vector> getEpetraColVector() const;
   Teuchos::RCP<Epetra_CrsMatrix> getEpetraMatrix() const;
   Teuchos::RCP<Epetra_CrsMatrix> getGhostedEpetraMatrix() const;

   void ghostToGlobalEpetraVector(const Epetra_Vector & in,Epetra_Vector & out,bool col=false) const;
   void ghostToGlobalEpetraMatrix(const Epetra_CrsMatrix & in,Epetra_CrsMatrix & out) const;
   void globalToGhostEpetraVector(const Epetra_Vector & in,Epetra_Vector & out,bool col=false) const;

   // get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> buildMap() const;
   virtual const Teuchos::RCP<Epetra_Map> buildColMap() const;
   virtual const Teuchos::RCP<Epetra_Map> buildGhostedMap() const;
   virtual const Teuchos::RCP<Epetra_Map> buildGhostedMap2() const;
   virtual const Teuchos::RCP<Epetra_Map> buildGhostedColMap() const;
   virtual const Teuchos::RCP<Epetra_Map> buildGhostedColMap2() const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGraph() const;
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGhostedGraph(bool optimizeStorage) const;

   // this method by defaults calls getGhostedGraph, however if the column global indexer
   // is filtered (i.e. a Filtered_UniqueGlobalIndexer) then the filtered columns are removed
   // in the version of the graph return by this. Note if the column UGI is filtered, then their
   // will be substantial computational cost including communication in this method.
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildFilteredGhostedGraph() const;

   // storage for Epetra graphs and maps
   Teuchos::RCP<const Epetra_Comm> comm_;
   mutable Teuchos::RCP<Epetra_Map> map_;
   mutable Teuchos::RCP<Epetra_Map> cMap_;
   mutable Teuchos::RCP<Epetra_Map> ghostedMap_;
   mutable Teuchos::RCP<Epetra_Map> cGhostedMap_;
   mutable Teuchos::RCP<Epetra_CrsGraph> graph_;
   mutable Teuchos::RCP<Epetra_CrsGraph> ghostedGraph_;

   // import/exporter storage
   mutable Teuchos::RCP<Epetra_Import> importer_, colImporter_;
   mutable Teuchos::RCP<Epetra_Export> exporter_, colExporter_;

   bool hasColProvider_;
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > gidProvider_;
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > colGidProvider_;
   Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm_;

   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace_;
 
   bool useDiscreteAdjoint_;
};

}

#else

#include "Panzer_BlockedEpetraLinearObjFactory.hpp"

namespace panzer {

template <typename Traits,typename LocalOrdinalT>
using EpetraLinearObjFactory = BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>;

}

#endif

#endif
