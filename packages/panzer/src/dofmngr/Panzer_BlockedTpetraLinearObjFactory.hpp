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

#ifndef __Panzer_BlockedTpetraLinearObjFactory_hpp__
#define __Panzer_BlockedTpetraLinearObjFactory_hpp__

#include <map>

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"

#include "Panzer_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_GatherOrientation.hpp"
#include "Panzer_GatherSolution_BlockedTpetra.hpp"
#include "Panzer_ScatterResidual_BlockedTpetra.hpp"
// #include "Panzer_ScatterDirichletResidual_BlockedEpetra.hpp"
// #include "Panzer_ScatterInitialCondition_BlockedEpetra.hpp"
#include "Panzer_ThyraObjFactory.hpp"

#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

namespace panzer {

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=Kokkos::DefaultNode::DefaultNodeType>
class BlockedTpetraLinearObjFactory : public LinearObjFactory<Traits>
                                    , public ThyraObjFactory<double> {
public:
   typedef BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> BTLOC;
   typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
   typedef Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> CrsMatrixType;
   typedef Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> OperatorType;
   typedef Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> CrsGraphType;
   typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
   typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;
   typedef Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> ExportType;

   typedef Thyra::TpetraVector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> ThyraVector;
   typedef Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> ThyraLinearOp;


   BlockedTpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                 const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,std::pair<int,GlobalOrdinalT> > > & blkProvider,
                                 const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > & gidProviders);
   BlockedTpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                 const Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT> > & gidProvider);

   virtual ~BlockedTpetraLinearObjFactory();

/*************** Linear object factory methods *******************/

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
                                             LinearObjContainer & ghostedObjs) const;

   Teuchos::MpiComm<int> getComm() const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   { return Teuchos::rcp(new ScatterResidual_BlockedTpetra<EvalT,Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(blockedDOFManager_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return Teuchos::rcp(new GatherSolution_BlockedTpetra<EvalT,Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(blockedDOFManager_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   { return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,std::pair<int,GlobalOrdinalT> >(blockProvider_)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   // { return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(blockedDOFManager_)); }
   { return Teuchos::null; }

   //! Use preconstructed initial condition scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterInitialCondition() const
   // { return Teuchos::rcp(new ScatterInitialCondition_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(blockedDOFManager_)); }
   { return Teuchos::null; }

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
   void initializeContainer(int mem,BTLOC & loc) const;

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
   void initializeGhostedContainer(int mem,BTLOC & loc) const;

/*************** Thyra based methods *******************/

   //! Get the domain vector space (x and dxdt)
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraDomainSpace() const;

   //! Get the range vector space (f)
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraRangeSpace() const;

   //! Get a domain vector
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > getThyraDomainVector() const;

   //! Get a range vector
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > getThyraRangeVector() const;

   //! Get a Thyra operator
   Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > getThyraMatrix() const;

   // and now the ghosted versions

   //! Get the domain vector space (x and dxdt)
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getGhostedThyraDomainSpace() const;

   //! Get the range vector space (f)
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getGhostedThyraRangeSpace() const;

   //! Get a domain vector
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > getGhostedThyraDomainVector() const;

   //! Get a range vector
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > getGhostedThyraRangeVector() const;

   //! Get a Thyra operator
   Teuchos::RCP<Thyra::BlockedLinearOpBase<ScalarT> > getGhostedThyraMatrix() const;

/*************** Tpetra based methods *******************/

   //! get the map from the matrix
   virtual Teuchos::RCP<const MapType> getMap(int i) const;

   //! get the ghosted map from the matrix
   virtual Teuchos::RCP<const MapType> getGhostedMap(int i) const;

   //! get the graph of the crs matrix
   virtual Teuchos::RCP<const CrsGraphType> getGraph(int i,int j) const;

   //! get the ghosted graph of the crs matrix
   virtual Teuchos::RCP<const CrsGraphType> getGhostedGraph(int i,int j) const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual Teuchos::RCP<const ImportType> getGhostedImport(int i) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual Teuchos::RCP<const ExportType> getGhostedExport(int j) const;

   Teuchos::RCP<CrsMatrixType> getTpetraMatrix(int i,int j) const;
   Teuchos::RCP<CrsMatrixType> getGhostedTpetraMatrix(int i,int j) const;

   //! how many block rows
   int getBlockRowCount() const;

   //! how many block columns
   int getBlockColCount() const;

   Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT> > getGlobalIndexer() const
   { return blockedDOFManager_; }

   //! exclude a block pair from the matrix
   void addExcludedPair(int rowBlock,int colBlock);

   //! exclude a vector of pairs from the matrix
   void addExcludedPairs(const std::vector<std::pair<int,int> > & exPairs);

protected:
/*************** Generic methods/members *******************/

   // Get the global indexer associated with a particular block
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > getGlobalIndexer(int i) const;

   //! Allocate the space in the std::vector objects so we can fill with appropriate Tpetra data
   void makeRoomForBlocks(std::size_t blockCnt);

   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,std::pair<int,GlobalOrdinalT> > > blockProvider_;
   Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT> > blockedDOFManager_;
   std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > gidProviders_;

   // which block entries are ignored
   boost::unordered_set<std::pair<int,int> > excludedPairs_;
  
/*************** Thyra based methods/members *******************/

   void ghostToGlobalThyraVector(const Teuchos::RCP<const Thyra::VectorBase<ScalarT> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & out) const;
   void ghostToGlobalThyraMatrix(const Thyra::LinearOpBase<ScalarT> & in,Thyra::LinearOpBase<ScalarT> & out) const;
   void globalToGhostThyraVector(const Teuchos::RCP<const Thyra::VectorBase<ScalarT> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & out) const;

   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > rangeSpace_;
   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > domainSpace_;

   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > ghostedRangeSpace_;
   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > ghostedDomainSpace_;

/*************** Tpetra based methods/members *******************/

   void adjustForDirichletConditions(const VectorType & local_bcs,
                                     const VectorType & global_bcs,
                                     const Teuchos::Ptr<VectorType> & f,
                                     const Teuchos::Ptr<CrsMatrixType> & A) const;

   void ghostToGlobalTpetraVector(int i,const VectorType & in,VectorType & out) const;
   void ghostToGlobalTpetraMatrix(int blockRow,const CrsMatrixType & in,CrsMatrixType & out) const;
   void globalToGhostTpetraVector(int i,const VectorType & in,VectorType & out) const;

   // get the map from the matrix
   virtual Teuchos::RCP<const MapType> buildTpetraMap(int i) const;
   virtual Teuchos::RCP<const MapType> buildTpetraGhostedMap(int i) const;

   // get the graph of the crs matrix
   virtual Teuchos::RCP<const CrsGraphType> buildTpetraGraph(int i,int j) const;
   virtual Teuchos::RCP<const CrsGraphType> buildTpetraGhostedGraph(int i,int j) const;

   // storage for Tpetra graphs and maps
   Teuchos::RCP<const Teuchos::MpiComm<int> > comm_;
   mutable std::vector<Teuchos::RCP<const MapType> > maps_;
   mutable std::vector<Teuchos::RCP<const MapType> > ghostedMaps_;
   mutable boost::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType> > graphs_ ;
   mutable boost::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType> > ghostedGraphs_;

   mutable std::vector<Teuchos::RCP<const ImportType> > importers_;
   mutable std::vector<Teuchos::RCP<const ExportType> > exporters_;
};

}

#include "Panzer_BlockedTpetraLinearObjFactory_impl.hpp"

#endif
