// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_BlockedTpetraLinearObjFactory_hpp__
#define __Panzer_BlockedTpetraLinearObjFactory_hpp__

#include <map>

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_HashUtils.hpp" // for pair_hash

#include "Panzer_GatherOrientation.hpp"
#include "Panzer_GatherSolution_BlockedTpetra.hpp"
#include "Panzer_GatherTangent_BlockedTpetra.hpp"
#include "Panzer_ScatterResidual_BlockedTpetra.hpp"
#include "Panzer_ScatterDirichletResidual_BlockedTpetra.hpp"
#include "Panzer_ThyraObjFactory.hpp"

#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

namespace panzer {

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=panzer::TpetraNodeType>
class BlockedTpetraLinearObjFactory : public LinearObjFactory<Traits>
                                    , public ThyraObjFactory<double> {
public:
   typedef TpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> TLOC;
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
                                 const Teuchos::RCP<const BlockedDOFManager> & gidProvider);

  /** \brief Ctor that takes a vector of DOFManagers instead of the
      BlockedDOFManager. Plan is to deprecate the BlockedDOFManager,
      but for now it is ingrained in all gather/scatter operators.
   */
   BlockedTpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                 const std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> & gidProviders);

   virtual ~BlockedTpetraLinearObjFactory();

/*************** Linear object factory methods *******************/

   virtual void readVector(const std::string & /* identifier */, LinearObjContainer & /* loc */, int /* id */) const
   { TEUCHOS_ASSERT(false); }

   virtual void writeVector(const std::string & /* identifier */, const LinearObjContainer & /* loc */, int /* id */) const
   { TEUCHOS_ASSERT(false); }

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
                                             bool zeroVectorRows=false, bool adjustX = false) const;

   /** Adjust a vector by replacing selected rows with the value of the evaluated
     * dirichlet conditions. This is handled through the standard container mechanism.
     */
   virtual void applyDirichletBCs(const LinearObjContainer & counter,
                                  LinearObjContainer & result) const;

   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> buildReadOnlyDomainContainer() const;

#ifdef PANZER_HAVE_EPETRA_STACK
   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<WriteVector_GlobalEvaluationData> buildWriteDomainContainer() const;
#endif

   Teuchos::MpiComm<int> getComm() const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   {
      if(gidProviders_.front()->containsBlockedDOFManager()) {
         return Teuchos::rcp(new ScatterResidual_BlockedTpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGlobalIndexer()));
      }
      return Teuchos::rcp(new ScatterResidual_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(getRangeGlobalIndexer()));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   {
      if(gidProviders_.front()->containsBlockedDOFManager()) {
         return Teuchos::rcp(new GatherSolution_BlockedTpetra<EvalT,Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGlobalIndexer()));
      }
      return Teuchos::rcp(new GatherSolution_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(getRangeGlobalIndexer()));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherTangent() const
   {
      if(gidProviders_.front()->containsBlockedDOFManager()) {
         return Teuchos::rcp(new GatherTangent_BlockedTpetra<EvalT,Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGlobalIndexer()));
      }
      return Teuchos::rcp(new GatherTangent_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(getRangeGlobalIndexer()));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherDomain() const
   {
      if(gidProviders_.back()->containsBlockedDOFManager()) {
         return Teuchos::rcp(new GatherSolution_BlockedTpetra<EvalT,Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getDomainBlockedIndexer()));
      }
      return Teuchos::rcp(new GatherSolution_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(getDomainGlobalIndexer()));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   {
      if(gidProviders_.front()->containsBlockedDOFManager()) {
         return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT>(nc2c_vector(getGlobalIndexer()->getFieldDOFManagers())));
      }
      return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT>(getRangeGlobalIndexer()));
   }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   {
      if(gidProviders_.front()->containsBlockedDOFManager()) {
         return Teuchos::rcp(new ScatterDirichletResidual_BlockedTpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGlobalIndexer()));
      }
      return Teuchos::rcp(new ScatterDirichletResidual_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(getRangeGlobalIndexer()));
   }

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

   //! Get the domain vector space (x and dxdt)
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraDomainSpace(int blk) const;

   //! Get the range vector space (f)
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraRangeSpace(int blk) const;

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
   virtual Teuchos::RCP<const MapType> getColMap(int i) const;

   //! get the ghosted map from the matrix
   virtual Teuchos::RCP<const MapType> getGhostedMap(int i) const;
   virtual Teuchos::RCP<const MapType> getColGhostedMap(int i) const;

   //! get the graph of the crs matrix
   virtual Teuchos::RCP<const CrsGraphType> getGraph(int i,int j) const;

   //! get the ghosted graph of the crs matrix
   virtual Teuchos::RCP<const CrsGraphType> getGhostedGraph(int i,int j) const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual Teuchos::RCP<const ImportType> getGhostedImport(int i) const;
   virtual Teuchos::RCP<const ImportType> getColGhostedImport(int i) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual Teuchos::RCP<const ExportType> getGhostedExport(int j) const;
   virtual Teuchos::RCP<const ExportType> getColGhostedExport(int j) const;

   Teuchos::RCP<CrsMatrixType> getTpetraMatrix(int i,int j) const;
   Teuchos::RCP<CrsMatrixType> getGhostedTpetraMatrix(int i,int j) const;

   Teuchos::RCP<VectorType> getTpetraDomainVector(int i) const;
   Teuchos::RCP<VectorType> getGhostedTpetraDomainVector(int i) const;

   Teuchos::RCP<VectorType> getTpetraRangeVector(int i) const;
   Teuchos::RCP<VectorType> getGhostedTpetraRangeVector(int i) const;

   //! how many block rows
   int getBlockRowCount() const;

   //! how many block columns
   int getBlockColCount() const;

   //! exclude a block pair from the matrix
   void addExcludedPair(int rowBlock,int colBlock);

   //! exclude a vector of pairs from the matrix
   void addExcludedPairs(const std::vector<std::pair<int,int> > & exPairs);

   virtual void beginFill(LinearObjContainer & loc) const;
   virtual void endFill(LinearObjContainer & loc) const;

   Teuchos::RCP<const panzer::BlockedDOFManager> getGlobalIndexer() const
   {
      TEUCHOS_ASSERT(!gidProviders_.empty());
      return gidProviders_.front()->getBlockedIndexer();
   }

   //! Get the domain unique global indexer this factory was created with.
   Teuchos::RCP<const panzer::GlobalIndexer> getDomainGlobalIndexer() const
   {
      TEUCHOS_ASSERT(!gidProviders_.empty());
      return gidProviders_.back()->getGlobalIndexer();
   }

   //! Get the domain unique global indexer this factory was created with.
   Teuchos::RCP<const panzer::BlockedDOFManager> getDomainBlockedIndexer() const
   {
      TEUCHOS_ASSERT(!gidProviders_.empty());
      return gidProviders_.back()->getBlockedIndexer();
   }

   //! Get the range unique global indexer this factory was created with.
   Teuchos::RCP<const panzer::GlobalIndexer> getRangeGlobalIndexer() const
   {
      TEUCHOS_ASSERT(!gidProviders_.empty());
      return gidProviders_.front()->getGlobalIndexer();
   }

protected:
   /*************** Utility class for handling blocked and nonblocked DOF managers *******************/

   /** This classes is mean to abstract away the different global indexer types and hide
    * if this is a blocked data structure or a unblocked one.
    */
   class DOFManagerContainer
   {
   public:
      DOFManagerContainer() {}
      DOFManagerContainer(const Teuchos::RCP<const GlobalIndexer> &ugi)
      {
         setGlobalIndexer(ugi);
      }

      void setGlobalIndexer(const Teuchos::RCP<const GlobalIndexer> &ugi)
      {
         using Teuchos::RCP;
         using Teuchos::rcp_dynamic_cast;

         auto blockedDOFManager = rcp_dynamic_cast<const BlockedDOFManager>(ugi);
         auto flatDOFManager = rcp_dynamic_cast<const GlobalIndexer>(ugi);

         if (blockedDOFManager != Teuchos::null)
         {
            // set BlockedDOFManager
            blockedDOFManager_ = blockedDOFManager;

            // get all GID providers
            auto dofManagers = blockedDOFManager_->getFieldDOFManagers();
            for (auto itr = dofManagers.begin(); itr != dofManagers.end(); ++itr)
               gidProviders_.push_back(*itr);
         }
         else if (flatDOFManager != Teuchos::null)
         {
            // for absolute clarity, nullify the blockedDOFManager_
            blockedDOFManager_ = Teuchos::null;

            // you have only a single GID provider
            gidProviders_.push_back(flatDOFManager);
         }
         else
         {
            TEUCHOS_ASSERT(false);
         }
      }

      //! Get the number of global indexers (not including the blocked one) contained.
      int getFieldBlocks() const
      {
         return Teuchos::as<int>(gidProviders_.size());
      }

      /** Return true if this contains a blocked DOFManager as opposed to only a single DOFManager.
       * If this returns true then <code>getGlobalIndexer</code> will return a <code>BlockedDOFManager<int,GO></code>,
       * other wise it will return a <code>GlobalIndexer<int,GO></code>.
       */
      bool containsBlockedDOFManager() const
      {
         return blockedDOFManager_ != Teuchos::null;
      }

      //! Get the "parent" global indexer (if <code>containsBlockedDOFManager()==false</code> this will throw)
      Teuchos::RCP<const BlockedDOFManager> getBlockedIndexer() const
      {
         TEUCHOS_ASSERT(containsBlockedDOFManager());
         return blockedDOFManager_;
      }

      //! Get the "parent" global indexer (if <code>getFieldBlocks()>1</code> this will be blocked, otherwise it may be either)
      Teuchos::RCP<const GlobalIndexer> getGlobalIndexer() const
      {
         if (blockedDOFManager_ != Teuchos::null)
            return blockedDOFManager_;

         TEUCHOS_ASSERT(gidProviders_.size() == 1);
         return gidProviders_[0];
      }

      //! Get DOFManagers associated with the blocks
      const std::vector<Teuchos::RCP<const GlobalIndexer>> &getFieldDOFManagers() const
      {
         return gidProviders_;
      }

   private:
      Teuchos::RCP<const BlockedDOFManager> blockedDOFManager_;
      std::vector<Teuchos::RCP<const GlobalIndexer>> gidProviders_;
   };

/*************** Generic methods/members *******************/

   // Get the global indexer associated with a particular block
   Teuchos::RCP<const GlobalIndexer> getGlobalIndexer(int i) const;
   Teuchos::RCP<const GlobalIndexer> getColGlobalIndexer(int i) const;

   //! Allocate the space in the std::vector objects so we can fill with appropriate Tpetra data
   void makeRoomForBlocks(std::size_t blockCnt, std::size_t colBlockCnt = 0);

   std::vector<Teuchos::RCP<const DOFManagerContainer> > gidProviders_;

   bool useColGidProviders_;

   // which block entries are ignored
  std::unordered_set<std::pair<int,int>,panzer::pair_hash> excludedPairs_;

/*************** Thyra based methods/members *******************/

   void ghostToGlobalThyraVector(const Teuchos::RCP<const Thyra::VectorBase<ScalarT> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & out, bool col) const;
   void ghostToGlobalThyraMatrix(const Thyra::LinearOpBase<ScalarT> & in,Thyra::LinearOpBase<ScalarT> & out) const;
   void globalToGhostThyraVector(const Teuchos::RCP<const Thyra::VectorBase<ScalarT> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & out, bool col) const;

   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > rangeSpace_;
   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > domainSpace_;

   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > ghostedRangeSpace_;
   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<ScalarT> > ghostedDomainSpace_;

/*************** Tpetra based methods/members *******************/

   void adjustForDirichletConditions(const VectorType & local_bcs,
                                     const VectorType & global_bcs,
                                     const Teuchos::Ptr<VectorType> & f,
                                     const Teuchos::Ptr<CrsMatrixType> & A,
                                     bool zeroVectorRows) const;

   void ghostToGlobalTpetraVector(int i,const VectorType & in,VectorType & out, bool col) const;
   void ghostToGlobalTpetraMatrix(int blockRow,const CrsMatrixType & in,CrsMatrixType & out) const;
   void globalToGhostTpetraVector(int i,const VectorType & in,VectorType & out, bool col) const;

   // get the map from the matrix
   virtual Teuchos::RCP<const MapType> buildTpetraMap(int i) const;
   virtual Teuchos::RCP<const MapType> buildColTpetraMap(int i) const;

   virtual Teuchos::RCP<const MapType> buildTpetraGhostedMap(int i) const;
   virtual Teuchos::RCP<const MapType> buildColTpetraGhostedMap(int i) const;

   // get the graph of the crs matrix
   virtual Teuchos::RCP<const CrsGraphType> buildTpetraGraph(int i,int j) const;
   virtual Teuchos::RCP<const CrsGraphType> buildTpetraGhostedGraph(int i,int j) const;

   // storage for Tpetra graphs and maps
   Teuchos::RCP<const Teuchos::MpiComm<int> > comm_;
   mutable std::vector<Teuchos::RCP<const MapType> > maps_;
   mutable std::vector<Teuchos::RCP<const MapType> > colMaps_;

   mutable std::vector<Teuchos::RCP<const MapType> > ghostedMaps_;
   mutable std::vector<Teuchos::RCP<const MapType> > colGhostedMaps_;

   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType>,panzer::pair_hash> graphs_ ;
   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType>,panzer::pair_hash> colGraphs_ ;

   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType>,panzer::pair_hash> ghostedGraphs_;
   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType>,panzer::pair_hash> colGhostedGraphs_;

   mutable std::vector<Teuchos::RCP<const ImportType> > importers_;
   mutable std::vector<Teuchos::RCP<const ImportType> > colImporters_;

   mutable std::vector<Teuchos::RCP<const ExportType> > exporters_;
   mutable std::vector<Teuchos::RCP<const ExportType> > colExporters_;
};

}

#endif
