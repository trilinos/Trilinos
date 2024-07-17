// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_BlockedEpetraLinearObjFactory_hpp__
#define   __Panzer_BlockedEpetraLinearObjFactory_hpp__

#include <map>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_GatherOrientation.hpp"
#include "Panzer_GatherSolution_BlockedEpetra.hpp"
#include "Panzer_GatherTangent_BlockedEpetra.hpp"
#include "Panzer_ScatterResidual_BlockedEpetra.hpp"
#include "Panzer_ScatterDirichletResidual_BlockedEpetra.hpp"
#include "Panzer_ScatterResidual_Epetra.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra.hpp"
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_GatherTangent_Epetra.hpp"
#include "Panzer_GatherOrientation.hpp"
#include "Panzer_ThyraObjFactory.hpp"
#include "Panzer_HashUtils.hpp"

#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

namespace panzer
{

template <typename Traits,typename LocalOrdinalT>
class BlockedEpetraLinearObjFactory : public LinearObjFactory<Traits>
                                    , public ThyraObjFactory<double> {
public:

   BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                 const Teuchos::RCP<const GlobalIndexer> & gidProvider,
                                 bool useDiscreteAdjoint=false);

   BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                 const Teuchos::RCP<const GlobalIndexer> & gidProvider,
                                 const Teuchos::RCP<const GlobalIndexer> & colGidProvider,
                                 bool useDiscreteAdjoint=false);

   virtual ~BlockedEpetraLinearObjFactory();

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
   virtual Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> buildReadOnlyDomainContainer() const;

   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<WriteVector_GlobalEvaluationData> buildWriteDomainContainer() const;

   virtual Teuchos::MpiComm<int> getComm() const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   {
     if(!colDOFManagerContainer_->containsBlockedDOFManager() &&
        !rowDOFManagerContainer_->containsBlockedDOFManager())
       return Teuchos::rcp(new ScatterResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers()[0],
                                                                                      colDOFManagerContainer_->getFieldDOFManagers()[0],
                                                                                      useDiscreteAdjoint_));

     return Teuchos::rcp(new ScatterResidual_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers(),
                                                                                           colDOFManagerContainer_->getFieldDOFManagers(),
                                                                                           useDiscreteAdjoint_));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   {
     if(!colDOFManagerContainer_->containsBlockedDOFManager() &&
        !rowDOFManagerContainer_->containsBlockedDOFManager())
       return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers()[0]));
     return Teuchos::rcp(new GatherSolution_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers()));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherTangent() const
   {
     if(!colDOFManagerContainer_->containsBlockedDOFManager() &&
        !rowDOFManagerContainer_->containsBlockedDOFManager())
       return Teuchos::rcp(new GatherTangent_Epetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers()[0]));
     return Teuchos::rcp(new GatherTangent_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers()));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherDomain() const
   {
     if(!colDOFManagerContainer_->containsBlockedDOFManager())
       return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits,LocalOrdinalT,int>(colDOFManagerContainer_->getFieldDOFManagers()[0]));
     return Teuchos::rcp(new GatherSolution_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(colDOFManagerContainer_->getFieldDOFManagers()));
   }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
  { return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,panzer::GlobalOrdinal>(rowDOFManagerContainer_->getFieldDOFManagers())); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   {
     if(!colDOFManagerContainer_->containsBlockedDOFManager() &&
        !rowDOFManagerContainer_->containsBlockedDOFManager())
       return Teuchos::rcp(new ScatterDirichletResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers()[0],
                                                                                               colDOFManagerContainer_->getFieldDOFManagers()[0]));
     return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(rowDOFManagerContainer_->getFieldDOFManagers(),
                                                                                                    colDOFManagerContainer_->getFieldDOFManagers()));
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
   void initializeGhostedContainer(int,LinearObjContainer & loc) const;

/*************** Thyra based methods *******************/

   //! Get the domain vector space (x and dxdt)
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraDomainSpace() const;

   //! Get the range vector space (f)
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraRangeSpace() const;

   //! Get a domain vector
   Teuchos::RCP<Thyra::VectorBase<double> > getThyraDomainVector() const;

   //! Get a range vector
   Teuchos::RCP<Thyra::VectorBase<double> > getThyraRangeVector() const;

   //! Get a Thyra operator
   Teuchos::RCP<Thyra::LinearOpBase<double> > getThyraMatrix() const;

   // and now the ghosted versions

   //! Get the domain vector space (x and dxdt)
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getGhostedThyraDomainSpace() const;

   /**
    *  \brief Get or create the ghosted `Thyra` domain space.
    *
    *  Get the vector space corresponding to the ghosted domain.  If it does
    *  not yet exist, create it from the ghosted column map(s).
    *
    *  \note This "version 2" routine works with non-overlapping owned and
    *        ghosted maps.
    *
    *  \returns The vector space corresponding to the ghosted domain.
    */
   Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
   getGhostedThyraDomainSpace2() const;

   //! Get the range vector space (f)
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getGhostedThyraRangeSpace() const;

   //! Get a domain vector
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedThyraDomainVector() const;

   //! Get a range vector
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedThyraRangeVector() const;

   //! Get a Thyra operator
   Teuchos::RCP<Thyra::LinearOpBase<double> > getGhostedThyraMatrix() const;

/*************** Epetra based methods *******************/

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getMap(int i) const;

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getColMap(int i) const;

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedMap(int i) const;

   /**
    *  \brief Get or create the `i`-th ghosted map.
    *
    *  \note This "version 2" routine works with non-overlapping owned and
    *        ghosted maps.
    *
    *  \param[in] i The index into the list of ghosted maps.
    *
    *  \returns The `i`-th ghosted map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   getGhostedMap2(
     int i) const;

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedColMap(int i) const;

   /**
    *  \brief Get or create the `i`-th ghosted column map.
    *
    *  \note This "version 2" routine works with non-overlapping owned and
    *        ghosted maps.
    *
    *  \param[in] i The index into the list of ghosted column maps.
    *
    *  \returns The `i`-th ghosted column map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   getGhostedColMap2(
     int i) const;

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGraph(int i,int j) const;

   //! get the ghosted graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGhostedGraph(int i,int j) const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getGhostedImport(int i) const;

   /**
    *  \brief Get or create the `i`-th ghosted importer corresponding to the
    *         `i`-th ghosted map.
    *
    *  \note This "version 2" routine works with non-overlapping owned and
    *        ghosted maps.
    *
    *  \param[in] i The index into the list of ghosted importers.
    *
    *  \returns The `i`-th ghosted importer.
    */
   virtual const Teuchos::RCP<Epetra_Import>
   getGhostedImport2(
     int i) const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getGhostedColImport(int i) const;

   /**
    *  \brief Get or create the `i`-th ghosted column importer corresponding to
    *         the `i`-th ghosted column map.
    *
    *  \note This "version 2" routine works with non-overlapping owned and
    *        ghosted maps.
    *
    *  \param[in] i The index into the list of ghosted column importers.
    *
    *  \returns The `i`-th ghosted column importer.
    */
   virtual const Teuchos::RCP<Epetra_Import>
   getGhostedColImport2(
     int i) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getGhostedExport(int j) const;

   /**
    *  \brief Get or create the `i`-th ghosted exporter corresponding to the
    *         `i`-th ghosted map.
    *
    *  \note This "version 2" routine works with non-overlapping owned and
    *        ghosted maps.
    *
    *  \param[in] i The index into the list of ghosted exporters.
    *
    *  \returns The `i`-th ghosted exporter.
    */
   virtual const Teuchos::RCP<Epetra_Export>
   getGhostedExport2(
     int i) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getGhostedColExport(int j) const;

   /**
    *  \brief Get or create the `i`-th ghosted column exporter corresponding to
    *         the `i`-th ghosted column map.
    *
    *  \note This "version 2" routine works with non-overlapping owned and
    *        ghosted maps.
    *
    *  \param[in] i The index into the list of ghosted column exporters.
    *
    *  \returns The `i`-th ghosted column exporter.
    */
   virtual const Teuchos::RCP<Epetra_Export>
   getGhostedColExport2(
     int i) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<const Epetra_Comm> getEpetraComm() const;

   Teuchos::RCP<Epetra_CrsMatrix> getEpetraMatrix(int i,int j) const;
   Teuchos::RCP<Epetra_CrsMatrix> getGhostedEpetraMatrix(int i,int j) const;

   //! how many block rows
   int getBlockRowCount() const;

   //! how many block columns
   int getBlockColCount() const;

   Teuchos::RCP<const panzer::BlockedDOFManager> getGlobalIndexer() const
   { return rowDOFManagerContainer_->getBlockedIndexer(); }

   Teuchos::RCP<const panzer::GlobalIndexer> getRangeGlobalIndexer() const
   { return rowDOFManagerContainer_->getGlobalIndexer(); }

   Teuchos::RCP<const panzer::GlobalIndexer> getDomainGlobalIndexer() const
   { return colDOFManagerContainer_->getGlobalIndexer(); }

   //! Get global indexers associated with the blocks
   const std::vector<Teuchos::RCP<const GlobalIndexer> > & getRangeGlobalIndexers() const
   { return rowDOFManagerContainer_->getFieldDOFManagers(); }

   //! Get global indexers associated with the blocks
   const std::vector<Teuchos::RCP<const GlobalIndexer> > & getDomainGlobalIndexers() const
   { return colDOFManagerContainer_->getFieldDOFManagers(); }

   //! exclude a block pair from the matrix
   void addExcludedPair(int rowBlock,int colBlock);

   //! exclude a vector of pairs from the matrix
   void addExcludedPairs(const std::vector<std::pair<int,int> > & exPairs);

protected:

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeContainer_internal(int mem,ThyraObjContainer<double> & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeGhostedContainer_internal(int mem,ThyraObjContainer<double> & loc) const;

/*************** Utility class for handling blocked and nonblocked DOF managers *******************/

  /** This classes is mean to abstract away the different global indexer types and hide
    * if this is a blocked data structure or a unblocked one.
    */
  class DOFManagerContainer {
  public:
    DOFManagerContainer() {}
    DOFManagerContainer(const Teuchos::RCP<const GlobalIndexer> & ugi)
    { setGlobalIndexer(ugi); }

    void setGlobalIndexer(const Teuchos::RCP<const GlobalIndexer> & ugi)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;

      auto blockedDOFManager = rcp_dynamic_cast<const BlockedDOFManager>(ugi);
      auto flatDOFManager    = rcp_dynamic_cast<const GlobalIndexer>(ugi);

      if(blockedDOFManager!=Teuchos::null) {
        // set BlockedDOFManager
        blockedDOFManager_ = blockedDOFManager;

        // get all GID providers
        auto dofManagers =  blockedDOFManager_->getFieldDOFManagers();
        for(auto itr=dofManagers.begin();itr!=dofManagers.end();++itr)
          gidProviders_.push_back(*itr);
      }
      else if(flatDOFManager!=Teuchos::null) {
        // for absolute clarity, nullify the blockedDOFManager_
        blockedDOFManager_ = Teuchos::null;

        // you have only a single GID provider
        gidProviders_.push_back(flatDOFManager);
      }
      else {
        TEUCHOS_ASSERT(false);
      }
    }

    //! Get the number of global indexers (not including the blocked one) contained.
    int getFieldBlocks() const
    { return Teuchos::as<int>(gidProviders_.size()); }

    /** Return true if this contains a blocked DOFManager as opposed to only a single DOFManager.
      * If this returns true then <code>getGlobalIndexer</code> will return a <code>BlockedDOFManager<int,GO></code>,
      * other wise it will return a <code>GlobalIndexer<int,GO></code>.
      */
    bool containsBlockedDOFManager() const
    { return blockedDOFManager_ !=Teuchos::null; }

    //! Get the "parent" global indexer (if <code>containsBlockedDOFManager()==false</code> this will throw)
    Teuchos::RCP<const BlockedDOFManager> getBlockedIndexer() const
    {
      TEUCHOS_ASSERT(containsBlockedDOFManager());
      return blockedDOFManager_;
    }

    //! Get the "parent" global indexer (if <code>getFieldBlocks()>1</code> this will be blocked, otherwise it may be either)
    Teuchos::RCP<const GlobalIndexer> getGlobalIndexer() const
    {
      if(blockedDOFManager_!=Teuchos::null)
        return blockedDOFManager_;

      TEUCHOS_ASSERT(gidProviders_.size()==1);
      return gidProviders_[0];
    }

    //! Get DOFManagers associated with the blocks
    const std::vector<Teuchos::RCP<const GlobalIndexer> > & getFieldDOFManagers() const
    { return gidProviders_; }

  private:
    Teuchos::RCP<const BlockedDOFManager> blockedDOFManager_;
    std::vector<Teuchos::RCP<const GlobalIndexer> > gidProviders_;
  };

/*************** Generic methods/members *******************/

   // Get the global indexer associated with a particular block
   Teuchos::RCP<const GlobalIndexer> getGlobalIndexer(int i) const;

   Teuchos::RCP<const GlobalIndexer> getColGlobalIndexer(int i) const;

   //! Allocate the space in the std::vector objects so we can fill with appropriate Epetra data
   void makeRoomForBlocks(std::size_t blockCnt,std::size_t colBlockCnt=0);

   Teuchos::RCP<const DOFManagerContainer> rowDOFManagerContainer_;
   Teuchos::RCP<const DOFManagerContainer> colDOFManagerContainer_;

   bool useColGidProviders_;

   // which block entries are ignored
   std::unordered_set<std::pair<int,int>,panzer::pair_hash> excludedPairs_;

/*************** Thyra based methods/members *******************/

   void ghostToGlobalThyraVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<double> > & out,bool col) const;
   void ghostToGlobalThyraMatrix(const Thyra::LinearOpBase<double> & in,Thyra::LinearOpBase<double> & out) const;
   void globalToGhostThyraVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<double> > & out,bool col) const;

   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace_;

   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ghostedRangeSpace_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ghostedDomainSpace_;

/*************** Epetra based methods/members *******************/

   void adjustForDirichletConditions(const Epetra_Vector & local_bcs,
                                     const Epetra_Vector & global_bcs,
                                     const Teuchos::Ptr<Epetra_Vector> & f,
                                     const Teuchos::Ptr<Epetra_CrsMatrix> & A,
                                     bool zeroVectorRows) const;

   void ghostToGlobalEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out,bool col) const;
   void globalToGhostEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out,bool col) const;
   void ghostToGlobalEpetraMatrix(int blockRow,const Epetra_CrsMatrix & in,Epetra_CrsMatrix & out) const;

   // get the map from the matrix

   /**
    *  \brief Build the `i`-th owned map from the owned indices of the `i`-th
    *         global indexer.
    *
    *  \param[in] i The index into the list of global indexers.
    *
    *  \returns The `i`-th owned map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   buildMap(
     int i) const;

   /**
    *  \brief Build the `i`-th ghosted map from the owned and ghosted indices
    *         of the `i`-th global indexer.
    *
    *  \param[in] i The index into the list of global indexers.
    *
    *  \returns The `i`-th owned and ghosted map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   buildGhostedMap(
     int i) const;

   /**
    *  \brief Build the `i`-th ghosted map from the ghosted indices of the
    *         `i`-th global indexer.
    *
    *  \param[in] i The index into the list of global indexers.
    *
    *  \returns The `i`-th ghosted map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   buildGhostedMap2(
     int i) const;

   // get the map from the matrix

   /**
    *  \brief Build the `i`-th owned column map from the owned indices of the
    *         `i`-th (column) global indexer.
    *
    *  \param[in] i The index into the list of (column) global indexers.
    *
    *  \returns The `i`-th owned column map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   buildColMap(
     int i) const;

   /**
    *  \brief Build the `i`-th ghosted column map from the owned and ghosted
    *         indices of the `i`-th (column) global indexer.
    *
    *  \param[in] i The index into the list of (column) global indexers.
    *
    *  \returns The `i`-th owned and ghosted column map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   buildColGhostedMap(
     int i) const;

   /**
    *  \brief Build the `i`-th ghosted column map from the ghosted indices of
    *         the `i`-th (column) global indexer.
    *
    *  \param[in] i The index into the list of (column) global indexers.
    *
    *  \returns The `i`-th ghosted column map.
    */
   virtual const Teuchos::RCP<Epetra_Map>
   buildColGhostedMap2(
     int i) const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGraph(int i,int j) const;
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGhostedGraph(int i,int j,bool optimizeStorage) const;
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildFilteredGhostedGraph(int i,int j) const;

   // storage for Epetra graphs and maps
   Teuchos::RCP<const Epetra_Comm> eComm_;
   Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm_;
   Teuchos::RCP<Teuchos::MpiComm<int> > tComm_;

   /**
    *  \brief The list of owned maps corresponding to the owned indices of the
    *         global indexers.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Map>> maps_;

   /**
    *  \brief The list of ghosted maps corresponding to the owned and ghosted
    *         indices of the global indexers.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Map>> ghostedMaps_;

   /**
    *  \brief The list of ghosted maps corresponding to the ghosted indices of
    *         the global indexers.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Map>> ghostedMaps2_;

   /**
    *  \brief The list of ghosted importers corresponding to `ghostedMaps_`.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Import>> importers_;

   /**
    *  \brief The list of ghosted importers corresponding to `ghostedMaps2_`.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Import>> importers2_;

   mutable std::vector<Teuchos::RCP<Epetra_Export>> exporters_;

   /**
    *  \brief The list of owned column maps corresponding to the owned indices
    *         of the (column) global indexers.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Map>> colMaps_;

   /**
    *  \brief The list of ghosted column maps corresponding to the owned and
    *         ghosted indices of the (column) global indexers.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Map>> colGhostedMaps_;

   /**
    *  \brief The list of ghosted column maps corresponding to the ghosted
    *         indices of the (column) global indexers.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Map>> colGhostedMaps2_;

   /**
    *  \brief The list of ghosted importers corresponding to `colGhostedMaps_`.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Import>> colImporters_;

   /**
    *  \brief The list of ghosted importers corresponding to
    *         `colGhostedMaps2_`.
    */
   mutable std::vector<Teuchos::RCP<Epetra_Import>> colImporters2_;

   mutable std::vector<Teuchos::RCP<Epetra_Export>> colExporters_;

   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph>,panzer::pair_hash> graphs_ ;
   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph>,panzer::pair_hash> ghostedGraphs_;

   bool useDiscreteAdjoint_;
};

} // end of namespace panzer

#endif // __Panzer_BlockedEpetraLinearObjFactory_hpp__