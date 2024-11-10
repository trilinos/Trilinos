// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_TpetraLinearObjFactory_decl_hpp__
#define __Panzer_TpetraLinearObjFactory_decl_hpp__

#include <map>

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_ScatterResidual_Tpetra.hpp"
#include "Panzer_ScatterDirichletResidual_Tpetra.hpp"
#include "Panzer_GatherSolution_Tpetra.hpp"
#include "Panzer_GatherTangent_Tpetra.hpp"
#include "Panzer_GatherOrientation.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_ThyraObjFactory.hpp"

#include"Panzer_NodeType.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

namespace panzer {

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=panzer::TpetraNodeType>
class TpetraLinearObjFactory : public LinearObjFactory<Traits>
                             , public ThyraObjFactory<ScalarT> {
public:
   typedef TpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> ContainerType;
   typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
   typedef Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> CrsMatrixType;
   typedef Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> CrsGraphType;
   typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
   typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;
   typedef Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> ExportType;

   TpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                          const Teuchos::RCP<const GlobalIndexer> & gidProvider);

   TpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                          const Teuchos::RCP<const GlobalIndexer> & rowProvider,
                          const Teuchos::RCP<const GlobalIndexer> & colProvider);

   virtual ~TpetraLinearObjFactory();

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

#ifdef PANZER_HAVE_EPETRA_STACK
   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<WriteVector_GlobalEvaluationData> buildWriteDomainContainer() const;
#endif

   /** Acess to the MPI Comm used in constructing this LOF.
     */
   virtual Teuchos::MpiComm<int> getComm() const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   { return Teuchos::rcp(new ScatterResidual_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return Teuchos::rcp(new GatherSolution_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherTangent() const
   { return Teuchos::rcp(new GatherTangent_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherDomain() const
   { return Teuchos::rcp(new GatherSolution_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   { return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT>(gidProvider_)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return Teuchos::rcp(new ScatterDirichletResidual_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(gidProvider_)); }

/*************** From ThyraObjFactory *******************/

   //! Get the domain space
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraDomainSpace() const;

   //! Get the range space
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraRangeSpace() const;

   //! Get a matrix operator
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > getThyraMatrix() const;

/*************** Tpetra construction functions *******************/

   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedTpetraVector() const;
   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedTpetraColVector() const;
   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getTpetraVector() const;
   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getTpetraColVector() const;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getTpetraMatrix() const;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedTpetraMatrix() const;

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
   void initializeContainer(int mem,TpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & loc) const;

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
   void initializeGhostedContainer(int mem,TpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & loc) const;

/*************** Tpetra based methods *******************/

   //! get the map from the matrix
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > getMap() const;
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > getColMap() const;

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedMap() const;
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedColMap() const;

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGraph() const;

   //! get the ghosted graph of the crs matrix
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedGraph() const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedImport() const;
   virtual const Teuchos::RCP<Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedColImport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedExport() const;
   virtual const Teuchos::RCP<Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedColExport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<const Teuchos::Comm<int> > getTeuchosComm() const;

   //! Get the domain global indexer this factory was created with.
   Teuchos::RCP<const panzer::GlobalIndexer> getDomainGlobalIndexer() const
   { return gidProvider_; }

   //! Get the domain global indexer this factory was created with.
   Teuchos::RCP<const panzer::GlobalIndexer> getRangeGlobalIndexer() const
   { return gidProvider_; }

   virtual void beginFill(LinearObjContainer & loc) const;
   virtual void endFill(LinearObjContainer & loc) const;

protected:

   void ghostToGlobalTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                                  Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out, bool col) const;
   void ghostToGlobalTpetraMatrix(const Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                                  Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out) const;
   void globalToGhostTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>& in,
                                  Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out, bool col) const;

   // get the map from the matrix
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildMap() const;
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildColMap() const;
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildGhostedMap() const;
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildGhostedColMap() const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildGraph() const;
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildGhostedGraph() const;

   // storage for Tpetra graphs and maps
   Teuchos::RCP<const Teuchos::Comm<int> > comm_;
   mutable Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > map_;
   mutable Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > cMap_;
   mutable Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > ghostedMap_;
   mutable Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > cGhostedMap_;
   mutable Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > graph_;
   mutable Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > ghostedGraph_;
   mutable Teuchos::RCP<ImportType> ghostedImporter_;
   mutable Teuchos::RCP<ImportType> ghostedColImporter_;
   mutable Teuchos::RCP<ExportType> ghostedExporter_;
   mutable Teuchos::RCP<ExportType> ghostedColExporter_;

   Teuchos::RCP<const GlobalIndexer> gidProvider_;
   Teuchos::RCP<const GlobalIndexer> colGidProvider_;

   bool hasColProvider_;

   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace_;
};

}

#endif
