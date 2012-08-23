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

#ifndef __Panzer_TpetraLinearObjFactory_decl_hpp__
#define __Panzer_TpetraLinearObjFactory_decl_hpp__

#include <map>

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"

#include "Panzer_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_ScatterResidual_Tpetra.hpp"
#include "Panzer_ScatterDirichletResidual_Tpetra.hpp"
#include "Panzer_ScatterInitialCondition_Tpetra.hpp"
#include "Panzer_GatherSolution_Tpetra.hpp"
#include "Panzer_GatherOrientation.hpp"
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_ThyraObjFactory.hpp"

#include"Kokkos_DefaultNode.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

namespace panzer {

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=Kokkos::DefaultNode::DefaultNodeType>
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
                          const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & gidProvider);

   virtual ~TpetraLinearObjFactory();

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
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   { return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT>(gidProvider_)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return Teuchos::rcp(new ScatterDirichletResidual_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(gidProvider_)); }

   //! Use preconstructed initial condition scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterInitialCondition() const
   { return Teuchos::rcp(new ScatterInitialCondition_Tpetra<EvalT,Traits,LocalOrdinalT,GlobalOrdinalT,NodeT>(gidProvider_)); }

/*************** From ThyraObjFactory *******************/

   //! Get the domain space
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraDomainSpace() const;

   //! Get the range space
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraRangeSpace() const;

   //! Get a matrix operator
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > getThyraMatrix() const;

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

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedMap() const;

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGraph() const;

   //! get the ghosted graph of the crs matrix
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedGraph() const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedImport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedExport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<const Teuchos::Comm<int> > getTeuchosComm() const;

protected:
   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedTpetraVector() const;
   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getTpetraVector() const;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getTpetraMatrix() const;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > getGhostedTpetraMatrix() const;

   void ghostToGlobalTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                                  Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out) const;
   void ghostToGlobalTpetraMatrix(const Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & in,
                                  Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out) const;
   void globalToGhostTpetraVector(const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>& in,
                                  Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> & out) const;

   // get the map from the matrix
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildMap() const;
   virtual const Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildGhostedMap() const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildGraph() const;
   virtual const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > buildGhostedGraph() const;

   // storage for Tpetra graphs and maps
   Teuchos::RCP<const Teuchos::Comm<int> > comm_;
   mutable Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > map_;
   mutable Teuchos::RCP<Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > ghostedMap_;
   mutable Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > graph_;
   mutable Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> > ghostedGraph_;

   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > gidProvider_;

   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace_;
};

}

#endif
