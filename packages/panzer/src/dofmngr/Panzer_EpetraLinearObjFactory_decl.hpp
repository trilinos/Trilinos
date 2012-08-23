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

#include <map>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "Panzer_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_ScatterResidual_Epetra.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra.hpp"
#include "Panzer_ScatterInitialCondition_Epetra.hpp"
#include "Panzer_GatherSolution_Epetra.hpp"
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

   EpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                          const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider);
   EpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                          const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > & gidProvider);

   virtual ~EpetraLinearObjFactory();

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

   //! Get the domain space
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraDomainSpace() const;

   //! Get the range space
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraRangeSpace() const;

   //! Get a matrix operator
   Teuchos::RCP<Thyra::LinearOpBase<double> > getThyraMatrix() const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   { return Teuchos::rcp(new ScatterResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return Teuchos::rcp(new GatherSolution_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   { return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return Teuchos::rcp(new ScatterDirichletResidual_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

   //! Use preconstructed initial condition scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterInitialCondition() const
   { return Teuchos::rcp(new ScatterInitialCondition_Epetra<EvalT,Traits,LocalOrdinalT,int>(gidProvider_)); }

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

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedMap() const;

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGraph() const;

   //! get the ghosted graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGhostedGraph() const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getGhostedImport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getGhostedExport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<const Epetra_Comm> getEpetraComm() const;

protected:
   Teuchos::RCP<Epetra_Vector> getGhostedEpetraVector() const;
   Teuchos::RCP<Epetra_Vector> getEpetraVector() const;
   Teuchos::RCP<Epetra_CrsMatrix> getEpetraMatrix() const;
   Teuchos::RCP<Epetra_CrsMatrix> getGhostedEpetraMatrix() const;

   void ghostToGlobalEpetraVector(const Epetra_Vector & in,Epetra_Vector & out) const;
   void ghostToGlobalEpetraMatrix(const Epetra_CrsMatrix & in,Epetra_CrsMatrix & out) const;
   void globalToGhostEpetraVector(const Epetra_Vector & in,Epetra_Vector & out) const;

   // get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> buildMap() const;
   virtual const Teuchos::RCP<Epetra_Map> buildGhostedMap() const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGraph() const;
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGhostedGraph() const;

   // storage for Epetra graphs and maps
   Teuchos::RCP<const Epetra_Comm> comm_;
   mutable Teuchos::RCP<Epetra_Map> map_;
   mutable Teuchos::RCP<Epetra_Map> ghostedMap_;
   mutable Teuchos::RCP<Epetra_CrsGraph> graph_;
   mutable Teuchos::RCP<Epetra_CrsGraph> ghostedGraph_;

   // import/exporter storage
   mutable Teuchos::RCP<Epetra_Import> importer_;
   mutable Teuchos::RCP<Epetra_Export> exporter_;

   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > gidProvider_;

   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace_;
};

}

#endif
