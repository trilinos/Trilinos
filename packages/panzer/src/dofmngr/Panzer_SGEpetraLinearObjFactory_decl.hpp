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

#ifndef __Panzer_SGEpetraLinearObjFactory_decl_hpp__
#define __Panzer_SGEpetraLinearObjFactory_decl_hpp__

#include "Panzer_config.hpp"
#ifdef HAVE_STOKHOS

#include <map>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "Panzer_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_SGEpetraLinearObjContainer.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraVectorOrthogPoly.hpp"

#include <vector>

namespace panzer {

/** Linear object factory for constructing Stochastic Galerkin epetra
  * linear object containers. Also handles some of the global to ghosting (and vice-versa)
  * communication.
  */
template <typename Traits,typename LocalOrdinalT>
class SGEpetraLinearObjFactory : public LinearObjFactory<Traits> {
public:

   SGEpetraLinearObjFactory(const Teuchos::RCP<EpetraLinearObjFactory<Traits,LocalOrdinalT> > & epetraFact,
                            const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
                            const Teuchos::RCP<const EpetraExt::MultiComm> & globalMultiComm);

   virtual ~SGEpetraLinearObjFactory();

/*************** Linear object factory methods *******************/

   virtual Teuchos::RCP<LinearObjContainer> buildLinearObjContainer() const;
   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveLinearObjContainer() const;

   virtual Teuchos::RCP<LinearObjContainer> buildGhostedLinearObjContainer() const;
   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveGhostedLinearObjContainer() const;

   virtual void globalToGhostContainer(const LinearObjContainer & container,
                                       LinearObjContainer & ghostContainer,int mem) const;
   virtual void ghostToGlobalContainer(const LinearObjContainer & ghostContainer,
                                       LinearObjContainer & container,int mem) const;

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
   { return epetraFact_->template buildScatter<EvalT>(); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return epetraFact_->template buildGather<EvalT>(); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   { return epetraFact_->template buildGatherOrientation<EvalT>(); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return epetraFact_->template buildScatterDirichlet<EvalT>(); }

   //! Use preconstructed initial condition scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterInitialCondition() const
   { return epetraFact_->template buildScatterInitialCondition<EvalT>(); }

/*************** Generic helper functions for container setup *******************/
   
   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeContainer(int mem,LinearObjContainer & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeGhostedContainer(int mem,LinearObjContainer & loc) const;

   /** Extract underlying epetra factory from SGEpetraLinearOpFactory
     */
   Teuchos::RCP<EpetraLinearObjFactory<Traits,LocalOrdinalT> > getEpetraFactory() const
   { return epetraFact_; }

   //! Accessor for the expansion object.
   Teuchos::RCP<const Stokhos::OrthogPolyExpansion<int,double> > getExpansion() const
   { return expansion_; }

   //! Accessor for the expansion object.
   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > getExpansion()
   { return expansion_; }

   //! Set orthog poly object, this serves as a template for converting vectors to block vectors
   Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> getVectorOrthogPoly() const;

   //! get the map from the matrix, this is the map for the solution vector
   Teuchos::RCP<const Epetra_Map> getMap();

   //! get the block map needed by Stokhos to describe the parallel layout of the SG unknowns
   Teuchos::RCP<const Epetra_Map> getSGBlockMap() const;

protected:

   Teuchos::RCP<EpetraLinearObjFactory<Traits,LocalOrdinalT> > epetraFact_;
   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion_;
   Teuchos::RCP<const EpetraExt::MultiComm> globalMultiComm_;

   mutable Teuchos::RCP<const Epetra_Map> sgBlockMap_; // constructed via lazy evaluation
};

}

#endif
#endif
