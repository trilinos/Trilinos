#ifndef __Panzer_SGEpetraLinearObjFactory_hpp__
#define __Panzer_SGEpetraLinearObjFactory_hpp__

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
                            const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion);

   virtual ~SGEpetraLinearObjFactory();

/*************** Linear object factory methods *******************/

   virtual Teuchos::RCP<LinearObjContainer> buildLinearObjContainer() const;

   virtual Teuchos::RCP<LinearObjContainer> buildGhostedLinearObjContainer() const;

   virtual void globalToGhostContainer(const LinearObjContainer & container,
                                       LinearObjContainer & ghostContainer) const;
   virtual void ghostToGlobalContainer(const LinearObjContainer & ghostContainer,
                                       LinearObjContainer & container) const;

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

protected:
   Teuchos::RCP<EpetraLinearObjFactory<Traits,LocalOrdinalT> > epetraFact_;
   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion_;
};

}

#include "Panzer_SGEpetraLinearObjFactoryT.hpp"

#endif
#endif
