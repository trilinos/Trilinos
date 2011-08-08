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

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include "Stokhos_OrthogPolyBasis.hpp"

#include <vector>

namespace panzer {

/** Linear object container for SG-Epetra objects.
  */
class SGEpetraLinearObjContainer : public LinearObjContainer {
public:
   typedef std::vector<Teuchos::RCP<EpetraLinearObjContainer> > CoeffVector;
   typedef CoeffVector::iterator iterator;
   typedef CoeffVector::const_iterator const_iterator;

   SGEpetraLinearObjContainer(const CoeffVector & coeffs,
                              const Teuchos::RCP<Stokhos::OrthogPolyBasis<int,double> > & basis);

   virtual void initialize();

   CoeffVector::iterator begin() { return coeffs_.begin(); }
   CoeffVector::iterator end() { return coeffs_.end(); }

   CoeffVector::const_iterator begin() const { return coeffs_.begin(); }
   CoeffVector::const_iterator end() const { return coeffs_.end(); }

   Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > getBasis() const
   { return basis_; }

private:
   CoeffVector coeffs_;
   Teuchos::RCP<Stokhos::OrthogPolyBasis<int,double> > basis_;
};

/** Linear object factory for constructing Stochastic Galerkin epetra
  * linear object containers. Also handles some of the global to ghosting (and vice-versa)
  * communication.
  */
template <typename Traits,typename LocalOrdinalT>
class SGEpetraLinearObjFactory : public LinearObjFactory<Traits> {
public:

   SGEpetraLinearObjFactory(const Teuchos::RCP<EpetraLinearObjFactory<Traits,LocalOrdinalT> > & epetraFact,
                            const Teuchos::RCP<Stokhos::OrthogPolyBasis<int,double> > & basis);

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

protected:
   Teuchos::RCP<EpetraLinearObjFactory<Traits,LocalOrdinalT> > epetraFact_;
   Teuchos::RCP<Stokhos::OrthogPolyBasis<int,double> > basis_;
};

}

#include "Panzer_SGEpetraLinearObjFactoryT.hpp"

#endif
#endif
