#ifndef __Panzer_Response_Functional_hpp__
#define __Panzer_Response_Functional_hpp__

#include <string>


#include <mpi.h> // need for comm

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

#include "Panzer_ResponseMESupport_Default.hpp"
#include "Panzer_GlobalEvaluationData.hpp"
#include "Panzer_ThyraObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_LinearObjFactory.hpp"


namespace panzer {

/** This class provides a response storage for
  * simple functionals of the solution (i.e. scalar
  * values).
  */
template <typename EvalT>
class Response_Functional : public ResponseMESupport_Default<EvalT> {
public:
   typedef typename EvalT::ScalarT ScalarT;

   Response_Functional(const std::string & responseName,MPI_Comm comm,const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & linObjFact=Teuchos::null)
     : ResponseMESupport_Default<EvalT>(responseName,comm), value(0.0), linObjFactory_(linObjFact)
   {
     if(linObjFactory_!=Teuchos::null) {
       // requires thyra object factory
       thyraObjFactory_ = Teuchos::rcp_dynamic_cast<const panzer::ThyraObjFactory<double> >(linObjFactory_,true);
       setSolnVectorSpace(thyraObjFactory_->getThyraDomainSpace());

       // build a ghosted container, with a solution vector
       uniqueContainer_ = linObjFactory_->buildLinearObjContainer();
       ghostedContainer_ = linObjFactory_->buildGhostedLinearObjContainer();

       // set ghosted container (work space for assembly)
       linObjFactory_->initializeGhostedContainer(panzer::LinearObjContainer::F,*ghostedContainer_);
     }
   }

   //! provide direct access, this thing is pretty simple
   ScalarT value;

   //! This simply does global summation, then shoves the result into a vector
   virtual void scatterResponse();

   virtual void initializeResponse()  
   { value = 0.0; }

   // from ResponseMESupport_Default

   //! What is the number of values you need locally
   virtual std::size_t localSizeRequired() const { return 1; }

   //! Is the vector distributed (or replicated)
   virtual bool vectorIsDistributed() const { return false; }

   //! Get ghosted responses (this will be filled by the evaluator)
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const
   { return Teuchos::rcp_dynamic_cast<const ThyraObjContainer<double> >(ghostedContainer_)->get_f_th(); }
    
private:
   //! Set solution vector space
   void setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs);

   // hide these methods
   Response_Functional();
   Response_Functional(const Response_Functional &);

   Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;
   Teuchos::RCP<const panzer::ThyraObjFactory<double> > thyraObjFactory_;

   Teuchos::RCP<LinearObjContainer> uniqueContainer_;
   Teuchos::RCP<LinearObjContainer> ghostedContainer_;
};

}

#include "Panzer_Response_Functional_impl.hpp"

#endif
