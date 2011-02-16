#ifndef PANZER_STK_MODEL_EVALUATOR_FACTORY_HPP
#define PANZER_STK_MODEL_EVALUATOR_FACTORY_HPP

#include "Panzer_ConfigDefs.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include <string>

namespace Thyra {
  template<typename ScalarT> class ModelEvaluator;
}

namespace panzer_stk {
  
  template<typename ScalarT>
  class ModelEvaluatorFactory_Epetra : public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    /** @name Overridden from ParameterListAcceptor */
    //@{
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    //@}

    void buildObjects(const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getPhysicsModelEvaluator();
    
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > getResponseOnlyModelEvaluator();
    
  private:

    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_physics_me;
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > m_rome_me;
  };

}

#include "Panzer_STK_ModelEvaluatorFactory_EpetraT.hpp"

#endif
