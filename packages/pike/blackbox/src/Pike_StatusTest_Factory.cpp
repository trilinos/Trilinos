#include "Pike_StatusTest_Factory.hpp"
#include "Teuchos_ParameterList.hpp"

// Status Tests
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_GlobalModelConvergence.hpp"
#include "Pike_StatusTest_LocalModelFailure.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

namespace pike {

  Teuchos::RCP<pike::StatusTest> 
  buildStatusTests(const Teuchos::RCP<Teuchos::ParameterList>& p)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!p->isType<std::string>("Type"),std::logic_error,
			       "The StatusTestFactory requires that each sublist declare its \"Type\".  This parameter does not exist for the sublist \"" << p->name() << "\".");
    std::string testType = p->get<std::string>("Type");

    Teuchos::RCP<pike::StatusTest> test;

    if ( (testType == "Composite AND") || (testType == "Composite OR") ) {
      Teuchos::RCP<pike::Composite> c = pike::composite();
      c->setParameterList(p);
      test = c;
    }
    else if (testType == "Maximum Iterations") {
      Teuchos::RCP<pike::MaxIterations> mi = Teuchos::rcp(new pike::MaxIterations);
      mi->setParameterList(p);
      test = mi;
    }
    else if (testType == "Global Model Convergence") {
      Teuchos::RCP<pike::GlobalModelConvergence> mc = Teuchos::rcp(new pike::GlobalModelConvergence);
      mc->setParameterList(p);
      test = mc;
    }
    else if (testType == "Local Model Failure") {
      Teuchos::RCP<pike::LocalModelFailure> mc = Teuchos::rcp(new pike::LocalModelFailure);
      mc->setParameterList(p);
      test = mc;
    }
    else if (testType == "Scalar Response Relative Tolerance") {
      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> rt = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance());
      rt->setParameterList(p);
      test = rt;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
				 "The status test factory failed to find a valid object type for the value \""
				 << testType <<"\"!");
    }
    
    return test;
  }

}
