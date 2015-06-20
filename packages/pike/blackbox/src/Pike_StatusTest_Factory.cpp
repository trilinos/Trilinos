#include "Pike_StatusTest_Factory.hpp"
#include "Teuchos_ParameterList.hpp"

// Status Tests
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_GlobalModelConvergence.hpp"
#include "Pike_StatusTest_LocalModelConvergence.hpp"
#include "Pike_StatusTest_LocalModelFailure.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

namespace pike {

  StatusTestFactory::StatusTestFactory()
  {
    myTypes_.push_back("Composite AND");
    myTypes_.push_back("Composite OR");
    myTypes_.push_back("Maximum Iterations");
    myTypes_.push_back("Global Model Convergence");
    myTypes_.push_back("Local Model Convergence");
    myTypes_.push_back("Local Model Failure");
    myTypes_.push_back("Scalar Response Relative Tolerance");
  }

  bool StatusTestFactory::supportsType(const std::string& type) const
  {
    return (std::find(myTypes_.begin(),myTypes_.end(),type) != myTypes_.end());
  }

  Teuchos::RCP<pike::StatusTest> 
  StatusTestFactory::buildStatusTests(const Teuchos::RCP<Teuchos::ParameterList>& p) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(!p->isType<std::string>("Type"),std::logic_error,
			       "The StatusTestFactory requires that each sublist declare its \"Type\".  This parameter does not exist for the sublist \"" 
			       << p->name() << "\".");

    std::string testType = p->get<std::string>("Type");

    Teuchos::RCP<pike::StatusTest> test;

    if ( (testType == "Composite AND") || (testType == "Composite OR") ) {
      Teuchos::RCP<pike::Composite> c = pike::composite();
      c->setParameterList(p,*this);
      test = c;
    }
    else if (testType == "Maximum Iterations") {
      Teuchos::RCP<pike::MaxIterations> mi = Teuchos::rcp(new pike::MaxIterations);
      mi->setParameterList(p);
      test = mi;
    }
    else if (testType == "Global Model Convergence") {
      Teuchos::RCP<pike::GlobalModelConvergence> gmc = Teuchos::rcp(new pike::GlobalModelConvergence);
      gmc->setParameterList(p);
      test = gmc;
    }
    else if (testType == "Local Model Convergence") {
      Teuchos::RCP<pike::LocalModelConvergence> lmc = Teuchos::rcp(new pike::LocalModelConvergence);
      lmc->setParameterList(p);
      test = lmc;
    }
    else if (testType == "Local Model Failure") {
      Teuchos::RCP<pike::LocalModelFailure> lmf = Teuchos::rcp(new pike::LocalModelFailure);
      lmf->setParameterList(p);
      test = lmf;
    }
    else if (testType == "Scalar Response Relative Tolerance") {
      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> rt = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance());
      rt->setParameterList(p);
      test = rt;
    }
    else {
      typedef std::vector<Teuchos::RCP<pike::StatusTestAbstractFactory> >::const_iterator it;
      for (it f=userFactories_.begin(); f != userFactories_.end(); ++f) {
	if ( (*f)->supportsType(testType) )
	  test = (*f)->buildStatusTests(p);
      }
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(test), std::logic_error,
			       "The status test factory failed to build a StatusTest of \"Type\" = \""
			       << testType <<"\"!");
    return test;
  }

  void StatusTestFactory::addFactory(const Teuchos::RCP<pike::StatusTestAbstractFactory>& f)
  {
    userFactories_.push_back(f);
  }

}
