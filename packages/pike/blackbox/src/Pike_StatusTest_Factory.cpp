#include "Pike_StatusTest_Factory.hpp"
#include "Teuchos_ParameterList.hpp"

// Status Tests
#include "Pike_StatusTest_Composite.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"

namespace pike {

  Teuchos::RCP<pike::StatusTest> 
  buildStatusTests(const Teuchos::RCP<Teuchos::ParameterList>& p)
  {
    TEUCHOS_ASSERT(p->numParams() == 1);
    TEUCHOS_TEST_FOR_EXCEPTION(!p->begin()->second.isList(),
			       std::logic_error,
			       "The ParameterList for building the Status Test with key \"" 
			       << p->begin()->first << "\" must be a sublist!");

    Teuchos::ParameterList& sublistRef = p->begin()->second.getValue(p.get());
    Teuchos::RCP<Teuchos::ParameterList> sublist = Teuchos::sublist(p,sublistRef.name(),true);
    
    Teuchos::RCP<pike::StatusTest> test;

    if (sublist->name() == "Composite") {
      Teuchos::RCP<pike::Composite> c = pike::composite();
      c->setParameterList(sublist);
      test = c;            
    }
    else if (sublist->name() == "Max Iterations") {
      Teuchos::RCP<pike::MaxIterations> mi = Teuchos::rcp(new pike::MaxIterations());
      mi->setParameterList(sublist);
      test = mi;
    }
    else if (sublist->name() == "Scalar Response Relative Tolerance") {
      Teuchos::RCP<pike::ScalarResponseRelativeTolerance> rt = 
	Teuchos::rcp(new pike::ScalarResponseRelativeTolerance());
      rt->setParameterList(sublist);
      test = rt;
    }
    else {
      
    }
    
    return test;
  }

}
