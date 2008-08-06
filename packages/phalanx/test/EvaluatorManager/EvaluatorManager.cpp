#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_Evaluator_Manager.hpp"
#include "Phalanx_DebugStrings.hpp"

// Evaluators
#include "evaluators/Evaluator_Constant.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"

// From test/Utilities directory
#include "Traits.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Data Container Testing
    // *********************************************************************
    {
      cout << "\nStarting EvaluatorManager Testing\n";

      
      cout << "\nConstructing EvaluatorManager...";
      EvaluatorManager<MyTraits> em;
      cout << "Passed!" << endl;

      RCP<DataLayout> nodes = rcp(new Generic("nodes",4));
      RCP<DataLayout> qp = rcp(new Generic("QP",4));

      RCP<FieldTag> den_n = rcp(new Tag<double>("Density", nodes));
      RCP<FieldTag> den_qp = rcp(new Tag<double>("Density", qp));
      
      
      cout << "\nTesting requireField()...";
      em.requireField(*den_n);
      em.requireField(*den_qp);
      cout << "Passed!" << endl;

      cout << "\nTesting registerEvaluator()...";
      { 
	ParameterList p;
	p.set<string>("Name", "Density");
	p.set<double>("Value", 2.0);
	p.set< RCP<DataLayout> >("Data Layout", nodes);
	Teuchos::RCP< PHX::Evaluator<MyTraits> > ptr = 
	  rcp(new Constant<MyTraits::Residual, MyTraits>(p));
	em.registerEvaluator(ptr);
      }
      { 
	ParameterList p;
	p.set<string>("Name", "Density");
	p.set<double>("Value", 2.0);
	p.set< RCP<DataLayout> >("Data Layout", qp);
	Teuchos::RCP< PHX::Evaluator<MyTraits> > ptr = 
	  rcp(new Constant<MyTraits::Residual, MyTraits>(p));
	em.registerEvaluator(ptr);
      }
      cout << "Passed!" << endl;
      
      cout << "\nTesting setEvaluationTypeName()...";
      em.setEvaluationTypeName(PHX::getTypeString<MyTraits::Residual,MyTraits>());
      cout << "Passed!" << endl;


      cout << "\nTesting sortAndOrderEvaluators()...";
      em.sortAndOrderEvaluators();
      cout << "Passed!" << endl;

      cout << "\nPrinting EvaluatorManager:\n";
      cout << em << endl;

    }
    
    // *********************************************************************
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
