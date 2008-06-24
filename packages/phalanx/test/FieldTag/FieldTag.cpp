#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"

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
    // Start of Field Tag Testing
    // *********************************************************************
    {

      // Dummy data layouts
      RCP<DataLayout> node4 = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_Nodes", 4));
      RCP<DataLayout> quad4 = 
	rcp(new Generic<MyTraits::MY_SCALAR>("Q1_QuadPoints", 4));
      RCP<DataLayout> gradQuad4 = 
	rcp(new Generic<MyTraits::MY_VECTOR>("Q1_QuadPoints", 4));
      
      // Tags with same name but different topology
      FieldTag nodal_density("density", node4);
      FieldTag qp_density("density", quad4);
      FieldTag grad_qp_density("density", gradQuad4);
      
      // test ostream
      cout << "Printing field tags" << endl;
      cout << nodal_density << endl;
      cout << qp_density << endl;
      cout << grad_qp_density << endl;
      cout << endl;
      
      // test operator ==
      cout << "Are nodal and qp fields equal (should be false)? = " 
	   << (nodal_density == qp_density) << endl;
      TEST_FOR_EXCEPTION(nodal_density == qp_density, std::logic_error,
			 "operator==() failed!");
      
      // New constructor that should be same as nodal_density
      FieldTag nodal_density_copy("density", node4);

      cout << "Are nodal and nodal copy fields equal (should be true)? = " 
	   << (nodal_density == nodal_density_copy) << endl;
      TEST_FOR_EXCEPTION(!(nodal_density == nodal_density_copy), 
			 std::logic_error,
			 "operator==() failed for unique copy comparison!");
      
      cout << "Are scalar and vector fields "
	   << "equal (should be false)? = " 
	   << (qp_density == grad_qp_density) << endl;
      TEST_FOR_EXCEPTION(qp_density == grad_qp_density, 
			 std::logic_error,
			 "operator==() failed for data layout comparison !");
      
      // test operator =
      FieldTag copy = nodal_density_copy;
      
      cout << "Compare copy of operator to itself (should be true)? = " 
	   << (copy == nodal_density_copy) << endl;
      TEST_FOR_EXCEPTION(!(copy == nodal_density_copy), 
			 std::logic_error,
			 "operator=() failed!");

      // name() accessor
      cout << "Testing name() accessor...";
      TEST_FOR_EXCEPTION( (nodal_density.name() != std::string("density") ), 
			 std::logic_error,
			 "name() accessor failed!");
      cout << "Passed." << endl;
      
      // dataLayout() accessor
      cout << "Testing dataLayout() accessor...";
      TEST_FOR_EXCEPTION(nodal_density.dataLayout() != node4, 
			 std::logic_error,
			 "dataLayout() accessor failed!");
      cout << "Passed." << endl;
      
    
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
