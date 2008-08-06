#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

/*! \brief Test to check performance of bracket operator

    The Field class uses a bracket operator to access data elements.
    The operation is inlined and compilers will optimize this function
    call away so that it should be as fast as raw access.  This test
    will allow a comparison against raw access to verify that your
    compiler is at an acceptable optimization level.

    This test shows that the fields are a factor of 2 slower if
    running at -O0 instead of -O3 on linux gnu 4.2.4 compilers.  For
    the high optimizaiton level, there is virtually no difference
    in runtimes between the field and the raw data.

*/
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    RCP<Time> arcp_time = TimeMonitor::getNewTimer("ArrayRCP<double> Time");
    RCP<Time> raw_time = TimeMonitor::getNewTimer("double* Time");
    
    TimeMonitor tm(*total_time);
    
    {
      
      const int num_loops = 1000;
      const int size = 200000;
      double value = 2.0;
      
      double* raw_density = new double[size];
      
      RCP<DataLayout> dl = 
	rcp(new Generic("cell_quantitiy", 1));
      Field<double> density("density", dl);
      ArrayRCP<double> a_density = arcp<double>(size);
      density.setFieldData(a_density);
      
      {
	TimeMonitor tm(*arcp_time);
	for (int i=0; i < num_loops; ++i)
	  for (int j=0; j < size; ++j)
	    density[j] = value;
      }
      
      {
	TimeMonitor tm(*raw_time);
	for (int i=0; i < num_loops; ++i)
	  for (int j=0; j < size; ++j)
	    raw_density[j] = value;
      }
      
      delete [] raw_density;
    }

    // *********************************************************************
    // *********************************************************************
    cout << "\nTest passed!\n" << endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const exception& e) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Exception Caught!" << endl;
    cout << "Error message is below\n " << e.what() << endl;
    cout << "************************************************" << endl;
  }
  catch (...) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Unknown Exception Caught!" << endl;
    cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
