#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

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
    // Start of debug naming testing
    // *********************************************************************
    {
      cout << "\nTesting Debug naming scheme.\n";
      // Default
      cout << "Next line should be a WARNING for an undefined type:" << endl;
      cout << getTypeString<int, MyTraits>() << endl;
      // Scalar Types
      cout << getTypeString<double, MyTraits>() << endl;
      cout << getTypeString<MyTraits::FadType, MyTraits>() << endl;
      // Data Types
      cout << getTypeString<double, MyTraits>() << endl;
      cout << getTypeString<MyTraits::FadType, MyTraits>() << endl;
      cout << getTypeString<MyVector<double>, MyTraits>() << endl;
      cout << getTypeString<MyVector<MyTraits::FadType>, MyTraits>() << endl;
      cout << getTypeString<MyTensor<double>, MyTraits>() << endl;
      cout << getTypeString<MyTensor<MyTraits::FadType>, MyTraits>() << endl;
    }
    
    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    cout << "\nRun has completed successfully!\n" << endl; 
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
