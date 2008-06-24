#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "AlgebraicTypes.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    {
      cout << "Vector Testing: a, b are vectors, c is scalar" << endl;
      MyVector<double> a;
      a.init(3.0);
      cout << "Printing a:\n" << a << endl;
      MyVector<double> b(2.0);
      cout << "Printing b:\n" << b << endl;
      double c = 4.0;
      cout << "Printing c: " << c << endl;
      cout << "Printing -a:\n" << -a << endl;
      cout << "Printing +a:\n" << +a << endl;
      cout << "Printing a+b:\n" << a+b << endl;
      cout << "Printing a-b:\n" << a-b << endl;
      cout << "Printing a*b:\n" << a*b << endl;
      cout << "Printing a/b:\n" << a/b << endl;
      cout << "Printing c*a:\n" << c*a << endl;;
      cout << "Printing a*c:\n" << a*c << endl;
    }

    // *********************************************************************
    // Finished all testing
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
