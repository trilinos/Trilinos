#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace Teuchos;
using std::string;

// declare Teuchos Time pointers here as global

RefCountPtr<Time> SetupTime;
RefCountPtr<Time> CompTime;

// functions declaration

void Setup();
void Comp1();
void Comp2();

// main driver

int main(int argc, char* argv[])
{
  try {
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

    // Create time monitors for each of the phases we want to
    // monitor. Each phase will have its own Time pointer
    SetupTime = TimeMonitor::getNewTimer("setup time");
    CompTime  = TimeMonitor::getNewTimer("comp time");

    Setup();

    Comp1();

    Comp2();
    
    // Summarize timings. This must be done before finalizing MPI  
    TimeMonitor::summarize();

  } catch(std::exception& e) {
    cerr << "caught exception " << e.what() << endl;
  }

  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}

// =======================================================================

void Setup()
{
  // construct a timer
  TimeMonitor LocalTimer(*SetupTime);

  // do something
  double sum = 0.0;
  
  for (int i=0; i<10000; i++) {
    sum += sqrt(sum);
  }

  // LocalTimer is destroyed here
  return;
}

// ========================================================================
void Comp1() 
{
  TimeMonitor LocalTimer(*CompTime);

  // do something
  double a = 0;
  for( int i=0 ; i<100000 ; ++i ) {
    a = pow(10,0.234);
  }
  
}

// ========================================================================

void Comp2() 
{
  TimeMonitor LocalTimer(*CompTime);

  // do something
  double a = 0;
  for( int i=0 ; i<100000 ; ++i ) {
    a = pow(10,0.444);
  }
  
}
