
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"

#include <iostream>

int main( int argc, char **argv )
{  
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  using std::cout;
  using std::endl;
  using Teuchos::RCP;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  if (comm->getRank() == 0) 
    cout << panzer_stk::version() << endl;

  cout << "Process " << comm->getRank() << " of " << comm->getSize() 
       << " is alive!" << endl; 

  return 0;
}
