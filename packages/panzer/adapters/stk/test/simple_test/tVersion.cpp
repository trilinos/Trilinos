#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"

namespace panzer_stk {

// triangle tests
TEUCHOS_UNIT_TEST(tVersion, version)
{
  using Teuchos::RCP;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  if (comm->getRank() == 0) 
    out << panzer_stk::version() << std::endl;

  out << "Process " << comm->getRank() << " of " << comm->getSize() 
       << " is alive!" << std::endl; 
}

}
