#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::Comm;
  using Tpetra::DefaultPlatform;
  using Teuchos::rcp;

  bool testMpi = true;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST( Platform, basic )
  {
    typedef DefaultPlatform::DefaultPlatformType::NodeType Node;
    // create a platform  
    DefaultPlatform::DefaultPlatformType &platform = DefaultPlatform::getDefaultPlatform();
    platform.setObjectLabel("not the default label");
    // get the comm for this platform
    RCP<Comm<int> > comm = platform.getComm();
    TEST_EQUALITY_CONST( comm != Teuchos::null, true );
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    out << "Default node type: " << Teuchos::TypeNameTraits<Node>::name() << std::endl;
    out << "Default platform type: " << Teuchos::typeName(platform) << std::endl;
    out << "numImages: " << numImages << std::endl;
    out << "myImageID: " << myImageID << std::endl;
    TEST_EQUALITY_CONST( myImageID < numImages, true );
    RCP<Node> node = platform.getNode();
    (void)node;
  }

}
