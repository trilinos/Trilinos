#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include "Tpetra_Map_def.hpp"
#include "Tpetra_Directory_def.hpp"
#endif

namespace {

  using Teuchos::RCP;
  using Teuchos::Array;
  using Tpetra::global_size_t;
  using Teuchos::outArg;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    return Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }

  RCP<Tpetra::DefaultPlatform::DefaultPlatformType::NodeType> getDefaultNode()
  {
    return Tpetra::DefaultPlatform::getDefaultPlatform().getNode();
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST( Map, Bug5401_NegativeBaseIndex )
  {
    // failure reading 1x4 matrix under MPI
    typedef int                       LO;
    typedef int                       GO;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
    typedef Tpetra::Map<LO,GO,Node>   Map;
    typedef Teuchos::Comm<int>        Comm;
    // create a comm  
    RCP<const Comm> comm = getDefaultComm();
    const int numImages = comm->getSize();
    TEUCHOS_TEST_FOR_EXCEPTION( numImages != 2, std::logic_error, "This test is appropriate only for MPI runs of rank 2.")
    RCP<Node>             node = getDefaultNode();

    const GO numElements = 78;
    const GO baseIndexIsNegOne = -1;
    const global_size_t GINV   = Teuchos::OrdinalTraits<global_size_t>::invalid();
    Array<int> elements(78);

    // first global element is -1
    for (int i = 0; i < elements.size(); ++i) elements[i] = i - 1;

    RCP<Map> map = rcp(new Map(GINV, elements(), baseIndexIsNegOne, comm));

    TEST_EQUALITY( map->getNodeNumElements(),   numElements );
    TEST_EQUALITY( map->getGlobalNumElements(), numElements*numImages );
    TEST_EQUALITY( map->getIndexBase(), -1 );
    TEST_EQUALITY( map->getMinGlobalIndex(),     -1 );
    TEST_EQUALITY( map->getMinAllGlobalIndex(),  -1 );

    // All procs fail if any proc fails 
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

}
