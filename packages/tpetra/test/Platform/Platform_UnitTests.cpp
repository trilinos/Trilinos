#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_SerialPlatform.hpp"
#ifdef HAVE_TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#endif

#include <Kokkos_SerialNode.hpp>

namespace {

  using Teuchos::RCP;
  using Teuchos::Comm;
  using Teuchos::rcp;
  using Tpetra::SerialPlatform;
#ifdef HAVE_TPETRA_MPI
  using Tpetra::MpiPlatform;
#endif
  using Kokkos::SerialNode;

  template <class PLAT>
  RCP<PLAT> getPlatform() {
    TEST_FOR_EXCEPTION(true,std::logic_error,"Platform type " << Teuchos::TypeNameTraits<PLAT>::name() << " not defined.");
  }

  RCP<SerialNode> snode;

  template <>
  RCP<SerialPlatform<SerialNode> > getPlatform() {
    if (snode == Teuchos::null) {
      Teuchos::ParameterList pl;
      snode = rcp(new SerialNode(pl));
    }
    return rcp(new SerialPlatform<SerialNode>(snode));
  }

#ifdef HAVE_TPETRA_MPI
  template <>
  RCP<MpiPlatform<SerialNode> > getPlatform() {
    if (snode == Teuchos::null) {
      Teuchos::ParameterList pl;
      snode = rcp(new SerialNode(pl));
    }
    return rcp(new MpiPlatform<SerialNode>(snode));
  }
#endif

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Platform, basic, PlatformType )
  {
    out << "Testing " << Teuchos::TypeNameTraits<PlatformType>::name() << std::endl;
    typedef typename PlatformType::NodeType            N;
    // create a platform  
    RCP<PlatformType> platform = getPlatform<PlatformType>();
    platform->setObjectLabel("not the default label");
    // get the comm for this platform
    RCP<const Comm<int> > comm = platform->getComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    TEST_EQUALITY( myImageID < numImages, true );
    TEST_EQUALITY_CONST( comm != Teuchos::null, true );
    RCP<N> node  = platform->getNode();
    (void)node;
  }


  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

typedef SerialPlatform<SerialNode> SP;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Platform, basic, SP)
#ifdef HAVE_TPETRA_MPI
typedef MpiPlatform<SerialNode> MP;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Platform, basic, MP )
#endif

}
