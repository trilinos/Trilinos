#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Cthulhu_ConfigDefs.hpp" //TODO
#include "Cthulhu_DefaultPlatform.hpp" //TODO
#include "Teuchos_as.hpp"

//#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_Operator.hpp"
#include "Cthulhu_CrsOperator.hpp"
#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraCrsMatrix.hpp" //TMP
#endif
namespace {
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  using Cthulhu::DefaultPlatform;
  using Cthulhu::Operator;
  using Cthulhu::CrsOperator;
#ifdef HAVE_CTHULHU_TPETRA
  using Cthulhu::TpetraCrsMatrix; //TMP
#endif
  using Cthulhu::Map;

  using Cthulhu::viewLabel_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  // 


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Operator, ViewSwitching, Scalar, LO, GO, Node ) //TODO: add template parameter <Node,...>
  {
#ifdef HAVE_CTHULHU_TPETRA
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Operator<Scalar, LO, GO, Node> Operator;
    typedef CrsOperator<Scalar, LO, GO, Node> CrsOperator;
    RCP<const Comm<int> > comm = getDefaultComm();

    const size_t numLocal = 10;
    const size_t INVALID = OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t
    RCP<const Map<LO,GO,Node> > map = Cthulhu::useTpetra::createContigMap<LO,GO>(INVALID,numLocal,comm);
     {
       TpetraCrsMatrix<Scalar, LO, GO, Node> t =  TpetraCrsMatrix<Scalar,LO,GO,Node>(map, numLocal);

       // Test of constructor
       CrsOperator op(map,1);
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.GetDefaultViewLabel());
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.SwitchToView(op.GetCurrentViewLabel()));

       // Test of CreateView
       TEST_THROW(op.CreateView(op.GetDefaultViewLabel(),op.getRowMap(),op.getColMap()), Cthulhu::Exceptions::RuntimeError); // a
       op.CreateView("newView",op.getRowMap(),op.getColMap());                                                               // b
       TEST_THROW(op.CreateView("newView",op.getRowMap(),op.getColMap()), Cthulhu::Exceptions::RuntimeError);                // c

       // Test of SwitchToView
       // a
       viewLabel_t viewLabel    = op.GetCurrentViewLabel();
       viewLabel_t oldViewLabel = op.SwitchToView("newView");
       TEST_EQUALITY_CONST(viewLabel, oldViewLabel);
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), "newView");
       // b
       TEST_THROW(op.SwitchToView("notAView"), Cthulhu::Exceptions::RuntimeError);
              
       // Test of SwitchToDefaultView()
       // a
       viewLabel    = op.GetCurrentViewLabel();
       oldViewLabel = op.SwitchToDefaultView();
       TEST_EQUALITY_CONST(viewLabel, oldViewLabel);
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.GetDefaultViewLabel());
       
       // Test of RemoveView()
       TEST_THROW(op.RemoveView(op.GetDefaultViewLabel()), Cthulhu::Exceptions::RuntimeError); // a
       TEST_THROW(op.RemoveView("notAView"), Cthulhu::Exceptions::RuntimeError);               // b
       op.RemoveView("newView");                                                               // c
       TEST_THROW(op.RemoveView("newView"), Cthulhu::Exceptions::RuntimeError);
       
       op.fillComplete();
     }
#endif
  }

  // 
  // INSTANTIATIONS
  //

#   define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Operator, ViewSwitching, SC, LO, GO, Node )
  
  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;
  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)

}
