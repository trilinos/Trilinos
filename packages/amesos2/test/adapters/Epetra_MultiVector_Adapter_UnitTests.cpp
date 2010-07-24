#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_DefaultPlatform.hpp>

#include "Amesos2_EpetraMultiVecAdapter.hpp"
#include "Amesos2_Util_is_same.hpp"

namespace {

using std::cout;
using std::endl;
using std::string;

using Teuchos::as;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::tuple;
using Teuchos::ScalarTraits;
using Teuchos::OrdinalTraits;
using Teuchos::FancyOStream;
using Teuchos::VerboseObjectBase;

using Amesos::MultiVecAdapter;

using Amesos::Util::is_same;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

bool testMpi = false;

// Where to look for input files
string filedir;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("filedir",&filedir,"Directory of matrix files.");
  clp.addOutputSetupOptions(true);
  clp.setOption("test-mpi", "test-serial", &testMpi,
    "Test Serial by default (for now) or force MPI test.  In a serial build,"
    " this option is ignored and a serial comm is always used." );
}

const Epetra_Comm& getDefaultComm()
{
#ifdef EPETRA_MPI
  return Epetra_MpiComm( MPI_COMM_WORLD );
#else
  return Epetra_SerialComm();
#endif
}

RCP<FancyOStream> getDefaultOStream()
{
  return( VerboseObjectBase::getDefaultOStream() );
}

/*
 * UNIT TESTS
 */

TEUCHOS_UNIT_TEST( MultiVecAdapter, Initialization )
{
  /* Test correct initialization of the MultiVecAdapter
   *
   * - All Constructors
   * - Correct initialization of class members
   * - Correct typedefs ( using Amesos::is_same<> )
   */
  typedef ScalarTraits<double> ST;
  typedef Epetra_MultiVector MV;
  typedef MultiVecAdapter<MV> ADAPT;

  Epetra_SerialComm comm;
  // create a Map
  const size_t numLocal = 10;
  Epetra_Map map(-1,numLocal,0,comm);

  RCP<MV> mv = rcp(new MV(map,11));
  mv->Random();

  RCP<ADAPT> adapter = rcp(new ADAPT(mv));
  // Test copy constructor
  RCP<ADAPT> adapter2 = rcp(new ADAPT(*adapter));

  // Check that the values remain the same (more comprehensive test of get1dView elsewhere...
  // TEST_EQUALITY( mv->get1dViewNonConst(),      adapter->get1dViewNonConst() );
  // TEST_EQUALITY( adapter->get1dViewNonConst(), adapter2->get1dViewNonConst() );

  // The following should all pass at compile time
  TEST_ASSERT( (is_same<double, ADAPT::scalar_type>::value) );
  TEST_ASSERT( (is_same<int,    ADAPT::local_ordinal_type>::value) );
  TEST_ASSERT( (is_same<int,    ADAPT::global_ordinal_type>::value) );
  TEST_ASSERT( (is_same<Node,   ADAPT::node_type>::value) );
  TEST_ASSERT( (is_same<size_t, ADAPT::global_size_type>::value) );
  TEST_ASSERT( (is_same<MV,     ADAPT::multivec_type>::value) );

}

TEUCHOS_UNIT_TEST( MultiVecAdapter, Dimensions )
{
  // Test that the dimensions reported by the adapter match those as reported
  // by the Tpetra::MultiVector
  typedef ScalarTraits<double> ST;
  typedef Epetra_MultiVector MV;
  typedef MultiVecAdapter<MV> ADAPT;

  Epetra_SerialComm comm;
  // create a Map
  const size_t numLocal = 10;
  Epetra_Map map(-1,numLocal,0,comm);

  RCP<MV> mv = rcp(new MV(map,11));
  mv->Random();

  RCP<ADAPT> adapter = rcp(new ADAPT(mv));

  TEST_EQUALITY( mv->MyLength(),     as<int>(adapter->getLocalLength())     );
  TEST_EQUALITY( mv->NumVectors(),   as<int>(adapter->getLocalNumVectors()) );
  TEST_EQUALITY( mv->NumVectors(),   as<int>(adapter->getGlobalNumVectors()));
  TEST_EQUALITY( mv->GlobalLength(), as<int>(adapter->getGlobalLength())    );
  TEST_EQUALITY( mv->Stride(),       as<int>(adapter->getStride())          );

}

TEUCHOS_UNIT_TEST( MultiVecAdapter, Copy )
{
  /* Test the get1dCopy() method of MultiVecAdapter.  We can check against a
   * known multivector and also check against what is returned by the
   * Tpetra::MultiVector.
   */
}

TEUCHOS_UNIT_TEST( MultiVecAdapter, View )
{
  /* Test the get1dViewNonConst() method of MultiVecAdapter.  We can check against a
   * known multivector and also check against what is returned by the
   * Tpetra::MultiVector.
   */
  typedef ScalarTraits<double> ST;
  typedef Epetra_MultiVector MV;
  typedef MultiVecAdapter<MV> ADAPT;

  Epetra_SerialComm comm;

  // create a Map
  const size_t numVectors = 7;
  const size_t numLocal = 13;
  Epetra_Map map(-1,numLocal,0,comm);

  RCP<MV> mv = rcp(new MV(map,numVectors));
  mv->Random();

  RCP<ADAPT> adapter = rcp(new ADAPT(mv));

  // Check that the values remain the same
  double* values;
  int lda;
  mv->ExtractView(&values,&lda);

  TEST_EQUALITY( Teuchos::arcp(values,0,numVectors*numLocal,false), adapter->get1dViewNonConst() );

  // Check that direct changes to the wrapped object are reflected in the adapter
  mv->Random();
  TEST_EQUALITY( Teuchos::arcp(values,0,numVectors*numLocal,false), adapter->get1dViewNonConst() );

  // check that 1dView and 1dCopy have the same values
  ArrayRCP<double> view;
  view = adapter->get1dViewNonConst();

  // clear view, ensure that mv is zero
  std::fill(view.begin(), view.end(), as<double>(0.0));
  view = Teuchos::null;
  Array<double> norms(numVectors), zeros(numVectors,as<double>(0.0));
  mv->Norm2(norms.getRawPtr());
  TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,as<double>(0.0));
}

/* Also test the updateValues() method, once it is implemented.  It should
 * take a representation from either the get1dCopy or get2dCopy methods and
 * place the values back into the matrix (either serial or distributed)
 */

} // end anonymous namespace
