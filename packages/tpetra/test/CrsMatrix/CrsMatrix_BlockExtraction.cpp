#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <TpetraExt_BlockExtraction.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <numeric>

namespace {

  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::ScalarTraits;
  using Teuchos::OrdinalTraits;
  using Teuchos::ArrayRCP;
  using Tpetra::CrsMatrix;
  using Tpetra::RowMatrix;
  using Tpetra::Map;
  using Tpetra::global_size_t;

  bool testMpi = true;
  string filedir;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected matrix files.");
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
    RCP<const Comm<int> > ret;
    if (testMpi) {
      ret = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( BlockExtraction, BlockDiagonalExtraction, LO, GO, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // set the block sizes
    Teuchos::ArrayView<LO> block_sizes( Teuchos::tuple<LO>(1,3,5,7,5,3,1) );
    // create a Map
    const size_t numLocal = std::accumulate( block_sizes.begin(), block_sizes.end(), (size_t)0 );
    RCP<const Map<LO,GO> > map = Tpetra::createContigMap<LO,GO>(INVALID,numLocal,comm);
    RCP<RowMatrix<Scalar,LO,GO> > mat;
    {
      RCP<MAT> mat_crs = rcp( new MAT(map,0,Tpetra::DynamicProfile) );
      mat_crs->fillComplete();
      mat = mat_crs;
    }
    Teuchos::ArrayRCP<Scalar> block_diagonals;
    Teuchos::ArrayRCP<LO>     block_offsets;
    Tpetra::Ext::extractBlockDiagonals<Scalar,LO,GO>( *mat, block_sizes, /*out*/ block_diagonals, /*out*/ block_offsets );
  }




  // 
  // INSTANTIATIONS
  //

#define UNIT_TEST_ORDINAL_SCALAR(LO, GO, SCALAR) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( BlockExtraction, BlockDiagonalExtraction, LO, GO, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
        UNIT_TEST_ORDINAL_SCALAR(LO, GO, double)

UNIT_TEST_GROUP_ORDINAL(int, int)

}

