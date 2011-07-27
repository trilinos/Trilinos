#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOS_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#endif

namespace {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;

  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::CrsMatrix;
  using Tpetra::createCrsMatrix;
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::MatrixMarket::Reader;
  using Tpetra::MatrixMatrix::Add;
  using Tpetra::RowMatrixTransposer;


  using Kokkos::SerialNode;
  RCP<SerialNode> snode;
#ifdef HAVE_KOKKOS_TBB
  using Kokkos::TBBNode;
  RCP<TBBNode> tbbnode;
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
  using Kokkos::TPINode;
  RCP<TPINode> tpinode;
#endif
#ifdef HAVE_KOKKOS_THRUST
  using Kokkos::ThrustGPUNode;
  RCP<ThrustGPUNode> thrustnode;
#endif
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

  RCP<const Comm<int> > getDefaultComm()
  {
    RCP<const Comm<int> > ret;
    if (testMpi) {
      ret = DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  template <class Node>
  RCP<Node> getNode() {
    assert(false);
  }

  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    if (snode == null) {
      Teuchos::ParameterList pl;
      snode = rcp(new SerialNode(pl));
    }
    return snode;
  }

#ifdef HAVE_KOKKOS_TBB
  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    if (tbbnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tbbnode = rcp(new TBBNode(pl));
    }
    return tbbnode;
  }
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
  template <>
  RCP<TPINode> getNode<TPINode>() {
    if (tpinode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tpinode = rcp(new TPINode(pl));
    }
    return tpinode;
  }
#endif

#ifdef HAVE_KOKKOS_THRUST
  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    if (thrustnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      pl.set<int>("Verbose",1);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }
#endif

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( RowMatrixTransposer, RectangularTranspose, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef OrdinalTraits<LO> LOT;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    int numProcs = comm->getSize();
    Tpetra::global_size_t numGlobal = 4*numProcs;

    RCP<const Map<LO,GO,Node> > rowMap = createUniformContigMapWithNode<LO,GO>(numGlobal,comm,node);

    RCP<MAT> matrix = Reader<MAT>::readSparseFile("a.mtx", comm, node);
    RCP<MAT> matrixT = Reader<MAT>::readSparseFile("atrans.mtx", comm, node);

    RowMatrixTransposer<Scalar, LO, GO, Node> at(*matrix);
    RCP<MAT> calculated = at.createTranspose();

    RCP<MAT> diffMatrix = rcp(new MAT(matrixT->getRowMap(), matrixT->getNodeMaxNumRowEntries()));

    Scalar sOne = ScalarTraits<Scalar>::one();
    Add(*calculated, false, -sOne, *matrixT, false, sOne, diffMatrix);
    diffMatrix->fillComplete(matrixT->getDomainMap(), matrixT->getRangeMap());

    Scalar diffNorm = diffMatrix->getFrobeniusNorm(); 
    Scalar realNorm = matrixT->getFrobeniusNorm(); 
    Scalar epsilon = diffNorm/realNorm;
    
    TEST_COMPARE(ScalarTraits<Scalar>::real(epsilon), <, 1e-10)
    TEST_COMPARE(ScalarTraits<Scalar>::imag(epsilon), <, 1e-10)

  }


//Not doing complex since we don't have a complex matrix to test on
//typedef std::complex<float>  ComplexFloat;
//typedef std::complex<double> ComplexDouble;

#define UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( RowMatrixTransposer, RectangularTranspose, LO, GO, SCALAR, NODE )  \





#define UNIT_TEST_SERIALNODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, SerialNode )

typedef Kokkos::DefaultNode::DefaultNodeType DefaultNode;
#define UNIT_TEST_DEFAULTNODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, DefaultNode )

#ifdef HAVE_KOKKOS_TBB
#define UNIT_TEST_TBBNODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, TBBNode )
#else
#define UNIT_TEST_TBBNODE(LO, GO, SCALAR)
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
#define UNIT_TEST_TPINODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, TPINode )
#else
#define UNIT_TEST_TPINODE(LO, GO, SCALAR)
#endif

// don't test Kokkos node for MPI builds, because we probably don't have multiple GPUs per node
#if defined(HAVE_KOKKOS_THRUST) && !defined(HAVE_TPETRA_MPI)
// float
#if defined(HAVE_KOKKOS_CUDA_FLOAT)
#  define UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, float, ThrustGPUNode )
#else
#  define UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO)
#endif
// double
#if defined(HAVE_KOKKOS_CUDA_DOUBLE)
#  define UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, double, ThrustGPUNode )
#else
#  define UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO)
#endif
#else
// none
# define UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO)
# define UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO)
# define UNIT_TEST_THRUSTGPUNODE_COMPLEX_FLOAT(LO, GO)
# define UNIT_TEST_THRUSTGPUNODE_COMPLEX_DOUBLE(LO, GO)
#endif

#define UNIT_TEST_ALLCPUNODES(LO, GO, SCALAR) \
    UNIT_TEST_DEFAULTNODE(LO, GO, SCALAR) 

#define UNIT_TEST_FLOAT(LO, GO) \
    UNIT_TEST_ALLCPUNODES(LO, GO, float) \
    UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO)

#define UNIT_TEST_DOUBLE(LO, GO) \
    UNIT_TEST_ALLCPUNODES(LO, GO, double) \
    UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO)


#if defined(HAVE_TPETRA_INST_DOUBLE)
  UNIT_TEST_DOUBLE(int, int)
#endif

}
