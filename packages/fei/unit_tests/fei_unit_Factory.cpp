
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_mpi.h>
#include <fei_CommUtils.hpp>
#include <fei_Factory_Trilinos.hpp>
#include <fei_Factory_Aztec.hpp>

#include <fei_Aztec_LinSysCore.hpp>
#include <snl_fei_Factory.hpp>

void init_VectorSpace_1(fei::VectorSpace& vspace, Teuchos::FancyOStream& out, bool& success)
{
  int fieldID = 1;
  int fieldSize = 2;
  int IDType = 1;

  vspace.defineFields(1, &fieldID, &fieldSize);
  vspace.defineIDTypes(1, &IDType);
}

void init_MatrixGraph_1(fei::MatrixGraph& mgraph, Teuchos::FancyOStream& out, bool& success)
{
  fei::SharedPtr<fei::VectorSpace> vspace = mgraph.getRowSpace();

  std::vector<int> IDTypes;
  vspace->getIDTypes(IDTypes);
  TEUCHOS_TEST_EQUALITY(IDTypes.size(), (size_t)1, out, success);

  std::vector<int> fieldIDs;
  vspace->getFields(fieldIDs);
  TEUCHOS_TEST_EQUALITY(fieldIDs.size(), (size_t)1, out, success);

  int idType = IDTypes[0];
  int fieldID = fieldIDs[0];
  int nodes_per_elem = 4;
  int num_elems = 1;

  int patternID = mgraph.definePattern(nodes_per_elem, idType, fieldID);
  int blockID = mgraph.initConnectivityBlock(num_elems, patternID);

  int elemID = 1;
  std::vector<int> elem_nodes(nodes_per_elem); 
  int thisProc = fei::localProc(MPI_COMM_WORLD);
  int firstNode = 2*thisProc;
  for(int i=0; i<nodes_per_elem; ++i) elem_nodes[i] = firstNode + i;

  int errorcode = mgraph.initConnectivity(blockID, elemID, &elem_nodes[0]);
  TEUCHOS_TEST_EQUALITY(errorcode, 0, out, success);

  errorcode = mgraph.initComplete();
  TEUCHOS_TEST_EQUALITY(errorcode, 0, out, success);
}

void test_factory_1(fei::Factory& factory, Teuchos::FancyOStream& out, bool& success)
{
  fei::SharedPtr<fei::VectorSpace> vspace = factory.createVectorSpace(MPI_COMM_WORLD, "test1_VecSpc");
  TEUCHOS_TEST_INEQUALITY(vspace.get(), (fei::VectorSpace*)NULL, out, success);

  fei::SharedPtr<fei::VectorSpace> nullvspace;
  fei::SharedPtr<fei::MatrixGraph> mgraph = factory.createMatrixGraph(vspace, nullvspace, "test1_MGrph");
  TEUCHOS_TEST_INEQUALITY(mgraph.get(), (fei::MatrixGraph*)NULL, out, success);
}

void test_Vector_1(fei::Vector& vec, Teuchos::FancyOStream& out, bool& success)
{
  fei::SharedPtr<fei::VectorSpace> vspc = vec.getVectorSpace();
  int thisProc = fei::localProc(vspc->getCommunicator());
  size_t num_owned = thisProc==0 ? 4 : 2;
  num_owned *= 2;//field-size is 2
  std::vector<int> indices;
  vspc->getIndices_Owned(indices);
  TEUCHOS_TEST_EQUALITY(indices.size(), num_owned, out, success);

  vec.putScalar(3.14);
  std::vector<double> check_coefs(num_owned, 3.14);
  std::vector<double> coefs(num_owned, 0);
  vec.copyOut(num_owned, &indices[0], &coefs[0]);
  TEST_COMPARE_ARRAYS(coefs, check_coefs);

  double coef = 3.14;
  vec.sumIn(1, &indices[0], &coef);
  check_coefs[0] += coef;
  vec.copyOut(1, &indices[0], &coef);

  TEUCHOS_TEST_EQUALITY(coef, check_coefs[0], out, success);
}

void test_Matrix_1(fei::Matrix& mat, Teuchos::FancyOStream& out, bool& success)
{
  mat.putScalar(3.14);

  fei::SharedPtr<fei::VectorSpace> vspc = mat.getMatrixGraph()->getRowSpace();
  int thisProc = fei::localProc(vspc->getCommunicator());
  size_t check_num_rows = thisProc==0 ? 4 : 2;
  check_num_rows *= 2;//field-size is 2
  size_t num_rows = mat.getLocalNumRows();
  TEUCHOS_TEST_EQUALITY(num_rows, check_num_rows, out, success);
}

void test_factory_2(fei::Factory& factory, Teuchos::FancyOStream& out, bool& success)
{
  fei::SharedPtr<fei::VectorSpace> vspace = factory.createVectorSpace(MPI_COMM_WORLD, "test1_VecSpc2");
  fei::SharedPtr<fei::VectorSpace> nullvspace;
  fei::SharedPtr<fei::MatrixGraph> mgraph = factory.createMatrixGraph(vspace, nullvspace, "test1_MGrph2");

  init_VectorSpace_1(*vspace, out, success);
  init_MatrixGraph_1(*mgraph, out, success);

  fei::SharedPtr<fei::Matrix> A = factory.createMatrix(mgraph);
  bool A_not_null = A.get() != NULL;
  TEUCHOS_TEST_EQUALITY(A_not_null, true, out, success);

  bool is_soln_vec = false;
  fei::SharedPtr<fei::Vector> b = factory.createVector(mgraph, is_soln_vec);
  bool b_not_null = b.get() != NULL;
  TEUCHOS_TEST_EQUALITY(b_not_null, true, out, success);

  test_Vector_1(*b, out, success);
  test_Matrix_1(*A, out, success);
}

TEUCHOS_UNIT_TEST(Factory, Trilinos1)
{
  Factory_Trilinos factory(MPI_COMM_WORLD);

  test_factory_1(factory, out, success);

  test_factory_2(factory, out, success);
}

TEUCHOS_UNIT_TEST(Factory, AZLSC1)
{
  fei::SharedPtr<LinearSystemCore> az_lsc(new fei_trilinos::Aztec_LinSysCore(MPI_COMM_WORLD));
  snl_fei::Factory factory(MPI_COMM_WORLD, az_lsc);

  test_factory_1(factory, out, success);

  test_factory_2(factory, out, success);
}

TEUCHOS_UNIT_TEST(Factory, Aztec1)
{
  Factory_Aztec factory(MPI_COMM_WORLD);

  test_factory_1(factory, out, success);

//  test_factory_2(factory, out, success);
}

