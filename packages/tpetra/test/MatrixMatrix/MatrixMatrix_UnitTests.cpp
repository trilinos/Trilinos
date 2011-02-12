#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <cmath>

namespace Tpetra{

static const double defaultEpsilon = 1e-10;
bool verbose = false;
std::string matnamesFile;
bool write_result_hb = false;

TEUCHOS_STATIC_SETUP(){
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("matnames-file", &matnamesFile, 
    "A file containing a list of matricies we'll import", true);
  clp.setOption("writeresults", "no-write-results", &write_result_hb, 
    "Whether or not to write the resutling matricies to hb files");
  clp.setOption("v", "not-verbose", &verbose, 
    "Whether or not to use verbose output");
}

template<class Ordinal>
int add_test(
    RCP<CrsMatrix<double,int> > A,
    RCP<CrsMatrix<double,int> > B,
    RCP<CrsMatrix<double,int> > C,
    bool AT,
    bool BT,
    double epsilon,
    RCP<const Comm<Ordinal> > comm,
    bool verbose)
{
  typedef Kokkos::DefaultNode::DefaultNodeType DNode;


  RCP<CrsMatrix<double,int> > computedC = null;
  RCP<const Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  computedC = rcp( new CrsMatrix<double,int>(rowmap, 1));

  MatrixMatrix::Add(*A, false, 1.0, *B, false, 1.0, computedC);

  computedC->fillComplete(C->getDomainMap(), C->getRangeMap());

  MatrixMatrix::Add(*C, false, -1.0, *computedC, 1.0);

  double calculated_euc_norm = computedC->getEuclideanNorm();
  double c_euc_norm = C->getEuclideanNorm();
  double resultVal1 = calculated_euc_norm/c_euc_norm;

  MatrixMatrix::Add(*A, false, 1.0, *B, 1.0);

  MatrixMatrix::Add(*C, false, -1.0, *B, 1.0);

  calculated_euc_norm = B->getEuclideanNorm();
  double resultVal2 = calculated_euc_norm/c_euc_norm;

  return (resultVal1 < epsilon && resultVal2 < epsilon) ? 0 : 1;

}

template<class Ordinal>
int multiply_test(
  RCP<CrsMatrix<double, int> > A,
  RCP<CrsMatrix<double, int> > B,
  bool AT,
  bool BT,
  double epsilon,
  RCP<const Comm<Ordinal> > comm,
  bool verbose,
  std::ostream& out)
{

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;
  typedef CrsMatrix<double,int> CrsMatrix_t;
  typedef Map<int, int> Map_t;
  typedef Vector<double> Vector_t;
  typedef CrsMatrixMultiplyOp<double> CrsMatrixMultiplyOp_t;
  RCP<const Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  RCP<CrsMatrix_t> computedC = rcp( new CrsMatrix_t(rowmap, 1));

  MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC);

  RCP<Vector_t> randomVector1 = rcp(new Vector_t(computedC->getGraph()->getImporter()->getSourceMap(), true));
  randomVector1->randomize();
  RCP<Vector_t> y1 = rcp(new Vector_t(computedC->getRowMap(), true));
  CrsMatrixMultiplyOp_t cTimesVOp(computedC);
  cTimesVOp.apply(*randomVector1, *y1);


  RCP<Vector_t> intermediatVector = rcp(new Vector_t(A->getGraph()->getImporter()->getSourceMap(), true));
  CrsMatrixMultiplyOp_t bTimesVOp(B);
  bTimesVOp.apply(*randomVector1,*intermediatVector);


  RCP<Vector_t> y2 = rcp(new Vector_t(A->getRowMap(), true));
  CrsMatrixMultiplyOp_t aTimesVOp(A);
  aTimesVOp.apply(*intermediatVector, *y2);



  double diff_result = y1->norm2()-y2->norm2();


  out << "Difference in norms: " << abs(diff_result);
  return abs(diff_result)<epsilon ? 0 : -1;

}


template<class Ordinal>
int run_test(RCP<const Comm<Ordinal> > comm,
             Teuchos::ParameterList matrixSystem,
             bool result_mtx_to_file,
             bool verbose,
             std::ostream& out)
{
  std::string A_file = matrixSystem.get<std::string>("A");
  std::string B_file = matrixSystem.get<std::string>("B");
  bool AT = matrixSystem.get<bool>("TransA");
  bool BT = matrixSystem.get<bool>("TransB");
  double epsilon = matrixSystem.get<double>("epsilon", defaultEpsilon);
  std::string op = matrixSystem.get<std::string>("op");


  RCP<CrsMatrix<double,int> > A = null;
  RCP<CrsMatrix<double,int> > B = null;

  Utils::readHBMatrix(A_file, comm, Kokkos::DefaultNode::getDefaultNode(), A);
  Utils::readHBMatrix(B_file, comm, Kokkos::DefaultNode::getDefaultNode(), B);

  if(op == "multiply"){
    out << "Running multiply test for " << matrixSystem.name() << std::endl;
    return multiply_test(A,B,AT,BT,epsilon,comm,verbose, out);
  }
  else if(op == "add"){
    std::string C_file = matrixSystem.get<std::string>("C");
    RCP<CrsMatrix<double,int> > C_check = null;
    Utils::readHBMatrix(C_file, comm, Kokkos::DefaultNode::getDefaultNode(), C_check);
    out << "Running add test for " << matrixSystem.name() << std::endl;
    return add_test(A,B,C_check,AT,BT,epsilon,comm,verbose);
  }
  else{
    out << "Unrecognize matrix operation: " << op << ".";
    return -1;
  }


}



TEUCHOS_UNIT_TEST(Tpetra_MatMat, test_find_rows){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  int numprocs = comm->getSize();
  int localproc = comm->getRank();
  int numlocalrows = 2;
  global_size_t numglobalrows = numprocs*numlocalrows;
  RCP<Map<int> > rowmap = 
    rcp(new Map<int>(numglobalrows, 0, comm));
  CrsMatrix<double, int> matrix(rowmap, numglobalrows);

  Array<int> cols(numglobalrows);
  Array<double> vals(numglobalrows);

  for(size_t j=0; j<numglobalrows; ++j) {
    cols[j] = j;
    vals[j] = 1.0;
  }

  RCP<Map<int> > colmap = 
    rcp(new Map<int>(-1, cols(), 0, comm));

  for(int i=0; i<numlocalrows; ++i) {
    Array<int> row(1,localproc*numlocalrows+i);
    matrix.insertGlobalValues(
      row[0], row.view(0,1),  vals.view(i,1) );
  }

  matrix.fillComplete();

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;
  typedef Kokkos::DefaultKernels<double, int, DNode>::SparseOps SpMatOps;

  RCP<const Map<int> > map_rows = 
    MMdetails::find_rows_containing_cols<double, int, int, DNode, SpMatOps>(matrix, colmap);

  TEST_EQUALITY(map_rows->getNodeNumElements(), numglobalrows);

}

TEUCHOS_UNIT_TEST(Tpetra_MatMat, operations_test){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  Teuchos::RCP<Teuchos::ParameterList> matrixSystems = 
    Teuchos::getParametersFromXmlFile(matnamesFile);
  for(
    Teuchos::ParameterList::ConstIterator it = matrixSystems->begin();
    it != matrixSystems->end();
    ++it)
  {
	  TEST_FOR_EXCEPTION(!it->second.isList(), std::runtime_error,
      "All top level items in the matrix "
	    "file names list must be ParameterLists! In otherwords, you always "
      "need to have matricies "
	    "encapsulated within a matrixsystem" << std::endl <<
      "Bad tag's name: " << it->first << 
      "Type name: " << it->second.getAny().typeName() << 
      std::endl << std::endl);
      
    TEST_EQUALITY(
    Tpetra::run_test<int>(
      comm, 
      matrixSystems->sublist(it->first), 
      write_result_hb, 
      verbose,
      out), 0);
  }
}


} //namespace Tpetra

