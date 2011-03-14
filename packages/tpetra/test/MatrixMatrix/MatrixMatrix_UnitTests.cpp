#include "Tpetra_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace Tpetra{

static const double defaultEpsilon = 1e-10;
bool verbose = false;
std::string matnamesFile;
bool write_result_hb = false;

typedef struct add_test_results_struct{
  double epsilon1;
//  double epsilon2;
  double y1Norm;
  double y2Norm;
  //double y3Norm;
} add_test_results;

TEUCHOS_STATIC_SETUP(){
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("matnames-file", &matnamesFile, 
    "A file containing a list of matricies we'll import", true);
  clp.setOption("writeresults", "no-write-results", &write_result_hb, 
    "Whether or not to write the resutling matricies to hb files");
  clp.setOption("v", "not-verbose", &verbose, 
    "Whether or not to use verbose output");
}
/*
template<class Ordinal>
add_test_results add_test(
    RCP<CrsMatrix<double,int> > A,
    RCP<CrsMatrix<double,int> > B,
    bool AT,
    bool BT,
    RCP<const Comm<Ordinal> > comm)
{
  typedef Kokkos::DefaultNode::DefaultNodeType DNode;
  typedef Vector<double> Vector_t;
  typedef CrsMatrixMultiplyOp<double> CrsMatrixMultiplyOp_t;
  add_test_results toReturn;


  RCP<CrsMatrix<double,int> > computedC = null;
  RCP<const Map<int> > rowmap = AT ? A->getDomainMap() : A->getRowMap();

  computedC = rcp( new CrsMatrix<double,int>(rowmap, 1));

  MatrixMatrix::Add(*A, false, 1.0, *B, false, 1.0, computedC);

  computedC->fillComplete(A->getDomainMap(), A->getRangeMap());


  Vector_t randomVector1(computedC->getDomainMap(), true);
  randomVector1.randomize();

  Vector_t y1(computedC->getRangeMap());
  CrsMatrixMultiplyOp_t cMultiOp(computedC);
  cMultiOp.apply(randomVector1, y1);

  Vector_t ax(A->getRangeMap());
  CrsMatrixMultiplyOp_t aMultiOp(A);
  aMultiOp.apply(randomVector1, ax);
  
  Vector_t bx(B->getRangeMap());
  CrsMatrixMultiplyOp_t bMultiOp(B);
  bMultiOp.apply(randomVector1, bx);

  Vector_t* y2 = &bx;
  y2->update(1,ax,1);

  toReturn.epsilon1 = fabs(y2->norm2()-y1.norm2());

  MatrixMatrix::Add(*A, false, 1.0, *B, 1.0);
  Vector_t randomVector3(B->getDomainMap());
  randomVector3.doImport(randomVector1, importer, REPLACE);
  Vector_t y3(B->getRangeMap());
  bMultiOp.apply(randomVector3, y3);
  toReturn.epsilon2 = fabs(y2->norm2()-y3.norm2());

  

  return toReturn;

}*/

template<class Ordinal>
double multiply_test(
  RCP<CrsMatrix<double, int> > A,
  RCP<CrsMatrix<double, int> > B,
  bool AT,
  bool BT,
  RCP<const Comm<Ordinal> > comm)
{

  typedef Kokkos::DefaultNode::DefaultNodeType DNode;
  typedef CrsMatrix<double,int> CrsMatrix_t;
  typedef Map<int, int> Map_t;
  typedef Vector<double> Vector_t;
  typedef CrsMatrixMultiplyOp<double> CrsMatrixMultiplyOp_t;
  RCP<const Map_t> rowmap = AT ? A->getDomainMap() : A->getRowMap();

  RCP<CrsMatrix_t> computedC = rcp( new CrsMatrix_t(rowmap, 1));

  MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC);

  Vector_t randomVector1(computedC->getDomainMap(), true);
  randomVector1.randomize();
  Vector_t y1(computedC->getRangeMap(), true);
  CrsMatrixMultiplyOp_t cTimesVOp(computedC);
  cTimesVOp.apply(randomVector1, y1);

  Vector_t intermediateVector(A->getDomainMap(), true);
  CrsMatrixMultiplyOp_t bTimesVOp(B);
  bTimesVOp.apply(randomVector1,intermediateVector);


  Vector_t y2(A->getRangeMap(), true);
  CrsMatrixMultiplyOp_t aTimesVOp(A);
  aTimesVOp.apply(intermediateVector, y2);

  return fabs(y1.norm2()-y2.norm2());

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

  RCP<const Map<int, int, DNode> > map_rows = 
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
    Teuchos::ParameterList currentSystem = matrixSystems->sublist(it->first); 
    std::string A_file = currentSystem.get<std::string>("A");
    std::string B_file = currentSystem.get<std::string>("B");
    bool AT = currentSystem.get<bool>("TransA");
    bool BT = currentSystem.get<bool>("TransB");
    double epsilon = currentSystem.get<double>("epsilon", defaultEpsilon);
    std::string op = currentSystem.get<std::string>("op");
  

    RCP<CrsMatrix<double,int> > A = null;
    RCP<CrsMatrix<double,int> > B = null;
  
    Utils::readHBMatrix(A_file, comm, Kokkos::DefaultNode::getDefaultNode(), A);
    Utils::readHBMatrix(B_file, comm, Kokkos::DefaultNode::getDefaultNode(), B);

    TEST_FOR_EXCEPTION(op != "multiply" && op != "add", std::runtime_error,
      "Unrecognized Matrix Operation: " << op << "!" << std::endl);
  
    if(op == "multiply"){
      if(verbose){
        out << "Running multiply test for " << currentSystem.name() << std::endl;
      }
      double result = multiply_test(A,B,AT,BT,comm);
      TEST_COMPARE(result, <, epsilon)
    }
    else if(op == "add"){
/*      if(verbose){
        out << "Running add test for " << currentSystem.name() << std::endl;
      }
      add_test_results results = add_test(A,B,AT,BT,comm);
      TEST_COMPARE(results.epsilon1, <, epsilon)
 //     TEST_COMPARE(results.epsilon2, <, epsilon)
      out << "Norm 1: " << results.y1Norm << std::endl;
      out << "Norm 2: " << results.y2Norm << std::endl;
//      out << "Norm 3: " << results.y3Norm << std::endl;*/
    }
  }   
}


} //namespace Tpetra

