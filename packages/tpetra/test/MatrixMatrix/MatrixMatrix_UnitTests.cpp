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
#include "MatrixMarket_Tpetra.hpp"
#include <cmath>

namespace {
  static const double defaultEpsilon = 1e-10;
  bool verbose = false;
  std::string matnamesFile;

  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Tpetra::global_size_t;
  using Teuchos::Comm;
  using Tpetra::CrsMatrix;
  using Tpetra::Map;
  using Teuchos::Array;
  using Tpetra::Vector;
  using Tpetra::CrsMatrixMultiplyOp;
  using Tpetra::DefaultPlatform;
  using Tpetra::MatrixMarket::Reader;
  using Teuchos::ArrayView;
  using Teuchos::FancyOStream;


TEUCHOS_STATIC_SETUP(){
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("matnames-file", &matnamesFile, 
    "A file containing a list of matricies we'll import", true);
  clp.setOption("v", "not-verbose", &verbose, 
    "Whether or not to use verbose output");
}



typedef struct add_test_results_struct{
  double epsilon1;
//  double epsilon2;
  double y1Norm;
  double y2Norm;
  //double y3Norm;
} add_test_results;

typedef struct mult_test_results_struct{
  double epsilon;
  double cNorm;
  double compNorm;
} mult_test_results;

template<class CrsMatrix_t> double getNorm(RCP<CrsMatrix_t> matrix){
  double mySum = 0;
  Array<int> inds(matrix->getNodeMaxNumRowEntries());
  Array<double> vals(matrix->getNodeMaxNumRowEntries());
  for(int i =0; ((size_t)i)<matrix->getNodeNumRows(); ++i){
    size_t numRowEnts = matrix->getNumEntriesInLocalRow(i);
    ArrayView<const int> indsView = inds();
    ArrayView<const double> valsView = vals();
    matrix->getLocalRowView(i, indsView, valsView);
    for(size_t j=0; ((size_t)j)<numRowEnts; ++j){
      mySum = valsView[j]*valsView[j];
    }
  }
  double totalSum = 0;
  Teuchos::reduceAll(*(matrix->getComm()), Teuchos::REDUCE_SUM, 1, &mySum, &totalSum);
  return sqrt(totalSum);

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
mult_test_results multiply_test(
  RCP<CrsMatrix<double, int> > A,
  RCP<CrsMatrix<double, int> > B,
  bool AT,
  bool BT,
  RCP<CrsMatrix<double, int> > C,
  RCP<const Comm<Ordinal> > comm,
  FancyOStream& out)
{

  typedef CrsMatrix<double,int> CrsMatrix_t;
  typedef Map<int, int> Map_t;
  RCP<const Map_t> map = A->getRowMap();

  RCP<CrsMatrix_t> computedC = rcp( new CrsMatrix_t(map, 1));

  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC);
  double cNorm = getNorm(C);
  Tpetra::MatrixMatrix::Add(*C, false, -1.0, *computedC, 1.0);
  double compNorm = getNorm(computedC);
  mult_test_results results;
  results.epsilon = compNorm/cNorm;
  results.cNorm = cNorm;
  results.compNorm = compNorm;
  return results;
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
    Tpetra::MMdetails::find_rows_containing_cols<double, int, int, DNode, SpMatOps>(matrix, colmap);

  TEST_EQUALITY(map_rows->getNodeNumElements(), numglobalrows);

}

TEUCHOS_UNIT_TEST(Tpetra_MatMat, operations_test){
  RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();
  RCP<Kokkos::DefaultNode::DefaultNodeType> node = 
    Kokkos::DefaultNode::getDefaultNode();
  Teuchos::RCP<Teuchos::ParameterList> matrixSystems = 
    Teuchos::getParametersFromXmlFile(matnamesFile);
  /*if(verbose){
    Tpetra::MMdebug::debug_stream = rcpFromRef(out);
    Tpetra::MMdebug::debug_level = Teuchos::VERB_HIGH;
  }*/
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
    std::string C_file = currentSystem.get<std::string>("C");
    bool AT = currentSystem.get<bool>("TransA");
    bool BT = currentSystem.get<bool>("TransB");
    double epsilon = currentSystem.get<double>("epsilon", defaultEpsilon);
    std::string op = currentSystem.get<std::string>("op");
  
    out << "A file: " << A_file << std::endl;

    RCP<CrsMatrix<double,int> > A = Reader<CrsMatrix<double, int> >::readFile(A_file, comm, node);
    RCP<CrsMatrix<double,int> > B = Reader<CrsMatrix<double, int> >::readFile(B_file, comm, node);
    RCP<CrsMatrix<double,int> > C = Reader<CrsMatrix<double, int> >::readFile(C_file, comm, node);
 

    TEST_FOR_EXCEPTION(op != "multiply" && op != "add", std::runtime_error,
      "Unrecognized Matrix Operation: " << op << "!" << std::endl);
  
    if(op == "multiply"){
      if(verbose){
        out << "Running multiply test for " << currentSystem.name() << std::endl;
      }
      mult_test_results results = multiply_test(A,B,AT,BT,C,comm, out);
      if(verbose){
        out << "Results:" <<std::endl;
        out << "\tEpsilon: " << results.epsilon << std::endl;
        out << "\tcNorm: " << results.cNorm << std::endl;
        out << "\tcompNorm: " << results.compNorm << std::endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon)
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

