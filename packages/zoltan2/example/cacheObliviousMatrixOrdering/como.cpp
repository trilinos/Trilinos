// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
//  TODO explain what's happening here
//   TODO - an example should not use RCPs.

#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_GraphModel.hpp>
//#include <Zoltan2_HypergraphModel.hpp>
//#include <Zoltan2_BasicMatrixInput.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_OrderingSolution.hpp>

#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <MatrixMarket_Tpetra.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::rcp_const_cast;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::ArrayView;

using std::string;
using std::vector;

#define STR_VALUE(path) #path
#define PATH_NAME(path) STR_VALUE(path)

// Data types are defined in Zoltan2_TestHelpers.hpp

typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  //////////////////////////////////////////////////////////
  // Use our test utilities to read in a matrix market file.
  //////////////////////////////////////////////////////////

#ifdef Z2_DATA_DIR
  string filePath(PATH_NAME(Z2_DATA_DIR));
#else
  string filePath(".");
#endif

  string mtxFileName(filePath + "/simple.mtx");

  RCP<UserInputForTests> uinput;

  try{
    uinput = rcp(new UserInputForTests(mtxFileName, comm));
  }
  catch(std::exception &e){
    if (rank==0)
      std::cerr << "Bad file name " << mtxFileName << std::endl;
    return 1;
  }

  RCP<tcrsMatrix_t> m = uinput->getTpetraCrsMatrix();

  RCP<const tcrsMatrix_t> matrix = rcp_const_cast<const tcrsMatrix_t>(m);

  //////////////////////////////////////////////////////////
  // Create a Matrix input adapter for Zoltan2
  //////////////////////////////////////////////////////////

  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> ia_t;

  RCP<ia_t> ia;

  try{
   ia = rcp(new ia_t(matrix));
  }
  catch(std::exception &e){
    if (rank==0)
      std::cerr << "Can not create input adapter" << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  // Create a partitioning problem.
  //////////////////////////////////////////////////////////

  typedef Zoltan2::PartitioningProblem<ia_t> parProb_t;

  RCP<parProb_t> parProb;

  ParameterList params;
  ParameterList &partitioningParams = params.sublist("partitioning");
  //partitioningParams.set("num_local_parts", size_t(1));
  partitioningParams.set("algorithm", "scotch");
  partitioningParams.set("approach", "partition");

  // If we have enough processes, let's do hierarchcial partitioning
  if (nprocs >= 8)
    partitioningParams.set("topology", "2,2,2");
  else if (nprocs >= 4)
    partitioningParams.set("topology", "2,2");

  try{
    parProb = rcp(new parProb_t(ia.getRawPtr(), &params));
  }
  catch(std::exception &e){
    if (rank==0)
      std::cerr << "Can not create problem " << e.what() << std::endl;
    return 1;
  }




  

  if (rank==0)
    std::cout << "PASS" << std::endl;
  return 0;
}

