// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// TODO: doxygen comments
// TODO Z2 could should throw errors and we should catch them
// TODO rewrite using Teuchos Unittest
//
// 3 cases:
//   Application GID is a Teuchos Global Ordinal type
//      GIDs are consecutive and increase with rank
//      GIDs are mixed up
//
//   Application GIDs can not be used as Teuchos Global Ordinals
//
// 2 cases:
//   Application supplies local IDs
//   Application does not supply local IDs
//
// Returns 0 on success, 1 on failure.

#include <string>
#include <ostream>
#include <iostream>
#include <exception>
#include <Teuchos_GlobalMPISession.hpp> // So we don't have to #ifdef HAVE_MPI
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Zoltan2_Partitioner.hpp>

using namespace std;

#define nobjects 10000

template< typename T1, typename T2>
 void show_result(std::string msg, std::vector<T1> &v1, std::vector<T2> &v2)
{
 std::cout << msg << std::endl;
 for (size_t i=0; i < v1.size(); i++){
    std::cout << v1[i] << "    " << v2[i] << std::endl;
 }
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int nprocs = session.getNProc();
  int rank = session.getRank();
  int errcode = 0;

  long numLocalObjects = nobjects /  nprocs;
  long leftOver = nobjects % nprocs;

  if (rank < leftOver) numLocalObjects++;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  if (!errcode) {

    // AppGID is long (a Teuchos Global Ordinal type).
    // AppGIDs are not consecutive
    // AppLID is an int.

    Teuchos::ArrayRCP<long> gids(new long [numLocalObjects], 
      0, numLocalObjects, true);
    Teuchos::ArrayRCP<int> lids(new int [numLocalObjects], 
      0, numLocalObjects, true);

    long base = nobjects * rank;   // nonconsecutive gids

    for (int i=0; i < numLocalObjects; i++){
      gids[i] = base + i;
      lids[i] = i;
    }

    // Template parameters: AppLID, AppGID, LNO, GNO

    std::string problem("IdentifierMap<int, long, int, long>");
    Z2::IdentifierMap<int, long, int, long> idmap;

    try {
      idmap.initialize(comm, gids, lids);
    }
    catch (std::exception &e){
      std::cerr << rank << ") initialize error: " << e.what();
      return 1;
    }

    if (idmap.gnosAreGids() != true){
      std::cerr << problem << " gnosAreGids" << std::endl;
      return 1;
    }

    Teuchos::Array<long> gnoArray1(numLocalObjects);
    Teuchos::Array<long> gnoArray2(numLocalObjects);

    Teuchos::ArrayView<long> gidArray = gids.view(0, numLocalObjects);

    try{
      idmap.gidTranslate(gidArray, gnoArray1, Z2::TRANSLATE_GID_TO_GNO);
    }
    catch (std::exception &e){
      std::cerr << rank << ") gidTranslate error: " << e.what();
      return 1;
    }

    Teuchos::ArrayView<int> lidArray = lids.view(0, numLocalObjects);

    try{
      idmap.lidTranslate(lidArray, gnoArray2, Z2::TRANSLATE_LID_TO_GNO);
    }
    catch (std::exception &e){
      std::cerr << rank << ") lidTranslate error: " << e.what();
      return 1;
    }

    for (int i=0; i < numLocalObjects; i++){
      if (gnoArray1[i] != gnoArray2[i]){
        errcode = 1;
        break;
      }
    }
  }

  if (errcode){
    std::cerr << rank << ") gno array is wrong" << std::endl;
  }

  return errcode;
}
