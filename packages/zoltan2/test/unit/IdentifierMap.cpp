// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER
//
// TODO: doxygen comments
// TODO Zoltan2 could should throw errors and we should catch them
// TODO rewrite using Teuchos Unittest
//     make this work if !HAVE_MPI
//
//  TODO we only test one case write more tests
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
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Zoltan2_Standards.hpp>    // for Zoltan2::default_node_t
#include <Zoltan2_InputTraits.hpp>

// user data structure 
template<typename AppLID, typename AppGID>
struct TestData{
  Teuchos::ArrayRCP<AppLID> lids;
  Teuchos::ArrayRCP<AppGID> gids;
};

// the InputTraits of our structure for Zoltan2
namespace Zoltan2{
template<>
template<typename AppLID, typename AppGID>
struct InputTraits<struct TestData<AppLID, AppGID> >
{
  typedef float scalar_t;
  typedef int lno_t;
  typedef long gno_t;
  typedef AppLID lid_t;
  typedef AppGID gid_t;
  typedef Zoltan2::default_node_t node_t;
};
}

#include <Zoltan2_IdentifierMap.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Comm;

template <typename AppLID, typename AppGID>
  int testIdentifierMap(
          RCP<const Comm<int> > &comm, RCP<Zoltan2::Environment> &envPtr,
          ArrayRCP<AppLID> &lids, ArrayRCP<AppGID> &gids,
          bool consecutiveGnosAreRequired, bool gnosShouldBeGids)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  typedef struct TestData<AppLID, AppGID> testdata_t;

  testdata_t test1;

  test1.gids = gids;
  test1.lids = lids;

  Zoltan2::IdentifierMap<AppLID, AppGID, AppLID, AppGID> idmap;

  try {
    idmap.initialize(comm, envPtr, test1.gids, test1.lids, consecutiveGnosAreRequired);
  }
  catch (std::exception &e){
    std::cerr << rank << ") initialize error: " << e.what();
    return 1;
  }

  if (idmap.gnosAreGids() != gnosShouldBeGids){
    std::cerr << " gnosAreGids" << std::endl;
    return 1;
  }

  int numLocalObjects = gids.size();

  Array<long> gnoArray1(numLocalObjects);
  Array<long> gnoArray2(numLocalObjects);

  ArrayView<AppGID> gidArray = test1.gids.view(0, numLocalObjects);

  try{
    idmap.gidTranslate(gidArray, gnoArray1, Zoltan2::TRANSLATE_APP_TO_LIB);
  }
  catch (std::exception &e){
    std::cerr << rank << ") gidTranslate error: " << e.what();
    return 1;
  }

  ArrayView<AppLID> lidArray = test1.lids.view(0, numLocalObjects);

  try{
    idmap.lidTranslate(lidArray, gnoArray2, Zoltan2::TRANSLATE_APP_TO_LIB);
  }
  catch (std::exception &e){
    std::cerr << rank << ") lidTranslate error: " << e.what();
    return 1;
  }

  for (int i=0; i < numLocalObjects; i++){
    if (gnoArray1[i] != gnoArray2[i]){
      std::cerr << rank << ") gnos don't match: " << std::endl;
      return 1;
    }
  }

  if (consecutiveGnosAreRequired){
    // TODO test htsi
  }
  return 0;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  int errcode = 0, globalcode=0;
  bool consecutiveGnosAreRequired=true; 
  bool gnosShouldBeGids=true;

  long numLocalObjects = 10000/  nprocs;
  long leftOver = 10000 % nprocs;

  if (rank < leftOver) numLocalObjects++;

  Teuchos::ParameterList params; 
  params.set(std::string("ERROR_CHECK_LEVEL"), 1);
  params.set(std::string("DEBUG_OSTREAM"), "std::cout");
  params.set(std::string("ERROR_OSTREAM"), "std::cerr");
  params.set(std::string("DEBUG_LEVEL"), 0);
    
  Teuchos::RCP<Zoltan2::Environment> envPtr = 
    Teuchos::rcp(new Zoltan2::Environment(params, comm));

  // Test GIDs are longs, but not consecutive, LIDs are ints.

  {
    // GIDs: long, LIDS: int, non consecutive global ids
    Teuchos::ArrayRCP<long> gids(new long [numLocalObjects], 0, numLocalObjects, true);
    Teuchos::ArrayRCP<int> lids(new int[numLocalObjects], 0, numLocalObjects, true);

    long base = 10000 * rank;   // nonconsecutive gids

    for (int i=0; i < numLocalObjects; i++){
      gids[i] = base + i;   
      lids[i] = i;
    }
    errcode = testIdentifierMap<int, long>(comm, envPtr, lids, gids,
                      !consecutiveGnosAreRequired, gnosShouldBeGids);

    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, 1, &errcode, &globalcode);

    if (rank == 0){
      if (globalcode){
        std::cout << "FAIL" << std::endl;
      }
      else{
        std::cout << "PASS" << std::endl;
      }
    }
    if (globalcode){
      return errcode;
    }
  }
}
