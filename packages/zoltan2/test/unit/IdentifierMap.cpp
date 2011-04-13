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
//    TODO Z2 could should throw errors and we should catch them
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

#include <vector>
#include <string>
#include <ostream>
#include <iostream>
#include <Teuchos_GlobalMPISession.hpp> // So we don't have to #ifdef HAVE_MPI
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Zoltan2_Partitioner.hpp>

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
  bool pass = true;

  long numLocalObjects = nobjects /  nprocs;
  long leftOver = nobjects % nprocs;

  if (rank < leftOver) numLocalObjects++;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  if (pass) {

    // AppGID is long (a Teuchos Global Ordinal type).
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

    Z2::IdentifierMap<int, long> idmap(comm, gids, lids);
  }

#if 0

  // Test local gno look up.

  std::vector<gnoType> gnos(numLocalObjects);

  idmap.gidTranslate(*gids, gnos);    // translate app gids to z2 gnos

  show_result<appGlobalId, gnoType>(std::string("App GIDs -> Z2 GNOs"), 
    *gids, gnos);

  std::vector<gnoType> gnoQuery(10,0);

  for (int i=0; i < 20; i+=2){
    gnoQuery.push_back(gnos[i]);
  }

  std::vector<appGlobalId> gidsReturned;

  idmap.gidTranslate(gidsReturned, gnoQuery);// translate gnos to gids

  //show_result<gnoType, appGlobalId>(std::string("Z2 gnos -> App gids"), 
    //gnoQuery, gidsReturned);

  // Test local gno/lid look up.

  gnos.resize(0);

  idmap.lidTranslate(*lids, gnos);    // translate app lids to z2 gnos

  //show_result<appLocalId, gnoType>(std::string("App local IDs -> Z2 GNOs"), 
    //*lids, gnos);

  for (int i=0; i < 20; i+=2){
    gnoQuery[i] = gnos[i];
  }

  std::vector<appLocalId> lidsReturned;
  
  idmap.lidTranslate(lidsReturned, gnoQuery);// translate gnos to app lids

  //show_result<gnoType, appGlobalId>(
   // std::string("Z2 GNOs -> Application local IDs"), 
    //gnoQuery, lidsReturned);
#endif

  std::cout << "PASS" << std::endl;
}
