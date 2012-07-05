// @HEADER
// ***********************************************************************
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
// ***********************************************************************
// @HEADER
//
// Test the GidLookupHelper class.

#include <Zoltan2_GidLookupHelper.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Comm;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int nprocs = comm->getSize();
  int rank = comm->getRank();
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);
  ArrayRCP<int> idx;

  if (rank == 0){
    int *buf = new int [10];
    for (int i=0; i < 10; i++)
      buf[i] = (i+1)*10;

    // ten distinct gids
    cout << endl << "10 distinct int gids" << endl;
    ArrayRCP<const int> intGids = Teuchos::arcp(buf, 0, 10, true);
    Zoltan2::GidLookupHelper<int, int> glhInt(env, intGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 10; i++){
      cout << "Gid: " << buf[i] << " idx: " << glhInt.lookup(buf[i]) << endl;
    }

    // ten gids with duplicates
    cout << endl << "10 int gids with duplicates" << endl;
    for (int i=0; i < 5; i++) 
      buf[i] = buf[i+5] = i;
    Zoltan2::GidLookupHelper<int, int> glhIntDup(env, intGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 10; i++){
      cout << "Gid: " << buf[i] << " idx: " << glhIntDup.lookup(buf[i]) << endl;
    }

    // string gids
    cout << endl << "4 string gids" << endl;
    string *sbuf = new string [4];
    sbuf[0] = string("gid number 0");
    sbuf[1] = string("gid number 1");
    sbuf[2] = string("gid number 2");
    sbuf[3] = string("gid number 3");
    ArrayRCP<const string> sGids = Teuchos::arcp(sbuf, 0, 4, true);
    Zoltan2::GidLookupHelper<string, int> glhString(env, sGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 4; i++){
      cout << "Gid: " << sbuf[i] << " idx: " << glhString.lookup(sbuf[i]) << endl;
    }
    
    // string gids with duplicates
    cout << endl << "4 string gids with a duplicate" << endl;
    sbuf[3] = string("gid number 0");
    Zoltan2::GidLookupHelper<string, int> glhStringDup(env, sGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 4; i++){
      cout << "Gid: " << sbuf[i] << " idx: " << glhStringDup.lookup(sbuf[i]) << endl;
    }

    cout << "PASS" << endl;
  }
}
