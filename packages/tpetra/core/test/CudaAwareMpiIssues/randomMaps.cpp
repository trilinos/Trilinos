#include "Tpetra_Map.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include <unordered_set>

typedef Tpetra::Map<>::global_ordinal_type GO_TYPE;

//////////////////////////////////////////////////////////////////////////////
//  Function to exercise the maps:  create a vector with each map, and 
//  import / export between them
void testImports(
  const Teuchos::RCP<const Tpetra::Map<> > &oneToOneMap,
  const Teuchos::RCP<const Tpetra::Map<> > &sharedMap,
  int me
)
{
  Tpetra::Import<> importer(oneToOneMap, sharedMap);

  Tpetra::Vector<double> oneToOneVec(oneToOneMap);
  oneToOneVec.putScalar(double(me));

  Tpetra::Vector<double> sharedVec(sharedMap);
  sharedVec.doImport(oneToOneVec, importer, Tpetra::REPLACE);
  oneToOneVec.doExport(sharedVec, importer, Tpetra::ADD);
}

//////////////////////////////////////////////////////////////////////////////
//  Create various shared and one-to-one maps related to a given one-to-own 
//  Map.  Then test the create maps with import / export.
void testMaps(
  const Teuchos::RCP<const Tpetra::Map<> > &ownedMap, 
  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
  int gsize,
  int n,
  int ntrials,
  int whichTest)
{
  int me = comm->getRank();
  int np = comm->getSize();

  std::cout << "Iter " << n << " of " << ntrials 
            << " gsize " << gsize 
            << " whichTest " << whichTest
            << std::endl;

  // Build "shared" map with up to lsize*2 elements
  int lsize = ownedMap->getNodeNumElements();

  std::unordered_set<GO_TYPE> shared(lsize*2);
  int cnt = 0;

  Teuchos::ArrayRCP<GO_TYPE> idx;

  if (whichTest == 0 || whichTest == 1) {
    idx = Teuchos::ArrayRCP<GO_TYPE>(lsize*2);

    if (whichTest == 0) {
      // Include owned entries in the shared map
      for (int i = 0; i < lsize; i++) {
        idx[i] = ownedMap->getNodeElementList()[i];
        shared.insert(idx[i]);
      }
      cnt = lsize;
    }

    // Generate random entries -- no duplicates
    for (int i = cnt; i < lsize*2; i++) {
      GO_TYPE tmp = (rand() % gsize);
      if (shared.find(tmp) == shared.end()) {
        shared.insert(tmp);
        idx[cnt++] = tmp;   
      }
    }

    std::cout << "SHARED MAP IS RANDOM " 
              << (whichTest == 0 ? "WITH" : "WITHOUT") << " OVERLAP "
              << sharedMap->getNodeNumElements() << std::endl;
  }
  else if (whichTest == 2) {
    // Deal out shared-map entries like playing cards (cyclic distribution)
    // among two subsets of processors (so that processors share entries)
    int halfnp = np / 2;
    int halfme = me / 2;
    idx = Teuchos::ArrayRCP<GO_TYPE>(2*gsize/np+1);
    for (int i = 0; i < gsize; i++) {
      if (i % halfnp == halfme) idx[cnt++] = i;
    }
    std::cout << "SHARED MAP IS STRIDED AND DUPLICATED "
              << sharedMap->getNodeNumElements() << std::endl;
  }

  Tpetra::global_size_t dummy = 
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const Tpetra::Map<> > sharedMap = 
           Teuchos::rcp(new Tpetra::Map<>(dummy, idx(0,cnt), 0, comm));

  // Make sure the maps work
  std::cout << "TEST IMPORTS OWNED" << std::endl;
  testImports(ownedMap, sharedMap, me);

  // Now create a one-to-one map from the shared map and do the same
  Teuchos::RCP<const Tpetra::Map<> > oneToOneMap = 
           Tpetra::createOneToOne(sharedMap);
  std::cout << "ONE2ONE MAP "
            << oneToOneMap->getNodeNumElements() << std::endl;

  std::cout << "TEST IMPORTS ONE2ONE" << std::endl;
  testImports(oneToOneMap, sharedMap, me);

}

//////////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[]) 
{

  Tpetra::ScopeGuard tpetraScope(&narg, &arg);
  auto comm = Tpetra::getDefaultComm();
  int np = comm->getSize();
  int me = comm->getRank();

  srand(comm->getRank());
  int maxsize = 10000;
  int ntrials = 100;
  if (narg > 1) ntrials = std::atoi(arg[1]);
  if (narg > 2) maxsize = std::atoi(arg[2]);
  
  for (int n = 0; n < ntrials; n++) {

    int lsize = (rand() % maxsize);
    int gsize;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &lsize, &gsize);

    // Build "owned" map with Trilinos default layout
    Teuchos::RCP<const Tpetra::Map<> > ownedMap;
    switch (n%3) {
    case 0:
      // lsize entries per processor, contiguously numbered
      ownedMap = Teuchos::rcp(new Tpetra::Map<>(gsize, lsize, 0, comm));
      std::cout << "OWNED MAP IS CONTIGUOUS NON-UNIFORM "
                <<  ownedMap->getNodeNumElements() << std::endl;
      break;
    case 1:
      // gsize/np per processor, contiguously numbered
      ownedMap = Teuchos::rcp(new Tpetra::Map<>(gsize, 0, comm));
      std::cout << "OWNED MAP IS CONTIGUOUS UNIFORM " 
                <<  ownedMap->getNodeNumElements() << std::endl;
      break;
    case 2:
      // processor has every np-th entry (cyclic distribution)
      Teuchos::ArrayRCP<GO_TYPE> mine(gsize / np + 1);
      int cnt = 0;
      for (int i = 0; i < gsize; i++)
        if (i % np == me) mine[cnt++] = i;

      Tpetra::global_size_t dummy = 
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      ownedMap = Teuchos::rcp(new Tpetra::Map<>(dummy, mine(0,cnt), 0, comm));
      std::cout << "OWNED MAP IS NONCONTIGUOUS " 
                <<  ownedMap->getNodeNumElements() << std::endl;
      break;
    }

    // test with overlapped entries (owned + some shared)
    testMaps(ownedMap, comm, gsize, n, ntrials, 0);  

    // test with random entries
    testMaps(ownedMap, comm, gsize, n, ntrials, 1); 

    // test with duplicate maps on some ranks
    testMaps(ownedMap, comm, gsize, n, ntrials, 2); 
  }

  return 0;
}
