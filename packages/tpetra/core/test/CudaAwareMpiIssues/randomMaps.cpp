#include "Tpetra_Map.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include <unordered_set>

typedef Tpetra::Map<>::global_ordinal_type GO_TYPE;

static int MagicIteration = -1;

//////////////////////////////////////////////////////////////////////////////
void printIds(const char *name, const Teuchos::ArrayView<GO_TYPE> &ids, int me)
{
  std::cout << "KDD PRINTING " << name << std::endl;
  char filename[180];
  sprintf(filename, "%sIDs.%02d", name, me);
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "%d\n", ids.size());
  for (size_t i = 0; i < ids.size(); i++)
   fprintf(fp, "%d\n", ids[i]);
  fclose(fp);
}
//////////////////////////////////////////////////////////////////////////////
void printMap(const char *name, Teuchos::RCP<const Tpetra::Map<> > &map)
{
  std::cout << "KDD PRINTING " << name << std::endl;
  char filename[180];
  sprintf(filename, "%s.%02d", name, map->getComm()->getRank());
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "%d\n", map->getNodeNumElements());
  for (size_t i = 0; i < map->getNodeNumElements(); i++)
   fprintf(fp, "%d\n", map->getNodeElementList()[i]);
  fclose(fp);
}

//////////////////////////////////////////////////////////////////////////////
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

  // Build "shared" map with up to lsize*2 elements
  int lsize = ownedMap->getNodeNumElements();

  std::unordered_set<GO_TYPE> shared(lsize*2);
  int cnt = 0;

  Teuchos::ArrayRCP<GO_TYPE> idx;

  if (whichTest == 0 || whichTest == 1) {
    idx = Teuchos::ArrayRCP<GO_TYPE>(lsize*2);

    if (whichTest == 0) {
      // Copy owned entries
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
    std::cout << "KDD SHARED MAP IS RANDOM " 
              << (whichTest == 0 ? " WITH " : " WITHOUT ") << " OVERLAP"
              << std::endl;
  }
  else if (whichTest == 2) {
    int halfnp = np / 2;
    int halfme = me / 2;
    idx = Teuchos::ArrayRCP<GO_TYPE>(2*gsize/np+1);
    for (int i = 0; i < gsize; i++) {
      if (i % halfnp == halfme) idx[cnt++] = i;
    }
    std::cout << "KDD SHARED MAP IS STRIDED AND DUPLICATED " << std::endl;
  }

  if (n == MagicIteration && whichTest == 1) printIds("sharedRandomWithout", idx(0,cnt), me);

  Tpetra::global_size_t dummy = 
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const Tpetra::Map<> > sharedMap = 
           Teuchos::rcp(new Tpetra::Map<>(dummy, idx(0,cnt), 0, comm));

  if (n == MagicIteration && whichTest == 1) printMap("sharedRandomWithout", sharedMap);

  // Make sure the maps work
  std::cout << "KDD IMPORTS ORIG " << std::endl;
  testImports(ownedMap, sharedMap, me);

  // Now create a one-to-one map from the shared map and do the same
  Teuchos::RCP<const Tpetra::Map<> > oneToOneMap = 
           Tpetra::createOneToOne(sharedMap);

  if (n == MagicIteration && whichTest == 1) printMap("oneToOne", oneToOneMap);

  std::cout << "KDD IMPORTS ONE2ONE " << std::endl;
  testImports(oneToOneMap, sharedMap, me);

  std::cout << "Iter " << n << " of " << ntrials 
            << " gsize " << gsize 
            << " owned " << ownedMap->getNodeNumElements()
            << " shared " << sharedMap->getNodeNumElements()
            << " one2one " << oneToOneMap->getNodeNumElements()
            << " whichTest " << whichTest
            << std::endl;

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
  if (narg > 3) MagicIteration = std::atoi(arg[3]);
  
  for (int n = 0; n < ntrials; n++) {

    int lsize = (rand() % maxsize);
    int gsize;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &lsize, &gsize);

    // Build "owned" map with Trilinos default layout
    Teuchos::RCP<const Tpetra::Map<> > ownedMap;
    switch (n%3) {
    case 0:
      // lsize per processor
      ownedMap = Teuchos::rcp(new Tpetra::Map<>(gsize, lsize, 0, comm));
      std::cout << "KDD OWNED MAP IS CONTIGUOUS NON-UNIFORM" << std::endl;
      break;
    case 1:
      // gsize/np per processor
      ownedMap = Teuchos::rcp(new Tpetra::Map<>(gsize, 0, comm));
      std::cout << "KDD OWNED MAP IS CONTIGUOUS UNIFORM" << std::endl;
      break;
    case 2:
      // every np-th entry
      Teuchos::ArrayRCP<GO_TYPE> mine(gsize / np + 1);
      int cnt = 0;
      for (int i = 0; i < gsize; i++)
        if (i % np == me) mine[cnt++] = i;

      Tpetra::global_size_t dummy = 
              Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
      ownedMap = Teuchos::rcp(new Tpetra::Map<>(dummy, mine(0,cnt), 0, comm));
      std::cout << "KDD OWNED MAP IS NONCONTIGUOUS" << std::endl;
      break;
    }

    if (n == MagicIteration) printMap("owned", ownedMap);

    // test with overlapped entries (owned + some shared)
    testMaps(ownedMap, comm, gsize, n, ntrials, 0);  

    // test with random entries
    testMaps(ownedMap, comm, gsize, n, ntrials, 1); 

    // test with duplicate maps on some ranks
    testMaps(ownedMap, comm, gsize, n, ntrials, 2); 
  }

  return 0;
}
