#include "Tpetra_Map.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include <unordered_set>

typedef Tpetra::Map<>::global_ordinal_type GO_TYPE;

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

int main(int narg, char *arg[]) 
{
  Tpetra::ScopeGuard tpetraScope(&narg, &arg);
  auto comm = Tpetra::getDefaultComm();
  int np = comm->getSize();
  int me = comm->getRank();

  int ntrials = 100;
  if (narg > 1) ntrials = std::atoi(arg[1]);

  size_t gsize = 3203254;
  if (narg > 2) gsize = std::atoi(arg[2]);
  
  Tpetra::global_size_t dummy;

  // Owned IDs:  every np-th entry --> non-contiguous
  Teuchos::ArrayRCP<GO_TYPE> ownedIds(gsize / np + 1);
  int nOwnedIds = 0;
  for (int i = 0; i < gsize; i++)
    if (i % np == me) ownedIds[nOwnedIds++] = i;

  // Shared IDs:  randomly selected; no local duplicates
  size_t maxShared = gsize / np * 2;
  Teuchos::ArrayRCP<GO_TYPE> sharedIds(maxShared);
  std::unordered_set<GO_TYPE> shared(maxShared);
  size_t nSharedIds = 0;
  for (int i = 0; i < maxShared; i++) {
    GO_TYPE tmp = (rand() % gsize);
    if (shared.find(tmp) == shared.end()) {
      shared.insert(tmp);
      sharedIds[nSharedIds++] = tmp;
    }
  }

  for (int n = 0; n < ntrials; n++) {

    // Build one-to-one owned map
    std::cout << "Iter " << n << " of " << ntrials << " BUILD OWNED MAP " 
              << std::endl;

    dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const Tpetra::Map<> > ownedMap = 
             Teuchos::rcp(new Tpetra::Map<>(dummy, ownedIds(0,nOwnedIds),
                                            0, comm));

    // Build overlapped shared map
    std::cout << "Iter " << n << " of " << ntrials << " BUILD SHARED MAP " 
              << std::endl;

    dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const Tpetra::Map<> > sharedMap = 
             Teuchos::rcp(new Tpetra::Map<>(dummy, sharedIds(0, nSharedIds),
                                            0, comm));

    // Create multivectors with owned and shared maps; then import
    std::cout << "Iter " << n << " of " << ntrials 
              << ": IMPORT OWNED + SHARED " 
              << ownedMap->getNodeNumElements() << " "
              << sharedMap->getNodeNumElements() << std::endl;

    testImports(ownedMap, sharedMap, me);

    // Create a one-to-one version of the shared map
    std::cout << "Iter " << n << " of " << ntrials << " BUILD ONETOONE MAP " 
              << std::endl;
    Teuchos::RCP<const Tpetra::Map<> > oneToOneMap =
             Tpetra::createOneToOne(sharedMap);

    // Create multivectors with created oneToOne and shared maps; then import
    std::cout << "Iter " << n << " of " << ntrials 
              << ": IMPORT ONETOONE + SHARED " 
              << oneToOneMap->getNodeNumElements() << " "
              << sharedMap->getNodeNumElements() << std::endl;

    testImports(oneToOneMap, sharedMap, me);
  }

  return 0;
}
