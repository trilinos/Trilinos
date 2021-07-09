///  Small test program showing how to take GIDs that may have
///  duplicates across processors (e.g., mesh vertices that are copied
///  at part boundaries in an element-based decomposition) and assign
///  unique owners to them.
///  Then, the test creates Vectors using the maps and transfers data
///  between them

#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"

#include <string>
#include <sstream>
#include <iostream>

/////////////////////////////////////////////////////////////////////

int main(int narg, char **arg)
{
  typedef Tpetra::Map<> map_t;
  typedef map_t::local_ordinal_type lno_t;
  typedef map_t::global_ordinal_type gno_t;
  typedef int scalar_t;

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  // Create a map with duplicated entries (mapWithCopies)
  // Each rank has 15 IDs, the last five of which overlap with the next rank.


  lno_t numLocalCoords = 15;
  lno_t offset = me * 10;

  Teuchos::Array<gno_t> gids(numLocalCoords);
  for (lno_t i = 0 ; i < numLocalCoords; i++)
    gids[i] = static_cast<gno_t> (offset + i);

  Tpetra::global_size_t numGlobalCoords =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> mapWithCopies =
          rcp(new map_t(numGlobalCoords, gids(), 0, comm));

  // Create a new map with IDs uniquely assigned to ranks (oneToOneMap)
  Teuchos::RCP<const map_t> oneToOneMap =
          Tpetra::createOneToOne<lno_t, gno_t>(mapWithCopies);

  // Create vectors with each map
  typedef Tpetra::Vector<scalar_t, lno_t, gno_t> vector_t;

  vector_t vecWithCopies(mapWithCopies);
  vector_t oneToOneVec(oneToOneMap);

  // Set values in oneToOneVec:  each entry == rank
  for (lno_t i = 0; i < lno_t(oneToOneMap->getNodeNumElements()); i++)
    oneToOneVec.replaceLocalValue(i, me);

  // Now import oneToOneVec's values back to vecWithCopies
  Teuchos::RCP<const Tpetra::Import<lno_t, gno_t> > importer =
      Tpetra::createImport<lno_t, gno_t>(oneToOneMap, mapWithCopies);
  vecWithCopies.doImport(oneToOneVec, *importer, Tpetra::REPLACE);

  // Print the entries of each vector
  std::cout << me << " ONE TO ONE VEC  ("
                  << oneToOneMap->getGlobalNumElements() << "):  ";
  lno_t nlocal = lno_t(oneToOneMap->getNodeNumElements());
  for (lno_t i = 0; i < nlocal; i++)
    std::cout << "[" << oneToOneMap->getGlobalElement(i) << " "
              << oneToOneVec.getData()[i] << "] ";
  std::cout << std::endl;

  // Should see copied vector values when print VEC WITH COPIES
  std::cout << me << " VEC WITH COPIES ("
                  << mapWithCopies->getGlobalNumElements() << "):  ";
  nlocal = lno_t(mapWithCopies->getNodeNumElements());
  for (lno_t i = 0; i < nlocal; i++)
    std::cout << "[" << mapWithCopies->getGlobalElement(i) << " "
              << vecWithCopies.getData()[i] << "] ";
  std::cout << std::endl;

  return 0;
}
