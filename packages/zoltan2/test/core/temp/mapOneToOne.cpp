// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

///  Small test program showing how to take GIDs that may have
///  duplicates across processors (e.g., mesh vertices that are copied
///  at part boundaries in an element-based decomposition) and assign
///  unique owners to them.

#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

#include <string>
#include <sstream>
#include <iostream>

/////////////////////////////////////////////////////////////////////

int main(int narg, char **arg)
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  typedef Tpetra::Map<> map_t;
  typedef map_t::local_ordinal_type lno_t;
  typedef map_t::global_ordinal_type gno_t;

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


  // Print the entries of each map
  std::cout << me << " MAP WITH COPIES ("
                  << mapWithCopies->getGlobalNumElements() << "):  ";
  lno_t nlocal = lno_t(mapWithCopies->getLocalNumElements());
  for (lno_t i = 0; i < nlocal; i++)
    std::cout << mapWithCopies->getGlobalElement(i) << " ";
  std::cout << std::endl;

  std::cout << me << " ONE TO ONE MAP  ("
                  << oneToOneMap->getGlobalNumElements() << "):  ";
  nlocal = lno_t(oneToOneMap->getLocalNumElements());
  for (lno_t i = 0; i < nlocal; i++)
    std::cout << oneToOneMap->getGlobalElement(i) << " ";
  std::cout << std::endl;

  return 0;
}
