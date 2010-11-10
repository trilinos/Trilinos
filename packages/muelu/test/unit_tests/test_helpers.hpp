#ifndef MUELU_UNITTEST_HELPERS_H
#define MUELU_UNITTEST_HELPERS_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"
#include <iostream>

namespace MueLu_UnitTest {

using Tpetra::global_size_t;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

inline
Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  return Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
}

//
// Function that creates a Tpetra map containing a specified number of local elements per process.
//
template<class LocalOrdinal,class GlobalOrdinal,class Node>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
create_tpetra_map(LocalOrdinal num_elements_per_proc)
{ 
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  
  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;
  
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tmap = Teuchos::rcp(new
Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_elements_per_proc, indexBase, comm));
  
  return tmap;
} //create_tpetra_map()

} //namespace MueLu_UnitTest

#endif //ifndef MUELU_UNITTEST_HELPERS_H
