#ifndef MUELU_UNITTEST_HELPERS_H
#define MUELU_UNITTEST_HELPERS_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_ConfigDefs.hpp" //TODO: use Cthulhu
#include "Tpetra_DefaultPlatform.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_TpetraMap.hpp"

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"
#include <iostream>

namespace MueLu_UnitTest {

  using Tpetra::global_size_t;
  using Teuchos::RCP;

  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  inline
  RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    return Tpetra::DefaultPlatform::getDefaultPlatform().getComm(); //TODO: use Cthulhu here
  }

  //
  // Function that creates a map containing a specified number of local elements per process.
  //
  template<class LocalOrdinal,class GlobalOrdinal,class Node>
  const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> >
  create_map(LocalOrdinal num_elements_per_proc)
  { 
    RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    const LocalOrdinal indexBase = 0;
  
    //TODO: use CthulhuMapFactory here
    return Teuchos::rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_elements_per_proc, indexBase, comm));
  
  } // create_map()

} // namespace MueLu_UnitTest

#endif // ifndef MUELU_UNITTEST_HELPERS_H
