#ifndef TIFPACK_UNITTESTHELPERS_HPP
#define TIFPACK_UNITTESTHELPERS_HPP

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsGraph.hpp>

namespace tif_utest {
using Tpetra::global_size_t;

typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

inline
Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  return Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > create_tridiag_graph(LocalOrdinal num_rows_per_proc)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_rows_per_proc, indexBase, comm));

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowmap, 0));

  Teuchos::Array<GlobalOrdinal> cols;
  for(GlobalOrdinal g_row = rowmap->getMinGlobalIndex(); g_row<=rowmap->getMaxGlobalIndex(); ++g_row) {
    if (g_row == rowmap->getMinAllGlobalIndex() ||
        g_row == rowmap->getMaxAllGlobalIndex()) cols.resize(2);
    else cols.resize(3);

    size_t coloffset = 0;
    if (g_row > rowmap->getMinAllGlobalIndex()) cols[coloffset++] = g_row-1;
    cols[coloffset++] = g_row;
    if (g_row < rowmap->getMaxAllGlobalIndex()) cols[coloffset++] = g_row+1;

    crsgraph->insertGlobalIndices(g_row, cols());
  }

  crsgraph->fillComplete();

  return crsgraph;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > create_test_graph(LocalOrdinal num_rows_per_proc)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_rows_per_proc, indexBase, comm));

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowmap, 0));

  size_t global_size = rowmap->getGlobalNumElements();
  Teuchos::Array<GlobalOrdinal> cols;
  for(GlobalOrdinal g_row = rowmap->getMinGlobalIndex(); g_row<=rowmap->getMaxGlobalIndex(); ++g_row) {

    if (g_row == rowmap->getMinAllGlobalIndex()) {
      cols.resize(global_size);
      for(size_t i=0; i<global_size; ++i) {
        GlobalOrdinal gcol = g_row + i;
        cols[i] = gcol;
      }
    }
    else {
      cols.resize(2);
      cols[0] = rowmap->getMinAllGlobalIndex();
      cols[1] = g_row;
    }

    crsgraph->insertGlobalIndices(g_row, cols());
  }

  crsgraph->fillComplete();

  return crsgraph;
}

}//namespace tif_utest

#endif

