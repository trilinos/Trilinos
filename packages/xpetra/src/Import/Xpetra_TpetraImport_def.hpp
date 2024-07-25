// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAIMPORT_DEF_HPP
#define XPETRA_TPETRAIMPORT_DEF_HPP
#include "Xpetra_TpetraConfigDefs.hpp"

#include "Xpetra_Import.hpp"
#include "Xpetra_TpetraImport_decl.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_TpetraMap.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Distributor.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target)
  : import_(Teuchos::rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target)))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target, const Teuchos::RCP<Teuchos::ParameterList> &plist)
  : import_(Teuchos::rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target), plist))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const Import<LocalOrdinal, GlobalOrdinal, Node> &import)
  : import_(Teuchos::rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(import)))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::~TpetraImport() {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::createRemoteOnlyImport(const Teuchos::RCP<const map_type> &remoteTarget) const {
  Teuchos::RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > newImport = import_->createRemoteOnlyImport(toTpetra(remoteTarget));
  return Teuchos::rcp(new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(newImport));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumSameIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumSameIDs");
  return import_->getNumSameIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumPermuteIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumPermuteIDs");
  return import_->getNumPermuteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getPermuteFromLIDs() const {
  XPETRA_MONITOR("TpetraImport::getPermuteFromLIDs");
  return import_->getPermuteFromLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getPermuteToLIDs() const {
  XPETRA_MONITOR("TpetraImport::getPermuteToLIDs");
  return import_->getPermuteToLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumRemoteIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumRemoteIDs");
  return import_->getNumRemoteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {
  XPETRA_MONITOR("TpetraImport::setDistributorParameters");
  import_->getDistributor().setParameterList(params);
  auto revDistor = import_->getDistributor().getReverse(false);
  if (!revDistor.is_null())
    revDistor->setParameterList(params);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getRemoteLIDs() const {
  XPETRA_MONITOR("TpetraImport::getRemoteLIDs");
  return import_->getRemoteLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumExportIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumExportIDs");
  return import_->getNumExportIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getExportLIDs() const {
  XPETRA_MONITOR("TpetraImport::getExportLIDs");
  return import_->getExportLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const int> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getExportPIDs() const {
  XPETRA_MONITOR("TpetraImport::getExportPIDs");
  return import_->getExportPIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getSourceMap() const {
  XPETRA_MONITOR("TpetraImport::getSourceMap");
  return toXpetra(import_->getSourceMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getTargetMap() const {
  XPETRA_MONITOR("TpetraImport::getTargetMap");
  return toXpetra(import_->getTargetMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::print(std::ostream &os) const {
  XPETRA_MONITOR("TpetraImport::print");
  import_->print(os);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import)
  : import_(import) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getTpetra_Import() const { return import_; }

#ifdef HAVE_XPETRA_EPETRA

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))

// stub implementation for GO=int and NO=EpetraNode
template <>
class TpetraImport<int, int, EpetraNode> : public Import<int, int, EpetraNode> {
 public:
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

  //! The specialization of Map used by this class.
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Construct an Import from the source and target Maps.
  TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "int", typeid(EpetraNode).name());
  }

  //! Constructor (with list of parameters).
  TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target, const Teuchos::RCP<Teuchos::ParameterList> &plist) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "int", typeid(EpetraNode).name());
  }

  //! Copy constructor.
  TpetraImport(const Import<LocalOrdinal, GlobalOrdinal, Node> &import) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "int", typeid(EpetraNode).name());
  }

  //! Destructor.
  ~TpetraImport() {}

  //! Special "constructor"
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  createRemoteOnlyImport(const Teuchos::RCP<const map_type> &remoteTarget) const {
    return Teuchos::null;
  }

  //@}

  //! @name Import Attribute Methods
  //@{

  //! Number of initial identical IDs.
  size_t getNumSameIDs() const { return 0; }

  //! Number of IDs to permute but not to communicate.
  size_t getNumPermuteIDs() const { return 0; }

  //! List of local IDs in the source Map that are permuted.
  ArrayView<const LocalOrdinal> getPermuteFromLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! List of local IDs in the target Map that are permuted.
  ArrayView<const LocalOrdinal> getPermuteToLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! Number of entries not on the calling process.
  size_t getNumRemoteIDs() const { return 0; }

  //! List of entries in the target Map to receive from other processes.
  ArrayView<const LocalOrdinal> getRemoteLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! Number of entries that must be sent by the calling process to other processes.
  size_t getNumExportIDs() const { return 0; }

  //! List of entries in the source Map that will be sent to other processes.
  ArrayView<const LocalOrdinal> getExportLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! List of processes to which entries will be sent.
  ArrayView<const int> getExportPIDs() const { return Teuchos::ArrayView<const int>(); }

  //! The Source Map used to construct this Import object.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getSourceMap() const { return Teuchos::null; }

  //! The Target Map used to construct this Import object.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getTargetMap() const { return Teuchos::null; }

  //! Set parameters on the underlying object
  void setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {}

  //@}

  //! @name I/O Methods
  //@{

  //! Print the Import's data to the given output stream.
  void print(std::ostream &os) const { /* noop */
  }

  //@}

  //! @name Xpetra specific
  //@{

  //! TpetraImport constructor to wrap a Tpetra::Import object
  TpetraImport(const RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "int", typeid(EpetraNode).name());
  }

  RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_Import() const { return Teuchos::null; }

  //@}

};  // TpetraImport class (stub implementation for GO=int, NO=EpetraNode)
#endif

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))

// stub implementation for GO=long long and NO=EpetraNode
template <>
class TpetraImport<int, long long, EpetraNode> : public Import<int, long long, EpetraNode> {
 public:
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

  //! The specialization of Map used by this class.
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Construct an Import from the source and target Maps.
  TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "long long", typeid(EpetraNode).name());
  }

  //! Constructor (with list of parameters).
  TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target, const Teuchos::RCP<Teuchos::ParameterList> &plist) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "long long", typeid(EpetraNode).name());
  }

  //! Copy constructor.
  TpetraImport(const Import<LocalOrdinal, GlobalOrdinal, Node> &import) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "long long", typeid(EpetraNode).name());
  }

  //! Destructor.
  ~TpetraImport() {}

  //! Special "constructor"
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  createRemoteOnlyImport(const Teuchos::RCP<const map_type> &remoteTarget) const {
    return Teuchos::null;
  }

  //@}

  //! @name Import Attribute Methods
  //@{

  //! Number of initial identical IDs.
  size_t getNumSameIDs() const { return 0; }

  //! Number of IDs to permute but not to communicate.
  size_t getNumPermuteIDs() const { return 0; }

  //! List of local IDs in the source Map that are permuted.
  ArrayView<const LocalOrdinal> getPermuteFromLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! List of local IDs in the target Map that are permuted.
  ArrayView<const LocalOrdinal> getPermuteToLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! Number of entries not on the calling process.
  size_t getNumRemoteIDs() const { return 0; }

  //! List of entries in the target Map to receive from other processes.
  ArrayView<const LocalOrdinal> getRemoteLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! Number of entries that must be sent by the calling process to other processes.
  size_t getNumExportIDs() const { return 0; }

  //! List of entries in the source Map that will be sent to other processes.
  ArrayView<const LocalOrdinal> getExportLIDs() const { return Teuchos::ArrayView<const LocalOrdinal>(); }

  //! List of processes to which entries will be sent.
  ArrayView<const int> getExportPIDs() const { return Teuchos::ArrayView<const int>(); }

  //! The Source Map used to construct this Import object.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getSourceMap() const { return Teuchos::null; }

  //! The Target Map used to construct this Import object.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getTargetMap() const { return Teuchos::null; }

  //! Set parameters on the underlying object
  void setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {}

  //@}

  //! @name I/O Methods
  //@{

  //! Print the Import's data to the given output stream.
  void print(std::ostream &os) const { /* noop */
  }

  //@}

  //! @name Xpetra specific
  //@{

  //! TpetraImport constructor to wrap a Tpetra::Import object
  TpetraImport(const RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), typeid(TpetraImport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(), "long long", typeid(EpetraNode).name());
  }

  RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_Import() const { return Teuchos::null; }

  //@}

};  // TpetraImport class (stub implementation for GO=long long, NO=EpetraNode)
#endif

#endif  // HAVE_XPETRA_EPETRA

}  // namespace Xpetra

#endif
