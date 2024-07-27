// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAEXPORT_DEF_HPP
#define XPETRA_TPETRAEXPORT_DEF_HPP

#include "Xpetra_TpetraExport_decl.hpp"
#include "Tpetra_Distributor.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(const Teuchos::RCP<const map_type>& source,
                 const Teuchos::RCP<const map_type>& target)
  : export_(Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target)))) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(const Teuchos::RCP<const map_type>& source,
                 const Teuchos::RCP<const map_type>& target,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : export_(Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target), plist))) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(const Export<LocalOrdinal, GlobalOrdinal, Node>& rhs)
  : export_(Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rhs)))) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    ~TpetraExport() {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumSameIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumSameIDs");
  return export_->getNumSameIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumPermuteIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumPermuteIDs");
  return export_->getNumPermuteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getPermuteFromLIDs() const {
  XPETRA_MONITOR("TpetraExport::getPermuteFromLIDs");
  return export_->getPermuteFromLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getPermuteToLIDs() const {
  XPETRA_MONITOR("TpetraExport::getPermuteToLIDs");
  return export_->getPermuteToLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumRemoteIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumRemoteIDs");
  return export_->getNumRemoteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteLIDs() const {
  XPETRA_MONITOR("TpetraExport::getRemoteLIDs");
  return export_->getRemoteLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumExportIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumExportIDs");
  return export_->getNumExportIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getExportLIDs() const {
  XPETRA_MONITOR("TpetraExport::getExportLIDs");
  return export_->getExportLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const int>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getExportPIDs() const {
  XPETRA_MONITOR("TpetraExport::getExportPIDs");
  return export_->getExportPIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getSourceMap() const {
  XPETRA_MONITOR("TpetraExport::getSourceMap");
  return toXpetra(export_->getSourceMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getTargetMap() const {
  XPETRA_MONITOR("TpetraExport::getTargetMap");
  return toXpetra(export_->getTargetMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {
  XPETRA_MONITOR("TpetraExport::setDistributorParameters");
  export_->getDistributor().setParameterList(params);
  auto revDistor = export_->getDistributor().getReverse(false);
  if (!revDistor.is_null())
    revDistor->setParameterList(params);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    print(std::ostream& os) const {
  XPETRA_MONITOR("TpetraExport::print");
  export_->print(os);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(
        const RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>>& exp)
  : export_(exp) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getTpetra_Export() const {
  return export_;
}

#ifdef HAVE_XPETRA_EPETRA

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))

// stub implementation for GO=int and NO=EpetraNode
template <>
class TpetraExport<int, int, EpetraNode>
  : public Export<int, int, EpetraNode> {
 public:
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

  //! The specialization of Map used by this class.
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Construct a Export object from the source and target Map.
  TpetraExport(const Teuchos::RCP<const map_type>& source, const Teuchos::RCP<const map_type>& target) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "int",
                                typeid(EpetraNode).name());
  }

  //! Constructor (with list of parameters).
  TpetraExport(const Teuchos::RCP<const map_type>& source,
               const Teuchos::RCP<const map_type>& target,
               const Teuchos::RCP<Teuchos::ParameterList>& plist) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "int",
                                typeid(EpetraNode).name());
  }

  //! Copy constructor.
  TpetraExport(const Export<LocalOrdinal, GlobalOrdinal, Node>& rhs) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "int",
                                typeid(EpetraNode).name());
  }

  //! Destructor.
  ~TpetraExport() {}

  //@}

  //! @name Export Attribute Methods
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

  //! The source Map used to construct this Export.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> getSourceMap() const { return Teuchos::null; }

  //! The target Map used to construct this Export.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> getTargetMap() const { return Teuchos::null; }

  //! Set parameters on the underlying object
  void setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {};

  //@}

  //! @name I/O Methods
  //@{

  //! Print the Export's data to the given output stream.
  void print(std::ostream& os) const { /* noop */
  }

  //@}

  //! @name Xpetra specific
  //@{

  //! TpetraExport constructor to wrap a Tpetra::Export object
  TpetraExport(const RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>>& exp) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "int",
                                typeid(EpetraNode).name());
  }

  RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>> getTpetra_Export() const { return Teuchos::null; }

  //@}

};      // TpetraExport class (specialization for LO=GO=int)
#endif  // #if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT)))

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))

// stub implementation for GO=long long and NO=EpetraNode
template <>
class TpetraExport<int, long long, EpetraNode>
  : public Export<int, long long, EpetraNode> {
 public:
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

  //! The specialization of Map used by this class.
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Construct a Export object from the source and target Map.
  TpetraExport(const Teuchos::RCP<const map_type>& source, const Teuchos::RCP<const map_type>& target) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "long long",
                                typeid(EpetraNode).name());
  }

  //! Constructor (with list of parameters).
  TpetraExport(const Teuchos::RCP<const map_type>& source,
               const Teuchos::RCP<const map_type>& target,
               const Teuchos::RCP<Teuchos::ParameterList>& plist) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "long long",
                                typeid(EpetraNode).name());
  }

  //! Copy constructor.
  TpetraExport(const Export<LocalOrdinal, GlobalOrdinal, Node>& rhs) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "long long",
                                typeid(EpetraNode).name());
  }

  //! Destructor.
  ~TpetraExport() {}

  //@}

  //! @name Export Attribute Methods
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

  //! The source Map used to construct this Export.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> getSourceMap() const { return Teuchos::null; }

  //! The target Map used to construct this Export.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> getTargetMap() const { return Teuchos::null; }

  //! Set parameters on the underlying object
  void setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {};

  //@}

  //! @name I/O Methods
  //@{

  //! Print the Export's data to the given output stream.
  void print(std::ostream& os) const { /* noop */
  }

  //@}

  //! @name Xpetra specific
  //@{

  //! TpetraExport constructor to wrap a Tpetra::Export object
  TpetraExport(const RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>>& exp) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraExport<LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "long long",
                                typeid(EpetraNode).name());
  }

  RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>> getTpetra_Export() const { return Teuchos::null; }

  //@}

};      // TpetraExport class (specialization for GO=long long, NO=EpetraNode)
#endif  // #if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG)))

#endif  // HAVE_XPETRA_EPETRA

}  // namespace Xpetra

#endif  // XPETRA_TPETRAEXPORT_DEF_HPP
