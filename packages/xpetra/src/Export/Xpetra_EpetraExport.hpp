// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EPETRAEXPORT_HPP
#define XPETRA_EPETRAEXPORT_HPP

#include "Xpetra_EpetraConfigDefs.hpp"

#include "Xpetra_Export.hpp"

#include "Xpetra_EpetraMap.hpp"  //TMP

#include "Epetra_Export.h"

#if defined(XPETRA_ENABLE_DEPRECATED_CODE)
#ifdef __GNUC__
#if defined(Xpetra_SHOW_DEPRECATED_WARNINGS)
#warning "The header file Trilinos/packages/xpetra/src/Export/Xpetra_EpetraExport.hpp is deprecated."
#endif
#endif
#else
#error "The header file Trilinos/packages/xpetra/src/Export/Xpetra_EpetraExport.hpp is deprecated."
#endif

// Note: 'export' is a reserved keyword in C++. Do not use 'export' as a variable name.

namespace Xpetra {

// TODO: move that elsewhere
template <class GlobalOrdinal, class Node>
XPETRA_DEPRECATED const Epetra_Export &toEpetra(const Export<int, GlobalOrdinal, Node> &);
template <class GlobalOrdinal, class Node>
XPETRA_DEPRECATED RCP<const Export<int, GlobalOrdinal, Node> > toXpetra(const Epetra_Export *exp);

template <class EpetraGlobalOrdinal, class Node>
class XPETRA_DEPRECATED EpetraExportT
  : public Export<int, EpetraGlobalOrdinal, Node> {
  typedef int LocalOrdinal;
  typedef EpetraGlobalOrdinal GlobalOrdinal;
  //! The specialization of Map used by this class.
  typedef Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Construct a Export object from the source and target Map.
  EpetraExportT(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target)
    : export_(rcp(new Epetra_Export(toEpetra<GlobalOrdinal, Node>(source), toEpetra<GlobalOrdinal, Node>(target)))) {}  // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

  ////! Constructor (with list of parameters).
  // Definition not in cpp, so comment out
  // EpetraExportT(const Teuchos::RCP< const map_type > &source, const Teuchos::RCP< const map_type > &target, const Teuchos::RCP< Teuchos::ParameterList > &plist);

  ////! Copy constructor.
  // Definition not in cpp, so comment out
  // EpetraExportT(const Export< LocalOrdinal, GlobalOrdinal, Node > &rhs);

  //! Destructor.
  ~EpetraExportT() {}

  //@}

  //! @name Export Attribute Methods
  //@{

  //! Number of initial identical IDs.
  size_t getNumSameIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumSameIDs");
    return export_->NumSameIDs();
  }

  //! Number of IDs to permute but not to communicate.
  size_t getNumPermuteIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumPermuteIDs");
    return export_->NumPermuteIDs();
  }

  //! List of local IDs in the source Map that are permuted.
  ArrayView<const LocalOrdinal> getPermuteFromLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getPermuteFromLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getExportImageIDs not implemented");
  }

  //! List of local IDs in the target Map that are permuted.
  ArrayView<const LocalOrdinal> getPermuteToLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getPermuteToLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getPermuteToLIDs not implemented");
  }

  //! Number of entries not on the calling process.
  size_t getNumRemoteIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumRemoteIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getNumRemoteIDs not implemented");
  }

  //! List of entries in the target Map to receive from other processes.
  ArrayView<const LocalOrdinal> getRemoteLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getRemoteLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getRemoteLIDs not implemented");
  }

  //! Number of entries that must be sent by the calling process to other processes.
  size_t getNumExportIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumExportIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getNumExportIDs not implemented");
  }

  //! List of entries in the source Map that will be sent to other processes.
  ArrayView<const LocalOrdinal> getExportLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getExportLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getExportLIDs not implemented");
  }

  //! List of processes to which entries will be sent.
  ArrayView<const int> getExportPIDs() const {
    XPETRA_MONITOR("EpetraExportT::getExportImageIDs");
    return ArrayView<const int>(export_->ExportPIDs(), export_->NumExportIDs());
  }

  //! The source Map used to construct this Export.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getSourceMap() const {
    XPETRA_MONITOR("EpetraExportT::getSourceMap");
    return toXpetra<GlobalOrdinal, Node>(export_->SourceMap());
  }

  //! The target Map used to construct this Export.
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getTargetMap() const {
    XPETRA_MONITOR("EpetraExportT::getTargetMap");
    return toXpetra<GlobalOrdinal, Node>(export_->TargetMap());
  }

  //! Set parameters on the underlying object
  void setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const { XPETRA_MONITOR("EpetraExportT::setDistributorParameters"); }

  //@}

  //! @name I/O Methods
  //@{

  //! Print the Export's data to the given output stream.
  void print(std::ostream &os) const {
    XPETRA_MONITOR("EpetraExportT::");
    export_->Print(os);
  }

  //@}

  //! @name Xpetra specific
  //@{

  //! EpetraExportT constructor to wrap a Epetra_Export object
  EpetraExportT(const RCP<const Epetra_Export> &exp)
    : export_(exp) {}

  //! Get the underlying Epetra export
  RCP<const Epetra_Export> getEpetra_Export() const { return export_; }

  //@}

 private:
  RCP<const Epetra_Export> export_;

};  // EpetraExportT class

}  // namespace Xpetra

#endif  // XPETRA_EPETRAEXPORT_HPP
