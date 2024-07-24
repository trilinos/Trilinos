// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_IMPORTEXPORTDATA_DECL_HPP
#define TPETRA_IMPORTEXPORTDATA_DECL_HPP

#include "Tpetra_ImportExportData_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_Distributor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
class ParameterList; // forward declaration
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

  /// \class ImportExportData
  /// \brief Implementation detail of Import and Export.
  ///
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  ///
  /// \warning This class is an implementation detail of Import and
  ///   Export.  It may change or disappear at any time.  Tpetra users
  ///   must not depend on this class.
  ///
  /// Import and Export both require the same data.  We use this class
  /// as a container for those data.  They include incoming ("remote")
  /// and outgoing ("export") local indices (LIDs), LIDs to permute on
  /// the source and target of the Import or Export, and process ranks
  /// ("image IDs") to which to send.
  template<class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  class ImportExportData {
  public:
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

    ImportExportData () = delete;

    /// \brief Constructor
    ///
    /// \param source [in] Source Map of the Import or Export
    /// \param target [in] Target Map of the Import or Export
    ImportExportData (const Teuchos::RCP<const map_type>& source,
                      const Teuchos::RCP<const map_type>& target);

    /// \brief Constructor with output stream.
    ///
    /// \param source [in] Source Map of the Import or Export
    /// \param target [in] Target Map of the Import or Export
    /// \param out [in/out] Output stream (for debugging output)
    ImportExportData (const Teuchos::RCP<const map_type>& source,
                      const Teuchos::RCP<const map_type>& target,
                      const Teuchos::RCP<Teuchos::FancyOStream>& out);

    /// \brief Constructor with ParameterList for Distributor
    ///
    /// \param source [in] Source Map of the Import or Export
    /// \param target [in] Target Map of the Import or Export
    /// \param plist [in/out] List of parameters for the Distributor
    ImportExportData (const Teuchos::RCP<const map_type>& source,
                      const Teuchos::RCP<const map_type>& target,
                      const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Constructor with output stream, and ParameterList for Distributor
    ///
    /// \param source [in] Source Map of the Import or Export
    /// \param target [in] Target Map of the Import or Export
    /// \param out [in/out] Output stream (for debugging output)
    /// \param plist [in/out] List of parameters for the Distributor
    ImportExportData (const Teuchos::RCP<const map_type>& source,
                      const Teuchos::RCP<const map_type>& target,
                      const Teuchos::RCP<Teuchos::FancyOStream>& out,
                      const Teuchos::RCP<Teuchos::ParameterList>& plist);
    //! Destructor
    ~ImportExportData () = default;

    /// \brief Copy the data, but reverse the direction of the
    ///   transfer as well as reversing the Distributor.
    ///
    /// "Reverse the direction of the transfer" means that an Import
    /// becomes an Export in the opposite direction, and vice versa.
    Teuchos::RCP<ImportExportData<LocalOrdinal, GlobalOrdinal, Node> > reverseClone();

    //! Source Map of the Import or Export
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > source_;

    //! Target Map of the Import or Export
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > target_;

    //! Output stream for verbose debugging output.
    Teuchos::RCP<Teuchos::FancyOStream> out_;

    //! Whether to print verbose debugging output.
    bool verbose_ = false;

    using execution_space = typename Node::device_type::execution_space;
    using memory_space =
      ::Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>;
    using device_type = Kokkos::Device<execution_space, memory_space>;

    /// \brief Index of target Map LIDs to which to permute.
    ///
    /// After the initial numSameIDs_ indices which are the same in
    /// both the source and target Map, zero or more global indices
    /// (GIDs) remain.  They exist in both the source and target Maps,
    /// but are in a different order.  Therefore, they may have
    /// different local indices (LIDs), and require permutation.
    ///
    /// For each remaining GIDs g in the target Map, if the source Map
    /// also owns g, then permuteToLIDs_ gets the corresponding LID in
    /// the target Map.
    Kokkos::DualView<LocalOrdinal*, device_type> permuteToLIDs_;

    /// \brief Index of source Map LIDs from which to permute.
    ///
    /// After the initial numSameIDs_ indices which are the same in
    /// both the source and target Map, zero or more global indices
    /// (GIDs) remain.  They exist in both the source and target Maps,
    /// but are in a different order.  Therefore, they may have
    /// different local indices (LIDs), and require permutation.
    ///
    /// For each remaining GID g in the target Map, if the source Map
    /// also owns g, then permuteFromLIDs_ gets the corresponding LID
    /// in the source Map.
    Kokkos::DualView<LocalOrdinal*, device_type> permuteFromLIDs_;

    /// \brief "Incoming" indices.
    ///
    /// This array holds the LIDs of the GIDs that are owned by the
    /// target Map, but not by the source Map.  The target object of
    /// the Import or Export will receive data for these LIDs from
    /// other processes.
    Kokkos::DualView<LocalOrdinal*, device_type> remoteLIDs_;

    //! Whether the remote LIDs are contiguous.
    bool remoteLIDsContiguous_ = false;

    //! Whether the export LIDs are contiguous.
    bool exportLIDsContiguous_ = false;

    /// \brief "Outgoing" local indices.
    ///
    /// This array holds the LIDs of the GIDs that are owned by the
    /// source Map, but not by the target Map.  The source object of
    /// the Import or Export will send data from these LIDs to other
    /// processes.
    Kokkos::DualView<LocalOrdinal*, device_type> exportLIDs_;

    //! Ranks of the processes to which the source object sends data.
    Teuchos::Array<int> exportPIDs_;

    /// \brief Number of initial identical indices.
    ///
    /// The number of initial indices (IDs) that are identical between
    /// the source and target Maps.  This count stops at the first
    /// different ID.
    ///
    /// Note that we didn't specify whether the IDs are global (GID)
    /// or local (LID).  That is because if the two Maps start with
    /// the same sequence of GIDs on the calling process, then those
    /// GIDs map to the same LIDs on the calling process.  Thus, when
    /// we say "ID" in the previous paragraph, we include both GID and
    /// LID.
    size_t numSameIDs_;

    /// \brief Object that actually distributes (sends and receives) data.
    ///
    /// The Import or Export object that controls this
    /// ImportExportData container is responsible for initializing the
    /// Distributor.  The Distributor's constructor just gives it the
    /// communicator; it does not complete initialization.
    Distributor distributor_;

    /// \brief Is this Export or Import locally complete?
    ///
    /// If this is an Export, then do all source Map indices on the
    /// calling process exist on at least one process (not necessarily
    /// this one) in the target Map?
    ///
    /// If this is an Import, then do all target Map indices on the
    /// calling process exist on at least one process (not necessarily
    /// this one) in the source Map?
    ///
    /// It's not necessarily an error for an Export or Import not to
    /// be locally complete on one or more processes.  For example,
    /// this may happen in the common use case of "restriction" --
    /// that is, taking a subset of a large object.  Nevertheless, you
    /// may find this predicate useful for figuring out whether you
    /// set up your Maps in the way that you expect.
    bool isLocallyComplete_;

  private:
    //! Copy constructor (declared but not defined, do not use)
    ImportExportData (const ImportExportData<LocalOrdinal,GlobalOrdinal,Node> &rhs);
    //! Assignment operator (declared but not defined, do not use)
    ImportExportData<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const ImportExportData<LocalOrdinal,GlobalOrdinal,Node> & rhs);
  }; // class ImportExportData

} // namespace Tpetra

#endif // TPETRA_IMPORTEXPORTDATA_DECL_HPP
