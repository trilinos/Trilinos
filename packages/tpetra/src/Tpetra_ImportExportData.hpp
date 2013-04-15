// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_IMPORTEXPORTDATA_HPP
#define TPETRA_IMPORTEXPORTDATA_HPP

#include "Tpetra_Distributor.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Import,Export, needed to prevent circular inclusions
  template<class LocalOrdinal, class GlobalOrdinal, class Node> class Import;
  template<class LocalOrdinal, class GlobalOrdinal, class Node> class Export;
#endif

  /// \class ImportExportData
  /// \brief Implementation detail of Import and Export.
  /// \tparam LocalOrdinal Same as the first template parameter of Map.
  /// \tparam GlobalOrdinal Same as the second template parameter of Map.
  /// \tparam Node Same as the third template parameter of Map.
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
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  struct ImportExportData {
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

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


    ~ImportExportData();

    //! Output stream for debug output.
    Teuchos::RCP<Teuchos::FancyOStream> out_;

    /// \brief Copys the import/export data, but reverse the direction of the transfer
    /// as well as reversing the distributor
    Teuchos::RCP<ImportExportData<LocalOrdinal, GlobalOrdinal, Node> > reverseClone();

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
    Teuchos::Array<LocalOrdinal> permuteToLIDs_;

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
    Teuchos::Array<LocalOrdinal> permuteFromLIDs_;

    /// \brief "Incoming" indices.
    ///
    /// This array holds the LIDs of the GIDs that are owned by the
    /// target Map, but not by the source Map.  The target object of
    /// the Import or Export will receive data for these LIDs from
    /// other processes.
    Teuchos::Array<LocalOrdinal> remoteLIDs_;

    /// \brief "Outgoing" local indices.
    ///
    /// This array holds the LIDs of the GIDs that are owned by the
    /// source Map, but not by the target Map.  The source object of
    /// the Import or Export will send data from these LIDs to other
    /// processes.
    Teuchos::Array<LocalOrdinal> exportLIDs_;

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

    //! Source Map of the Import or Export
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > source_;

    //! Target Map of the Import or Export
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > target_;

    //! Communicator over which the source and target objects are distributed.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

    /// \brief Object that actually distributes (sends and receives) data.
    ///
    /// The Import or Export object that controls this
    /// ImportExportData container is responsible for initializing the
    /// Distributor.  The Distributor's constructor just gives it the
    /// communicator; it does not complete initialization.
    Distributor distributor_;

  private:
    //! Copy constructor (declared but not defined, do not use)
    ImportExportData (const ImportExportData<LocalOrdinal,GlobalOrdinal,Node> &rhs);
    //! Assignment operator (declared but not defined, do not use)
    ImportExportData<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const ImportExportData<LocalOrdinal,GlobalOrdinal,Node> & rhs);
  }; // class ImportExportData

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::
  ImportExportData (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target)
    : out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)))
    , numSameIDs_ (0)
    , source_ (source)
    , target_ (target)
    , comm_ (source->getComm())
    , distributor_ (comm_, out_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::
  ImportExportData (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out)
    : out_ (out)
    , numSameIDs_ (0)
    , source_ (source)
    , target_ (target)
    , comm_ (source->getComm())
    , distributor_ (comm_, out_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::
  ImportExportData (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
                    const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : out_ (Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)))
    , numSameIDs_ (0)
    , source_ (source)
    , target_ (target)
    , comm_ (source->getComm())
    , distributor_ (comm_, out_, plist)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::
  ImportExportData (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
                    const Teuchos::RCP<Teuchos::FancyOStream>& out,
                    const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : out_ (out)
    , numSameIDs_ (0)
    , source_ (source)
    , target_ (target)
    , comm_ (source->getComm())
    , distributor_ (comm_, out_, plist)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<ImportExportData<LocalOrdinal, GlobalOrdinal, Node> >
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::reverseClone()
  {
    Teuchos::RCP<ImportExportData<LocalOrdinal,GlobalOrdinal,Node> > tData = rcp(new ImportExportData<LocalOrdinal,GlobalOrdinal,Node>(target_,source_));

    // Things that stay the same
    tData->comm_             = comm_;
    tData->numSameIDs_       = numSameIDs_;

    // Things that reverse
    tData->distributor_      = *distributor_.getReverse();
    tData->permuteToLIDs_    = permuteFromLIDs_;
    tData->permuteFromLIDs_  = permuteToLIDs_;

    // Remotes / exports (easy part)
    tData->exportLIDs_       = remoteLIDs_;
    tData->remoteLIDs_       = exportLIDs_;
    tData->exportPIDs_.resize(tData->exportLIDs_.size());

    // Remotes / exports (hard part) - extract the exportPIDs from the remotes of my distributor
    size_t NumReceives                  = distributor_.getNumReceives();
    ArrayView<const int> ProcsFrom      = distributor_.getImagesFrom();
    ArrayView<const size_t> LengthsFrom = distributor_.getLengthsFrom();

    for(size_t i=0,j=0;i<NumReceives;i++) {
      int pid=ProcsFrom[i];
      for(size_t k=0;k<LengthsFrom[i];k++){
	tData->exportPIDs_[j]=pid;
	j++;
      }    
    }
    return tData;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ImportExportData<LocalOrdinal,GlobalOrdinal,Node>::~ImportExportData() 
  {}

} // namespace Tpetra

#endif // TPETRA_IMPORTEXPORTDATA_HPP
