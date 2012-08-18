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

#ifndef TPETRA_IMPORT_HPP
#define TPETRA_IMPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ImportExportData.hpp"
#include "Tpetra_Distributor.hpp"
#include <iterator>

namespace Tpetra {

  /// \brief Communication plan for data redistribution from a uniquely-owned to a (possibly) multiply-owned distribution.
  ///
  /// Tpetra users should use this class to construct a communication
  /// plan between two data distributions (i.e., two \c Map objects).
  /// The plan can be called repeatedly by computational classes to
  /// perform communication according to the same pattern.
  /// Constructing the plan may be expensive, but it can be reused
  /// inexpensively.
  ///
  /// Tpetra has two classes for data redistribution: \c Import and \c
  /// Export.  \c Import is for redistributing data from a
  /// uniquely-owned distribution to a possibly multiply-owned
  /// distribution.  \c Export is for redistributing data from a
  /// possibly multiply-owned distribution to a uniquely-owned
  /// distribution.
  ///
  /// One use case of Import is bringing in remote source vector data
  /// for a distributed sparse matrix-vector multiply.  The source
  /// vector itself is uniquely owned, but must be brought in into an
  /// overlapping distribution so that each process can compute its
  /// part of the target vector without further communication.
  ///
  /// Epetra separated \c Import and \c Export for performance
  /// reasons.  The implementation is different, depending on which
  /// direction is the uniquely-owned Map.  Tpetra retains this
  /// convention.
  ///
  /// This class is templated on the same template arguments as \c
  /// Map: the local ordinal type (\c LocalOrdinal), the global
  /// ordinal type (\c GlobalOrdinal), and the Kokkos Node type (\c
  /// Node).
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Import: public Teuchos::Describable {

  public:
    //! The specialization of Map used by this class.
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Construct an Import from the source and target Maps.
    ///
    /// \param source [in] The source distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    Import (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target);

    /// \brief Constructor (with list of parameters)
    ///
    /// \param source [in] The source distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param plist [in/out] List of parameters.  Currently passed
    ///   directly to the Distributor that implements communication.
    Import (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Copy constructor.
    ///
    /// \note Currently this only makes a shallow copy of the Import's
    ///   underlying data.
    Import(const Import<LocalOrdinal,GlobalOrdinal,Node> & import);

    //! Destructor.
    ~Import();

    //@}

    //! @name Import Attribute Methods
    //@{

    /// \brief Number of initial identical IDs.
    ///
    /// The number of IDs that are identical between the source and
    /// target Maps, up to the first different ID.
    inline size_t getNumSameIDs() const;

    /// \brief Number of IDs to permute but not to communicate.
    ///
    /// The number of IDs that are local to the calling process, but
    /// not part of the first \c getNumSameIDs() entries.  The Import
    /// will permute these entries locally (without distributed-memory
    /// communication).
    inline size_t getNumPermuteIDs() const;

    //! List of local IDs in the source Map that are permuted.
    inline ArrayView<const LocalOrdinal> getPermuteFromLIDs() const;

    //! List of local IDs in the target Map that are permuted.
    inline ArrayView<const LocalOrdinal> getPermuteToLIDs() const;

    //! Number of entries not on the calling process.
    inline size_t getNumRemoteIDs() const;

    //! List of entries in the target Map to receive from other processes.
    inline ArrayView<const LocalOrdinal> getRemoteLIDs() const;

    //! Number of entries that must be sent by the calling process to other processes.
    inline size_t getNumExportIDs() const;

    //! List of entries in the source Map that will be sent to other processes.
    inline ArrayView<const LocalOrdinal> getExportLIDs() const;

    /// \brief List of processes to which entries will be sent.
    ///
    /// The entry with Local ID \c getExportLIDs()[i] will be sent to
    /// process getExportImageIDs()[i].
    inline ArrayView<const int> getExportImageIDs() const;

    //! The Source Map used to construct this Import object.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getSourceMap() const;

    //! The Target Map used to construct this Import object.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getTargetMap() const;

    //! The \c Distributor that this \c Import object uses to move data.
    inline Distributor & getDistributor() const;

    //! Assignment operator.
    Import<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const Import<LocalOrdinal,GlobalOrdinal,Node>& Source);

    //@}

    //! @name I/O Methods
    //@{

    //! Print method
    virtual void print(std::ostream& os) const;

    //@}

  private:

    RCP<ImportExportData<LocalOrdinal,GlobalOrdinal,Node> > ImportData_;
    RCP<Array<GlobalOrdinal> > remoteGIDs_;

    //! @name Initialization helper functions (called by the constructor)
    //@{

    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and numRemoteIDs_
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, and remoteLIDs_

    /// \brief Compute the necessary receives for the Import.
    ///
    /// This routine fills in the following fields of ImportData_:
    ///
    ///   - numSameIDs_ (the number of consecutive initial GIDs owned
    ///     by both the source and target Maps)
    ///   - permuteToLIDs_ (for each of the remaining GIDs g in the
    ///     target Map, if the source Map also owns g, then
    ///     permuteToLIDs_ gets the corresponding LID in the target,
    ///     and permuteFromLIDs_ gets the corresponding LID in the
    ///     source)
    ///   - permuteFromLIDs_ (see permuteToLIDs_)
    ///   - remoteLIDs_ (the LID of each GID that are owned by the
    ///     target Map but not by the source Map)
    ///
    /// It also fills in the temporary remoteGIDs_ array with the GIDs
    /// that are owned by the target Map but not by the source Map.
    ///
    /// The name for this routine comes from what it does.  It first
    /// finds the GIDs that are the same (representing elements which
    /// require neither communication nor permutation).  Then it finds
    /// permutation IDs (which require permutation, but no
    /// communication, because they are in a possibly different order
    /// in the source and target Maps, but owned by the same process)
    /// and remote IDs (which require communication, because they are
    /// owned by the target Map but not by the source Map).
    ///
    /// This routine does not communicate, except perhaps for the
    /// TPETRA_ABUSE_WARNING (that is only triggered if there are
    /// remote IDs but the source is not distributed).
    void setupSamePermuteRemote();

    /// \brief Compute the send communication plan from the receives.
    ///
    /// This routine is called after \c setupSamePermuteRemote(), if
    /// the source Map is distributed.  It uses the remoteGIDs_
    /// temporary array that was allocated by that routine.  After
    /// this routine completes, the remoteGIDs_ array is no longer
    /// needed.
    ///
    /// Algorithm:
    ///
    /// 1. Identify which GIDs are in the target Map but not in the
    ///    source Map.  These correspond to required receives.  Store
    ///    them for now in remoteGIDs_.  Find the process IDs of the
    ///    remote GIDs to receive.
    ///
    /// 2. Invoke Distributor's createFromRecvs() using the above
    ///    remote GIDs and remote process IDs as input.  This sets up
    ///    the Distributor and computes the send GIDs and process IDs.
    ///
    /// 3. Use the source Map to compute the send LIDs from the send
    ///    GIDs.
    ///
    /// This routine fills in the following fields of ImportData_:
    ///
    ///   - remoteLIDs_
    void setupExport();

    //@}
  };

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::Import(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & source,
                                                  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & target) {
    ImportData_ = rcp(new ImportExportData<LocalOrdinal,GlobalOrdinal,Node>(source, target));
    // call subfunctions
    setupSamePermuteRemote();
    if (source->isDistributed()) {
      setupExport();
    }
    // don't need remoteGIDs_ anymore
    remoteGIDs_ = null;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using Teuchos::rcp;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    ImportData_ = rcp (new data_type (source, target, plist));
    setupSamePermuteRemote();
    if (source->isDistributed()) {
      setupExport();
    }
    remoteGIDs_ = null; // Don't need this anymore
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::Import(const Import<LocalOrdinal,GlobalOrdinal,Node> & import)
  : ImportData_(import.ImportData_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::~Import()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumSameIDs() const {
    return ImportData_->numSameIDs_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumPermuteIDs() const {
    return ImportData_->permuteFromLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteFromLIDs() const {
    return ImportData_->permuteFromLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteToLIDs() const {
    return ImportData_->permuteToLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumRemoteIDs() const {
    return ImportData_->remoteLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getRemoteLIDs() const {
    return ImportData_->remoteLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumExportIDs() const {
    return ImportData_->exportLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportLIDs() const {
    return ImportData_->exportLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const int>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportImageIDs() const {
    return ImportData_->exportImageIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getSourceMap() const {
    return ImportData_->source_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getTargetMap() const {
    return ImportData_->target_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Distributor &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getDistributor() const {
    return ImportData_->distributor_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>&
  Import<LocalOrdinal,GlobalOrdinal,Node>::operator=(const Import<LocalOrdinal,GlobalOrdinal,Node> & source) {
    ImportData_ = source.ImportData_;
    return *this;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::print(std::ostream& os) const {
    using Teuchos::getFancyOStream;
    using Teuchos::rcpFromRef;
    using std::endl;

    ArrayView<const LocalOrdinal> av;
    ArrayView<const int> avi;
    const RCP<const Comm<int> > & comm = getSourceMap()->getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) {
        os << endl;
        if (myImageID == 0) { // this is the root node (only output this info once)
          os << "Import Data Members:" << endl;
        }
        os << "Image ID       : " << myImageID << endl;

        os << "permuteFromLIDs: "; os << toString (getPermuteFromLIDs()) << endl;
        //av = getPermuteFromLIDs(); std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << "}" << endl;

        os << "permuteToLIDs  : ";
        os << toString (getPermuteToLIDs()) << endl;
        //av = getPermuteToLIDs();   std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << "}" << endl;

        os << "remoteLIDs     : ";
        os << toString (getRemoteLIDs()) << endl;
        //av = getRemoteLIDs();      std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << "}" << endl;

        os << "exportLIDs     : ";
        os << toString (getExportLIDs()) << endl;
        //av = getExportLIDs();      std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << "}" << endl;

        os << "exportImageIDs : ";
        os << toString (getExportImageIDs()) << endl;
        //avi = getExportImageIDs();  std::copy(avi.begin(),avi.end(),std::ostream_iterator<int>(os," ")); os << "}" << endl;

        os << "numSameIDs     : " << getNumSameIDs() << endl;
        os << "numPermuteIDs  : " << getNumPermuteIDs() << endl;
        os << "numRemoteIDs   : " << getNumRemoteIDs() << endl;
        os << "numExportIDs   : " << getNumExportIDs() << endl;
      }

      // A few global barriers give I/O a chance to complete.
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }

    const bool printMaps = false;
    if (printMaps) {
      if (myImageID == 0) {
        os << endl << endl << "Source Map:" << endl << std::flush;
      }
      comm->barrier();
      os << *getSourceMap();
      comm->barrier();

      if (myImageID == 0) {
        os << endl << endl << "Target Map:" << endl << std::flush;
      }
      comm->barrier();
      os << *getTargetMap();
      comm->barrier();
    }

    // It's also helpful for debugging to print the Distributor
    // object.  Epetra_Import::Print() does this, so we can do a
    // side-by-side comparison.
    if (myImageID == 0) {
      os << endl << endl << "Distributor:" << endl << std::flush;
    }
    comm->barrier();
    getDistributor().describe (*(getFancyOStream (rcpFromRef (os))),
                               Teuchos::VERB_EXTREME);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::setupSamePermuteRemote() {
    const Map<LocalOrdinal,GlobalOrdinal,Node> & source = *getSourceMap();
    const Map<LocalOrdinal,GlobalOrdinal,Node> & target = *getTargetMap();
    ArrayView<const GlobalOrdinal> sourceGIDs = source.getNodeElementList();
    ArrayView<const GlobalOrdinal> targetGIDs = target.getNodeElementList();

    // Compute numSameIDs_:
    //
    // Iterate through the source and target GID lists.  If the i-th
    // GID of both is the same, increment numSameIDs_ and try the
    // next.  As soon as you come to a nonmatching pair, give up.
    //
    // The point of numSameIDs_ is for the common case of an Import
    // where all the overlapping GIDs are at the end of the source
    // Map, but otherwise the source and target Maps are the same.
    // This allows a fast contiguous copy for the initial "same IDs."
    typename ArrayView<const GlobalOrdinal>::iterator sourceIter = sourceGIDs.begin(),
                                                      targetIter = targetGIDs.begin();
    while (sourceIter != sourceGIDs.end() && targetIter != targetGIDs.end() && *sourceIter == *targetIter) {
      ++ImportData_->numSameIDs_;
      ++sourceIter;
      ++targetIter;
    }
    // targetIter should now point either to the GID of the first
    // non-same entry in targetGIDs, or to the end of targetGIDs (if
    // all the entries were the same).

    // Compute IDs to be permuted, vs. IDs to be received ("imported";
    // called "remote" IDs).
    //
    // IDs to permute are in both the source and target Maps, which
    // means we don't have to send or receive them, but we do have to
    // rearrange (permute) them in general.  (We've already identified
    // an initial stretch of IDs which can be copied without
    // rearrangement; the iterator targetIter is past that point.)
    // IDs to receive are in the target Map, but not the source Map.
    //
    // How do the following code and its equivalent in Export differ?
    //
    // 1. Export uses sourceIter, whereas Import uses targetIter.
    //
    // 2. Import collects remoteGIDs_ (target Map GIDs that are not in
    //    the source Map), which is a separate array.  Export can use
    //    exportGIDs_, which is an array belonging to its
    //    ImportExportData object.
    remoteGIDs_ = rcp( new Array<GlobalOrdinal>() );
    for (; targetIter != targetGIDs.end(); ++targetIter) {
      const GlobalOrdinal curTargetGID = *targetIter;
      if (source.isNodeGlobalElement (curTargetGID)) {
        // The current process owns this GID, for both the source and
        // the target Maps.  Determine the LIDs for this GID on both
        // Maps and add them to the permutation lists.
        ImportData_->permuteToLIDs_.push_back (target.getLocalElement (curTargetGID));
        ImportData_->permuteFromLIDs_.push_back (source.getLocalElement (curTargetGID));
      }
      else {
        // The current GID is owned by this process in the target Map,
        // but is not owned by this process in the source Map.  That
        // means the Import operation has to receive it from another
        // process.  Store it in the "remote" (incoming) list, along
        // with its destination LID on this process.
        //
        // remoteLIDs_ is the list of this process' LIDs that it has
        // to receive from other processes.  Since this is an Import,
        // and therefore the source Map is nonoverlapping, we know
        // that each remote LID can receive from only one process.
        remoteGIDs_->push_back (curTargetGID);
        ImportData_->remoteLIDs_.push_back (target.getLocalElement (curTargetGID));
      }
    }

    TPETRA_ABUSE_WARNING(
      getNumRemoteIDs() > 0 && ! source.isDistributed(),
      std::runtime_error,
      "::setupSamePermuteRemote(): Target has remote LIDs but Source is not "
      "distributed globally." << std::endl
      << "Importing to a submap of the target map.");
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::setupExport() {
    typedef typename Array<int>::difference_type size_type;
    const Map<LocalOrdinal,GlobalOrdinal,Node> & source = *getSourceMap();

    // For each entry remoteGIDs[i], remoteImageIDs[i] will contain
    // the process ID of the process that owns that GID.
    ArrayView<GlobalOrdinal> remoteGIDs = (*remoteGIDs_)();
    Array<int> remoteImageIDs(remoteGIDs.size());
    // lookup == IDNotPresent means that the source Map wasn't able to
    // figure out to which processes one or more of the GIDs in the
    // given list of remoteGIDs belong.
    //
    // The previous abuse warning said "The target Map has GIDs not
    // found in the source Map."  This statement could be confusing,
    // because it doesn't refer to ownership by the current process,
    // but rather to ownership by _any_ process participating in the
    // Map.  (It could not possibly refer to ownership by the current
    // process, since remoteGIDs is exactly the list of GIDs owned by
    // the target Map but not owned by the source Map.  It was
    // constructed that way by setupSamePermuteRemote().)
    //
    // What this statement means is that the source and target Maps
    // don't contain the same set of GIDs globally (over all
    // processes).  That is, there is at least one GID owned by some
    // process in the target Map, which is not owned by _any_ process
    // in the source Map.
    const LookupStatus lookup = source.getRemoteIndexList(remoteGIDs, remoteImageIDs());
    TPETRA_ABUSE_WARNING( lookup == IDNotPresent, std::runtime_error,
      "::setupExport(): the source Map wasn't able to figure out which process "
      "owns one or more of the GIDs in the list of remote GIDs.  This probably "
      "means that there is at least one GID owned by some process in the target"
      " Map which is not owned by any process in the source Map.  (That is, the"
      " source and target Maps do not contain the same set of GIDs globally.)");

    // Ignore remote GIDs that aren't owned by any process in the
    // source Map.  getRemoteIndexList() gives each of these a process
    // ID of -1.
    if ( lookup == IDNotPresent ) {
      const size_type numInvalidRemote = std::count_if( remoteImageIDs.begin(), remoteImageIDs.end(), std::bind1st(std::equal_to<int>(),-1) );
      // if all of them are invalid, we can delete the whole array
      const size_type totalNumRemote = getNumRemoteIDs();
      if (numInvalidRemote == totalNumRemote) {
        // all remotes are invalid; we have no remotes; we can delete the remotes
        remoteImageIDs.clear();
        (*remoteGIDs_).clear();
        ImportData_->remoteLIDs_.clear();
      }
      else {
        // Some remotes are valid; we need to keep the valid ones.
        // Pack and resize remoteImageIDs, remoteGIDs_, and
        // remoteLIDs_.
        size_type numValidRemote = 0;
        for (size_type r = 0; r < totalNumRemote; ++r) {
          if (remoteImageIDs[r] != -1) {
            remoteImageIDs[numValidRemote] = remoteImageIDs[r];
            (*remoteGIDs_)[numValidRemote] = (*remoteGIDs_)[r];
            ImportData_->remoteLIDs_[numValidRemote] = ImportData_->remoteLIDs_[r];
            ++numValidRemote;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          numValidRemote != totalNumRemote - numInvalidRemote, std::logic_error,
          "Tpetra::Import::setupExport(): After removing invalid remote GIDs and"
          " packing the valid remote GIDs, numValidRemote = " << numValidRemote
          << " != totalNumRemote - numInvalidRemote = "
          << totalNumRemote - numInvalidRemote
          << ".  Please report this bug to the Tpetra developers.");

        remoteImageIDs.resize(numValidRemote);
        (*remoteGIDs_).resize(numValidRemote);
        ImportData_->remoteLIDs_.resize(numValidRemote);
      }
      remoteGIDs = (*remoteGIDs_)();
    }

    // Sort remoteImageIDs in ascending order, and apply the resulting
    // permutation to remoteGIDs_ and remoteLIDs_.  This ensures that
    // remoteImageIDs[i], remoteGIDs_[i], and remoteLIDs_[i] all refer
    // to the same thing.
    sort3 (remoteImageIDs.begin(),
           remoteImageIDs.end(),
           remoteGIDs.begin(),
           ImportData_->remoteLIDs_.begin());

    // Call the Distributor's createFromRecvs() method to turn the
    // remote GIDs and their owning processes into a send-and-receive
    // communication plan.  remoteGIDs and remoteImageIDs_ are input;
    // exportGIDs and exportImageIDs_ are output arrays which are
    // allocated by createFromRecvs().
    ArrayRCP<GlobalOrdinal> exportGIDs;
    ImportData_->distributor_.createFromRecvs(remoteGIDs().getConst(), remoteImageIDs, exportGIDs, ImportData_->exportImageIDs_);

    // Find the LIDs corresponding to the (outgoing) GIDs in
    // exportGIDs.  For sparse matrix-vector multiply, this tells the
    // calling process how to index into the source vector to get the
    // elements which it needs to send.
    if (exportGIDs != null) {
      ImportData_->exportLIDs_ = arcp<LocalOrdinal>(exportGIDs.size());
    }
    typename ArrayRCP<LocalOrdinal>::iterator dst = ImportData_->exportLIDs_.begin();
    typename ArrayRCP<GlobalOrdinal>::const_iterator src = exportGIDs.begin();
    while (src != exportGIDs.end()) {
      (*dst++) = source.getLocalElement(*src++);
    }
  }

  /** \brief Non-member constructor for Import objects.

      Creates a Import object from the given source and target maps.
      \pre <tt>src != null</tt>
      \pre <tt>tgt != null</tt>
      \return Returns the Import object. If <tt>src == tgt</tt>, returns \c null. (Debug mode: throws std::runtime_error if one of \c src or \c tgt is \c null.)

      \relatesalso Import
    */
  template <class LO, class GO, class Node>
  RCP< const Import<LO,GO,Node> >
  createImport( const RCP<const Map<LO,GO,Node> > & src,
                const RCP<const Map<LO,GO,Node> > & tgt )
  {
    if (src == tgt) return null;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(src == null || tgt == null, std::runtime_error,
        "Tpetra::createImport(): neither source nor target map may be null:\nsource: " << src << "\ntarget: " << tgt << "\n");
#endif
    return rcp(new Import<LO,GO,Node>(src,tgt));
  }

} // namespace Tpetra

#endif // TPETRA_IMPORT_HPP
