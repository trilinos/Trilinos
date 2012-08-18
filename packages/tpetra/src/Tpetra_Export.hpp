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

#ifndef TPETRA_EXPORT_HPP
#define TPETRA_EXPORT_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ImportExportData.hpp"
#include <iterator>

namespace Tpetra {

  /// \brief Communication plan for data redistribution from a (possibly) multiply-owned to a uniquely-owned distribution.
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
  /// One use case of Export is finite element assembly.  For example,
  /// one way to compute a distributed forcing term vector is to use
  /// an overlapping distribution for the basis functions' domains.
  /// An Export with the SUM combine mode combines each process'
  /// contribution to the integration into a single nonoverlapping
  /// distribution.
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
  class Export: public Teuchos::Describable {

  public:
    //! The specialization of Map used by this class.
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Construct a Export object from the source and target Map.
    ///
    /// \param source [in] The source distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    Export (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target);

    /// \brief Constructor (with list of parameters)
    ///
    /// \param source [in] The source distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param plist [in/out] List of parameters.  Currently passed
    ///   directly to the Distributor that implements communication.
    Export (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Copy constructor.
    ///
    /// \note Currently this only makes a shallow copy of the Export's
    ///   underlying data.
    Export (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs);

    //! Destructor.
    ~Export();

    //@}

    //! @name Export Attribute Methods
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
    /// The entry with local ID getExportLIDs()[i] will be sent to
    /// process getExportImageIDs()[i].
    inline ArrayView<const int> getExportImageIDs() const;

    //! The source \c Map used to construct this Export.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getSourceMap() const;

    //! The target \c Map used to construct this Export.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getTargetMap() const;

    //! The Distributor that this \c Export object uses to move data.
    inline Distributor & getDistributor() const;

    //! Assignment operator
    Export<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs);

    //@}

    //! @name I/O Methods
    //@{

    /// \brief Print the Export's data to the given output stream.
    ///
    /// This method assumes that the given output stream can be
    /// written on all process(es) in the Export's communicator.  The
    /// resulting output is useful mainly for debugging.
    ///
    /// \note This method tries its best (by using barriers at the end
    ///   of each iteration of a for loop over all communicator ranks)
    ///   to ensure ordered deterministic output.  However, the
    ///   assumption that all processes can write to the stream means
    ///   that there are no ordering guarantees other than what the
    ///   operating and run-time system provide.  (MPI synchronization
    ///   may be separate from output stream synchronization, so the
    ///   barriers only improve the chances that output can complete
    ///   before the next process starts writing.)
    virtual void print (std::ostream& os) const;

    //@}

  private:

    RCP<ImportExportData<LocalOrdinal,GlobalOrdinal,Node> > ExportData_;

    //! @name Initialization helper functions (called by the constructor)
    //@{

    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and the export IDs
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, exportGIDs_, and exportLIDs_
    void setupSamePermuteExport();
    void setupRemote();

    //@}
  };


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target)
  {
    using Teuchos::rcp;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    ExportData_ = rcp (new data_type (source, target));
    setupSamePermuteExport();
    if (source->isDistributed()) {
      setupRemote();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::
  Export (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    using Teuchos::rcp;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    ExportData_ = rcp (new data_type (source, target, plist));
    setupSamePermuteExport();
    if (source->isDistributed()) {
      setupRemote();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::Export(const Export<LocalOrdinal,GlobalOrdinal,Node> & rhs)
  : ExportData_(rhs.ExportData_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>::~Export()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumSameIDs() const {
    return ExportData_->numSameIDs_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumPermuteIDs() const {
    return ExportData_->permuteFromLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getPermuteFromLIDs() const {
    return ExportData_->permuteFromLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getPermuteToLIDs() const {
    return ExportData_->permuteToLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumRemoteIDs() const {
    return ExportData_->remoteLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getRemoteLIDs() const {
    return ExportData_->remoteLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Export<LocalOrdinal,GlobalOrdinal,Node>::getNumExportIDs() const {
    return ExportData_->exportLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getExportLIDs() const {
    return ExportData_->exportLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const int>
  Export<LocalOrdinal,GlobalOrdinal,Node>::getExportImageIDs() const {
    return ExportData_->exportImageIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Export<LocalOrdinal,GlobalOrdinal,Node>::getSourceMap() const {
    return ExportData_->source_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Export<LocalOrdinal,GlobalOrdinal,Node>::getTargetMap() const {
    return ExportData_->target_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Distributor &
  Export<LocalOrdinal,GlobalOrdinal,Node>::getDistributor() const {
    return ExportData_->distributor_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Export<LocalOrdinal,GlobalOrdinal,Node>&
  Export<LocalOrdinal,GlobalOrdinal,Node>::operator=(const Export<LocalOrdinal,GlobalOrdinal,Node> & rhs) {
    if (&rhs != this) {
      // It's bad form to clobber your own data in a self-assignment.
      // This can result in dangling pointers if some member data are
      // raw pointers that the class deallocates in the constructor.
      // It doesn't matter in this case, because ExportData_ is an
      // RCP, which defines self-assignment sensibly.  Nevertheless,
      // we include the check for self-assignment, because it's good
      // form and not expensive (just a raw pointer comparison).
      ExportData_ = rhs.ExportData_;
    }
    return *this;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Export<LocalOrdinal,GlobalOrdinal,Node>::print(std::ostream& os) const {
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
          os << "Export Data Members:" << endl;
        }
        os << "Image ID       : " << myImageID << endl;
        os << "permuteFromLIDs: {"; av = getPermuteFromLIDs(); std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << " }" << endl;
        os << "permuteToLIDs  : {"; av = getPermuteToLIDs();   std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << " }" << endl;
        os << "remoteLIDs     : {"; av = getRemoteLIDs();      std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << " }" << endl;
        os << "exportLIDs     : {"; av = getExportLIDs();      std::copy(av.begin(),av.end(),std::ostream_iterator<LocalOrdinal>(os," ")); os << " }" << endl;
        os << "exportImageIDs : {"; avi = getExportImageIDs();  std::copy(avi.begin(),avi.end(),std::ostream_iterator<int>(os," ")); os << " }" << endl;
        os << "numSameIDs     : " << getNumSameIDs() << endl;
        os << "numPermuteIDs  : " << getNumPermuteIDs() << endl;
        os << "numRemoteIDs   : " << getNumRemoteIDs() << endl;
        os << "numExportIDs   : " << getNumExportIDs() << endl;
      }
      // Do a few global ops to give I/O a chance to complete
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }
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

    // It's also helpful for debugging to print the Distributor
    // object.  Epetra_Import::Print() does this (or _should_ do this,
    // but doesn't, as of 05 Jan 2012), so we can do a side-by-side
    // comparison.
    if (myImageID == 0) {
      os << endl << endl << "Distributor:" << endl << std::flush;
    }
    comm->barrier();
    getDistributor().describe (*(getFancyOStream (rcpFromRef (os))),
                               Teuchos::VERB_EXTREME);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Export<LocalOrdinal,GlobalOrdinal,Node>::setupSamePermuteExport()
  {
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
    // The point of numSameIDs_ is for the common case of an Export
    // where all the overlapping GIDs are at the end of the target
    // Map, but otherwise the source and target Maps are the same.
    // This allows a fast contiguous copy for the initial "same IDs."
    typename ArrayView<const GlobalOrdinal>::iterator sourceIter = sourceGIDs.begin(),
                                                      targetIter = targetGIDs.begin();
    while (sourceIter != sourceGIDs.end() && targetIter != targetGIDs.end() && *sourceIter == *targetIter) {
      ++ExportData_->numSameIDs_;
      ++sourceIter;
      ++targetIter;
    }
    // sourceIter should now point either to the GID of the first
    // non-same entry in sourceGIDs, or to the end of sourceGIDs (if
    // all the entries were the same).

    // Compute IDs to be permuted, vs. IDs to be sent out ("exported";
    // called "export" IDs).
    //
    // IDs to permute are in both the source and target Maps, which
    // means we don't have to send or receive them, but we do have to
    // rearrange (permute) them in general.  (We've already identified
    // an initial stretch of IDs which can be copied without
    // rearrangement; the iterator sourceIter is past that point.)
    // IDs to send out are in the source Map, but not the target Map.
    for (; sourceIter != sourceGIDs.end(); ++sourceIter) {
      const GlobalOrdinal curSourceGID = *sourceIter;
      if (target.isNodeGlobalElement (curSourceGID)) {
        // The current process owns this GID, for both the source and
        // the target Maps.  Determine the LIDs for this GID on both
        // Maps and add them to the permutation lists.
        ExportData_->permuteToLIDs_.push_back (target.getLocalElement (curSourceGID));
        ExportData_->permuteFromLIDs_.push_back (source.getLocalElement (curSourceGID));
      }
      else {
        // The current GID is owned by this process in the source Map,
        // but is not owned by this process in the target Map.  That
        // means the Export operation has to send it to another
        // process.  Store such GIDs in the "export" (outgoing) list.
        //
        // QUESTION (mfh 18 Aug 2012) Import at this point computes
        // remoteLIDs_.  Would it makes sense to compute them here,
        // instead of passing over exportGIDs_ again below?  That
        // would make the implementations of Export and Import look
        // more alike.
        ExportData_->exportGIDs_.push_back (curSourceGID);
      }
    }

    // Above, we filled exportGIDs_ with all the "outgoing" GIDs (that
    // is, the GIDs which we own in the source Map, but not in the
    // target Map).  Now allocate exportLIDs_, and fill it with the
    // LIDs (from the source Map) corresponding to those GIDs.
    //
    // exportLIDs_ is the list of this process' LIDs that it has to
    // send out.  Since this is an Export, and therefore the target
    // Map is nonoverlapping, we know that each export LID only needs
    // to be sent to one process.  However, the source Map may be
    // overlapping, so multiple processes might send to the same LID
    // on a receiving process.
    if (ExportData_->exportGIDs_.size()) {
      ExportData_->exportLIDs_ = arcp<LocalOrdinal>(ExportData_->exportGIDs_.size());
    }
    {
      typename ArrayRCP<LocalOrdinal>::iterator liditer = ExportData_->exportLIDs_.begin();
      typename Array<GlobalOrdinal>::iterator   giditer = ExportData_->exportGIDs_.begin();
      for (; giditer != ExportData_->exportGIDs_.end(); ++liditer, ++giditer) {
        *liditer = source.getLocalElement(*giditer);
      }
    }

    TPETRA_ABUSE_WARNING(
      getNumExportIDs() > 0 && ! source.isDistributed(),
      std::runtime_error,
      "::setupSamePermuteExport(): Source has export LIDs but Source is not "
      "distributed globally." << std::endl
      << "Importing to a submap of the target map.");

    // Compute exportImageIDs_ ("outgoing" process IDs).
    //
    // For each GID in exportGIDs_ (GIDs to which this process must
    // send), find its corresponding owning process (a.k.a. "image")
    // ID in the target Map.  Store these process IDs in
    // exportImageIDs_.  These are the process IDs to which the Export
    // needs to send data.
    //
    // We only need to do this if the source Map is distributed;
    // otherwise, the Export doesn't have to perform any
    // communication.
    if (source.isDistributed()) {
      ExportData_->exportImageIDs_ = arcp<int>(ExportData_->exportGIDs_.size());
      // This call will assign any GID in the target Map with no
      // corresponding process ID a fake process ID of -1.  We'll use
      // this below to remove exports for processses that don't exist.
      const LookupStatus lookup = target.getRemoteIndexList(ExportData_->exportGIDs_(), ExportData_->exportImageIDs_());
      TPETRA_ABUSE_WARNING( lookup == IDNotPresent, std::runtime_error,
        "::setupSamePermuteExport(): The source Map has GIDs not found in the target Map.");

      // Get rid of process IDs not in the target Map.  This prevents
      // exporting to GIDs which don't belong to any process in the
      // target Map.
      typedef typename ArrayRCP<int>::difference_type size_type;
      if (lookup == IDNotPresent) {
        const size_type numInvalidExports =
          std::count_if (ExportData_->exportImageIDs_().begin(),
                         ExportData_->exportImageIDs_().end(),
                         std::bind1st (std::equal_to<int>(), -1));

        // count number of valid and total number of exports
        const size_type totalNumExports = ExportData_->exportImageIDs_.size();
        if (numInvalidExports == totalNumExports) {
          // all exports are invalid; we have no exports; we can delete all exports
          ExportData_->exportGIDs_.resize(0);
          ExportData_->exportLIDs_ = null;
          ExportData_->exportImageIDs_ = null;
        }
        else {
          // some exports are valid; we need to keep the valid exports
          // pack and resize
          size_type numValidExports = 0;
          for (size_type e = 0; e < totalNumExports; ++e) {
            if (ExportData_->exportImageIDs_[e] != -1) {
              ExportData_->exportGIDs_[numValidExports]     = ExportData_->exportGIDs_[e];
              ExportData_->exportLIDs_[numValidExports]     = ExportData_->exportLIDs_[e];
              ExportData_->exportImageIDs_[numValidExports] = ExportData_->exportImageIDs_[e];
              ++numValidExports;
            }
          }
          ExportData_->exportGIDs_.resize(numValidExports);
          ExportData_->exportLIDs_.resize(numValidExports);
          ExportData_->exportImageIDs_.resize(numValidExports);
        }
      }
    }
  } // setupSamePermuteExport()


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Export<LocalOrdinal,GlobalOrdinal,Node>::setupRemote()
  {
    const Map<LocalOrdinal,GlobalOrdinal,Node>& target = *getTargetMap();

    // Sort exportImageIDs_ in ascending order, and apply the same
    // permutation to exportGIDs_ and exportLIDs_.  This ensures that
    // exportImageIDs_[i], exportGIDs_[i], and remoteLIDs_[i] all
    // refer to the same thing.
    sort3 (ExportData_->exportImageIDs_.begin(),
           ExportData_->exportImageIDs_.end(),
           ExportData_->exportGIDs_.begin(),
           ExportData_->exportLIDs_.begin());

    // Construct the list of entries that calling image needs to send
    // as a result of everyone asking for what it needs to receive.
    //
    // mfh 05 Jan 2012: I understand the above comment as follows:
    // Construct the communication plan from the list of image IDs to
    // which we need to send.
    size_t numRemoteIDs;
    numRemoteIDs = ExportData_->distributor_.createFromSends (ExportData_->exportImageIDs_());

    // Use the communication plan with ExportGIDs to find out who is
    // sending to us and get the proper ordering of GIDs for incoming
    // remote entries (these will be converted to LIDs when done).
    Array<GlobalOrdinal> remoteGIDs(numRemoteIDs);
    ExportData_->distributor_.doPostsAndWaits (ExportData_->exportGIDs_().getConst(), 1, remoteGIDs());

    // Remote (incoming) IDs come in as GIDs; convert to LIDs.  LIDs
    // tell this process where to store the incoming remote data.
    ExportData_->remoteLIDs_.resize(numRemoteIDs);
    {
      typename Array<GlobalOrdinal>::const_iterator i = remoteGIDs.begin();
      typename Array<LocalOrdinal>::iterator        j = ExportData_->remoteLIDs_.begin();
      while (i != remoteGIDs.end()) {
        *j++ = target.getLocalElement(*i++);
      }
    }
  }

  /** \brief Non-member constructor for Export objects.

      Creates a Export object from the given source and target maps.
      \pre <tt>src != null</tt>
      \pre <tt>tgt != null</tt>
      \return Returns the Export object. If <tt>src == tgt</tt>, returns \c null. (Debug mode: throws std::runtime_error if one of \c src or \c tgt is \c null.)

      \relatesalso Export
    */
  template <class LO, class GO, class Node>
  RCP< const Export<LO,GO,Node> >
  createExport( const RCP<const Map<LO,GO,Node> > & src,
                const RCP<const Map<LO,GO,Node> > & tgt )
  {
    if (src == tgt) return null;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(src == null || tgt == null, std::runtime_error,
        "Tpetra::createExport(): neither source nor target map may be null:\nsource: " << src << "\ntarget: " << tgt << "\n");
#endif
    return rcp(new Export<LO,GO,Node>(src,tgt));
  }

} // namespace Tpetra

#endif // TPETRA_EXPORT_HPP
