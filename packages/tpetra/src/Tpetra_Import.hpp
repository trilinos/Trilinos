// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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

  //! \brief This class builds an object containing information necesary for efficiently importing off-processor entries.
  /*! Import is used to construct a communication plan that can be called repeatedly by computational
      classes to efficiently import entries from other nodes.
      For example, an exporter is used when we start out with a multiple-ownership distribution,
      and we want to merge that into a uniquely-owned distribution.

      This class currently has one constructor, taking two Map objects
      specifying the distributions of the distributed objects on which the Export class will operate.

      This class is templated on \c LocalOrdinal and \c GlobalOrdinal. 
      The \c GlobalOrdinal type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Import: public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Import object from the source and target Maps.
    Import(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & source, 
           const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & target);

    //! copy constructor. 
    Import(const Import<LocalOrdinal,GlobalOrdinal,Node> & import);

    //! destructor.
    ~Import();

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    inline size_t getNumSameIDs() const;

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    inline size_t getNumPermuteIDs() const;

    //! List of entries in the source Map that are permuted. (non-persisting view)
    inline ArrayView<const LocalOrdinal> getPermuteFromLIDs() const;

    //! List of entries in the target Map that are permuted. (non-persisting view)
    inline ArrayView<const LocalOrdinal> getPermuteToLIDs() const;

    //! Returns the number of entries that are not on the calling image.
    inline size_t getNumRemoteIDs() const;

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    inline ArrayView<const LocalOrdinal> getRemoteLIDs() const;

    //! Returns the number of entries that must be sent by the calling image to other images.
    inline size_t getNumExportIDs() const;

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    inline ArrayView<const LocalOrdinal> getExportLIDs() const;

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    inline ArrayView<const int> getExportImageIDs() const;

    //! Returns the Source Map used to construct this importer.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getSourceMap() const;

    //! Returns the Target Map used to construct this importer.
    inline const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getTargetMap() const;

    inline Distributor & getDistributor() const;

    //! Assignment operator
    Import<LocalOrdinal,GlobalOrdinal,Node>& operator = (const Import<LocalOrdinal,GlobalOrdinal,Node> & Source);

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method
    virtual void print(std::ostream& os) const;

    //@}

  private:

    RCP<ImportExportData<LocalOrdinal,GlobalOrdinal,Node> > ImportData_;
    RCP<Array<GlobalOrdinal> > remoteGIDs_;

    // subfunctions used by constructor
    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and numRemoteIDs_
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, and remoteLIDs_
    void setupSamePermuteRemote();
    void setupExport();
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
    using std::endl;
    ArrayView<const LocalOrdinal> av;
    ArrayView<const int> avi;
    const RCP<const Comm<int> > & comm = getSourceMap()->getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) 
      {
        os << endl;
        if(myImageID == 0) { // this is the root node (only output this info once)
          os << "Import Data Members:" << endl;
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
      os << "\nSource Map: " << endl; 
    }
    os << *getSourceMap();
    if (myImageID == 0) {
      os << "\nTarget Map: " << endl; 
    }
    os << *getTargetMap();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::setupSamePermuteRemote() {
    const Map<LocalOrdinal,GlobalOrdinal,Node> & source = *getSourceMap();
    const Map<LocalOrdinal,GlobalOrdinal,Node> & target = *getTargetMap();
    ArrayView<const GlobalOrdinal> sourceGIDs = source.getNodeElementList();
    ArrayView<const GlobalOrdinal> targetGIDs = target.getNodeElementList();

    // -- compute numSameIDs_ ---
    // go through GID lists of source and target. if the ith GID on both is the same, 
    // increment numSameIDs_ and try the next. as soon as you come to a pair that don't
    // match, give up.
    typename ArrayView<const GlobalOrdinal>::iterator sourceIter = sourceGIDs.begin(),
                                                      targetIter = targetGIDs.begin();
    while( sourceIter != sourceGIDs.end() && targetIter != targetGIDs.end() && *sourceIter == *targetIter )
    {
      ++ImportData_->numSameIDs_;
      ++sourceIter;
      ++targetIter;
    }
    // targetIter should now point to the GID of the first non-same entry or the end of targetGIDs

    // -- compute numPermuteIDs and numRemoteIDs --
    // -- fill permuteToLIDs_, permuteFromLIDs_, remoteGIDs_, and remoteLIDs_ --
    // go through remaining entries in targetGIDs. if source owns that GID, 
    // increment numPermuteIDs_, and add entries to permuteToLIDs_ and permuteFromLIDs_.
    // otherwise increment numRemoteIDs_ and add entries to remoteLIDs_ and remoteGIDs_.
    remoteGIDs_ = rcp( new Array<GlobalOrdinal>() );
    for (; targetIter != targetGIDs.end(); ++targetIter) {
      if (source.isNodeGlobalElement(*targetIter)) {
        // both source and target list this GID (*targetIter)
        // determine the LIDs for this GID on both Maps and add them to the permutation lists
        ImportData_->permuteToLIDs_.push_back(target.getLocalElement(*targetIter));
        ImportData_->permuteFromLIDs_.push_back(source.getLocalElement(*targetIter));
      }
      else {
        // this GID is on another processor; store it, along with its destination LID on this processor
        remoteGIDs_->push_back(*targetIter);
        ImportData_->remoteLIDs_.push_back(target.getLocalElement(*targetIter));
      }
    }

    TPETRA_ABUSE_WARNING((getNumRemoteIDs() > 0) && !source.isDistributed(),std::runtime_error,
        "::setupSamePermuteRemote(): Target has remote LIDs but Source is not distributed globally."
        << std::endl << "Importing to a submap of the target map.");
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::setupExport() {
    typedef typename Array<int>::difference_type size_type;
    const Map<LocalOrdinal,GlobalOrdinal,Node> & source = *getSourceMap();

    // create remoteImageID list: for each entry remoteGIDs[i],
    // remoteImageIDs[i] will contain the ImageID of the image that owns that GID.
    // check for GIDs that exist in target but not in source: we see this if getRemoteIndexList returns true
    ArrayView<GlobalOrdinal> remoteGIDs = (*remoteGIDs_)();
    Array<int> remoteImageIDs(remoteGIDs.size());
    const LookupStatus lookup = source.getRemoteIndexList(remoteGIDs, remoteImageIDs());
    TPETRA_ABUSE_WARNING( lookup == IDNotPresent, std::runtime_error, "::setupExport(): Target has GIDs not found in Source." );

    // get rid of ids that don't exist in the source map
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
        // some remotes are valid; we need to keep the valid ones
        // pack and resize
        size_type numValidRemote = 0;
        for (size_type r = 0; r < totalNumRemote; ++r) {
          if (remoteImageIDs[r] != -1) {
            remoteImageIDs[numValidRemote] = remoteImageIDs[r];
            (*remoteGIDs_)[numValidRemote] = (*remoteGIDs_)[r];
            ImportData_->remoteLIDs_[numValidRemote] = ImportData_->remoteLIDs_[r];
            ++numValidRemote;
          }
        }
        TEST_FOR_EXCEPTION( numValidRemote != totalNumRemote - numInvalidRemote, std::logic_error,
            typeName(*this) << "::setupExport(): internal logic error. Please contact Tpetra team.")
        remoteImageIDs.resize(numValidRemote);
        (*remoteGIDs_).resize(numValidRemote);
        ImportData_->remoteLIDs_.resize(numValidRemote);
      }
      remoteGIDs = (*remoteGIDs_)();
    }

    // sort remoteImageIDs in ascending order
    // apply same permutation to remoteGIDs_
    sort2(remoteImageIDs.begin(), remoteImageIDs.end(), remoteGIDs.begin());

    // call Distributor.createFromRecvs()
    // takes in remoteGIDs and remoteImageIDs_
    // returns exportLIDs_, exportImageIDs_ 
    ArrayRCP<GlobalOrdinal> exportGIDs;
    ImportData_->distributor_.createFromRecvs(remoteGIDs().getConst(), remoteImageIDs, exportGIDs, ImportData_->exportImageIDs_);
    // -- exportGIDs and exportImageIDs_ allocated by createFromRecvs (the former contains GIDs, we will convert to LIDs below) --

    // convert exportGIDs from GIDs to LIDs
    if (exportGIDs != null) {
      ImportData_->exportLIDs_ = arcp<LocalOrdinal>(exportGIDs.size());
    }
    typename ArrayRCP<LocalOrdinal>::iterator dst = ImportData_->exportLIDs_.begin();
    typename ArrayRCP<GlobalOrdinal>::const_iterator src = exportGIDs.begin();
    while (src != exportGIDs.end())
    {
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
    TEST_FOR_EXCEPTION(src == null || tgt == null, std::runtime_error,
        "Tpetra::createImport(): neither source nor target map may be null:\nsource: " << src << "\ntarget: " << tgt << "\n");
#endif
    return rcp(new Import<LO,GO,Node>(src,tgt));
  }

  /** \brief Deprecated. Use createImport().
    */
  template <class LO, class GO, class Node> 
  RCP< const Import<LO,GO,Node> >
  TPETRA_DEPRECATED makeImport( const RCP<const Map<LO,GO,Node> > & src, 
                                const RCP<const Map<LO,GO,Node> > & tgt )
  {
    return createImport<LO,GO,Node>(src,tgt);
  }

} // namespace Tpetra

#endif // TPETRA_IMPORT_HPP
