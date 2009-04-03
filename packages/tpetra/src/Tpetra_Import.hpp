// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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

#include <Teuchos_Describable.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ImportExportData.hpp"
#include <iterator>

namespace Tpetra {

  //! Tpetra::Import: This class builds an import object for efficiently importing off-processor entries.

  /*! Import is used to construct a communication plan that can be called repeatedly by computational
      classes such the Tpetra CisMatrix and Vector classes to efficiently import entries from other
      images. An importer is used when we start out with a uniquely-owned distribution,
      and want to distribute that into a multiple-ownership distribution.

      This class currently has one constructor, taking two Map objects.
      The first Map specifies the distribution we have now. The second 
      Map specifies the distribution we want to have after importing.

      NOTE: Behavior is undefined if the source Map is not uniquely-owned.
  */

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal>
  class Import: public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Import object from the source and target Maps.
    Import(const Map<LocalOrdinal,GlobalOrdinal> & source, const Map<LocalOrdinal,GlobalOrdinal> & target);

    //! copy constructor. 
    Import(const Import<LocalOrdinal,GlobalOrdinal> & import);

    //! destructor.
    ~Import();

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    Teuchos_Ordinal getNumSameIDs() const;

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    Teuchos_Ordinal getNumPermuteIDs() const;

    //! List of entries in the source Map that are permuted. (non-persisting view)
    Teuchos::ArrayView<const LocalOrdinal> getPermuteFromLIDs() const;

    //! List of entries in the target Map that are permuted. (non-persisting view)
    Teuchos::ArrayView<const LocalOrdinal> getPermuteToLIDs() const;

    //! Returns the number of entries that are not on the calling image.
    Teuchos_Ordinal getNumRemoteIDs() const;

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    Teuchos::ArrayView<const LocalOrdinal> getRemoteLIDs() const;

    //! Returns the number of entries that must be sent by the calling image to other images.
    Teuchos_Ordinal getNumExportIDs() const;

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    Teuchos::ArrayView<const LocalOrdinal> getExportLIDs() const;

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    Teuchos::ArrayView<const int> getExportImageIDs() const;

    //! Returns the Source Map used to construct this importer.
    const Map<LocalOrdinal,GlobalOrdinal> & getSourceMap() const;

    //! Returns the Target Map used to construct this importer.
    const Map<LocalOrdinal,GlobalOrdinal> & getTargetMap() const;

    Distributor & getDistributor() const;

    //! Assignment operator
    Import<LocalOrdinal,GlobalOrdinal>& operator = (const Import<LocalOrdinal,GlobalOrdinal> & Source);

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method 
    virtual void print(std::ostream& os) const;

    //@}

  private:

    Teuchos::RCP<ImportExportData<LocalOrdinal,GlobalOrdinal> > ImportData_;
    Teuchos::RCP<Teuchos::Array<GlobalOrdinal> > remoteGIDs_;

    // subfunctions used by constructor
    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and numRemoteIDs_
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, and remoteLIDs_
    void setupSamePermuteRemote();
    void setupExport();
  };

  template <class LocalOrdinal, class GlobalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal>::Import(const Map<LocalOrdinal,GlobalOrdinal> & source, const Map<LocalOrdinal,GlobalOrdinal> & target)
  : ImportData_()
  {
    ImportData_ = Teuchos::rcp(new ImportExportData<LocalOrdinal,GlobalOrdinal>(source, target));
    // call subfunctions
    setupSamePermuteRemote();
    if(source.isDistributed()) {
      setupExport();
    }
    // don't need remoteGIDs_ anymore
    remoteGIDs_ = Teuchos::null;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal>::Import(const Import<LocalOrdinal,GlobalOrdinal> & import)
  : ImportData_(import.ImportData_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal>::~Import() 
  {}

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal Import<LocalOrdinal,GlobalOrdinal>::getNumSameIDs() const {
    return ImportData_->numSameIDs_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal Import<LocalOrdinal,GlobalOrdinal>::getNumPermuteIDs() const {
    return ImportData_->permuteFromLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::ArrayView<const LocalOrdinal> 
  Import<LocalOrdinal,GlobalOrdinal>::getPermuteFromLIDs() const {
    return ImportData_->permuteFromLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal>::getPermuteToLIDs() const {
    return ImportData_->permuteToLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal Import<LocalOrdinal,GlobalOrdinal>::getNumRemoteIDs() const {
    return ImportData_->remoteLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::ArrayView<const LocalOrdinal> 
  Import<LocalOrdinal,GlobalOrdinal>::getRemoteLIDs() const {
    return ImportData_->remoteLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos_Ordinal Import<LocalOrdinal,GlobalOrdinal>::getNumExportIDs() const {
    return ImportData_->exportLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::ArrayView<const LocalOrdinal> 
  Import<LocalOrdinal,GlobalOrdinal>::getExportLIDs() const {
    return ImportData_->exportLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Teuchos::ArrayView<const int> 
  Import<LocalOrdinal,GlobalOrdinal>::getExportImageIDs() const {
    return ImportData_->exportImageIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  Import<LocalOrdinal,GlobalOrdinal>::getSourceMap() const {
    return ImportData_->source_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & 
  Import<LocalOrdinal,GlobalOrdinal>::getTargetMap() const {
    return ImportData_->target_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Distributor & 
  Import<LocalOrdinal,GlobalOrdinal>::getDistributor() const {
    return ImportData_->distributor_;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal>& 
  Import<LocalOrdinal,GlobalOrdinal>::operator=(const Import<LocalOrdinal,GlobalOrdinal> & Source) 
  {
    ImportData_ = Source.ImportData_;
    return *this;
  }

  template <class LocalOrdinal, class GlobalOrdinal>
  void Import<LocalOrdinal,GlobalOrdinal>::print(std::ostream& os) const 
  {
    using std::endl;
    Teuchos::ArrayView<const LocalOrdinal> av;
    Teuchos::ArrayView<const int> avi;
    int myImageID = getSourceMap().getComm()->getRank();
    int numImages = getSourceMap().getComm()->getSize();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if(myImageID == imageCtr) {
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
      getSourceMap().getComm()->barrier();
      getSourceMap().getComm()->barrier();
      getSourceMap().getComm()->barrier();
    }
    if (myImageID == 0) {
      os << "\nSource Map: " << endl; 
    }
    os << getSourceMap();
    if (myImageID == 0) {
      os << "\nTarget Map: " << endl; 
    }
    os << getTargetMap();
  }


  template <class LocalOrdinal, class GlobalOrdinal>
  void Import<LocalOrdinal,GlobalOrdinal>::setupSamePermuteRemote() 
  {
    const Map<LocalOrdinal,GlobalOrdinal> & source = getSourceMap();
    const Map<LocalOrdinal,GlobalOrdinal> & target = getTargetMap();
    Teuchos::ArrayView<const GlobalOrdinal> sourceGIDs = source.getMyGlobalEntries();
    Teuchos::ArrayView<const GlobalOrdinal> targetGIDs = target.getMyGlobalEntries();

    // -- compute numSameIDs_ ---
    // go through GID lists of source and target. if the ith GID on both is the same, 
    // increment numSameIDs_ and try the next. as soon as you come to a pair that don't
    // match, give up.
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator sourceIter = sourceGIDs.begin(),
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
    remoteGIDs_ = Teuchos::rcp( new Teuchos::Array<GlobalOrdinal>() );
    for (; targetIter != targetGIDs.end(); ++targetIter) {
      if (source.isMyGlobalIndex(*targetIter)) {
        // both source and target list this GID (*targetIter)
        // determine the LIDs for this GID on both Maps and add them to the permutation lists
        ImportData_->permuteToLIDs_.push_back(target.getLocalIndex(*targetIter));
        ImportData_->permuteFromLIDs_.push_back(source.getLocalIndex(*targetIter));
      }
      else {
        // this GID is on another processor; store it, along with its destination LID on this processor
        remoteGIDs_->push_back(*targetIter);
        ImportData_->remoteLIDs_.push_back(target.getLocalIndex(*targetIter));
      }
    }

    TEST_FOR_EXCEPTION( (getNumRemoteIDs() > 0) && !source.isDistributed(), std::runtime_error, 
        Teuchos::typeName(*this) << "::setupSamePermuteRemote(): Target has remote LIDs but Source is not distributed globally.");
  }


  template <class LocalOrdinal, class GlobalOrdinal>
  void Import<LocalOrdinal,GlobalOrdinal>::setupExport()
  {
    const Map<LocalOrdinal,GlobalOrdinal> & source = getSourceMap();

    // create remoteImageID list: for each entry remoteGIDs[i],
    // remoteImageIDs[i] will contain the ImageID of the image that owns that GID.
    // check for GIDs that exist in target but not in source: we see this if getRemoteIndexList returns true
    Teuchos::ArrayView<GlobalOrdinal> remoteGIDs = (*remoteGIDs_)();
    Teuchos::Array<int> remoteImageIDs(remoteGIDs.size());
    TEST_FOR_EXCEPTION( source.getRemoteIndexList(remoteGIDs, remoteImageIDs()) == true, std::runtime_error,
        "Tpetra::Import::setupExport(): Target has GIDs not found in Source.");

    // sort remoteImageIDs in ascending order
    // apply same permutation to remoteGIDs_
    sort2(remoteImageIDs.begin(), remoteImageIDs.end(), remoteGIDs.begin());

    // call Distributor.createFromRecvs()
    // takes in remoteGIDs and remoteImageIDs_
    // returns exportLIDs_, exportImageIDs_ 
    Teuchos::ArrayRCP<GlobalOrdinal> exportGIDs;
    ImportData_->distributor_.createFromRecvs(remoteGIDs().getConst(), remoteImageIDs, exportGIDs, ImportData_->exportImageIDs_);
    // -- exportGIDs and exportImageIDs_ allocated by createFromRecvs (the former contains GIDs, we will convert to LIDs below) --

    // convert exportGIDs from GIDs to LIDs
    if (exportGIDs != Teuchos::null) {
      ImportData_->exportLIDs_ = Teuchos::arcp<LocalOrdinal>(exportGIDs.size());
    }
    typename Teuchos::ArrayRCP<LocalOrdinal>::const_iterator dst = ImportData_->exportLIDs_.begin();
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator src = exportGIDs.begin();
    while (src != exportGIDs.end())
    {
      (*dst++) = source.getLocalIndex(*src++);
    }
  }

} // namespace Tpetra

#endif // TPETRA_IMPORT_HPP
