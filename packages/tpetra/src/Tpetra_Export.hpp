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

#ifndef TPETRA_EXPORT_HPP
#define TPETRA_EXPORT_HPP

#include <Teuchos_Object.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ImportExportData.hpp"

namespace Tpetra {

  //! Tpetra::Export: This class builds an export object for efficiently exporting entries off-processor.

  /*! Export is used to construct a communication plan that can be called repeatedly by computational
      classes such the Tpetra CisMatrix and Vector classes to efficiently export entries to other
      images. An exporter is used when we start out with a multiple-ownership distribution,
      and we want to merge that into a uniquely-owned distribution.

      This class currently has one constructor, taking two Map objects.
      The first Map specifies the distribution we have now. The second 
      Map specifies the distribution we want to have after exporting.

      NOTE: Behavior is undefined if the destination Map is not uniquely-owned.
   */

  template <typename Ordinal>
  class Export: public Teuchos::Object {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructs a Export object from the source and target Map.
    Export(const Map<Ordinal> & source, const Map<Ordinal> & target);

    //! copy constructor. 
    Export(Export<Ordinal> const& rhs);

    //! destructor.
    ~Export();

    //@}

    //! @name Export Attribute Methods
    //@{ 

    //! Returns the number of entries that are identical between the source and target Maps, up to the first different ID.
    Teuchos_Ordinal getNumSameIDs() const;

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    Teuchos_Ordinal getNumPermuteIDs() const;

    //! List of entries in the source Map that are permuted. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getPermuteFromLIDs() const;

    //! List of entries in the target Map that are permuted. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getPermuteToLIDs() const;

    //! Returns the number of entries that are not on the calling image.
    Teuchos_Ordinal getNumRemoteIDs() const;

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getRemoteLIDs() const;

    //! Returns the number of entries that must be sent by the calling image to other images.
    Teuchos_Ordinal getNumExportIDs() const;

    //! List of entries in the source Map that will be sent to other images. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getExportLIDs() const;

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i]. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getExportImageIDs() const;

    //! Returns the Source Map used to construct this exporter.
    const Map<Ordinal> & getSourceMap() const;

    //! Returns the Target Map used to construct this exporter.
    const Map<Ordinal> & getTargetMap() const;

    Distributor<Ordinal>& getDistributor() const;

    //! Assignment operator
    Export<Ordinal>& operator = (const Export<Ordinal> & Source);

    //@}

    //! @name I/O Methods
    //@{ 

    //! Print method inherited from Teuchos::Object
    virtual void print(std::ostream& os) const;

    //@}

  private:

    Teuchos::RCP<ImportExportData<Ordinal> > ExportData_;

    // subfunctions used by constructor
    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and numExportIDs_
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, exportGIDs_, and exportLIDs_
    void setupSamePermuteExport();
    void setupRemote();
  };

  template <typename Ordinal>
    Export<Ordinal>::Export(const Map<Ordinal> & source, const Map<Ordinal> & target)
    : Teuchos::Object("Tpetra::Export")
    , ExportData_()
  {
    ExportData_ = Teuchos::rcp(new ImportExportData<Ordinal>(source, target));
    // call subfunctions
    setupSamePermuteExport();
    if(source.isDistributed()) {
      setupRemote();
    }
  }

  template <typename Ordinal>
  Export<Ordinal>::Export(const Export<Ordinal> & rhs)
  : Teuchos::Object(rhs.label())
  , ExportData_(rhs.ExportData_)
  {}

  template <typename Ordinal>
  Export<Ordinal>::~Export() 
  {}

  template <typename Ordinal>
  Teuchos_Ordinal Export<Ordinal>::getNumSameIDs() const {
    return ExportData_->numSameIDs_;
  }

  template <typename Ordinal>
  Teuchos_Ordinal Export<Ordinal>::getNumPermuteIDs() const {
    return ExportData_->permuteFromLIDs_.size();
  }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> 
  Export<Ordinal>::getPermuteFromLIDs() const {
    return ExportData_->permuteFromLIDs_();
  }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal>
  Export<Ordinal>::getPermuteToLIDs() const {
    return ExportData_->permuteToLIDs_();
  }

  template <typename Ordinal>
  Teuchos_Ordinal Export<Ordinal>::getNumRemoteIDs() const {
    return ExportData_->remoteLIDs_.size();
  }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> 
  Export<Ordinal>::getRemoteLIDs() const {
    return ExportData_->remoteLIDs_();
  }

  template <typename Ordinal>
  Teuchos_Ordinal Export<Ordinal>::getNumExportIDs() const {
    return ExportData_->exportLIDs_.size();
  }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> 
  Export<Ordinal>::getExportLIDs() const {
    return ExportData_->exportLIDs_();
  }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> 
  Export<Ordinal>::getExportImageIDs() const {
    return ExportData_->exportImageIDs_();
  }

  template <typename Ordinal>
  const Map<Ordinal> & 
  Export<Ordinal>::getSourceMap() const {
    return ExportData_->source_;
  }

  template <typename Ordinal>
  const Map<Ordinal> & 
  Export<Ordinal>::getTargetMap() const {
    return ExportData_->target_;
  }

  template <typename Ordinal>
  Distributor<Ordinal>& 
  Export<Ordinal>::getDistributor() const {
    return ExportData_->distributor_;
  }

  template <typename Ordinal>
  Export<Ordinal>& 
  Export<Ordinal>::operator=(const Export<Ordinal> & Source) 
  {
    ExportData_ = Source.ExportData_;
    return *this;
  }

  template <typename Ordinal>
  void Export<Ordinal>::print(std::ostream& os) const 
  {
    using std::endl;
    Teuchos::ArrayView<const Ordinal> av;
    int myImageID = getSourceMap().getComm()->getRank();
    int numImages = getSourceMap().getComm()->getSize();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if(myImageID == imageCtr) {
        os << endl;
        if(myImageID == 0) { // this is the root node (only output this info once)
          os << "Export Data Members:" << endl;
        }
        os << "Image ID       : " << myImageID << endl;
        os << "permuteFromLIDs: {"; av = getPermuteFromLIDs(); std::copy(av.begin(),av.end(),std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
        os << "permuteToLIDs  : {"; av = getPermuteToLIDs();   std::copy(av.begin(),av.end(),std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
        os << "remoteLIDs     : {"; av = getRemoteLIDs();      std::copy(av.begin(),av.end(),std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
        os << "exportLIDs     : {"; av = getExportLIDs();      std::copy(av.begin(),av.end(),std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
        os << "exportImageIDs : {"; av = getExportImageIDs();  std::copy(av.begin(),av.end(),std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
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


  template <typename Ordinal>
  void Export<Ordinal>::setupSamePermuteExport() 
  {
    const Map<Ordinal> & source = getSourceMap();
    const Map<Ordinal> & target = getTargetMap();
    Teuchos::ArrayView<const Ordinal> sourceGIDs = source.getMyGlobalEntries();
    Teuchos::ArrayView<const Ordinal> targetGIDs = target.getMyGlobalEntries();
    const Ordinal ZERO = Teuchos::OrdinalTraits<Ordinal>::zero();

    // -- compute numSameIDs_ ---
    // go through GID lists of source and target. if the ith GID on both is the same, 
    // increment numSameIDs_ and try the next. as soon as you come to a pair that don't
    // match, give up.
    typename Teuchos::ArrayView<const Ordinal>::iterator sourceIter = sourceGIDs.begin(),
                                                         targetIter = targetGIDs.begin();
    while( sourceIter != sourceGIDs.end() && targetIter != targetGIDs.end() && *sourceIter == *targetIter )
    {
      ++ExportData_->numSameIDs_;
      ++sourceIter;
      ++targetIter;
    }
    // sourceIter should now point to the GID of the first non-same entry or at the end of targetGIDs

    // -- compute numPermuteIDs --
    // -- fill permuteToLIDs_, permuteFromLIDs_ --
    // go through remaining entries in sourceGIDs. if target owns that GID, 
    // increment numPermuteIDs_, and add entries to permuteToLIDs_ and permuteFromLIDs_.
    // otherwise increment numExportIDs_ and add entries to exportLIDs_ and exportGIDs_.
    for(; sourceIter != sourceGIDs.end(); ++sourceIter) {
      if(target.isMyGlobalIndex(*sourceIter)) {
        // both source and target list this GID (*targetIter)
        // determine the LIDs for this GID on both Maps and add them to the permutation lists
        ExportData_->permuteToLIDs_.push_back(  target.getLocalIndex(*sourceIter));
        ExportData_->permuteFromLIDs_.push_back(source.getLocalIndex(*sourceIter));
      }
      else {
        ExportData_->exportGIDs_.push_back(*sourceIter);
      }
    }

    TEST_FOR_EXCEPTION( (getNumExportIDs() > ZERO) && (!source.isDistributed()), std::runtime_error, 
        "Tpetra::Export<" << Teuchos::OrdinalTraits<Ordinal>::name() 
        << ">::setupSamePermuteExport(): Source has export LIDs but Source is not distributed globally.");

    // -- compute exportImageIDs_ --
    // get list of images that own the GIDs in exportGIDs_ (in the target Map)
    ExportData_->exportImageIDs_ = Teuchos::arcp<Ordinal>(ExportData_->exportGIDs_.size());
    ExportData_->exportLIDs_     = Teuchos::arcp<Ordinal>(ExportData_->exportGIDs_.size());
    {
      typename Teuchos::ArrayRCP<Ordinal>::iterator liditer = ExportData_->exportLIDs_.begin();
      typename Teuchos::Array<Ordinal>::iterator    giditer = ExportData_->exportGIDs_.begin();
      for (; giditer != ExportData_->exportGIDs_.end(); ++liditer, ++giditer) {
        *liditer = source.getLocalIndex(*giditer);
      }
    }
    TEST_FOR_EXCEPTION( target.getRemoteIndexList(ExportData_->exportGIDs_(), ExportData_->exportImageIDs_()) == true,
        std::runtime_error, "Tpetra::Export::setupSamePermuteExport(): Source has GIDs not found in Target.");
  }


  template <typename Ordinal>
  void Export<Ordinal>::setupRemote() 
  {
    const Map<Ordinal> & target = getTargetMap();
    const Ordinal ONE  = Teuchos::OrdinalTraits<Ordinal>::one();

    // make sure export IDs are ordered by image
    // sort exportImageIDs_ in ascending order,
    // and apply the same permutation to exportGIDs_ and exportLIDs_.
    sortArrays(ExportData_->exportImageIDs_(), ExportData_->exportGIDs_(), ExportData_->exportLIDs_());

    // Construct list of entries that calling image needs to send as a result
    // of everyone asking for what it needs to receive.
    Teuchos_Ordinal numRemoteIDs;
    ExportData_->distributor_.createFromSends(ExportData_->exportImageIDs_(), numRemoteIDs);

    // Use comm plan with ExportGIDs to find out who is sending to us and
    // get proper ordering of GIDs for remote entries 
    // (that we will convert to LIDs when done).
    Teuchos::RCP< Teuchos::Comm<Ordinal> > comm = ExportData_->platform_->createComm();
    Teuchos::Array<Ordinal> remoteGIDs(numRemoteIDs);
    ExportData_->distributor_.doPostsAndWaits(ExportData_->exportGIDs_().getConst(),ONE,remoteGIDs());

    // Remote IDs come in as GIDs, convert to LIDs
    ExportData_->remoteLIDs_.resize(numRemoteIDs);
    {
      typename Teuchos::Array<Ordinal>::const_iterator i = remoteGIDs.begin();
      typename Teuchos::Array<Ordinal>::iterator       j = ExportData_->remoteLIDs_.begin();
      while (i != remoteGIDs.end()) 
      {
        *j++ = target.getLocalIndex(*i++);
      }
    }
  }

} // namespace Tpetra

#endif // TPETRA_EXPORT_HPP
