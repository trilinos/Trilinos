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

  //! Tpetra::Export: This class builds an export object for efficient exporting of off-processor entries.

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
    Ordinal getNumSameIDs() const;

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    Ordinal getNumPermuteIDs() const;

    //! List of entries in the source Map that are permuted. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getPermuteFromLIDs() const;

    //! List of entries in the target Map that are permuted. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getPermuteToLIDs() const;

    //! Returns the number of entries that are not on the calling image.
    Ordinal getNumRemoteIDs() const;

    //! List of entries in the target Map that are coming from other images. (non-persisting view)
    Teuchos::ArrayView<const Ordinal> getRemoteLIDs() const;

    //! Returns the number of entries that must be sent by the calling image to other images.
    Ordinal getNumExportIDs() const;

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
  Export<Ordinal>::Export(Export<Ordinal> const& rhs)
  : Teuchos::Object(rhs.label())
  , ExportData_(rhs.ExportData_)
  {}

  template <typename Ordinal>
  Export<Ordinal>::~Export() 
  {}

  template <typename Ordinal>
  Ordinal Export<Ordinal>::getNumSameIDs() const {
    return ExportData_->numSameIDs_;
  }

  template <typename Ordinal>
  Ordinal Export<Ordinal>::getNumPermuteIDs() const {
    return ExportData_->numPermuteIDs_;
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
  Ordinal Export<Ordinal>::getNumRemoteIDs() const {
    return ExportData_->numRemoteIDs_;
  }

  template <typename Ordinal>
  Teuchos::ArrayView<const Ordinal> 
  Export<Ordinal>::getRemoteLIDs() const {
    return ExportData_->remoteLIDs_();
  }

  template <typename Ordinal>
  Ordinal Export<Ordinal>::getNumExportIDs() const {
    return ExportData_->numExportIDs_;
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
    os << "Export Data Members:" << endl;
    os << "permuteToLIDs:   {"; std::copy(getPermuteToLIDs().begin()  ,getPermuteToLIDs().end()  ,std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
    os << "permuteFromLIDs: {"; std::copy(getPermuteFromLIDs().begin(),getPermuteFromLIDs().end(),std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
    os << "remoteLIDs:      {"; std::copy(getRemoteLIDs().begin()     ,getRemoteLIDs().end()     ,std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
    os << "exportLIDs:      {"; std::copy(getExportLIDs().begin()     ,getExportLIDs().end()     ,std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
    os << "exportImageIDs:  {"; std::copy(getExportImageIDs().begin() ,getExportImageIDs().end() ,std::ostream_iterator<Ordinal>(os," ")); os << " }" << endl;
    os << "numSameIDs   : " << getNumSameIDs() << endl;
    os << "numPermuteIDs: " << getNumPermuteIDs() << endl;
    os << "numRemoteIDs : " << getNumRemoteIDs() << endl;
    os << "numExportIDs : " << getNumExportIDs() << endl;
    os << "\nsource: " << endl << getSourceMap();
    os << "\ntarget: " << endl << getTargetMap();
  }


  template <typename Ordinal>
  void Export<Ordinal>::setupSamePermuteExport() 
  {
    Ordinal const zero = Teuchos::OrdinalTraits<Ordinal>::zero();
    Ordinal const one = Teuchos::OrdinalTraits<Ordinal>::one();
    Ordinal const negOne = zero - one;
    Map<Ordinal> const& source = getSourceMap();
    Map<Ordinal> const& target = getTargetMap();
    std::vector<Ordinal> const& sourceGIDs = source.getMyGlobalMap();
    std::vector<Ordinal> const& targetGIDs = target.getMyGlobalMap();

    // -- compute numSameIDs_ ---
    // go through GID lists of source and target. if the ith GID on both is the same, 
    // increment numSameIDs_ and try the next. as soon as you come to a pair that don't
    // match, give up.
    typename std::vector<Ordinal>::const_iterator sourceIter = sourceGIDs.begin();
    typename std::vector<Ordinal>::const_iterator targetIter = targetGIDs.begin();
    while((sourceIter != sourceGIDs.end()) && 
        (targetIter != targetGIDs.end()) && 
        (*targetIter == *sourceIter)) {
      ExportData_->numSameIDs_++;
      sourceIter++;
      targetIter++;
    }
    // sourceIter should now point to the GID of the first non-same entry

    // -- compute numPermuteIDs and numRemoteIDs --
    // -- fill permuteToLIDs_, permuteFromLIDs_, remoteGIDs_, and remoteLIDs_ --
    // go through remaining entries in sourceGIDs. if target owns that GID, 
    // increment numPermuteIDs_, and add entries to permuteToLIDs_ and permuteFromLIDs_.
    // otherwise increment numExportIDs_ and add entries to exportLIDs_ and exportGIDs_.
    for(; sourceIter != sourceGIDs.end(); sourceIter++) {
      if(target.isMyGID(*sourceIter)) {
        ExportData_->numPermuteIDs_++;
        ExportData_->permuteToLIDs_.push_back(target.getLID(*sourceIter));
        ExportData_->permuteFromLIDs_.push_back(source.getLID(*sourceIter));
      }
      else {
        ExportData_->numExportIDs_++;
        ExportData_->exportLIDs_.push_back(source.getLID(*sourceIter));
        ExportData_->exportGIDs_.push_back(*sourceIter);
      }
    }

    if((ExportData_->numExportIDs_ > zero) && (!source.isGlobal())) {
      throw reportError("Source has export LIDs but is not distributed globally.", 1); 
      //*** what do we do here??? ***
    }

    // -- compute exportImageIDs_ --
    // get list of images that own the GIDs in exportGIDs_ (in the target Map)
    // check exportImageIDs_ for any -1 entries (nobody owns that GID in the target Map)
    target.getRemoteIDList(ExportData_->exportGIDs_, ExportData_->exportImageIDs_);
    Ordinal count = std::count(ExportData_->exportImageIDs_.begin(), ExportData_->exportImageIDs_.end(), negOne);
    if(count > zero) {
      throw reportError("Source has GIDs not found in Target.", 2);
    }
  }

  template <typename Ordinal>
  void Export<Ordinal>::setupRemote() 
  {
    Ordinal const zero = Teuchos::OrdinalTraits<Ordinal>::zero();
    Ordinal const one = Teuchos::OrdinalTraits<Ordinal>::one();
    Map<Ordinal> const& target = getTargetMap();

    // make sure export IDs are ordered by image
    // sort exportImageIDs_ in ascending order,
    // and apply the same permutation to exportGIDs_ and exportLIDs_.
    sortArrays(ExportData_->exportImageIDs_, ExportData_->exportGIDs_, ExportData_->exportLIDs_);

    // Construct list of entries that calling image needs to send as a result
    // of everyone asking for what it needs to receive.
    ExportData_->distributor_.createFromSends(ExportData_->numExportIDs_, ExportData_->exportImageIDs_, true, ExportData_->numRemoteIDs_);
    // -- numRemoteIDs_ is now defined --

    // Use comm plan with ExportGIDs to find out who is sending to us and
    // get proper ordering of GIDs for remote entries 
    // (that we will convert to LIDs when done).
    Teuchos::RCP< Teuchos::Comm<Ordinal> > comm = ExportData_->platform_->createOrdinalComm();
    comm->doPostsAndWaits(ExportData_->distributor_, ExportData_->exportGIDs_, one, ExportData_->remoteGIDs_);
    // -- remoteGIDs_ is now defined --

    // Remote IDs come in as GIDs, convert to LIDs
    for(Ordinal i = zero; i < ExportData_->numRemoteIDs_; i++) {
      ExportData_->remoteLIDs_.push_back(target.getLID(ExportData_->remoteGIDs_[i]));
    }
  }

} // namespace Tpetra

#endif // TPETRA_EXPORT_HPP
