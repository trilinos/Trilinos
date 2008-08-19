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

#include <Teuchos_Object.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of ImportData, needed to prevent circular inclusions
  // actual #include statement at the end of this file
  template<typename Ordinal> class ImportData;
#endif

  //! Tpetra::Import: This class builds an import object for efficient importing of off-processor entries.

  /*! Import is used to construct a communication plan that can be called repeatedly by computational
      classes such the Tpetra CisMatrix and Vector classes to efficiently import entries from other
      images. An importer is used when we start out with a uniquely-owned distribution,
      and want to distribute that into a multiple-ownership distribution.

      This class currently has one constructor, taking two Map objects.
      The first Map specifies the distribution we have now. The second 
      Map specifies the distribution we want to have after importing.

      NOTE: Behavior is undefined if the source Map is not uniquely-owned.
  */

  template <typename Ordinal>
  class Import: public Teuchos::Object {

  public:

    //@{ \name Constructor/Destructor Methods

    //! Constructs a Import object from the source and target Maps.
    Import(const Map<Ordinal> & source, const Map<Ordinal> & target);

    //! copy constructor. 
    Import(const Import<Ordinal> & import);

    //! destructor.
    ~Import();

    //@}

    //@{ \name Export Attribute Methods

    //! Returns the number of entries that are identical between the source and target maps, up to the first different ID.
    Ordinal getNumSameIDs() const;

    //! Returns the number of entries that are local to the calling image, but not part of the first getNumSameIDs() entries.
    Ordinal getNumPermuteIDs() const;

    //! List of entries in the source map that are permuted.
    const std::vector<Ordinal> & getPermuteFromLIDs() const;

    //! List of entries in the target map that are permuted.
    const std::vector<Ordinal> & getPermuteToLIDs() const;

    //! Returns the number of entries that are not on the calling image.
    Ordinal getNumRemoteIDs() const;

    //! List of entries in the target map that are coming from other images.
    const std::vector<Ordinal> & getRemoteLIDs() const;

    //! Returns the number of entries that must be sent by the calling image to other images.
    Ordinal getNumExportIDs() const;

    //! List of entries in the source map that will be sent to other images.
    const std::vector<Ordinal> & getExportLIDs() const;

    //! List of images to which entries will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i].
    const std::vector<Ordinal> & getExportImageIDs() const;

    //! Returns the Source Map used to construct this importer.
    const Map<Ordinal> & getSourceMap() const;

    //! Returns the Target Map used to construct this importer.
    const Map<Ordinal> & getTargetMap() const;

    const Distributor<Ordinal> & getDistributor() const;

    //! Assignment operator
    Import<Ordinal>& operator = (const Import<Ordinal> & Source);

    //@}

    //@{ \name I/O Methods

    //! print method inherited from Teuchos::Object
    virtual void print(std::ostream& os) const;

    //@}

  private:

    Teuchos::RCP< ImportData<Ordinal> > ImportData_;

    // subfunctions used by constructor
    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and numRemoteIDs_
    // these variables are already initialized to 0 by the ImportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, and remoteLIDs_
    void setupSamePermuteRemote();

    //==============================================================================
    void setupExport();

    //==============================================================================
  };

  template <typename Ordinal>
  Import<Ordinal>::Import(const Map<Ordinal> & source, const Map<Ordinal> & target)
  : Teuchos::Object("Tpetra::Import")
  , ImportData_()
  {
    ImportData_ = Teuchos::rcp(new ImportData<Ordinal>(source, target));
    // call subfunctions
    setupSamePermuteRemote();
    if(source.isDistributed()) {
      setupExport();
    }
  }

  template <typename Ordinal>
  Import<Ordinal>::Import(const Import<Ordinal> & import)
  : Teuchos::Object(import.label())
  , ImportData_(import.ImportData_)
  {}

  template <typename Ordinal>
  Import<Ordinal>::~Import() 
  {}

  template <typename Ordinal>
  Ordinal Import<Ordinal>::getNumSameIDs() const {
    return ImportData_->numSameIDs_;
  }

  template <typename Ordinal>
  Ordinal Import<Ordinal>::getNumPermuteIDs() const {
    return ImportData_->numPermuteIDs_;
  }

  template <typename Ordinal>
  const std::vector<Ordinal> & Import<Ordinal>::getPermuteFromLIDs() const {
    return ImportData_->permuteFromLIDs_;
  }

  template <typename Ordinal>
  const std::vector<Ordinal> & Import<Ordinal>::getPermuteToLIDs() const {
    return ImportData_->permuteToLIDs_;
  }

  template <typename Ordinal>
  Ordinal Import<Ordinal>::getNumRemoteIDs() const {
    return ImportData_->numRemoteIDs_;
  }

  template <typename Ordinal>
  const std::vector<Ordinal> & Import<Ordinal>::getRemoteLIDs() const {
    return ImportData_->remoteLIDs_;
  }

  template <typename Ordinal>
  Ordinal Import<Ordinal>::getNumExportIDs() const {
    return ImportData_->numExportIDs_;
  }

  template <typename Ordinal>
  const std::vector<Ordinal> & Import<Ordinal>::getExportLIDs() const {
    return ImportData_->exportLIDs_;
  }

  template <typename Ordinal>
  const std::vector<Ordinal> & Import<Ordinal>::getExportImageIDs() const {
    return ImportData_->exportImageIDs_;
  }

  template <typename Ordinal>
  const Map<Ordinal> & Import<Ordinal>::getSourceMap() const {
    return ImportData_->source_;
  }

  template <typename Ordinal>
  const Map<Ordinal> & Import<Ordinal>::getTargetMap() const {
    return ImportData_->target_;
  }

  template <typename Ordinal>
  const Distributor<Ordinal>& Import<Ordinal>::getDistributor() const {
    return ImportData_->distributor_;
  }

  template <typename Ordinal>
  Import<Ordinal>& Import<Ordinal>::operator=(const Import<Ordinal> & Source) {
    ImportData_ = Source.ImportData_;
    return *this;
  }

  template <typename Ordinal>
  void Import<Ordinal>::print(std::ostream& os) const {
    using std::endl;
    os << "Import Data Members:" << endl;
    os << "permuteToLIDs_: " << toString(getPermuteToLIDs()) << endl;;
    os << "permuteFromLIDs_: " << toString(getPermuteFromLIDs()) << endl;
    os << "remoteLIDs_: " << toString(getRemoteLIDs()) << endl;
    os << "exportLIDs_: " << toString(getExportLIDs()) << endl;
    os << "exportImageIDs_: " << toString(getExportImageIDs()) << endl;
    os << "numSameIDs_: " << getNumSameIDs() << endl;
    os << "numPermuteIDs_: " << getNumPermuteIDs() << endl;
    os << "numRemoteIDs_: " << getNumRemoteIDs() << endl;
    os << "numExportIDs_: " << getNumExportIDs() << endl;
    os << "\nsource_: " << endl << getSourceMap();
    os << "\ntarget_: " << endl << getTargetMap();
  }


  template <typename Ordinal>
  void Import<Ordinal>::setupSamePermuteRemote() 
  {
    const Map<Ordinal> & source = getSourceMap();
    const Map<Ordinal> & target = getTargetMap();
    const std::vector<Ordinal> & sourceGIDs = source.getMyGlobalEntries();
    const std::vector<Ordinal> & targetGIDs = target.getMyGlobalEntries();

    // -- compute numSameIDs_ ---
    // go through GID lists of source and target. if the ith GID on both is the same, 
    // increment numSameIDs_ and try the next. as soon as you come to a pair that don't
    // match, give up.
    typename std::vector<Ordinal>::const_iterator sourceIter = sourceGIDs.begin(),
                                                  targetIter = targetGIDs.begin();
    while( sourceIter != sourceGIDs.end() && targetIter != targetGIDs.end() && 
        *sourceIter == *targetIter )
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
    for (; targetIter != targetGIDs.end(); ++targetIter) {
      if (source.isMyGlobalIndex(*targetIter)) {
        // both source and target list this GID (*targetIter)
        // determine the LIDs for this GID on both Maps and add them to the permutation lists
        ImportData_->permuteToLIDs_.push_back(target.getLocalIndex(*targetIter));
        ImportData_->permuteFromLIDs_.push_back(source.getLocalIndex(*targetIter));
      }
      else {
        ++ImportData_->numRemoteIDs_;
        ImportData_->remoteLIDs_.push_back(target.getLocalIndex(*targetIter));
        ImportData_->remoteGIDs_.push_back(*targetIter);
      }
    }
    ImportData_->numPermuteIDs_ = Teuchos::as<Ordinal>(ImportData_->permuteToLIDs_.size());
    ImportData_->numRemoteIDs_  = Teuchos::as<Ordinal>(ImportData_->remoteLIDs_.size());

    TEST_FOR_EXCEPTION( (ImportData_->numRemoteIDs_ > 0) && !source.isDistributed(), std::runtime_error, 
        "Tpetra::Import<" << Teuchos::OrdinalTraits<Ordinal>::name() 
        << ">::setupSamePermuteRemote(): Target has remote LIDs but Source is not distributed globally.");
  }


  template <typename Ordinal>
  void Import<Ordinal>::setupExport()
  {
    const Map<Ordinal> & source = getSourceMap();

    // create remoteImageID list: for each entry remoteGIDs[i],
    // remoteImageIDs[i] will contain the ImageID of the image that owns that GID.
    std::vector<Ordinal> remoteImageIDs(ImportData_->remoteGIDs_.size());
    source.getRemoteIndexList(ImportData_->remoteGIDs_, remoteImageIDs);

    // check for GIDs that exist in target but not in source
    // getRemoteIndexList will return -1 for the ImageID for any GIDs where this is the case
    if(ImportData_->numRemoteIDs_ > Teuchos::OrdinalTraits<Ordinal>::zero()) {
      const Ordinal negOne = Teuchos::OrdinalTraits<Ordinal>::zero() - Teuchos::OrdinalTraits<Ordinal>::one();
#ifdef HAVE_STD_NEW_COUNT_SYNTAX
      Ordinal count = std::count(remoteImageIDs.begin(), remoteImageIDs.end(), negOne);
#else
      Ordinal count = 0;
      std::count(remoteImageIDs.begin(), remoteImageIDs.end(), negOne, count);
#endif
      if(count > Teuchos::OrdinalTraits<Ordinal>::zero()) {
        throw reportError("Target has GIDs not found in Source.", 1);
      }
    }

    // sort remoteImageIDs in ascending order
    // apply same permutation to remoteGIDs_
    sortArrays(remoteImageIDs, ImportData_->remoteGIDs_);

    // call Distributor.createFromRecvs()
    // takes in numRemoteIDs_, remoteGIDs_, and remoteImageIDs_
    // returns numExportIDs_, exportLIDs_, and exportImageIDs_
    TEST_FOR_EXCEPT(ImportData_->remoteGIDs_.size() != (unsigned int)ImportData_->numRemoteIDs_); // FINISH: remove this
    TEST_FOR_EXCEPT(ImportData_->exportLIDs_.size() != (unsigned int)ImportData_->numExportIDs_ ); // FINISH: remove this
    ImportData_->distributor_.createFromRecvs(ImportData_->remoteGIDs_, remoteImageIDs, ImportData_->exportLIDs_, ImportData_->exportImageIDs_);

    // convert exportLIDs_ from GIDs to their LIDs in target
    for(typename std::vector<Ordinal>::iterator i = ImportData_->exportLIDs_.begin(); 
        i != ImportData_->exportLIDs_.end(); ++i) 
    {
      *i = source.getLocalIndex(*i);
    }
  }

} // namespace Tpetra

#include "Tpetra_ImportData.hpp"

#endif // TPETRA_IMPORT_HPP
