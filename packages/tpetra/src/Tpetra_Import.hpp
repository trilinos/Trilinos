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

#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra {
	
  // forward declaration of ImportData, needed to prevent circular inclusions
  // actual #include statement at the end of this file
  template<typename OrdinalType> class ImportData;
	
  //! Tpetra::Import: This class builds an import object for efficient importing of off-processor elements.
  
  /*! Import is used to construct a communication plan that can be called repeatedly by computational
      classes such the Tpetra CisMatrix and Vector classes to efficiently obtain off-processor
      elements.
    
      This class currently has one constructor, taking two VectorSpace objects.
      The first VectorSpace specifies the distribution we have now. The second 
      VectorSpace specifies the distribution we want to have after importing.
  */
  
  template <typename OrdinalType>
  class Import: public Object {
    
  public:
    
    //! Constructs a Import object from the source and target ElementSpaces.
    Import(ElementSpace<OrdinalType> const& source, ElementSpace<OrdinalType> const& target)
      : Object("Tpetra::Import")
      , ImportData_()
    {
      ImportData_ = Teuchos::rcp(new ImportData<OrdinalType>(source, target));
      
      // call subfunctions
      setupSamePermuteRemote();
      if(source.isGlobal())
        setupExport();
    }
    
    //! copy constructor. 
    Import(Import<OrdinalType> const& import)
      : Object(import.label())
      , ImportData_(import.ImportData_)
    {}
    
  	//! destructor.
  	~Import() {};
    
  	//! Returns the number of elements that are identical between the source and target spaces, up to the first different ID.
  	OrdinalType getNumSameIDs() const {return(data().numSameIDs_);};
    
  	//! Returns the number of elements that are local to the calling image, but not part of the first getNumSameIDs() elements.
  	OrdinalType getNumPermuteIDs() const {return(data().numPermuteIDs_);};

  	//! List of elements in the source space that are permuted.
  	std::vector<OrdinalType> const& getPermuteFromLIDs() const {return(data().permuteFromLIDs_);};

  	//! List of elements in the target space that are permuted.
  	std::vector<OrdinalType> const& getPermuteToLIDs() const {return(data().permuteToLIDs_);};

  	//! Returns the number of elements that are not on the calling image.
  	OrdinalType getNumRemoteIDs() const {return(data().numRemoteIDs_);};
  
  	//! List of elements in the target space that are coming from other images.
  	std::vector<OrdinalType> const& getRemoteLIDs() const {return(data().remoteLIDs_);};

  	//! Returns the number of elements that must be sent by the calling image to other images.
  	OrdinalType getNumExportIDs() const {return(data().numExportIDs_);};

  	//! List of elements in the source space that will be sent to other images.
  	std::vector<OrdinalType> const& getExportLIDs() const {return(data().exportLIDs_);};

  	//! List of images to which elements will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i].
  	std::vector<OrdinalType> const& getExportImageIDs() const {return(data().exportImageIDs_);};

  	//! Returns the Source ElementSpace used to construct this importer.
  	ElementSpace<OrdinalType> const& getSourceSpace() const {return(data().source_);};

    //! Returns the Target ElementSpace used to construct this importer.
    ElementSpace<OrdinalType> const& getTargetSpace() const {return(data().target_);};
    
    //Distributor<ScalarType, OrdinalType>const& getDistributor() const {return(data().distributor_);}; // ST is PT
  	
  	//! Assignment operator
    Import<OrdinalType>& operator = (Import<OrdinalType> const& Source) {
      ImportData_ = Source.ImportData_;
      return(*this);
    }
    
  	//@{ \name I/O Methods
    //! print method inherited from Object
  	virtual void print(ostream& os) const {
      os << "Import Data Members:" << endl;
      os << "permuteToLIDs_: " << getPermuteToLIDs() << endl;;
      os << "permuteFromLIDs_: " << getPermuteFromLIDs() << endl;
      os << "remoteLIDs_: " << getRemoteLIDs() << endl;
      //os << "exportLIDs_: N/A" << endl;
      //os << "exportImageIDs_: N/A" << endl;
      os << "numSameIDs_: " << getNumSameIDs() << endl;
      os << "numPermuteIDs_: " << getNumPermuteIDs() << endl;
      os << "numRemoteIDs_: " << getNumRemoteIDs() << endl;
      //os << "numExportIDs_: N/A" << endl;
      os << "\nsource_: " << endl << getSourceSpace();
      os << "\ntarget_: " << endl << getTargetSpace();
      //os << "distributor_: N/A" << endl;
    };
  	//@}
    
  private:
 
    Teuchos::RefCountPtr< ImportData<OrdinalType> > ImportData_;
    
  	// convenience functions for returning inner data class, both const and nonconst versions.
  	ImportData<OrdinalType>& data() {return(*ImportData_);}
  	ImportData<OrdinalType> const& data() const {return(*ImportData_);}

    // subfunctions used by constructor
    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and numRemoteIDs_
    // these variables are already initialized to 0 by the ImportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, and remoteLIDs_
    void setupSamePermuteRemote() {
      ElementSpace<OrdinalType> const& source = getSourceSpace();
      ElementSpace<OrdinalType> const& target = getTargetSpace();
      std::vector<OrdinalType> const& sourceGIDs = source.getMyGlobalElements();
      std::vector<OrdinalType> const& targetGIDs = target.getMyGlobalElements();
      
      // -- compute numSameIDs_ ---
      // go through GID lists of source and target. if the ith GID on both is the same, 
      // increment numSameIDs_ and try the next. as soon as you come to a pair that don't
      // match, give up.
      typename std::vector<OrdinalType>::const_iterator sourceIter = sourceGIDs.begin();
      typename std::vector<OrdinalType>::const_iterator targetIter = targetGIDs.begin();
      while((sourceIter != sourceGIDs.end()) && 
            (targetIter != targetGIDs.end()) && 
            (*sourceIter == *targetIter)) {
        data().numSameIDs_++;
        sourceIter++;
        targetIter++;
      }
      // targetIter should now point to the GID of the first non-same entry
      
      // -- compute numPermuteIDs and numRemoteIDs --
      // -- fill permuteToLIDs_, permuteFromLIDs_, remoteGIDs_, and remoteLIDs_ --
      // go through remaining entries in targetGIDs. if source owns that GID, 
      // increment numPermuteIDs_, and add entries to permuteToLIDs_ and permuteFromLIDs_.
      // otherwise increment numRemoteIDs_ and add entries to remoteLIDs_ and remoteGIDs_.
      for(; targetIter != targetGIDs.end(); targetIter++) {
        if(source.isMyGID(*targetIter)) {
          data().numPermuteIDs_++;
          data().permuteToLIDs_.push_back(target.getLID(*targetIter));
          data().permuteFromLIDs_.push_back(source.getLID(*targetIter));
        }
        else {
          data().numRemoteIDs_++;
          data().remoteLIDs_.push_back(target.getLID(*targetIter));
          data().remoteGIDs_.push_back(*targetIter);
        }
      }

      if((data().numRemoteIDs_ > 0) && (!source.isGlobal()))
        throw reportError("Target has remote LIDs but Source is not distributed globally.", 1); //*** what do we do here??? ***
    };

    //==============================================================================
    void setupExport() {
      ElementSpace<OrdinalType> const& source = getSourceSpace();

      // create remoteImageID list: for each entry remoteGIDs[i],
      // remoteImageIDs[i] will contain the ImageID of the image that owns that GID.
      std::vector<OrdinalType> remoteImageIDs(data().remoteGIDs_.size());
      source.getRemoteIDList(data().remoteGIDs_, remoteImageIDs);
      
      // check for GIDs that exist in target but not in source
      // getRemoteIDList will return -1 for the ImageID for any GIDs where this is the case
      if(data().numRemoteIDs_ > Teuchos::OrdinalTraits<OrdinalType>::zero()) {
        OrdinalType const negOne = Teuchos::OrdinalTraits<OrdinalType>::zero() - Teuchos::OrdinalTraits<OrdinalType>::one();
        OrdinalType count = std::count(remoteImageIDs.begin(), remoteImageIDs.end(), negOne); // count = number of occurances of -1 in remoteImageIDs
        if(count > Teuchos::OrdinalTraits<OrdinalType>::zero()) {
          throw reportError("Target has GIDs not found in Source.", 1);
        }
      }

      // sort remoteImageIDs in ascending order
      // apply same permutation to remoteGIDs_
      sortArrays(remoteImageIDs, data().remoteGIDs_);

      // create Distributor instance
      // call Distributor.createFromRecvs()
      // takes in numRemoteIDs_, remoteGIDs_, and remoteImageIDs_
      // returns numExportIDs_, exportLIDs_, and exportImageIDs_
      data().distributor_->createFromRecvs(data().numRemoteIDs_, data().remoteGIDs_, remoteImageIDs, true, data().numExportIDs_, data().exportLIDs_, data().exportImageIDs_);

      // convert exportLIDs_ from GIDs to their LIDs in target
      for(typename std::vector<OrdinalType>::iterator i = data().exportLIDs_.begin(); i != data().exportLIDs_.end(); i++)
        *i = source.getLID(*i);
      };
    
    //==============================================================================
  };
  
} // namespace Tpetra

#include "Tpetra_ImportData.hpp"

#endif // TPETRA_IMPORT_HPP
