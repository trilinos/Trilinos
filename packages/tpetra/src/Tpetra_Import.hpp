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
      setupIDLists();
      setupLIDs();
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
  	std::vector<OrdinalType> const& getNumExportIDs() const {return(&data().numExportIDs_);};

  	//! List of elements that will be sent to other images.
  	std::vector<OrdinalType> const& getExportLIDs() const {return(&data().exportLIDs_);};

  	//! List of images to which elements will be sent, getExportLIDs() [i] will be sent to image getExportImageIDs() [i].
  	std::vector<OrdinalType> const& getExportImageIDs() const {return(&data().exportImageIDs_);};

  	//! Total number of elements to be sent.
  	OrdinalType getNumSend() const {return(data().numSend_);};

  	//! Total number of elements to be received.
  	OrdinalType getNumRecv() const {return(data().numRecv_);};

  	//! Returns the Source ElementSpace used to construct this importer.
  	ElementSpace<OrdinalType> const& getSourceSpace() const {return(data().source_);};

    //! Returns the Target ElementSpace used to construct this importer.
    ElementSpace<OrdinalType> const& getTargetSpace() const {return(data().target_);};
    
    //Distributor<ScalarType, OrdinalType>const& getDistributor() const {return(*Distributor_);}; // ST is PT
  	
  	//! Assignment operator
    Import<OrdinalType>& operator = (Import<OrdinalType> const& Source) {
      ImportData_ = Source.ImportData_;
      return(*this);
    }
    
  	//@{ \name I/O Methods
    //! print method inherited from Object
  	virtual void print(ostream& os) const {
      os << "Import Data Members:" << endl;
      os << "permuteToLIDs_: "; printVector(os, data().permuteToLIDs_);
      os << "permuteFromLIDs_: "; printVector(os, getPermuteFromLIDs());
      os << "remoteLIDs_: "; printVector(os, getRemoteLIDs());
      //os << "exportLIDs_: " << endl;
      //os << "exportImageIDs_: " << endl;
      os << "numSameIDs_: " << getNumSameIDs() << endl;
      os << "numPermuteIDs_: " << getNumPermuteIDs() << endl;
      os << "numRemoteIDs_: " << getNumRemoteIDs() << endl;
      //os << "numExportIDs_: " << endl;
      //os << "numSend_: " << endl;
      //os << "numRecv_: " << endl;
      os << "\nsource_: " << endl << getSourceSpace();
      os << "\ntarget_: " << endl << getTargetSpace();
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
    void setupIDLists() {
      ElementSpace<OrdinalType> const& source = getSourceSpace();
      ElementSpace<OrdinalType> const& target = getTargetSpace();
      OrdinalType const sourceLength = source.getNumMyElements();
      OrdinalType const targetLength = target.getNumMyElements();
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
      // go through remaining entries in targetGIDs. if source owns that GID, 
      // increment numPermuteIDs_, otherwise increment numRemoteIDs_.
      for(; targetIter != targetGIDs.end(); targetIter++) {
        if(source.isMyGID(*targetIter))
          data().numPermuteIDs_++;
        else
          data().numRemoteIDs_++;
      }
    };

    //==============================================================================
    // sets up permuteToLIDs_, permuteFromLIDs_, remoteLIDs_, and numRecv_
    // numRecv_ has been initialized to 0, others are uninitialized
    // ifGlobal section = ???
    void setupLIDs() {
      ElementSpace<OrdinalType> const& source = getSourceSpace();
      ElementSpace<OrdinalType> const& target = getTargetSpace();
      OrdinalType const targetLength = target.getNumMyElements();

      // go through non-same elements. if that GID is owned by source, add entries to
      // permuteToLIDs_ and permuteFromLIDs_. otherwise add entries to remoteGIDs and remoteLIDs_.
      std::vector<OrdinalType> remoteGIDs(data().numRemoteIDs_);
      for(OrdinalType i = data().numSameIDs_; i < targetLength; i++) {
        if(source.isMyGID(target.getGID(i))) {
          data().permuteToLIDs_.push_back(i);
          data().permuteFromLIDs_.push_back(source.getLID(target.getGID(i)));
        }
        else {
          data().numRecv_++;
          remoteGIDs.push_back(target.getGID(i));
          data().remoteLIDs_.push_back(i);
        }
      }

      if((data().numRemoteIDs_ > 0) && (!source.isGlobal()))
        throw reportError("Serial Import has remote IDs. (Importing to Subset of Target ElementSpace)", 1);
    };

    //==============================================================================
    /*void setupIfGlobal() {
      if(source.isGlobal()) {
        // create remoteImageID list: for each entry remoteGIDs[i],
        // remoteImageIDs[i] will contain the ImageID of the image that owns that GID.
        std::vector<OrdinalType> remoteImageIDs(remoteGIDs.size());
        source.getRemoteIDList(data().numRemoteIDs_, remoteGIDs, remoteImageIDs);
        
        // remove IDs that don't exist in source
        if(data().numRemoteIDs_ > 0) {
          OrdinalType const negOne = Teuchos::OrdinalTraits<OrdinalType>::zero() - Teuchos::OrdinalTraits<OrdinalType>::one();
          OrdinalType cnt = std::count(remoteImageIDs.begin(), remoteImageIDs.end(), negOne); // cnt = number of occurances of -1 in remoteImageIDs
          if(cnt > 0 && cnt < data().numRemoteIDs_) {
            std::vector<OrdinalType> newRemoteGIDs(data().numRemoteIDs_ - cnt);
            std::vector<OrdinalType> newRemoteImageIDs(data().numRemoteIDs_ - cnt);
            cnt = 0;
            for(OrdinalType i = Teuchos::OrdinalTraits<OrdinalType>::zero(); i < data().numRemoteIDs_; i++) {
              if(remoteImageIDs[i] == negOne) {
                newRemoteGIDs[cnt] = remoteGIDs[i];
                newRemoteImageIDs[cnt] = remoteImageIDs[i];
                cnt++;
              }
            }
            data().numRemoteIDs_ = cnt;
            remoteGIDs.assign(newRemoteGIDs.begin(), newRemoteGIDs.end());
            remoteImageIDs.assign(newRemoteImageIDs.begin(), newRemoteImageIDs.end());
            throw reportError("Target IDs not found in Source ElementSpace (Do you want to import to subset of Target ElementSpace?", 1);
          }
          if(cnt >= data().numRemoteIDs_) {
            data().numRemoteIDs_ = Teuchos::OrdinalTraits<OrdinalType>::zero();
            remoteGIDs.clear();
            remoteImageIDs.clear();
          }
        } 

        // sort remote IDs by ImageID
      }
    };*/
    
    //==============================================================================
  };
  
} // namespace Tpetra

#include "Tpetra_ImportData.hpp"

#endif // TPETRA_IMPORT_HPP
