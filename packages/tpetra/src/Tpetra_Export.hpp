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

#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_Util.hpp"

namespace Tpetra {
	
	// forward declaration of ExportData, needed to prevent circular inclusions
	// actual #include statement at the end of this file
	template<typename OrdinalType> class ExportData;
	
	//! Tpetra::Export: This class builds an export object for efficient exporting of off-processor elements.
  
	/*! Export is used to construct a communication plan that can be called repeatedly by computational
        classes such the Tpetra CisMatrix and Vector classes to efficiently export elements.
    
		This class currently has one constructor, taking two ElementSpace objects.
		The first ElementSpace specifies the distribution we have now. The second 
		ElementSpace specifies the distribution we want to have after exporting.
	*/
  
	template <typename OrdinalType>
	class Export: public Object {
    
	public:
    
		//! Constructs a Export object from the source and target ElementSpaces.
		Export(ElementSpace<OrdinalType> const& source, ElementSpace<OrdinalType> const& target)
			: Object("Tpetra::Export")
			, ExportData_()
		{
			ExportData_ = Teuchos::rcp(new ExportData<OrdinalType>(source, target));
      
			// call subfunctions
			setupSamePermuteRemote();
			if(source.isGlobal())
				setupExport();
		}
    
		//! copy constructor. 
		Export(Export<OrdinalType> const& export)
			: Object(export.label())
			, ExportData_(export.ExportData_)
		{}
    
		//! destructor.
		~Export() {};
    
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

		//! Returns the Source ElementSpace used to construct this exporter.
		ElementSpace<OrdinalType> const& getSourceSpace() const {return(data().source_);};

		//! Returns the Target ElementSpace used to construct this exporter.
		ElementSpace<OrdinalType> const& getTargetSpace() const {return(data().target_);};
    
		//Distributor<ScalarType, OrdinalType>const& getDistributor() const {return(data().distributor_);}; // ST is PT
  	
		//! Assignment operator
		Export<OrdinalType>& operator = (Export<OrdinalType> const& Source) {
			ExportData_ = Source.ExportData_;
			return(*this);
		}
    
		//@{ \name I/O Methods
		//! print method inherited from Object
		virtual void print(ostream& os) const {
			os << "Export Data Members:" << endl;
			os << "permuteToLIDs_: " << toString(getPermuteToLIDs()) << endl;;
			os << "permuteFromLIDs_: " << toString(getPermuteFromLIDs()) << endl;
			os << "remoteLIDs_: " << toString(getRemoteLIDs()) << endl;
			os << "exportLIDs_: " << toString(getExportLIDs()) << endl;
			os << "exportImageIDs_: " << toString(getExportImageIDs()) << endl;
			os << "numSameIDs_: " << getNumSameIDs() << endl;
			os << "numPermuteIDs_: " << getNumPermuteIDs() << endl;
			os << "numRemoteIDs_: " << getNumRemoteIDs() << endl;
			os << "numExportIDs_: " << getNumExportIDs() << endl;
			os << "\nsource_: " << endl << getSourceSpace();
			os << "\ntarget_: " << endl << getTargetSpace();
			//os << "distributor_: N/A" << endl;
		};
		//@}
    
	private:
 
		Teuchos::RefCountPtr< ExportData<OrdinalType> > ExportData_;
    
		// convenience functions for returning inner data class, both const and nonconst versions.
		ExportData<OrdinalType>& data() {return(*ExportData_);}
		ExportData<OrdinalType> const& data() const {return(*ExportData_);}

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
#ifdef HAVE_STD_NEW_COUNT_SYNTAX
				OrdinalType count = std::count(remoteImageIDs.begin(), remoteImageIDs.end(), negOne);
#else
				OrdinalType count = 0;
				std::count(remoteImageIDs.begin(), remoteImageIDs.end(), negOne, count);
#endif
				if(count > Teuchos::OrdinalTraits<OrdinalType>::zero()) {
					throw reportError("Target has GIDs not found in Source.", 1);
				}
			}

			// sort remoteImageIDs in ascending order
			// apply same permutation to remoteGIDs_
			sortArrays(remoteImageIDs, data().remoteGIDs_);

			// create Distributor instance
			data().distributor_ = data().platform_->createDistributor();

			// call Distributor.createFromRecvs()
			// takes in numRemoteIDs_, remoteGIDs_, and remoteImageIDs_
			// returns numExportIDs_, exportLIDs_, and exportImageIDs_
			data().distributor_->createFromRecvs(data().numRemoteIDs_, data().remoteGIDs_, 
												 remoteImageIDs, true, data().numExportIDs_, 
												 data().exportLIDs_, data().exportImageIDs_);

			// convert exportLIDs_ from GIDs to their LIDs in target
			for(typename std::vector<OrdinalType>::iterator i = data().exportLIDs_.begin(); i != data().exportLIDs_.end(); i++)
				*i = source.getLID(*i);
		};
    
		//==============================================================================
	};
  
} // namespace Tpetra

#include "Tpetra_ExportData.hpp"

#endif // TPETRA_EXPORT_HPP
