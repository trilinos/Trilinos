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

#ifndef TPETRA_DISTOBJECT_HPP
#define TPETRA_DISTOBJECT_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Comm.hpp"

namespace Tpetra {

	//! Tpetra::DistObject: A class for constructing and using dense multi-vectors, vectors and matrices in parallel.

	/*! The DistObject is a base class for all Tpetra distributed global objects.  It provides the basic
	    mechanisms and interface specifications for importing and exporting operations using Tpetra::Import and
		Tpetra::Export objects.
		
		<b> Distributed Global vs. Replicated Local.</b>
		
		<ul>
		<li> Distributed Global objects - In most instances, a distributed object will be partitioned
		across multiple memory images associated with multiple processors.  In this case, there is 
		a unique copy of each element and elements are spread across all images specified by 
		the Tpetra::Platform object.
		<li> Replicated Local Objects - Some algorithms use objects that are too small to
		be distributed across all processors, the Hessenberg matrix in a GMRES
		computation.  In other cases, such as with block iterative methods,  block dot product 
		functions produce small dense matrices that are required by all images.  
		Replicated local objects handle these types of situation.
		</ul>

		DistObject error codes (positive for non-fatal, negative for fatal):
		<ol>
		<li> -1  importer's target ElementSpace does not match my ElementSpace
		<li> -2  importer's source ElementSpace does not match sourceObj's ElementSpace
		<li> -99 Internal DistObject error. Contact developer.
		</ol>
		
	*/

	template<typename OrdinalType, typename ScalarType>
	class DistObject: public Object {

	public:

		//@{ \name Constructor/Destructor Methods

		//! constructor
		DistObject(ElementSpace<OrdinalType> const elementspace, 
				   Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> > comm)
			: Object("Tpetra::DistObject")
			, ElementSpace_(elementspace)
			, Comm_(comm)
			, imports_()
			, exports_()
			, sizes_()
		{}

		//! constructor, taking label
		DistObject(ElementSpace<OrdinalType> const elementspace, 
				   Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> > comm,
				   std::string const& label)
			: Object(label)
			, ElementSpace_(elementspace)
			, Comm_(comm)
			, imports_()
			, exports_()
			, sizes_()
		{}

		//! copy constructor
		DistObject(DistObject<OrdinalType, ScalarType> const& DistObject)
			: Object(DistObject.label())
			, ElementSpace_(DistObject.ElementSpace_)
			, Comm_(DistObject.Comm_)
			, imports_(DistObject.imports_)
			, exports_(DistObject.exports_)
			, sizes_(DistObject.sizes_)
		{}

		//! destructor
		virtual ~DistObject() {}

		//@}

		//@{ \name Import/Export Methods

		//! Import
		int doImport(DistObject<OrdinalType, ScalarType> const& sourceObj, 
					 Import<OrdinalType> const& importer, CombineMode CM) 
		{
			// throw exception -1 if my ElementSpace != importer.getTargetSpace()
			if(elementspace() != importer.getTargetSpace())
				throw reportError("Target ElementSpaces don't match", -1);
			// throw exception -2 if sourceObj's ElementSpace != importer.getSourceSpace()
			if(sourceObj.elementspace() != importer.getSourceSpace())
				throw reportError("Source ElementSpaces don't match", -2);

			// copy variables from importer
			OrdinalType numSameIDs = importer.getNumSameIDs();
			OrdinalType numPermuteIDs = importer.getNumPermuteIDs();
			OrdinalType numRemoteIDs = importer.getNumRemoteIDs();
			OrdinalType numExportIDs = importer.getNumExportIDs();
			std::vector<OrdinalType> const& exportLIDs = importer.getExportLIDs();
			std::vector<OrdinalType> const& remoteLIDs = importer.getRemoteLIDs();
			std::vector<OrdinalType> const& permuteToLIDs = importer.getPermuteToLIDs();
			std::vector<OrdinalType> const& permuteFromLIDs = importer.getPermuteFromLIDs();

			// call doTransfer
			doTransfer(sourceObj, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
					   permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
					   exports_, imports_, importer.getDistributor(), false);

			return(0);
		}

		//! Export
		int doExport(DistObject<OrdinalType, ScalarType> const& sourceObj, 
					 Export<OrdinalType> const& exporter, CombineMode CM) 
		{
			// throw exception -1 if my ElementSpace != exporter.getTargetSpace()
			if(elementspace() != exporter.getTargetSpace())
				throw reportError("Target ElementSpaces don't match", -1);
			// throw exception -2 if sourceObj's ElementSpace != exporter.getSourceSpace()
			if(sourceObj.elementspace() != exporter.getSourceSpace())
				throw reportError("Source ElementSpaces don't match", -2);

			// copy variables from exporter
			OrdinalType numSameIDs = exporter.getNumSameIDs();
			OrdinalType numPermuteIDs = exporter.getNumPermuteIDs();
			OrdinalType numRemoteIDs = exporter.getNumRemoteIDs();
			OrdinalType numExportIDs = exporter.getNumExportIDs();
			std::vector<OrdinalType> const& exportLIDs = exporter.getExportLIDs();
			std::vector<OrdinalType> const& remoteLIDs = exporter.getRemoteLIDs();
			std::vector<OrdinalType> const& permuteToLIDs = exporter.getPermuteToLIDs();
			std::vector<OrdinalType> const& permuteFromLIDs = exporter.getPermuteFromLIDs();

			// call doTransfer
			doTransfer(sourceObj, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
					   permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
					   exports_, imports_, exporter.getDistributor(), false);

			return(0);
		}

		//! Import (using an Exporter)
		int doImport(DistObject<OrdinalType, ScalarType> const& sourceObj,
					 Export<OrdinalType> const& exporter, CombineMode CM)
		{
			// throw exception -1 if my ElementSpace != exporter.getSourceSpace()
			if(elementspace() != exporter.getSourceSpace())
				throw reportError("Target ElementSpaces don't match", -1);
			// throw exception -2 if sourceObj's ElementSpace != exporter.getTargetSpace()
			if(sourceObj.elementspace() != exporter.getTargetSpace())
				throw reportError("Source ElementSpaces don't match", -2);

			// copy variables from exporter
			OrdinalType numSameIDs = exporter.getNumSameIDs();
			OrdinalType numPermuteIDs = exporter.getNumPermuteIDs();
			OrdinalType numRemoteIDs = exporter.getNumRemoteIDs();
			OrdinalType numExportIDs = exporter.getNumExportIDs();
			std::vector<OrdinalType> const& exportLIDs = exporter.getExportLIDs();
			std::vector<OrdinalType> const& remoteLIDs = exporter.getRemoteLIDs();
			std::vector<OrdinalType> const& permuteToLIDs = exporter.getPermuteToLIDs();
			std::vector<OrdinalType> const& permuteFromLIDs = exporter.getPermuteFromLIDs();

			// call doTransfer
			doTransfer(sourceObj, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
					   permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
					   imports_, exports_, exporter.getDistributor(), true);

			return(0);
		}

		//! Export (using an Importer)
		int doExport(DistObject<OrdinalType, ScalarType> const& sourceObj,
					 Import<OrdinalType> const& importer, CombineMode CM)
		{
			// throw exception -1 if my ElementSpace != importer.getSourceSpace()
			if(elementspace() != importer.getSourceSpace())
				throw reportError("Target ElementSpaces don't match", -1);
			// throw exception -2 if sourceObj's ElementSpace != importer.getTargetSpace()
			if(sourceObj.elementspace() != importer.getTargetSpace())
				throw reportError("Source ElementSpaces don't match", -2);

			// copy variables from importer
			OrdinalType numSameIDs = importer.getNumSameIDs();
			OrdinalType numPermuteIDs = importer.getNumPermuteIDs();
			OrdinalType numRemoteIDs = importer.getNumRemoteIDs();
			OrdinalType numExportIDs = importer.getNumExportIDs();
			std::vector<OrdinalType> const& exportLIDs = importer.getExportLIDs();
			std::vector<OrdinalType> const& remoteLIDs = importer.getRemoteLIDs();
			std::vector<OrdinalType> const& permuteToLIDs = importer.getPermuteToLIDs();
			std::vector<OrdinalType> const& permuteFromLIDs = importer.getPermuteFromLIDs();

			// call doTransfer
			doTransfer(sourceObj, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
					   permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
					   imports_, exports_, importer.getDistributor(), true);

			return(0);
		}

		//@}

		//@{ \name I/O methods

		//! print method
		virtual void print(ostream& os) const {
			// ...
		}

		//@}

		//@{ \name Attribute Accessor Methods

		//! Accessor for whether or not this is a global object
		bool isGlobal() const {
			return(ElementSpace_.isGlobal());
		}

		//! Access function for the Tpetra::ElementSpace this DistObject was constructed with.
		ElementSpace<OrdinalType> const& elementspace() const {
			return(ElementSpace_);
		}

		//@}

	protected:
 
		//! Perform actual transfer (redistribution) of data across memory images.
		virtual int doTransfer(DistObject<OrdinalType, ScalarType> const& sourceObj,
							   CombineMode const CM,
							   OrdinalType const numSameIDs,
							   OrdinalType const numPermuteIDs,
							   OrdinalType const numRemoteIDs,
							   OrdinalType const numExportIDs,
							   std::vector<OrdinalType> const& permuteToLIDs,
							   std::vector<OrdinalType> const& permuteFromLIDs,
							   std::vector<OrdinalType> const& remoteLIDs,
							   std::vector<OrdinalType> const& exportLIDs,
							   std::vector<ScalarType>& exports,
							   std::vector<ScalarType>& imports,
							   Distributor<OrdinalType> const& distor,
							   bool doReverse) 
		{
			OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();

			checkSizes(sourceObj);

			if(numSameIDs + numPermuteIDs > zero)
				copyAndPermute(sourceObj, numSameIDs, numPermuteIDs, permuteToLIDs, permuteFromLIDs);

			// we don't have a "Zero" CombineMode like Epetra does, so we don't have to check for that

			OrdinalType packetSize = zero; // dummy value
			bool varSizes = false;
			if((!sizes_.empty()) && (numExportIDs > zero))
				sizes_.resize(numExportIDs);
			packAndPrepare(sourceObj, numExportIDs, exportLIDs, exports, packetSize, distor);

			if((isGlobal() && doReverse) || (sourceObj.elementspace().isGlobal() && !doReverse)) {
				// call one of the doPostsAndWaits functions
				if(doReverse) {
					if(varSizes)
						throw reportError("var-sized doReversePostsAndWaits not implemented yet", -99);
					else
						Comm_->doReversePostsAndWaits(distor, exports, packetSize, imports);
				}
				else {
					if(varSizes)
						throw reportError("var-sized doPostsAndWaits not implemented yet", -99);
					else
						Comm_->doPostsAndWaits(distor, exports, packetSize, imports);
				}
				unpackAndCombine(numRemoteIDs, remoteLIDs, imports, distor, CM);
			}

			return(0);
		}

		// The following four methods must be implemented by the derived class

		//! Allows the source and target (\e this) objects to be compared for compatibility.
		/*! Return true if they are compatible, return false if they aren't. */ 
		virtual bool checkSizes(DistObject<OrdinalType, ScalarType> const& sourceObj) = 0;

		//! Perform copies and permutations that are local to this image.
		/*!
		  \param sourceObj In
		         On entry, the DistObject that we are importing from.
		  \param numSameIDs In
		         On entry, the number of elements that are the same on the source and dest objects.
				 (i.e. The element is owned by the same image in both source and dest, 
				 and no permutation occurs.)
		  \param numPermuteIDs In
		         On entry, the number of elements that are locally permuted between source and dest objects.
		  \param permuteToLIDs In
		         On entry, contains a list of the elements that are permuted. (Listed by their LID in the
				 destination DistObject.)
		  \param permuteFromLIDs In
		         On entry, contains a list of the elements that are permuted. (Listed by their LID in the
				 source DistObject.)
		*/
		virtual int copyAndPermute(DistObject<OrdinalType, ScalarType> const& sourceObj,
								   OrdinalType const numSameIDs,
								   OrdinalType const numPermuteIDs,
								   std::vector<OrdinalType> const& permuteToLIDs,
								   std::vector<OrdinalType> const& permuteFromLIDs) = 0;

		//! Perform any packing or preparation required for communication.
		/*!
		  \param sourceObj In
		         On entry, the DistObject that we are importing from.
		  \param numExportIDs In
		         On entry, the number of elements we will be sending to other images.
		  \param exportLIDs In
		         On entry, a list of the elements we will be sending to other images.
				 (Listed by their LID in the source DistObject.)
		  \param exports Out
		         On exit, buffer for data we will be sending out.
		  \param packetSize Out
		         On exit, will contain the number of ScalarType variables used to pack
				 a single element.
		  \param distor In
		         On entry, contains the Distributor object we are using.				 
		*/
		virtual int packAndPrepare(DistObject<OrdinalType, ScalarType> const& sourceObj,
								   OrdinalType const numExportIDs,
								   std::vector<OrdinalType> const& exportLIDs,
								   std::vector<ScalarType>& exports,
								   OrdinalType& packetSize,
								   Distributor<OrdinalType> const& distor) = 0;
  
		//! Perform any unpacking and combining after communication.
		/*!
		  \param numImportIDs In
		         The number of elements we received from other images.
		  \param importLIDs In
		         On entry, a list of the elements we received from other images.
				 (Listed by their LID in the target DistObject.)
		  \param imports In
		         Buffer for data we received.
		  \param distor In
		         The Distributor object we are using.
		  \param CM In
		         The Tpetra::CombineMode to use when combining the imported
				 entries with existing entries.
		*/
		virtual int unpackAndCombine(OrdinalType const numImportIDs,
									 std::vector<OrdinalType> const& importLIDs,
									 std::vector<ScalarType> const& imports,
									 Distributor<OrdinalType> const& distor,
									 CombineMode const CM) = 0;

	private:
		
		ElementSpace<OrdinalType> const ElementSpace_;
		Teuchos::RefCountPtr< Comm<ScalarType, OrdinalType> > Comm_;
		std::vector<ScalarType> imports_;
		std::vector<ScalarType> exports_;
		std::vector<OrdinalType> sizes_;

	}; // class DistObject

} // namespace Tpetra

#endif /* TPETRA_DISTOBJECT_HPP */
