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

#include "Tpetra_ConfigDefs.h"
#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_CombineMode.hpp"

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
		
	*/

	template<typename OrdinalType, typename ScalarType>
	class DistObject: public Object {

	public:

		//! constructor
		DistObject(ElementSpace<OrdinalType> const elementspace)
			: Object("Tpetra::DistObject")
			, ElementSpace_(elementspace)
		{}

		//! constructor, taking label
		DistObject(ElementSpace<OrdinalType> const elementspace, std::string const& label)
			: Object(label)
			, ElementSpace_(elementspace)
		{}

		//! copy constructor
		DistObject(DistObject<OrdinalType, ScalarType> const& DistObject)
			: Object(DistObject.label())
			, ElementSpace_(DistObject.ElementSpace_)
		{}

		//! destructor
		virtual ~DistObject() {}

		//! Import
		int doImport(DistObject<OrdinalType, ScalarType> const& sourceObj, 
					 Import<OrdinalType> const& importer, 
					 CombineMode CM) 
		{
			// throw exception if my ElementSpace != importer.getTargetSpace()
			// throw exception if sourceObj's ElementSpace != importer.getSourceSpace()

			// copy variables from importer
			OrdinalType numSameIDs = importer.getNumSameIDs();
			OrdinalType numPermuteIDs = importer.getNumPermuteIDs();
			OrdinalType numRemoteIDs = importer.getNumRemoteIDs();
			OrdinalType numExportIDs = importer.getNumExportIDs();
			std::vector<OrdinalType> exportLIDs = importer.getExportLIDs();
			std::vector<OrdinalType> remoteLIDs = importer.getRemoteLIDs();
			std::vector<OrdinalType> permuteToLIDs = importer.getPermuteToLIDs();
			std::vector<OrdinalType> permuteFromLIDs = importer.getPermuteFromLIDs();

			// call doTransfer
			doTransfer(sourceObj, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
					   permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
					   ordinalExports_, scalarExports_, ordinalImports_, scalarImports_, 
					   importer.getDistributor(), false);
		}

		// ** TO DO: Add the three other forms of doImport/doExport **

		//! print method
		virtual void print(ostream& os) const {
			// ...
		}

		//! Accessor for whether or not this is a global object
		bool isGlobal() const {
			return(ElementSpace_.isGlobal());
		}

		//! Access function for the Tpetra::ElementSpace this DistObject was constructed with.
		ElementSpace<OrdinalType> const& elementspace() const {
			return(ElementSpace_);
		}

	protected:
 
		//! Perform actual transfer (redistribution) of data across memory images.
		virtual int doTransfer(DistObject<OrdinalType, ScalarType> const& sourceObj,
							   CombineMode CM,
							   OrdinalType numSameIDs,
							   OrdinalType numPermuteIDs,
							   OrdinalType numRemoteIDs,
							   OrdinalType numExportIDs,
							   std::vector<OrdinalType> permuteToLIDs,
							   std::vector<OrdinalType> permuteFromLIDs,
							   std::vector<OrdinalType> remoteLIDs,
							   std::vector<OrdinalType> exportLIDs,
							   std::vector<OrdinalType> ordinalExports,
							   std::vector<ScalarType> scalarExports,
							   std::vector<OrdinalType> ordinalImports,
							   std::vector<ScalarType> scalarImports,
							   Distributor<OrdinalType>& distor,
							   bool doReverse) {
			// ...
		}

		// The following four methods must be implemented by the derived class

		//! Allows the source and target (\e this) objects to be compared for compatibility, return false if not.
		virtual bool checkCompatibility(DistObject<OrdinalType, ScalarType> const& sourceObj) = 0;

		//! Perform copies and permutations that are local to this image.
		virtual int copyAndPermute(DistObject<OrdinalType, ScalarType> const& sourceObj,
								   OrdinalType numImportIDs,
								   OrdinalType numPermuteIDs,
								   std::vector<OrdinalType> permuteToLIDs,
								   std::vector<OrdinalType> permuteFromLIDs) = 0;

		//! Perform any packing or preparation required for call to doTransfer().
		virtual int packAndPrepare(DistObject<OrdinalType, ScalarType> const& sourceObj,
								   OrdinalType numExportIDs,
								   std::vector<OrdinalType> exportLIDs,
								   std::vector<OrdinalType> ordinalExports,
								   std::vector<ScalarType> scalarExports,
								   Distributor<OrdinalType>& distor) = 0;
  
		//! Perform any unpacking and combining after call to doTransfer().
		virtual int unpackAndCombine(DistObject<OrdinalType, ScalarType> const& sourceObj,
									 OrdinalType numImportIDs,
									 std::vector<OrdinalType> importLIDs,
									 std::vector<OrdinalType> ordinalImports,
									 std::vector<ScalarType> scalarImports,
									 Distributor<OrdinalType>& distor,
									 CombineMode CM) = 0;

	private:
		
		ElementSpace<OrdinalType> const ElementSpace_;
		std::vector<OrdinalType> ordinalImports_;
		std::vector<ScalarType> scalarImports_;
		std::vector<OrdinalType> ordinalExports_;
		std::vector<ScalarType> scalarExports_;

	}; // class DistObject

} // namespace Tpetra

#endif /* TPETRA_DISTOBJECT_HPP */
