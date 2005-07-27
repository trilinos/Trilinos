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

#ifndef TPETRA_EXPORTDATA_HPP
#define TPETRA_EXPORTDATA_HPP

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_Export.hpp"

namespace Tpetra {
  
	template<typename OrdinalType>
	class ExportData : public Object {
		friend class Export<OrdinalType>;
	public:
		// default constructor
		ExportData(ElementSpace<OrdinalType> const& source, ElementSpace<OrdinalType> const& target)
			: Object("Tpetra::ExportData")
			// permuteToLIDs_, permuteFromLIDs_, remoteLIDs_, remoteGIDs_, exportLIDs_, exportGIDs_, 
			// and exportImageIDs_ don't need/use a member initialization. (These are all STL vectors.)
			, numSameIDs_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, numPermuteIDs_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, numRemoteIDs_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, numExportIDs_(Teuchos::OrdinalTraits<OrdinalType>::zero())
			, source_(source)
			, target_(target)
			, platform_(source.platform().clone())
			, distributor_(source.platform().createOrdinalComm())
		{};
    
		// destructor. no heap-data, so no need to override
		~ExportData() {};
    
	protected:
		// OT vectors
		Teuchos::Array<OrdinalType> permuteToLIDs_;
		Teuchos::Array<OrdinalType> permuteFromLIDs_;
		Teuchos::Array<OrdinalType> remoteLIDs_;
		Teuchos::Array<OrdinalType> remoteGIDs_;
		Teuchos::Array<OrdinalType> exportLIDs_;
		Teuchos::Array<OrdinalType> exportGIDs_;
		Teuchos::Array<OrdinalType> exportImageIDs_;
    
		// OTs
		OrdinalType numSameIDs_;
		OrdinalType numPermuteIDs_;
		OrdinalType numRemoteIDs_;
		OrdinalType numExportIDs_;
    
		// ElementSpaces
		ElementSpace<OrdinalType> const source_;
		ElementSpace<OrdinalType> const target_;
    
		// Platform, Comm, Distributor, etc.
		Teuchos::RefCountPtr< Platform<OrdinalType, OrdinalType> const > platform_;
		Distributor<OrdinalType> distributor_;
    
	private:
		//! Copy constructor (declared but not defined, do not use)
		ExportData(ExportData<OrdinalType> const& rhs);
		//! Assignment operator (declared but not defined, do not use)
		ExportData<OrdinalType>& operator = (ExportData<OrdinalType> const& rhs);
    
	}; // class ExportData
  
} // namespace Tpetra

#endif // TPETRA_EXPORTDATA_HPP
