// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#include "Teuchos_StringToIntMap.hpp"

namespace Teuchos {

StringToIntMap::StringToIntMap(
  const std::string& defaultGroupName_in, int n, const char* strings[]
  ) : defaultGroupName_(defaultGroupName_in)
{
	typedef map_t::value_type val_t;
	for( int i = 0; i < n; ++i ) {
		const bool unique = map_.insert( val_t( strings[i], i ) ).second;
		TEST_FOR_EXCEPTION(
			!unique, AlreadyExists
			,"Teuchos::StringToIntMap::StringToIntMap(...): "
			<< "Error, the std::string \"" << strings[i] << "\" is a duplicate for "
			<< defaultGroupName_ );
	}
}

int StringToIntMap::get( const std::string& option, const std::string& groupName ) const
{
	map_t::const_iterator itr = map_.find( option );
	TEST_FOR_EXCEPTION(
		itr == map_.end(), DoesNotExist
		,"Teuchos::StringToIntMap:::get(\""<<option<<"\",...): "
		<< "Error, the std::string \"" << option << "\" is not recongnised for "
		<< ( groupName.length() ? groupName : defaultGroupName_ )
    << "; valid selections include " << this->validSelections() << "."
    );
	return (*itr).second;	
}

// private

std::string StringToIntMap::validSelections() const
{
  std::ostringstream oss;
  oss << "{";
  map_t::const_iterator itr = map_.begin();
  for( int i = 0; itr != map_.end(); ++itr, ++i ) {
    if(i > 0)
      oss << ",";
    oss << "\""<<itr->first<<"\":"<<itr->second;
  }
  oss << "}";
  return oss.str();
}

} // end namespace Teuchos
