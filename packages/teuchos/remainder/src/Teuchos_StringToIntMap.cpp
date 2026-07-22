// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
		TEUCHOS_TEST_FOR_EXCEPTION(
			!unique, AlreadyExists
			,"Teuchos::StringToIntMap::StringToIntMap(...): "
			<< "Error, the std::string \"" << strings[i] << "\" is a duplicate for "
			<< defaultGroupName_ );
	}
}

int StringToIntMap::get( const std::string& option, const std::string& groupName ) const
{
	map_t::const_iterator itr = map_.find( option );
	TEUCHOS_TEST_FOR_EXCEPTION(
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
