// ////////////////////////////////////////////////////
// RefCountPtr.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

void PrivateUtilityPack::assert_not_null(const void *ptr)
{
	TEST_FOR_EXCEPTION(
		!ptr, std::logic_error
		,"RefCountPtr<...>::assert_not_null() : You can not "
		" call operator->() or operator*() if get() == 0" );
}

namespace PrivateUtilityPack {

int RefCountPtr_node::set_extra_data( const any &extra_data, int ctx )
{
	TEST_FOR_EXCEPTION(
		(extra_data_array_==NULL) && (ctx >= 0), std::invalid_argument
		,"Error, not existing extra data is set yet (set ctx to -1)!" );
	TEST_FOR_EXCEPTION(
		(extra_data_array_!=NULL) && (ctx > (int)extra_data_array_->size()-1), std::invalid_argument
		,"Error, the value of ctx = " << ctx << " > 0 is outside of the range [0,"
		<< (extra_data_array_->size()-1) << "] and therefore could not have been returned"
		" from a previous call to set_extra_data(...)!" );
	if(extra_data_array_==NULL) extra_data_array_ = new extra_data_array_t;
	if( ctx > 0 ) {
		(*extra_data_array_)[ctx] = extra_data;
	}
	else {
		extra_data_array_->push_back( extra_data );
		ctx = extra_data_array_->size()-1;
	}
	return ctx;
}

any& RefCountPtr_node::get_extra_data( int ctx )
{
	TEST_FOR_EXCEPTION(
		extra_data_array_==NULL, std::invalid_argument
		,"Error, no extra data has been set yet!" );
	TEST_FOR_EXCEPTION(
		(ctx < 0) || (extra_data_array_->size()-1) < ctx, std::invalid_argument
		,"Error, the value of ctx = " << ctx << " is outside of the range [0,"
		<< (extra_data_array_->size()-1) << "] and therefore could not have been returned"
		" from a previous call to set_extra_data(...)!" );
	return (*extra_data_array_)[ctx];
}

} // namespace PrivateUtilityPack
} // namespace Teuchos

