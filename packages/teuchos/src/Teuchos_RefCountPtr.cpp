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

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

void PrivateUtilityPack::throw_null( const std::string &type_name )
{
	TEST_FOR_EXCEPTION(
		true, std::logic_error
		,"RefCountPtr<"<<type_name<<">::assert_not_null() : You can not "
		" call operator->() or operator*() if get()==NULL!" );
}

namespace PrivateUtilityPack {

void RefCountPtr_node::set_extra_data( const any &extra_data, const std::string& name, bool force_unique )
{
	if(extra_data_map_==NULL) {
		extra_data_map_ = new extra_data_map_t;
	}
	const std::string type_and_name( extra_data.type().name() + std::string(":") + name );
	if( !extra_data_map_->empty() && force_unique ) {
		extra_data_map_t::iterator itr = extra_data_map_->find(type_and_name);
		TEST_FOR_EXCEPTION(
			itr != extra_data_map_->end(), std::invalid_argument
			,"Error, the type:name pair \'" << type_and_name << "\' already exists and force_unique==true!" );
	}
	(*extra_data_map_)[type_and_name] = extra_data; // This may add or replace!
}

any& RefCountPtr_node::get_extra_data( const std::string& type_name, const std::string& name )
{
	TEST_FOR_EXCEPTION(
		extra_data_map_==NULL, std::invalid_argument
		,"Error, no extra data has been set yet!" );
	const std::string type_and_name( type_name + std::string(":") + name );
	extra_data_map_t::iterator itr = extra_data_map_->find(type_and_name);
	TEST_FOR_EXCEPTION(
		itr == extra_data_map_->end(), std::invalid_argument
		,"Error, the type:name pair \'" << type_and_name << "\' is not found!" );
	return (*itr).second;
}

} // namespace PrivateUtilityPack
} // namespace Teuchos

