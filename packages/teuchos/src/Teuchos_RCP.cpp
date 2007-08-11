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

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace {


struct InfoAndCallNumber {
  InfoAndCallNumber()
    :call_number(-1)
  {}
  InfoAndCallNumber(
		    const std::string &info_in,
		    const int &call_number_in
		    )
    :info(info_in), call_number(call_number_in)

  {}
  std::string info;
  int call_number;
};


typedef std::map<Teuchos::PrivateUtilityPack::RCP_node*,InfoAndCallNumber>
rcp_node_list_t;


// Here we must let the PrintActiveRCPNodes constructor and destructor handle
// the creation and destruction of this map object.  This will ensure that
// this map object will be valid when any global/static RCP objects are
// destroyed!
rcp_node_list_t *rcp_node_list = 0;


} // namespace


namespace Teuchos {


void PrivateUtilityPack::throw_null( const std::string &type_name )
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"RCP<"<<type_name<<">::assert_not_null() : You can not"
    " call operator->() or operator*() if get()==NULL!" );
}


namespace PrivateUtilityPack {


void RCP_node::set_extra_data(
  const any &extra_data, const std::string& name
  ,EPrePostDestruction destroy_when
  ,bool force_unique
  )
{
  if(extra_data_map_==NULL) {
    extra_data_map_ = new extra_data_map_t;
  }
  const std::string type_and_name( extra_data.typeName() + std::string(":") + name );
  if( !extra_data_map_->empty() && force_unique ) {
    extra_data_map_t::iterator itr = extra_data_map_->find(type_and_name);
    TEST_FOR_EXCEPTION(
      itr != extra_data_map_->end(), std::invalid_argument
      ,"Error, the type:name pair \'" << type_and_name << "\' already exists and force_unique==true!" );
  }
  (*extra_data_map_)[type_and_name] = extra_data_entry_t(extra_data,destroy_when); // This may add or replace!
}


any& RCP_node::get_extra_data( const std::string& type_name, const std::string& name )
{
  TEST_FOR_EXCEPTION(
    extra_data_map_==NULL, std::invalid_argument
    ,"Error, no extra data has been set yet!" );
  any *extra_data = get_optional_extra_data(type_name,name);
  if(extra_data) return *extra_data;
  const std::string type_and_name( type_name + std::string(":") + name );
  TEST_FOR_EXCEPTION(
    extra_data == NULL, std::invalid_argument
    ,"Error, the type:name pair \'" << type_and_name << "\' is not found!" );
  return *extra_data; // Will never be executed!
}


any* RCP_node::get_optional_extra_data( const std::string& type_name, const std::string& name )
{
  if( extra_data_map_ == NULL ) return NULL;
  const std::string type_and_name( type_name + std::string(":") + name );
  extra_data_map_t::iterator itr = extra_data_map_->find(type_and_name);
  if(itr != extra_data_map_->end())
    return &(*itr).second.extra_data;
  return NULL;
}


void RCP_node::impl_pre_delete_extra_data()
{
  for( extra_data_map_t::iterator itr = extra_data_map_->begin(); itr != extra_data_map_->end(); ++itr ) {
    extra_data_map_t::value_type &entry = *itr;
    if(entry.second.destroy_when == PRE_DESTROY)
      entry.second.extra_data = any();
  }
}


} // namespace PrivateUtilityPack


// Define this macro here locally and rebuild just this *.cpp file and update
// the Teuchos library and you will get node tracing turned on when debugging
// support is enabled!  Note that you also have to TEUCHOS_DEBUG defined as
// well (using --enable-teuchos-debug at configure time).

//#define TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES


void PrivateUtilityPack::add_new_RCP_node( RCP_node* rcp_node, const std::string &info )
{
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
  TEST_FOR_EXCEPT(0==rcp_node_list);
  static int call_number = 0;
  (*rcp_node_list)[rcp_node] = InfoAndCallNumber(info,call_number);
  ++call_number;
#endif
}


void PrivateUtilityPack::remove_RCP_node( RCP_node* rcp_node )
{
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
  TEST_FOR_EXCEPT(0==rcp_node_list);
  const rcp_node_list_t::iterator itr = rcp_node_list->find(rcp_node);
  TEST_FOR_EXCEPT_PRINT(itr==rcp_node_list->end(),&std::cerr);
  rcp_node_list->erase(itr);
#endif
}


namespace PrivateUtilityPack {


PrintActiveRCPNodes::PrintActiveRCPNodes()
{
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
  std::cerr << "\nCalled PrintActiveRCPNodes::PrintActiveRCPNodes() : count = " << count_ << "\n";
#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
  if (!rcp_node_list)
    rcp_node_list = new rcp_node_list_t;
  ++count_;
}


PrintActiveRCPNodes::~PrintActiveRCPNodes()
{
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
  std::cerr << "\nCalled PrintActiveRCPNodes::~PrintActiveRCPNodes() : count = " << count_ << "\n";
#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
  if(--count_ == 0 ) {
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
    std::cerr << "\nPrint active nodes!\n";
#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
    std::cout << std::flush;
    print_active_RCP_nodes(std::cerr);
    TEST_FOR_EXCEPT(0==rcp_node_list);
    delete rcp_node_list;
  }
}


void PrintActiveRCPNodes::foo()
{
  int dummy = count_;
  ++dummy; // Avoid unused variable warning (bug 2664)
}


int PrintActiveRCPNodes::count_ = 0;


} // namespace PrivateUtilityPack


} // namespace Teuchos
 

void Teuchos::print_active_RCP_nodes(std::ostream &out)
{
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
  out
    << "\nCalled PrivateUtilityPack::print_active_RCP_nodes() :"
    << " rcp_node_list.size() = " << rcp_node_list.size() << "\n";
#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODE_TRACE
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
  TEST_FOR_EXCEPT(0==rcp_node_list);
  rcp_node_list_t::const_iterator itr = rcp_node_list->begin();
  if(itr != rcp_node_list->end()) {
    out
      << "\n***"
      << "\n*** Warning! The following Teucho::RCP_node objects were created but have"
      << "\n*** not been destoryed yet.  This may be an indication that these objects may"
      << "\n*** be involved in a circular dependency!  A memory checking tool may complain"
      << "\n*** that these objects are not destoryed correctly."
      << "\n***\n";
    while( itr != rcp_node_list->end() ) {
      const rcp_node_list_t::value_type &entry = *itr;
      out
        << "\n  RCP_node address = \'" << entry.first << "\',"
        << " information = " << entry.second.info << ","
        << " call number = " << entry.second.call_number;
      ++itr;
    }
    out << "\n";
  }
#endif // TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
}
