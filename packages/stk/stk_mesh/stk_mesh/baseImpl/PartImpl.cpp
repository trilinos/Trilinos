/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


//----------------------------------------------------------------------
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/util/string_case_compare.hpp>
#include <stk_mesh/base/PartRelation.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/baseImpl/PartImpl.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

namespace impl {


void PartImpl::add_part_to_subset( Part & part) 
{ 
  insert( m_subsets, part ); 
}


void PartImpl::add_part_to_superset( Part & part ) 
{ 
  insert( m_supersets, part ); 
}

void PartImpl::add_relation( PartRelation relation )
{
  m_relations.push_back(relation);
}

PartImpl::~PartImpl()
{
}


// Subset part constructor:
PartImpl::PartImpl( MetaData          * arg_meta_data ,
            const std::string & arg_name ,
            EntityRank            arg_rank ,
            size_t              arg_ordinal )
  : m_name( arg_name ),
    m_attribute(),
    m_subsets() , m_supersets() , m_intersect() , m_relations() ,
    m_mesh_meta_data( arg_meta_data ),
    m_universe_ordinal( arg_ordinal ),
    m_entity_rank( arg_rank )
{}



//----------------------------------------------------------------------



//----------------------------------------------------------------------

} // namespace impl 

} // namespace mesh
} // namespace stk


