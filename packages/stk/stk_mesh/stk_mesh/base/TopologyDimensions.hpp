// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_mesh_TopologyDimensions_hpp
#define stk_mesh_TopologyDimensions_hpp

#include <Shards_Array.hpp>             // for ArrayDimTag
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <string>                       // for string
#include "stk_topology/topology.hpp"    // for topology, etc



namespace stk {
namespace mesh {

/**
 * The file contains Field typed and ArrayDimTags that are useful for
 * setting up various types of field relations.
 */

//----------------------------------------------------------------------
/** \ingroup stk_mesh_field_dimension_tags
 *  \brief  Define an array dimension of the number of nodes per element.
 */
class ElementNode : public shards::ArrayDimTag {
public:
  const char * name() const ;
  static const ElementNode & tag(); ///< \brief Singleton
private:
  ElementNode() {}
  ElementNode( const ElementNode & );
  ElementNode & operator = ( const ElementNode & );
};

/** \ingroup stk_mesh_relation_stencil
 *  \brief  An element Field defining an array of values, one value per
 *          node of the element.
 */
typedef Field<double,ElementNode> ElementNodeField ;

/** \ingroup stk_mesh_relation_stencil
 *  \brief  A Field defining an array of pointers
 *          to an element's nodal field data.
 */

/** \ingroup stk_mesh_relation_stencil
 *  \brief  Declare an element-node field.
 */
inline
ElementNodeField &
declare_element_node_field( MetaData & md , const std::string & s )
{

  ElementNodeField & f =
    md.declare_field< ElementNodeField >( stk::topology::ELEMENT_RANK, s, 1 /* 1 state */ );

  return f ;
}

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_field_dimension_tags
 *  \{
 */

struct QuadratureTag : public shards::ArrayDimTag {
  const char * name() const ;
  static const QuadratureTag & tag(); ///< \brief Singleton
private:
  QuadratureTag() {}
  QuadratureTag( const QuadratureTag & );
  QuadratureTag & operator = ( const QuadratureTag & );
};

struct BasisTag : public shards::ArrayDimTag {
  const char * name() const ;
  static const BasisTag & tag(); ///< \brief Singleton
private:
  BasisTag() {}
  BasisTag( const BasisTag & );
  BasisTag & operator = ( const BasisTag & );
};


/** \} */
//----------------------------------------------------------------------

}//namespace mesh
}//namespace stk

#endif

