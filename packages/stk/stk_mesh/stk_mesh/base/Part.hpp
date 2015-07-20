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

#ifndef stk_mesh_Part_hpp
#define stk_mesh_Part_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for int64_t
#include <algorithm>                    // for sort, unique
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for PartVector, OrdinalVector, etc
#include <stk_mesh/baseImpl/PartImpl.hpp>  // for PartImpl
#include <string>                       // for string, basic_string
#include <vector>                       // for vector, vector<>::iterator
#include "stk_topology/topology.hpp"    // for topology
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { namespace impl { class PartRepository; } } }


//----------------------------------------------------------------------

namespace stk {
namespace mesh {

namespace impl {
} // namespace impl

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  An application-defined subset of a problem domain.
 *
 *  An application may define parts corresponding to geometric subdomains,
 *  modeling subdomains, material subdomains, parallel distributed subdomains,
 *  entities of the same type of discretization
 *  (e.g., hexahedrons, tetrahedrons, ...),
 *  or any other application need for subsetting.
 *
 *  A Part is created, owned, and modified by a
 *  \ref stk::mesh::PartRepository "Part manager".
 */
class Part {
public:
    enum { INVALID_ID = -1 };

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns the PartRepository which created this part.
   */
  MetaData & mesh_meta_data() const { return m_partImpl.mesh_meta_data(); }

  BulkData & mesh_bulk_data() const;

  /** \brief  The primary entity type for this part.
   *
   *   For example, the primary purpose of an Element part
   *   is to define a collection of elements.  However, the
   *   nodes of those elements are also members of an element part.
   *   Return InvalidEntityRank if no primary entity type.
   */
  EntityRank primary_entity_rank() const { return m_partImpl.primary_entity_rank(); }

  stk::topology topology() const { return m_partImpl.topology(); }

  /** \brief  Application-defined text name of this part, must be unique within the set of parts owned by a MetaData*/
  const std::string & name() const { return m_partImpl.name(); }

  bool force_no_induce() const { return m_partImpl.force_no_induce(); }

  /** \brief Should we induce this part
   *
   * We need the from_rank argument because part induction does not chain.
   * IE, induced parts are not transitive. from_rank represents the rank
   * of the entity that is the src entity of the connectivity.
   */
  bool should_induce(EntityRank from_rank) const
  { return primary_entity_rank() == from_rank && !force_no_induce(); }

  /** \brief Could an entity of a certain rank be induced into part
   */
  bool was_induced(EntityRank rank) const
  {
    return primary_entity_rank() != InvalidEntityRank &&
           rank < primary_entity_rank() &&
           !force_no_induce();
  }

  //whether an entity's membership in this part is required to be the same on each processor that shares/ghosts the entity.
  bool entity_membership_is_parallel_consistent() const { return m_partImpl.entity_membership_is_parallel_consistent(); }
  void entity_membership_is_parallel_consistent(bool trueOrFalse) { m_partImpl.entity_membership_is_parallel_consistent(trueOrFalse); }

  int64_t id() const { return m_partImpl.id(); }

  /** \brief  Internally generated ordinal of this part that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_partImpl.mesh_meta_data_ordinal(); }

  /** \brief  Parts that are supersets of this part. */
  const PartVector & supersets() const { return m_partImpl.supersets(); }

  /** \brief  Parts that are subsets of this part. */
  const PartVector & subsets() const { return m_partImpl.subsets(); }

  /** \brief  Check if argument is subset of this */
  bool contains(const Part& part) const;

  /** \brief  Equality comparison */
  bool operator == ( const Part & rhs ) const { return this == & rhs ; }

  /** \brief  Inequality comparison */
  bool operator != ( const Part & rhs ) const { return this != & rhs ; }

  /** \brief  Query attribute that has been attached to this part */
  template<class A>
  const A * attribute() const { return m_partImpl.attribute<A>(); }

private:

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns the PartRepository which created this part.
   */
  MetaData & meta_data() const { return m_partImpl.mesh_meta_data(); }


  impl::PartImpl m_partImpl;

  /* \brief  A part is owned by a PartRepository, as such only the owning
   *         PartRepository can create, delete, or modify a part.
   *         The owner-modifies rule is enforced by the implementation being
   *         a private data object on the Part and the PartRepository is a
   *         friend.
   */
  friend class ::stk::mesh::impl::PartRepository ;
  friend class ::stk::mesh::MetaData ;

#ifndef DOXYGEN_COMPILE

  /** Construct a subset part within a given mesh.
   *  Is used internally by the two 'declare_part' methods on PartRepository.
   */
  Part( MetaData * arg_meta_data , const std::string & arg_name, EntityRank arg_rank, size_t arg_ordinal, bool arg_force_no_induce = false)
    : m_partImpl(arg_meta_data, arg_name, arg_rank, arg_ordinal, arg_force_no_induce)
  { }

  ~Part() {}
  Part();
  Part( const Part & );
  Part & operator = ( const Part & );

#endif /* DOXYGEN_COMPILE */
};

//----------------------------------------------------------------------
/** \brief  Ordering operator for parts. */
struct PartLess {

  inline bool operator()( const Part & lhs , const Part & rhs ) const
    { return lhs.mesh_meta_data_ordinal() < rhs.mesh_meta_data_ordinal(); }

  inline bool operator()( const Part & lhs , const Part * rhs ) const
    { return lhs.mesh_meta_data_ordinal() < rhs->mesh_meta_data_ordinal(); }

  inline bool operator()( const Part * lhs , const Part & rhs ) const
    { return lhs->mesh_meta_data_ordinal() < rhs.mesh_meta_data_ordinal(); }

  inline bool operator()( const Part * lhs , const Part * rhs ) const
    { return lhs->mesh_meta_data_ordinal() < rhs->mesh_meta_data_ordinal(); }
};

/** \brief  Order a collection of parts: invoke sort and then unique */
void sort_and_unique( PartVector &partVector );

/** \brief  Insert a part into a properly ordered collection of parts.
 *          Returns true if this is a new insertion.
 */
bool insert( ConstPartVector & , const Part & );
bool insert( PartVector & , Part & );

inline
bool insert_ordinal( OrdinalVector & v , unsigned part_ordinal )
{
  for(OrdinalVector::iterator i=v.begin(), e=v.end(); i!=e; ++i) {
    if (*i == part_ordinal) return false;
    if (*i > part_ordinal) {
      v.insert(i, part_ordinal);
      return true;
    }
  }

  v.push_back(part_ordinal);
  return true ;
}

void get_part_and_all_subsets(const Part& part, ConstPartVector& part_and_all_subsets);

/** \brief  Remove a part from a properly ordered collection of parts. */
void remove( PartVector & , Part & );

/** \brief  Find a part by name in a collection of parts. */
Part * find( const PartVector & , const std::string & );

/** \brief  Query containment within properly ordered PartVector */
bool contain( const ConstPartVector & , const Part & );
bool contain( const PartVector & , const Part & );

template<class Iterator>
inline
bool contains_ordinal( Iterator beg, Iterator end, unsigned part_ordinal )
{
  for(Iterator i=beg; i!=end; ++i) {
    if (*i == part_ordinal) return true;
  }

  return false;
}

inline
bool contains_ordinal_part( PartVector::const_iterator beg, PartVector::const_iterator end, unsigned part_ordinal )
{
  for(PartVector::const_iterator i=beg; i!=end; ++i) {
    if ((*i)->mesh_meta_data_ordinal() == part_ordinal) return true;
  }

  return false;
}

/** \brief  Query containment for two properly ordered PartVector */
bool contain( const PartVector & , const PartVector & );

/** \brief  Query cardinality of intersection of two PartVectors */
size_t intersect( const PartVector & , const PartVector & );

/** \brief  Generate the intersection of two PartVectors */
size_t intersect( const PartVector & , const PartVector & , PartVector & );

/** \brief  Query if two parts intersect; i.e.,
 *          if one is a subset of the other or they share a common subset
 */
bool intersect( const Part & , const Part & );

//----------------------------------------------------------------------
/** \brief  Print a part's information including supersets, subsets, and
 *          intersection.  Each line starts with the given leader string.
 */
std::ostream & print( std::ostream & , const char * const , const Part & );


/** \} */


} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

