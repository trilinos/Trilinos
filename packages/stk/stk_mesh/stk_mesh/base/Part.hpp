// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include <stk_mesh/base/Types.hpp>      // for PartVector, OrdinalVector, etc
#include <stk_topology/topology.hpp>    // for topology
#include <stk_util/util/CSet.hpp>
#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for int64_t
#include <iosfwd>                       // for ostream
#include <string>                       // for string, basic_string
#include <vector>                       // for vector, vector<>::iterator
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { namespace impl { class PartRepository; } } }


//----------------------------------------------------------------------

namespace stk {
namespace mesh {

namespace impl {
stk::CSet & get_attributes(stk::mesh::Part & part);
}

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
  MetaData & mesh_meta_data() const { return *m_mesh_meta_data; }

  BulkData & mesh_bulk_data() const;

  /** \brief  The primary entity type for this part.
   *
   *   For example, the primary purpose of an Element part
   *   is to define a collection of elements.  However, the
   *   nodes of those elements are also members of an element part.
   *   Return InvalidEntityRank if no primary entity type.
   */
  EntityRank primary_entity_rank() const { return m_entity_rank; }

  stk::topology topology() const { return m_topology; }

  /** \brief  Application-defined text name of this part, must be unique within the set of parts owned by a MetaData*/
  const std::string & name() const { return m_name; }

  bool force_no_induce() const { return m_force_no_induce; }

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
  bool entity_membership_is_parallel_consistent() const { return m_entity_membership_is_parallel_consistent; }
  void entity_membership_is_parallel_consistent(bool trueOrFalse) { m_entity_membership_is_parallel_consistent = trueOrFalse; }

  int64_t id() const { return m_id; }

  /** \brief  Internally generated ordinal of this part that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_ordinal; }

  /** \brief  Parts that are supersets of this part. */
  const PartVector & supersets() const { return m_supersets; }

  /** \brief  Parts that are subsets of this part. */
  const PartVector & subsets() const { return m_subsets; }

  /** \brief  Check if argument is subset of this */
  bool contains(const Part& part) const;

  /** \brief  Equality comparison */
  bool operator == ( const Part & rhs ) const { return this == & rhs ; }

  /** \brief  Inequality comparison */
  bool operator != ( const Part & rhs ) const { return this != & rhs ; }

  /** \brief  Query attribute that has been attached to this part */
  template<class A>
  const A * attribute() const { return m_attribute.template get<A>(); }

private:

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns the PartRepository which created this part.
   */
  MetaData & meta_data() const { return *m_mesh_meta_data; }

  void set_id(int64_t inputId) { m_id = inputId; }

  CSet & get_attributes() { return m_attribute; }

  void set_name(const std::string& newName);
  void set_primary_entity_rank( EntityRank entity_rank );

  void set_topology( stk::topology topo )
  {
    if ( topo == stk::topology::INVALID_TOPOLOGY || topo == m_topology ) return;
  
    STK_ThrowErrorMsgIf( m_topology != stk::topology::INVALID_TOPOLOGY && m_topology != topo,
        "Error set_topology: part " << name()
        << " already defined with " << m_topology
        << " conflicts with  " << topo);
  
    STK_ThrowErrorMsgIf(m_entity_rank != stk::topology::INVALID_RANK && m_entity_rank != topo.rank(),
        "Error set_topology: part " << name()
        << " already defined with " << m_entity_rank
        << " conflicts with  " << topo << " which has rank " << topo.rank());
  
    m_entity_rank = topo.rank();
    m_topology = topo;
  }

  bool add_part_to_subset(Part & part);
  bool add_part_to_superset(Part & part);

  void set_force_no_induce(bool input) { m_force_no_induce = input; }

  template<class T>
  const T * declare_attribute_with_delete( const T *);
  template<class T>
  const T * declare_attribute_no_delete( const T *);
  template<class T>
  bool remove_attribute( const T *);

  MetaData* const   m_mesh_meta_data;
  std::string m_name;
  EntityRank        m_entity_rank;
  stk::topology     m_topology;
  int64_t           m_id;
  CSet              m_attribute;
  PartVector        m_subsets;
  bool              m_subsetsEmpty;
  PartVector        m_supersets;
  const unsigned    m_ordinal;
  bool              m_force_no_induce;
  bool              m_entity_membership_is_parallel_consistent;

  /* \brief  A part is owned by a PartRepository, as such only the owning
   *         PartRepository can create, delete, or modify a part.
   *         The owner-modifies rule is enforced by the implementation being
   *         a private data object on the Part and the PartRepository is a
   *         friend.
   */
  friend class ::stk::mesh::impl::PartRepository;
  friend class ::stk::mesh::MetaData;
  friend CSet & impl::get_attributes(stk::mesh::Part & part);

#ifndef DOXYGEN_COMPILE

  /** Construct a subset part within a given mesh.
   *  Is used internally by the two 'declare_part' methods on PartRepository.
   */
  Part(MetaData * arg_meta_data,
       const std::string & arg_name,
       EntityRank arg_rank,
       size_t arg_ordinal,
       bool arg_force_no_induce = false)
    : m_mesh_meta_data(arg_meta_data),
      m_name(arg_name),
      m_entity_rank(arg_rank),
      m_topology(stk::topology::INVALID_TOPOLOGY),
      m_id(Part::INVALID_ID),
      m_attribute(),
      m_subsets(),
      m_subsetsEmpty(true),
      m_supersets(),
      m_ordinal(arg_ordinal),
      m_force_no_induce(arg_force_no_induce),
      m_entity_membership_is_parallel_consistent(true)
  { }

  ~Part() {}
  Part();
  Part( const Part & );
  Part & operator = ( const Part & );

#endif /* DOXYGEN_COMPILE */
};

template<class T>
inline
const T *
Part::declare_attribute_with_delete( const T * a )
{ 
  return m_attribute.template insert_with_delete<T>( a );
}

template<class T>
inline
const T *
Part::declare_attribute_no_delete( const T * a )
{
  return m_attribute.template insert_no_delete<T>( a );
}

template<class T>
inline
bool
Part::remove_attribute( const T * a )
{
  return m_attribute.template remove<T>( a );
}

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

  bool operator()(const stk::mesh::Part *lhs, unsigned rhs) const {
    return lhs->mesh_meta_data_ordinal() < rhs;
  }
};

struct PartLessById {

  inline bool operator()( const Part & lhs , const Part & rhs ) const
    { return lhs.id() < rhs.id(); }

  inline bool operator()( const Part & lhs , const Part * rhs ) const
    { return lhs.id() < rhs->id(); }

  inline bool operator()( const Part * lhs , const Part & rhs ) const
    { return lhs->id() < rhs.id(); }

  inline bool operator()( const Part * lhs , const Part * rhs ) const
    { return lhs->id() < rhs->id(); }

  bool operator()(const stk::mesh::Part *lhs, unsigned rhs) const {
    return lhs->id() < rhs;
  }
};

/** \brief  Insert a part into a properly ordered collection of parts.
 *          Returns true if this is a new insertion.
 */
bool insert( ConstPartVector & , const Part & );
bool insert( PartVector & , Part & );
bool contains( const PartVector & v , const Part & part );

void get_part_and_all_subsets(const Part& part, ConstPartVector& part_and_all_subsets);

/** \brief  Remove a part from a properly ordered collection of parts. */
void remove( PartVector & , Part & );

/** \brief  Find a part by name in a collection of parts. */
Part * find( const PartVector & , const std::string & );

/** \brief  Query containment within properly ordered PartVector */
template<class PARTVECTOR>
bool contain( const PARTVECTOR & v , const Part & part )
{
  typename PARTVECTOR::const_iterator e = v.end();
  typename PARTVECTOR::const_iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  return i != e && *i == & part ;
}

template<class Iterator>
inline
bool contains_ordinal( Iterator beg, Iterator end, unsigned part_ordinal )
{
  while(beg!=end) {
    if (*beg++ == part_ordinal) return true;
  }

  return false;
}

inline
bool contains_ordinal(const OrdinalVector& ordinals, Ordinal ord)
{
  OrdinalVector::const_iterator iter = std::lower_bound(ordinals.begin(), ordinals.end(), ord);
  return iter != ordinals.end() && *iter == ord;
}

template<class Iterator>
inline
bool contains_ordinal_part( Iterator beg, Iterator end, unsigned part_ordinal )
{
  for(Iterator i=beg; i!=end; ++i) {
    if ((*i)->mesh_meta_data_ordinal() == part_ordinal) return true;
  }

  return false;
}

/** \brief  Query cardinality of intersection of two PartVectors */
template<class PARTVECTOR>
size_t intersect( const PARTVECTOR & v , const PARTVECTOR & p )
{
  // Both lists must be sorted, assume v.size() > p.size()

  typename PARTVECTOR::const_iterator ev = v.end();
  typename PARTVECTOR::const_iterator iv = v.begin();

  typename PARTVECTOR::const_iterator ep = p.end();
  typename PARTVECTOR::const_iterator ip = p.begin();

  size_t count = 0 ;

  for ( ; ip != ep && iv != ev ; ++ip ) {
    Part * const q = *ip ;
    iv = std::lower_bound( iv , ev , q , PartLess() );
    if ( iv != ev && *iv == q ) { ++count ; }
  }

  return count ;
}

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

