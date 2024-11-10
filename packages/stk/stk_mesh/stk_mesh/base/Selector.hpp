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

#ifndef stk_mesh_Selector_hpp
#define stk_mesh_Selector_hpp

#include <stddef.h>                     // for NULL
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <stk_mesh/base/Part.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowAssert
#include <vector>                       // for vector, operator!=, etc
#include <utility>                       // for vector, operator!=, etc
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class FieldBase; } }

namespace stk { namespace mesh {

// Specify what can be in the expression tree of a selector
struct SelectorNodeType
{
  enum node_type {
    PART,
    PART_UNION,
    PART_INTERSECTION,
    UNION,
    INTERSECTION,
    DIFFERENCE,
    COMPLEMENT,
    FIELD
  };
};

namespace impl {

// A node in the expression tree of a selector
struct SelectorNode
{
  SelectorNode(Part const* arg_part = nullptr)
   : m_type(SelectorNodeType::PART),
     m_partOrds(),
     m_subsetPartOrds(),
     m_field_ptr(nullptr),
     left_offset(0)
  {
    m_partOrd = (arg_part==nullptr) ? InvalidPartOrdinal : arg_part->mesh_meta_data_ordinal();
    m_partIsRanked = arg_part==nullptr ? false : arg_part->primary_entity_rank()!=stk::topology::INVALID_RANK;
  }

  SelectorNode(FieldBase const* arg_field)
   : m_type(SelectorNodeType::FIELD),
     m_partOrd(InvalidPartOrdinal),
     m_partOrds(),
     m_subsetPartOrds(),
     m_partIsRanked(false),
     left_offset(0)
  {
    m_field_ptr = arg_field;
  }

  SelectorNode(SelectorNodeType::node_type selectorNodeType, const OrdinalVector& partOrdinals, const OrdinalVector& subsetPartOrdinals)
   : m_type(selectorNodeType),
     m_partOrd(InvalidPartOrdinal),
     m_partOrds(partOrdinals),
     m_subsetPartOrds(subsetPartOrdinals),
     m_partIsRanked(false),
     m_field_ptr(nullptr),
     left_offset(0)
  {
  }

  //
  //  Data structure design:
  //
  //  Selector is orgainzed such that the right most (back) entry of the selector is the root node.
  //
  //  Examples: 
  //
  //    left_offset    = L
  //    right_offset   = R
  //    uranary_offset = U 
  //
  //  Unary selector:  !(S1)
  //
  //     +--U--+<--(Root Node)
  //     |     |
  //     V     |
  //    S1     ! 
  //
  //  Complex selector: (S1 & S2) & (!S3)
  //
  //               +------L------+<--(Root Node)
  //               |             |
  //               V        +--R-+
  //     +----L----+        |    |
  //     |         |        V    |
  //     |     +-R-+    +-U-+    |
  //     |     |   |    |   |    |
  //     V     V   |    V   |    |
  //    S1    S2   &   S3   !    &
  //

  // Either leaf (partOrd) or unary (no data) or binary (offset from current pos to rhs)
  enum { right_offset = 1 };
  enum { unary_offset = 1 };

  SelectorNode const* lhs() const
  {
    STK_ThrowAssert(m_type == SelectorNodeType::UNION || m_type == SelectorNodeType::INTERSECTION || m_type == SelectorNodeType::DIFFERENCE);
    return this - left_offset;
  }

  SelectorNode const* rhs() const
  {
    STK_ThrowAssert(m_type == SelectorNodeType::UNION || m_type == SelectorNodeType::INTERSECTION || m_type == SelectorNodeType::DIFFERENCE);
    return this - right_offset;
  }

  SelectorNode const* unary() const
  {
    STK_ThrowAssert(m_type == SelectorNodeType::COMPLEMENT);
    return this - unary_offset;
  }

  PartOrdinal part() const
  {
    STK_ThrowAssert(m_type == SelectorNodeType::PART);
    return m_partOrd;
  }

  bool part_is_ranked() const
  {
    return m_partIsRanked;
  }

  FieldBase const* field() const
  {
    STK_ThrowAssert(m_type == SelectorNodeType::FIELD);
    return m_field_ptr;
  }

  SelectorNodeType::node_type node_type() const
  {
      return m_type;
  }

  bool operator==(SelectorNode const& arg_rhs) const;

  SelectorNodeType::node_type  m_type;
  Ordinal                      m_partOrd;
  OrdinalVector                m_partOrds;
  OrdinalVector                m_subsetPartOrds;
  bool                         m_partIsRanked;
  FieldBase const*             m_field_ptr;
  unsigned left_offset; // for binary op
};

} // namespace impl

/**
 * Selects subsets of the mesh. Allows for creating set expressions from
 * parts and set operators (union |, intersection &, difference -, complement !).
 */
class Selector {
public:

  Selector()
    : m_expr(1) // default Selector is null part (selects nothing)
    , m_meta(nullptr)
  {}

  /** \brief  A part that is required */
  Selector(const Part & part)
    : m_expr(1, impl::SelectorNode(&part))
    , m_meta(((m_expr[0].part()==InvalidPartOrdinal) ? nullptr : &part.mesh_meta_data()))
  {
    STK_ThrowAssertMsg(m_meta!=nullptr, "constructing Selector with bad part reference, probably from dereferencing a null part pointer.");
  }

  /** \brief  Bucket has field */
  Selector(const FieldBase & field);

  Selector(SelectorNodeType::node_type selectorNodeType,
           const MetaData* metaPtr,
           const OrdinalVector& partOrdinals,
           const OrdinalVector& subsetPartOrdinals);

  bool operator == (const Selector & rhs) const
  { return m_expr == rhs.m_expr; }

  bool operator != (const Selector & rhs) const
  { return m_expr != rhs.m_expr; }

  /** \brief  Intersection: this = this INTERSECT ( expression ) */
  Selector & operator &= ( const Selector & selector)
  { return add_binary_op(SelectorNodeType::INTERSECTION, selector); }

  Selector & operator &= ( const Part& part)
  { return add_binary_op(SelectorNodeType::INTERSECTION, Selector(part)); }

  Selector & operator |= ( const Part& part )
  { return operator|=(Selector(part)); }

  /** \brief  Union: this = this UNION ( expression ) */
  Selector & operator |= ( const Selector & selector)
  {
    if (is_null()) {
      m_expr = selector.m_expr;
      m_meta = selector.m_meta;
    }
    else {
      add_binary_op(SelectorNodeType::UNION, selector);
    }
    return *this;
  }

  /** \brief  Difference: this = this - ( expression ) */
  Selector & operator -= ( const Selector & selector)
  { return add_binary_op(SelectorNodeType::DIFFERENCE, selector); }

  bool operator<(const Selector& rhs) const;

  bool operator<=(const Selector& rhs) const {
    return *this < rhs || *this == rhs;
  }

  bool operator>(const Selector& rhs) const {
    return rhs < *this;
  }

  bool operator>=(const Selector& rhs) const {
    return *this > rhs || *this == rhs;
  }

  /** \brief  Complement: this = !(this)
   * Complements this selector in-place
   * */
  Selector & complement()
  {
    impl::SelectorNode root;
    root.m_type = SelectorNodeType::COMPLEMENT;

    m_expr.push_back(root);

    return *this;
  }

  /** \brief Complement:  return !(this) */
  Selector operator ! () const
  {
    Selector tmp( *this );
    return tmp.complement();
  }

  bool is_all_unions() const;

  template<typename PARTVECTOR>
  void get_parts(PARTVECTOR& parts) const;

  BucketVector const& get_buckets(EntityRank entity_rank) const;

  bool is_empty(EntityRank entity_rank) const;

  /** \brief  Is this part a member of the
   *          set defined by the selector expression.
   */
  bool operator()( const Part & part ) const;

  /** \brief  Is this part a member of the
   *          set defined by the selector expression.
   */
  bool operator()( const Part * part ) const;

  /** \brief  Is this bucket a subset of the
   *          set defined by the selector expression.
   */
  bool operator()( const Bucket & bucket ) const;

  /** \brief  Is this bucket a subset of the
   *          set defined by the selector expression.
   */
  bool operator()( const Bucket * bucket ) const;


  bool operator()(const PartVector& parts) const;


  stk::mesh::Selector clone_for_different_mesh(const stk::mesh::MetaData &differentMeta) const;

  friend std::ostream & operator << ( std::ostream & out, const Selector & selector);

private:
  bool select_part_impl(Part const& part, impl::SelectorNode const* root) const;
  bool select_bucket_impl(Bucket const& bucket, impl::SelectorNode const* root) const;

  const BulkData* find_mesh() const;

  bool is_null() const {
    if(m_expr.size() > 1) return false;
    if(m_expr.back().m_type == SelectorNodeType::PART  && m_expr.back().part()  == InvalidPartOrdinal) {
      return true;
    } else if(m_expr.back().m_type == SelectorNodeType::FIELD && m_expr.back().field() == nullptr) {
      return true;
    }
    return false;
  }

  Selector& add_binary_op(SelectorNodeType::node_type type, const Selector& rhs);

  std::vector<impl::SelectorNode> m_expr;
  mutable const MetaData* m_meta;
};

inline
Selector operator & ( const Part & A , const Part & B )
{
  Selector S( A );
  S &= Selector( B );
  return S;
}

inline
Selector operator & ( const Part & A , const Selector & B )
{
  Selector S( A );
  S &= B;
  return S;
}

inline
Selector operator & ( const Selector & A, const Part & B )
{
  Selector S( A );
  S &= Selector(B);
  return S;
}

inline
Selector operator & ( const Selector & A, const Selector & B )
{
  Selector S( A );
  S &= Selector(B);
  return S;
}

inline
Selector operator | ( const Part & A , const Part & B )
{
  Selector S( A );
  S |= Selector( B );
  return S;
}

inline
Selector operator | ( const Part & A , const Selector & B )
{
  Selector S( A );
  S |= B;
  return S;
}

inline
Selector operator | ( const Selector & A, const Part & B  )
{
  Selector S( A );
  S |= Selector(B);
  return S;
}

inline
Selector operator | ( const Selector & A, const Selector & B  )
{
  Selector S( A );
  S |= Selector(B);
  return S;
}

inline
Selector operator - ( const Part & A , const Part & B )
{
  Selector S( A );
  S -= Selector( B );
  return S;
}

inline
Selector operator - ( const Part & A , const Selector & B )
{
  Selector S( A );
  S -= B;
  return S;
}

inline
Selector operator - ( const Selector & A, const Part & B  )
{
  Selector S( A );
  S -= Selector(B);
  return S;
}

inline
Selector operator - ( const Selector & A, const Selector & B  )
{
  Selector S( A );
  S -= Selector(B);
  return S;
}

inline
Selector operator ! ( const Part & A )
{
  Selector S(A);
  return S.complement();
}

/** \brief .
 * \relates Selector
 * */
template <typename PartVectorType>
Selector selectUnion( const PartVectorType & union_part_vector );

/** \brief .
 * \relates Selector
 * */
Selector selectIntersection( const PartVector& intersection_part_vector );
Selector selectIntersection( const ConstPartVector& intersection_part_vector );

/** \brief Return a selector for the union of the parts where field exists.
 * \relates Selector
 * */
Selector selectField( const FieldBase& field );

/** \brief Is lhs a subset of rhs, only works for simple selectors (parts and unions)
 */
bool is_subset(Selector const& lhs, Selector const& rhs);
bool is_subset(Selector const& lhs, Selector const& rhs, PartVector& scratchLhs, PartVector& scratchRhs);
bool is_subset(const PartVector& lhsParts, const PartVector& rhsParts);
bool is_subset(const OrdinalVector& lhsParts, const OrdinalVector& rhsParts);

/** \} */

}} // namespace stk::mesh

#endif // stk_mesh_Selector_hpp

