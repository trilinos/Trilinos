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

#ifndef stk_mesh_Selector_hpp
#define stk_mesh_Selector_hpp

#include <stddef.h>                     // for NULL
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowAssert
#include <vector>                       // for vector, operator!=, etc
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class Part; } }




namespace stk { namespace mesh {

// Specify what can be in the expression tree of a selector
struct SelectorNodeType
{
  enum node_type {
    PART,
    UNION,
    INTERSECTION,
    DIFFERENCE,
    COMPLEMENT
  };
};

namespace impl {

// A node in the expression tree of a selector
struct SelectorNode
{
  SelectorNode(Part const* arg_part = NULL) : m_type(SelectorNodeType::PART)
  {
    m_value.part_ptr = arg_part;
  }

  // Either leaf (part_ptr) or unary (no data) or binary (offset from current pos to rhs)
  union value_type
  {
    enum { left_offset = 1 };
    enum { unary_offset = 1 };

    Part const* part_ptr;
    unsigned right_offset; // for binary op
    // no storage required for unary op
  };

  SelectorNode const* lhs() const
  {
    ThrowAssert(m_type == SelectorNodeType::UNION || m_type == SelectorNodeType::INTERSECTION || m_type == SelectorNodeType::DIFFERENCE);
    return this + m_value.left_offset;
  }

  SelectorNode const* rhs() const
  {
    ThrowAssert(m_type == SelectorNodeType::UNION || m_type == SelectorNodeType::INTERSECTION || m_type == SelectorNodeType::DIFFERENCE);
    return this + m_value.right_offset;
  }

  SelectorNode const* unary() const
  {
    ThrowAssert(m_type == SelectorNodeType::COMPLEMENT);
    return this + m_value.unary_offset;
  }

  Part const* part() const
  {
    ThrowAssert(m_type == SelectorNodeType::PART);
    return m_value.part_ptr;
  }

  SelectorNodeType::node_type node_type() const
  {
      return m_type;
  }

  bool operator==(SelectorNode const& arg_rhs) const
  {
    if (m_type != arg_rhs.m_type) {
      return false;
    }
    else if (m_type == SelectorNodeType::COMPLEMENT) {
      // there's no rhs for a SelectorNode of type Complement
      return true;
    }
    else if (m_type == SelectorNodeType::PART) {
      return m_value.part_ptr == arg_rhs.m_value.part_ptr;
    }
    else {
      return m_value.right_offset == arg_rhs.m_value.right_offset;
    }
  }

  SelectorNodeType::node_type  m_type;
  value_type                   m_value;
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
  {}

  /** \brief  A part that is required */
  Selector(const Part & part)
    : m_expr(1, impl::SelectorNode(&part))
  {}

  bool operator == (const Selector & rhs) const
  { return m_expr == rhs.m_expr; }

  bool operator != (const Selector & rhs) const
  { return m_expr != rhs.m_expr; }

  /** \brief  Intersection: this = this INTERSECT ( expression ) */
  Selector & operator &= ( const Selector & selector)
  { return add_binary_op(SelectorNodeType::INTERSECTION, selector); }

// Remove once Fmwk::MeshPart no longer exists
#ifdef SIERRA_MIGRATION
  Selector & operator &= ( const Part& part)
  { return add_binary_op(SelectorNodeType::INTERSECTION, Selector(part)); }

  Selector & operator |= ( const Part& part )
  { return operator|=(Selector(part)); }
#endif

  /** \brief  Union: this = this UNION ( expression ) */
  Selector & operator |= ( const Selector & selector)
  {
    if (is_null()) {
      m_expr = selector.m_expr;
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

    m_expr.insert(m_expr.begin(), root);
    return *this;
  }

  /** \brief Complement:  return !(this) */
  Selector operator ! () const
  {
    Selector tmp( *this );
    return tmp.complement();
  }

  bool is_all_unions() const;

  void get_parts(PartVector& parts) const;

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

  friend std::ostream & operator << ( std::ostream & out, const Selector & selector);

private:

  BulkData* find_mesh() const;

  bool is_null() const {
    return m_expr.size() == 1 && m_expr[0].m_type == SelectorNodeType::PART && m_expr[0].m_value.part_ptr == NULL;
  }

  Selector& add_binary_op(SelectorNodeType::node_type type, const Selector& rhs)
  {
    impl::SelectorNode root;
    root.m_type = type;
    root.m_value.right_offset = 1 + m_expr.size();

    m_expr.insert(m_expr.begin(), root);
    m_expr.insert(m_expr.end(), rhs.m_expr.begin(), rhs.m_expr.end());

    return *this;
  }

  std::vector<impl::SelectorNode> m_expr;
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

/** \} */

}} // namespace stk::mesh

#endif // stk_mesh_Selector_hpp

