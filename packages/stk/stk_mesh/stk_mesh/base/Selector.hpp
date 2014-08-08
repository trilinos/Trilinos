/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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

  size_t size(EntityRank entity_rank) const;

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
Selector selectUnion( const PartVector& union_part_vector );

/** \brief .
 * \relates Selector
 * */
Selector selectIntersection( const PartVector& intersection_part_vector );

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

