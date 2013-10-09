/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <sstream>
#include <iostream>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {

namespace {

struct print_selector_impl :
  boost::static_visitor<std::ostream &>
{
  std::ostream & m_out;

  print_selector_impl(std::ostream & out)
    : m_out(out)
  {}

  std::ostream& operator() (const Part * p) const
  {
    if (p != NULL) {
      m_out <<  p->name();
    } else {
      m_out << "NOTHING";
    }
    return m_out;
  }

  template <typename Operation>
  std::ostream& operator() (impl::binary_op_<Operation> const& expr) const
  {
    m_out << "(";
    boost::apply_visitor( print_selector_impl(m_out), expr.m_lhs );
    m_out << Operation();
    boost::apply_visitor( print_selector_impl(m_out), expr.m_rhs );
    m_out << ")";
    return m_out;
  }

  std::ostream& operator() (impl::complement_ const& expr) const
  {
    m_out << "!(";
    boost::apply_visitor( print_selector_impl(m_out), expr.m_selector );
    m_out << ")";
    return m_out;
  }
};


struct select_bucket_impl : public boost::static_visitor<bool>
{
  Bucket const& m_bucket;

  select_bucket_impl(Bucket const& bucket)
    : m_bucket(bucket)
  {}

  bool operator() (const Part * p) const
  { return (p != NULL)? has_superset(m_bucket, *p) : false; }

  bool operator() (impl::binary_op_<impl::union_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) || boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::binary_op_<impl::intersect_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::binary_op_<impl::difference_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && !boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::complement_ const& set) const
  { return !boost::apply_visitor(*this, set.m_selector); }
};

struct select_part_impl : public boost::static_visitor<bool>
{
  const Part * m_part;

  select_part_impl(Part const& part)
    : m_part(&part)
  {}

  bool operator() (const Part * p) const
  {
    if (p == NULL) {
      return false;
    }
    return p->contains(*m_part);
  }

  bool operator() (impl::binary_op_<impl::union_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) || boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::binary_op_<impl::intersect_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::binary_op_<impl::difference_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && !boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::complement_ const& set) const
  { return !boost::apply_visitor(*this, set.m_selector); }
};

struct is_all_union_visitor : public boost::static_visitor<bool>
{
  bool operator() (const Part * p) const
  { return p != NULL; }

  bool operator() (impl::binary_op_<impl::union_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::binary_op_<impl::intersect_> const& set) const
  { return false; }

  bool operator() (impl::binary_op_<impl::difference_> const& set) const
  { return false; }

  bool operator() (impl::complement_ const& set) const
  { return false; }
};

struct gather_parts : public boost::static_visitor<void>
{
  PartVector& m_parts;

  gather_parts(PartVector & parts)
    : m_parts(parts)
  {}

  void operator() (const Part * p) const
  {
    if (p != NULL) {
      m_parts.push_back(const_cast<Part*>(p));
    }
  }

  void operator() (impl::binary_op_<impl::union_> const& set) const
  { boost::apply_visitor(*this, set.m_lhs);
    boost::apply_visitor(*this, set.m_rhs);
  }

  void operator() (impl::binary_op_<impl::intersect_> const& set) const
  {
    // HACK: Only first part (picks up Context part)
    boost::apply_visitor(*this, set.m_lhs);
  }

  void operator() (impl::binary_op_<impl::difference_> const& set) const
  {
    ThrowRequireMsg(false, "Cannot get_parts from a selector with differences");
  }

  void operator() (impl::complement_ const& set) const
  {
    ThrowRequireMsg(false, "Cannot get_parts from a selector with complements");
  }
};

struct select_part_vector_impl : public boost::static_visitor<bool>
{
  const PartVector& m_parts;

  select_part_vector_impl(PartVector const& parts)
    : m_parts(parts)
  {}

  bool operator() (const Part * p) const
  {
    if (p == NULL) {
      return false;
    }

    for (size_t i = 0, ie = m_parts.size(); i < ie; ++i) {
      if (m_parts[i] == p) {
        return true;
      }
    }

    return false;
  }

  bool operator() (impl::binary_op_<impl::union_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) || boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::binary_op_<impl::intersect_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::binary_op_<impl::difference_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && !boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (impl::complement_ const& set) const
  { return !boost::apply_visitor(*this, set.m_selector); }
};

struct compare_less_impl : public boost::static_visitor<bool>
{
  bool operator()(const Part* a, const Part* b) const
  {
    if (a != NULL && b != NULL) {
      return a->mesh_meta_data_ordinal() < b->mesh_meta_data_ordinal();
    }
    else if (a == NULL && b != NULL) {
      return false;
    }
    else if (a != NULL && b == NULL) {
      return true;
    }
    return false;
  }

  template <typename Expr>
  bool operator()(const Part*, const Expr&) const
  {
    return true;
  }

  template <typename Expr>
  bool operator()(const Expr&, const Part*) const
  {
    return false;
  }

  bool operator()(impl::binary_op_<impl::union_> const& a, impl::binary_op_<impl::union_> const& b) const
  {
    return boost::apply_visitor(*this, a.m_lhs, b.m_lhs) && boost::apply_visitor(*this, a.m_rhs, b.m_rhs);
  }

  bool operator()(impl::binary_op_<impl::intersect_> const& a, impl::binary_op_<impl::intersect_> const& b) const
  {
    return boost::apply_visitor(*this, a.m_lhs, b.m_lhs) && boost::apply_visitor(*this, a.m_rhs, b.m_rhs);
  }

  bool operator()(impl::binary_op_<impl::difference_> const& a, impl::binary_op_<impl::difference_> const& b) const
  {
    return boost::apply_visitor(*this, a.m_lhs, b.m_lhs) && boost::apply_visitor(*this, a.m_rhs, b.m_rhs);
  }

  bool operator()(impl::complement_ const& a, impl::complement_ const& b) const
  {
    return boost::apply_visitor(*this, a.m_selector, b.m_selector);
  }

  //

  bool operator()(impl::binary_op_<impl::union_> const&, impl::binary_op_<impl::intersect_> const&) const
  {
    return true;
  }

  bool operator()(impl::binary_op_<impl::union_> const&, impl::binary_op_<impl::difference_> const&) const
  {
    return true;
  }

  bool operator()(impl::binary_op_<impl::union_> const&, impl::complement_ const&) const
  {
    return true;
  }


  bool operator()(impl::binary_op_<impl::intersect_> const&, impl::binary_op_<impl::union_> const&) const
  {
    return false;
  }

  bool operator()(impl::binary_op_<impl::difference_> const&, impl::binary_op_<impl::union_> const&) const
  {
    return false;
  }

  bool operator()(impl::complement_ const&, impl::binary_op_<impl::union_> const&) const
  {
    return false;
  }

  //

  bool operator()(impl::binary_op_<impl::intersect_> const&, impl::binary_op_<impl::difference_> const&) const
  {
    return true;
  }

  bool operator()(impl::binary_op_<impl::difference_> const&, impl::binary_op_<impl::intersect_> const&) const
  {
    return false;
  }

  bool operator()(impl::binary_op_<impl::intersect_> const&, impl::complement_ const&) const
  {
    return true;
  }

  bool operator()(impl::complement_ const&, impl::binary_op_<impl::intersect_> const&) const
  {
    return false;
  }

  bool operator()(impl::binary_op_<impl::difference_> const&, impl::complement_ const&) const
  {
    return true;
  }

  bool operator()(impl::complement_ const&, impl::binary_op_<impl::difference_> const&) const
  {
    return false;
  }
};

}// namespace


std::ostream & operator<<( std::ostream & out, const Selector & selector)
{
  return boost::apply_visitor(print_selector_impl(out),selector.m_selector);
}


bool Selector::operator()( const Part & part ) const
{
  return boost::apply_visitor(select_part_impl(part),m_selector);
}

bool Selector::operator()( const Part * part ) const
{
  if (part == NULL) return false;
  return boost::apply_visitor(select_part_impl(*part),m_selector);
}

bool Selector::operator()( const Bucket & bucket ) const
{
  return boost::apply_visitor(select_bucket_impl(bucket),m_selector);
}

bool Selector::operator()( const Bucket * bucket ) const
{
  if (bucket != NULL)
    return boost::apply_visitor(select_bucket_impl(*bucket),m_selector);
  return false;
}

bool Selector::operator()(const PartVector& parts) const
{
  return boost::apply_visitor(select_part_vector_impl(parts),m_selector);
}

bool Selector::operator<(const Selector& rhs) const
{
  return boost::apply_visitor(compare_less_impl(), m_selector, rhs.m_selector);
}

void Selector::get_parts(PartVector& parts) const
{
  boost::apply_visitor(gather_parts(parts), m_selector);
}

bool Selector::is_all_unions() const
{
  return boost::apply_visitor(is_all_union_visitor(), m_selector);
}

Selector operator & ( const Part & A , const Part & B )
{
  Selector S( A );
  S &= Selector( B );
  return S;
}

Selector operator & ( const Part & A , const Selector & B )
{
  Selector S( A );
  S &= B;
  return S;
}

Selector operator & ( const Selector & A, const Part & B )
{
  Selector S( A );
  S &= Selector(B);
  return S;
}

Selector operator & ( const Selector & A, const Selector & B )
{
  Selector S( A );
  S &= Selector(B);
  return S;
}

Selector operator | ( const Part & A , const Part & B )
{
  Selector S( A );
  S |= Selector( B );
  return S;
}


Selector operator | ( const Part & A , const Selector & B )
{
  Selector S( A );
  S |= B;
  return S;
}

Selector operator | ( const Selector & A, const Part & B  )
{
  Selector S( A );
  S |= Selector(B);
  return S;
}

Selector operator | ( const Selector & A, const Selector & B  )
{
  Selector S( A );
  S |= Selector(B);
  return S;
}


Selector operator - ( const Part & A , const Part & B )
{
  Selector S( A );
  S -= Selector( B );
  return S;
}


Selector operator - ( const Part & A , const Selector & B )
{
  Selector S( A );
  S -= B;
  return S;
}

Selector operator - ( const Selector & A, const Part & B  )
{
  Selector S( A );
  S -= Selector(B);
  return S;
}

Selector operator - ( const Selector & A, const Selector & B  )
{
  Selector S( A );
  S -= Selector(B);
  return S;
}


Selector operator ! ( const Part & A )
{
  Selector S(A);
  return S.complement();
}


Selector selectUnion( const PartVector& union_part_vector )
{
  Selector selector;
  if (union_part_vector.size() > 0) {
    selector = *union_part_vector[0];
    for (unsigned i = 1 ; i < union_part_vector.size() ; ++i) {
      selector |= *union_part_vector[i];
    }
  }
  return selector;
}


Selector selectIntersection( const PartVector& intersection_part_vector )
{
  Selector selector;
  if (intersection_part_vector.size() > 0) {
    selector = *intersection_part_vector[0];
    for (unsigned i = 1 ; i < intersection_part_vector.size() ; ++i) {
      selector &= *intersection_part_vector[i];
    }
  }
  return selector;
}

Selector selectField( const FieldBase& field )
{
  Selector selector;
  const FieldRestrictionVector& sel_rvec = field.restrictions();
  for(size_t i=0; i<sel_rvec.size(); ++i) {
    selector |= sel_rvec[i].selector();
  }

  return selector;
}

bool is_subset(Selector const& lhs, Selector const& rhs)
{
  if (lhs.is_all_unions() && rhs.is_all_unions()) {
    PartVector lhs_parts, rhs_parts;
    lhs.get_parts(lhs_parts);
    rhs.get_parts(rhs_parts);
    for (size_t l = 0, le = lhs_parts.size(); l < le; ++l) {
      Part const& lhs_part = *lhs_parts[l];
      bool found = false;
      for (size_t r = 0, re = rhs_parts.size(); !found && r < re; ++r) {
        Part const& rhs_part = *rhs_parts[r];
        found = rhs_part.contains(lhs_part);
      }
      if (!found) {
        return false;
      }
    }
    return true;
  }
  else {
    // If either selector has complements or intersections, it becomes
    // much harder to determine if one is a subset of the other
    return false;
  }
}

} // namespace mesh
} // namespace stk
