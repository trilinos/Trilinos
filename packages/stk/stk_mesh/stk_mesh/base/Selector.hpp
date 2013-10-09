/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_Selector_hpp
#define stk_mesh_Selector_hpp

#include <stk_mesh/base/Types.hpp>

#include <boost/type_traits.hpp>
#include <boost/variant.hpp>

#include <iosfwd>

namespace stk { namespace mesh {

namespace impl {

struct union_
{
  friend inline std::ostream& operator<<(std::ostream& out, union_)
  { return out << " | "; }
};

struct intersect_
{
  friend inline std::ostream& operator<<(std::ostream& out, intersect_)
  { return out << " & "; }
};

struct difference_
{
  friend inline std::ostream& operator<<(std::ostream& out, difference_)
  { return out << " - "; }
};

template <typename Op>
struct binary_op_;

struct complement_;


typedef
  boost::variant<
     const Part *
   , boost::recursive_wrapper< binary_op_<union_> >
   , boost::recursive_wrapper< binary_op_<intersect_> >
   , boost::recursive_wrapper< binary_op_<difference_> >
   , boost::recursive_wrapper< complement_ >
  > selector_impl_;

template <typename Op>
struct binary_op_
{
  typedef Op op_type;

  selector_impl_ m_lhs;
  selector_impl_ m_rhs;

  binary_op_()
    : m_lhs()
    , m_rhs()
  {}

  binary_op_(selector_impl_ const& lhs, selector_impl_ const& rhs)
    : m_lhs(lhs)
    , m_rhs(rhs)
  {}

  bool operator==(const binary_op_<Op> & b) const
  { return (m_lhs == b.m_lhs) && (m_rhs == b.m_rhs); }

  bool operator!=(const binary_op_<Op> & b) const
  { return !(*this == b); }
};

struct complement_
{
  selector_impl_ m_selector;

  complement_()
    : m_selector()
  {}

  complement_(selector_impl_ const& select)
    : m_selector(select)
  {}

  bool operator==(const complement_ & b) const
  { return m_selector == b.m_selector; }

  bool operator!=(const complement_ & b) const
  { return !(*this == b); }
};

struct is_empty_selector : public boost::static_visitor<bool>
{
  bool operator() (const Part * p) const
  {
    return (p == NULL);
  }

  template <typename Op>
  bool operator() (Op const& expr) const
  { return false; }
};

} // namespace impl


class Selector {
public:

  Selector()
    : m_selector((Part*)NULL)
  {}

  /** \brief  A part that is required */
  Selector(const Part & part)
    : m_selector(&part)
  {}

  bool operator == (const Selector & rhs) const
  { return m_selector == rhs.m_selector; }

  bool operator != (const Selector & rhs) const
  { return !(m_selector == rhs.m_selector); }


  /** \brief  Intersection: this = this INTERSECT ( expression ) */
  Selector & operator &= ( const Selector & selector)
  {
    m_selector = impl::binary_op_<impl::intersect_>(m_selector,selector.m_selector);
    return *this;
  }

  /** \brief  Union: this = this UNION ( expression ) */
  Selector & operator |= ( const Selector & selector)
  {
    if (boost::apply_visitor(impl::is_empty_selector(),m_selector)) {
      m_selector = selector.m_selector;
    }
    else {
      m_selector = impl::binary_op_<impl::union_>(m_selector,selector.m_selector);
    }
    return *this;
  }

  /** \brief  Difference: this = this - ( expression ) */
  Selector & operator -= ( const Selector & selector)
  {
    m_selector = impl::binary_op_<impl::difference_>(m_selector,selector.m_selector);
    return *this;
  }

  bool operator<(const Selector& rhs) const;

  /** \brief  Complement: this = !(this)
   * Postcondition:  this is a compound expression
   * */
  Selector & complement()
  {
    m_selector = impl::complement_(m_selector);
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

  impl::selector_impl_ m_selector;

};


std::ostream & operator<<( std::ostream & out, const Selector & selector);

/** \brief .
 * \relates Selector
 * */
Selector operator & ( const Part & A , const Part & B );
Selector operator & ( const Part & A , const Selector & B );
Selector operator & ( const Selector & A, const Part & B );
Selector operator & ( const Selector & A, const Selector & B );

/** \brief .
 * \relates Selector
 * */
Selector operator | ( const Part & A , const Part & B );
Selector operator | ( const Part & A , const Selector & B );
Selector operator | ( const Selector & A, const Part & B  );
Selector operator | ( const Selector & A , const Selector & B );

/** \brief .
 * \relates Selector
 * */
Selector operator - ( const Part & A , const Part & B );
Selector operator - ( const Part & A , const Selector & B );
Selector operator - ( const Selector & A, const Part & B  );
Selector operator - ( const Selector & A , const Selector & B );

/** \brief .
 * \relates Selector
 * */
Selector operator ! ( const Part & A );


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

bool is_subset(Selector const& lhs, Selector const& rhs);

/** \} */

}} // namespace stk::mesh

#endif // stk_mesh_Selector_hpp

