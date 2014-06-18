#ifndef SAMBA_SAMBA_SET_EXPRESSION_OPERATORS_HPP
#define SAMBA_SAMBA_SET_EXPRESSION_OPERATORS_HPP

#include <samba/set_expression/set_expression_impl.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/logical.hpp>

namespace samba {

//*******************************************************************
//intersect
//*******************************************************************

//term & term
template <typename Left, typename Right>
inline typename boost::enable_if<
  boost::mpl::and_<
    detail::convertable_to_entity_part<Left>
   ,detail::convertable_to_entity_part<Right>
  >
 ,set_expression >::type
operator & (Left const& l, Right const& r)
{ return set_expression::intersect_(entity_part(l),(r)); }

//expr & term
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression >::type
operator & (set_expression const& l, T const& r)
{ return set_expression::intersect_(l,entity_part(r)); }

//term & expr
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression >::type
operator & (T const& l, set_expression const& r)
{ return set_expression::intersect_(entity_part(l),r); }

//expr & expr
inline
set_expression
operator & (set_expression const& l, set_expression const& r)
{ return set_expression::intersect_(l,r); }

//expr &= term
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression& >::type
operator &= (set_expression& l, T const& r)
{ return l = set_expression::intersect_(l,entity_part(r)); }

//expr &= expr
inline
set_expression &
operator &= (set_expression& l, set_expression const& r)
{ return l = set_expression::intersect_(l,r); }


//*******************************************************************
//union
//*******************************************************************

//term | term
template <typename Left, typename Right>
inline typename boost::enable_if<
  boost::mpl::and_<
    detail::convertable_to_entity_part<Left>
   ,detail::convertable_to_entity_part<Right>
  >
 ,set_expression >::type
operator | (Left const& l, Right const& r)
{ return set_expression::union_(entity_part(l),(r)); }

//expr | term
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression >::type
operator | (set_expression const& l, T const& r)
{ return set_expression::union_(l,entity_part(r)); }

//term | expr
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression >::type
operator | (T const& l, set_expression const& r)
{ return set_expression::union_(entity_part(l),r); }

//expr | expr
inline
set_expression
operator | (set_expression const& l, set_expression const& r)
{ return set_expression::union_(l,r); }

//expr |= term
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression& >::type
operator |= (set_expression& l, T const& r)
{ return l = set_expression::union_(l,entity_part(r)); }

//expr |= expr
inline
set_expression &
operator |= (set_expression& l, set_expression const& r)
{ return l = set_expression::union_(l,r); }

//*******************************************************************
//difference
//*******************************************************************

//term - term
template <typename Left, typename Right>
inline typename boost::enable_if<
  boost::mpl::and_<
    detail::convertable_to_entity_part<Left>
   ,detail::convertable_to_entity_part<Right>
  >
 ,set_expression >::type
operator - (Left const& l, Right const& r)
{ return set_expression::difference_(entity_part(l),(r)); }

//expr - term
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression >::type
operator - (set_expression const& l, T const& r)
{ return set_expression::difference_(l,entity_part(r)); }

//term - expr
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression >::type
operator - (T const& l, set_expression const& r)
{ return set_expression::difference_(entity_part(l),r); }

//expr - expr
inline
set_expression
operator - (set_expression const& l, set_expression const& r)
{ return set_expression::difference_(l,r); }

//expr -= term
template <typename T>
inline typename boost::enable_if<
  detail::convertable_to_entity_part<T>
 ,set_expression& >::type
operator -= (set_expression& l, T const& r)
{ return l = set_expression::difference_(l,entity_part(r)); }

//expr -= expr
inline
set_expression &
operator -= (set_expression& l, set_expression const& r)
{ return l = set_expression::difference_(l,r); }

} // namespace samba

#endif //SAMBA_SAMBA_SET_EXPRESSION_OPERATORS_HPP
