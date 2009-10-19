#ifndef stk_util_util_CSet_hpp
#define stk_util_util_CSet_hpp

#include <typeinfo>
#include <vector>

namespace stk {

//----------------------------------------------------------------------
/** \ingroup util_module
 * \class CSet
 * \brief Set of entities of arbitrary types.
 *
 * @todo REFACTOR Use smart pointers to avoid destruction issue.
 *
 *  Example usage of the three methods:
 *
 * <PRE>
 *  class A { ... };
 *  class B { ... };
 *
 *  CSet cset ;
 *
 *  // Insert pointers to objects:
 *
 *  cset.insert<A>( new A ); // Do not delete on destruction
 *  cset.insert<B>( new B , true ); // Delete on destruction
 *
 *  // Query the collection of objects of a given type:
 *
 *  const A * sa = cset.get<A>();
 *  const B * sb = cset.get<B>();
 *
 *  // Remove a member:
 *
 *  {
 *    B * b = ... ;
 *    cset.remove<B>( b ); // Remove never deletes
 *    delete b ;
 *  }
 * </PRE>
 */
class CSet {
public:

  /** Get member conforming to the given type. */
  template<class T> const T * get() const ;

  /** Insert a new member.  Invoke 'delete' upon destruction.
   *  If already exists then return existing member, insert fails.
   */
  template<class T> const T * insert_with_delete( const T *);

  /** Insert a new member.  Do nothing to it upon destruction.
   *  If already exists then return existing member, insert fails.
   */
  template<class T> const T * insert_no_delete( const T * );

  /** Erase a member without deleting.
   *  Return if the remove operation was successful.
   */
  template<class T> bool remove( const T * );

  //--------------------------------

  ~CSet();
  CSet();

private:

  typedef void (*DeleteFunction)(void *);

  typedef std::pair< const std::type_info * , DeleteFunction > Manager ;

  const void * p_get( const std::type_info & ) const ;

  const void * p_insert( const Manager & , const void * );

  bool p_remove( const std::type_info & , const void * );

  std::vector< Manager > m_manager ;
  std::vector< const void * > m_value ;

  CSet( const CSet & );
  CSet & operator = ( const CSet & );
};

} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Inlined template methods have casting.

#ifndef DOXYGEN_COMPILE

namespace stk {

namespace {
template<class T>
void cset_member_delete( void * v ) { delete reinterpret_cast<T*>( v ); }
}

template<class T>
inline
const T * CSet::get() const
{ return (const T*) p_get( typeid(T) ); }

template<class T>
inline
const T * CSet::insert_with_delete( const T * arg_value)
{
  Manager m ;
  m.first = & typeid(T);
  m.second = & cset_member_delete<T> ;

  return (const T *) p_insert( m , arg_value );
}

template<class T>
inline
const T * CSet::insert_no_delete( const T * arg_value)
{
  Manager m ;
  m.first = & typeid(T);
  m.second = NULL ;

  return (const T *) p_insert( m , arg_value );
}

template<class T>
inline
bool CSet::remove( const T * arg_value )
{ return p_remove( typeid(T) , arg_value ); }

} // namespace stk

#endif /* DOXYGEN_COMPILE */

#endif // stk_util_util_CSet_hpp
