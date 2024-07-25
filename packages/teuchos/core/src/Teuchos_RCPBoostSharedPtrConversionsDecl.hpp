// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_RCP_BOOST_SHAREDPTR_CONVERSIONS_DECL_HPP
#define TEUCHOS_RCP_BOOST_SHAREDPTR_CONVERSIONS_DECL_HPP

#include "Teuchos_RCPDecl.hpp"
#include "boost/shared_ptr.hpp"


namespace Teuchos {


/** \defgroup Teuchos_RCPBoostSharedPtrConversions_grp Conversion utilities for going between Teuchos::RCP and boost::shared_ptr.

The smart pointer classes <tt>Teuchos::RCP</tt> and <tt>boost::shared_ptr</tt>
are easily compatible.  The two templated conversion functions
<tt>Teuchos::rcp( const boost::shared_ptr<T> & )</tt> and
<tt>Teuchos::shared_pointer( const RCP<T> & )</tt> have been created for
converting back and forth (see the related non-member functions <tt>rcp()</tt>
and <tt>shared_pointer()</tt> the <tt>RCP</tt> classes' documentation).

\ingroup teuchos_mem_mng_grp

*/


/** \brief <tt>Teuchos::RCP</tt> Deallocator class that wraps a
 * <tt>boost::shared_ptr</tt>
 *
 * \ingroup Teuchos_RCPBoostSharedPtrConversions_grp
 */
template<class T>
class DeallocBoostSharedPtr
{
public:
  /** \brief. */
  DeallocBoostSharedPtr( const boost::shared_ptr<T> &sptr ) : sptr_(sptr) {}
  /** \brief. */
	typedef T ptr_t;
  /** \brief. */
	void free( T* ptr_in ) const { sptr_.reset(); }
  /** \brief. */
  const boost::shared_ptr<T>& ptr() const { return sptr_; }
private:
  mutable boost::shared_ptr<T> sptr_;
  DeallocBoostSharedPtr(); // Not defined and not to be called!
};


/** \brief <tt>boost::shared_ptr</tt> deleter class that wraps a
 * <tt>Teuchos::RCP</tt>.
 *
 * \ingroup Teuchos_RCPBoostSharedPtrConversions_grp
 */
template<class T>
class RCPDeleter
{
public:
  /** \brief. */
  RCPDeleter( const RCP<T> &rcp ) : rcp_(rcp) {}
  /** \brief. */
  typedef void result_type;
  /** \brief. */
  typedef T * argument_type;
  /** \brief. */
  void operator()(T * x) const { rcp_ = null; }
  /** \brief. */
  const RCP<T>& ptr() const { return rcp_; }
private:
  mutable RCP<T> rcp_;
  RCPDeleter(); // Not defined and not to be called!
};


/** \brief Conversion function that takes in a <tt>boost::shared_ptr</tt>
 * object and spits out a <tt>Teuchos::RCP</tt> object.
 *
 * If the input <tt>boost::shared_ptr</tt> already wraps a <tt>Teuchos::RCP</tt>
 * object, then that <tt>Teuchos::RCP</tt> object will be copied and returned.
 *
 * This function is not complicated, just look at its defintion below.
 *
 * \relates RCP
 */
template<class T>
RCP<T> rcp( const boost::shared_ptr<T> &sptr );


/** \brief Conversion function that takes in a <tt>Teuchos::RCP</tt>
 * object and spits out a <tt>boost::shared_ptr</tt> object.
 *
 * If the input <tt>Teuchos::RCP</tt> already wraps a
 * <tt>boost::shared_ptr</tt> object, then that <tt>boost::shared_ptr</tt>
 * object will be copied and returned.
 *
 * This function is not complicated, just look at its defintion below.
 *
 * \relates RCP
 */
template<class T>
boost::shared_ptr<T> shared_pointer( const RCP<T> &rcp );


} // namespace Teuchos


namespace boost {


/** \brief Returns true if <tt>p.get()==NULL</tt>.
 *
 * \ingroup Teuchos_RCPBoostSharedPtrConversions_grp
 */
template<class T> inline
bool is_null( const boost::shared_ptr<T> &p )
{
  return p.get() == 0;
}


/** \brief Returns true if <tt>p.get()!=NULL</tt>.
 *
 * \ingroup Teuchos_RCPBoostSharedPtrConversions_grp
 */
template<class T> inline
bool nonnull( const boost::shared_ptr<T> &p )
{
  return p.get() != 0;
}


} // namespace boost


#endif	// TEUCHOS_RCP_BOOST_SHAREDPTR_CONVERSIONS_DECL_HPP
