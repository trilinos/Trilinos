// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_REFCOUNTPTR_SHAREDPTR_CONVERSIONS_DECL_H
#define TEUCHOS_REFCOUNTPTR_SHAREDPTR_CONVERSIONS_DECL_H

#include "Teuchos_RefCountPtrDecl.hpp"
#include "boost/shared_ptr.hpp"

namespace Teuchos {

/** \defgroup Teuchos_RefCountPtrBoostSharedPtrConversions_grp Conversion utilities for going between Teuchos::RefCountPtr and boost::shared_ptr.

The smart pointer classes <tt>Teuchos::RefCountPtr</tt> and
<tt>boost::shared_ptr</tt> are easily compatible.  The two templated
conversion functions <tt>Teuchos::rcp( const boost::shared_ptr<T> & )</tt> and
<tt>Teuchos::shared_pointer( const RefCountPtr<T> & )</tt> have been created
for converting back and forth.

The following code shows how to convert back and forth between these two smart
pointer types:

\ingroup teuchos_mem_mng_grp

*/
//@{

/** \brief <tt>Teuchos::RefCountPtr</tt> Deallocator class that wraps a
 * <tt>boost::shared_ptr</tt> */
template<class T>
class DeallocBoostSharedPtr
{
public:
  /** \brief. */
  DeallocBoostSharedPtr( const boost::shared_ptr<T> &sptr ) : sptr_(sptr) {}
  /** \brief. */
	typedef T ptr_t;
  /** \brief. */
	void free( T* ptr ) const { sptr_.reset(); }
  /** \brief. */
  const boost::shared_ptr<T>& ptr() const { return sptr_; }
private:
  mutable boost::shared_ptr<T> sptr_;
  DeallocBoostSharedPtr(); // Not defined and not to be called!
};

/** \brief <tt>boost::shared_ptr</tt> deleter class that wraps a
 * <tt>Teuchos::RefCountPtr</tt>.
 */
template<class T>
class RefCountPtrDeleter
{
public:
  /** \brief. */
  RefCountPtrDeleter( const RefCountPtr<T> &rcp ) : rcp_(rcp) {}
  /** \brief. */
  typedef void result_type;
  /** \brief. */
  typedef T * argument_type;
  /** \brief. */
  void operator()(T * x) const { rcp_ = null; }
  /** \brief. */
  const RefCountPtr<T>& ptr() const { return rcp_; }
private:
  mutable RefCountPtr<T> rcp_;
  RefCountPtrDeleter(); // Not defined and not to be called!
};

/** \brief Conversion function that takes in a <tt>boost::shared_ptr</tt>
 * object and spits out a <tt>Teuchos::RefCountPtr</tt> object.
 *
 * If the input <tt>boost::shared_ptr</tt> already wraps a <tt>Teuchos::RefCountPtr</tt>
 * object, then that <tt>Teuchos::RefCountPtr</tt> object will be copied and returned.
 *
 * This function is not complicated, just look at its defintion below.
 */
template<class T>
RefCountPtr<T> rcp( const boost::shared_ptr<T> &sptr );

/** \brief Conversion function that takes in a <tt>Teuchos::RefCountPtr</tt>
 * object and spits out a <tt>boost::shared_ptr</tt> object.
 *
 * If the input <tt>Teuchos::RefCountPtr</tt> already wraps a
 * <tt>boost::shared_ptr</tt> object, then that <tt>boost::shared_ptr</tt>
 * object will be copied and returned.
 *
 * This function is not complicated, just look at its defintion below.
 */
template<class T>
boost::shared_ptr<T> shared_pointer( const RefCountPtr<T> &rcp );

//@}

} // namespace Teuchos

#endif	// TEUCHOS_REFCOUNTPTR_SHAREDPTR_CONVERSIONS_DECL_H
