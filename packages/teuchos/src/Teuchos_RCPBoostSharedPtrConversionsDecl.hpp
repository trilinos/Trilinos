// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_RCP_SHAREDPTR_CONVERSIONS_DECL_HPP
#define TEUCHOS_RCP_SHAREDPTR_CONVERSIONS_DECL_HPP

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


#endif	// TEUCHOS_RCP_SHAREDPTR_CONVERSIONS_DECL_HPP
