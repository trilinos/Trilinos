/*
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
*/

#ifndef TEUCHOS_SHARED_ACCESS_POLICY_HPP
#define TEUCHOS_SHARED_ACCESS_POLICY_HPP


#include "Teuchos_ConfigDefs.hpp"


//
// WARNING: This current file is just for iterating on the design of thread
// safety in Teuchos and is not working code yet!
//


namespace Teuchos {


/** \brief Basic portable thread lock primative class.
 *
 * This is a heavy-weight appraoch to lock a thread.
 *
 * This class may or may not be exposed in the Teuchos namespace.  It might
 * just be exposed as a typedef in SharedAccessPolicy.
 */
class ThreadLock {
public:
  /** \brief . */
  ThreadLock();
  /** \brief . */
  ~ThreadLock();
  /** \brief . */
  bool try_lock();
  /** \brief . */
  void lock();
  /** \brief . */
  void unlock();
};


/** \brief Stack-based object for locking a thread.
 *
 * This class is used to lock a function so that no matter how the function is
 * existed the lock will be released.
 *
 * This class may or may not be exposed in the Teuchos namespace.  It might
 * just be exposed as a typedef in SharedAccessPolicy.
 */
template<class T>
class ScopedThreadLock {
public:
  /** \brief . */
  explicit ScopedThreadLock(ThreadLock &lock);
  /** \brief . */
  ~ScopedThreadLock();
};


/** \brief Single policy class defining an approach for sharing an integral
 * object across threads as well as a general heavy-weight locking policy.
 *
 * This policy class provides an primative integral type
 * (atomic_integral_type) and a set of atomic primitives for
 * incrementing/decrementing and/or setting/fetching the object.
 *
 * This policy class also provides typdefs for generic heavier-weight thread
 * locking objects that can be used for more general locking.
 *
 * This class is designed, in purpose, to provide sufficient functionality to
 * make the Teuchos MM classes thread-safe.  However, the will have other uses
 * in Teuchos and other C++ code also.
 *
 * There will likely be a typedef called DefaultSharedAccessPolicy in the
 * Teuchos namespace that will select the default global policy (i.e. the one
 * used by the Teuchos MM classes by default).  Also, there well as concrete
 * implementations NonthreadedSharedAccessPolicy, TbbSharedAccessPolicy,
 * PtheadsSharedAccessPolicy, and others (when support is enabled).
 */
class SharedAccessPolicy {
public:
  /** \brief . */
  typedef ThreadLock lock_type;
  /** \brief . */
  typedef ScopedThreadLock scoped_lock_type;
  /** \brief Supported type for shared integral objects. */
  typedef int atomic_integral_type;
  /** \brief Atomic setting a shared integral object. */
  inline static void atomic_set( atomic_integral_type * p,
    const atomic_integral_type v );
  /** \brief Atomic fetch a shared integral object. */
  inline static const atomic_integral_type
  atomic_fetch( const atomic_integral_type * p );
  /** \brief Atomic increment of a shared integral object. */
  inline static void atomic_increment( atomic_integral_type * p );
  /** \brief Atomic decrement of a shared integral object. */
  inline static void atomic_decrement( atomic_integral_type * p );
  // ToDo: Define some other basic fetch/increment primatives needed for
  // better performance, for example, for Teuchos::RCPNode.
};


} // namespace Teuchos


#endif /* TEUCHOS_SHARED_ACCESS_POLICY_HPP */
