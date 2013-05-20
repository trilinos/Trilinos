// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PIRO_PROVIDER_H
#define PIRO_PROVIDER_H

#include "Piro_ProviderBase.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include <functional>

namespace Piro {

//! \cond DETAILS

template <typename T, typename Functor>
class ProviderImpl : public ProviderBase<T> {
public:
  ProviderImpl() :
    functor_()
  {}

  explicit ProviderImpl(const Functor &functor) :
    functor_(functor)
  {}

  const Functor &functor() const { return functor_; }
  Functor &functor() { return functor_; }

  virtual Teuchos::RCP<T> getInstance(const Teuchos::RCP<Teuchos::ParameterList> &params);

private:
  Functor functor_;
};

template <typename T, typename Functor>
Teuchos::RCP<T>
ProviderImpl<T, Functor>::getInstance(const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  return functor_(params);
}

template <typename T>
struct ProviderFunctorBase :
  public std::unary_function<const Teuchos::RCP<Teuchos::ParameterList> &, Teuchos::RCP<T> > {};


template <typename T>
struct NullProviderFunctor : public ProviderFunctorBase<T> {
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &/*params*/) const
  {
    return Teuchos::null;
  }
};


template <typename T>
class SharingProviderFunctor : public ProviderFunctorBase<T> {
public:
  SharingProviderFunctor() :
    instance_(Teuchos::null)
  {}

  explicit SharingProviderFunctor(Teuchos::ENull) :
    instance_(Teuchos::null)
  {}

  explicit SharingProviderFunctor(const Teuchos::RCP<T> &instance) :
    instance_(instance)
  {}

  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &/*params*/) const {
    return instance_;
  }

private:
  Teuchos::RCP<T> instance_;
};

template <typename T>
inline
SharingProviderFunctor<T>
makeSharingProviderFunctor(const Teuchos::RCP<T> &instance)
{
  return SharingProviderFunctor<T>(instance);
}

//! \endcond

/*! \brief Handle for auxiliary object factories
 *
 *  A Provider is a smart handle that offers functor (function object) and value semantics
 *  as well as several implicit conversions as a convenience to manipulate objects that implement the ProviderBase interface.
 *
 *  It is essentially a thin wrapper around <tt>Teuchos::RCP<ProviderBase<T>></tt>.
 *
 *  The solver factory Piro::Epetra::SolverFactory uses \link Provider Providers\endlink
 *  as sources of auxiliary objects for the different %Piro solvers.
 */
template <typename T>
class Provider : public ProviderFunctorBase<T> {
public:
  //! \name Constructors
  //@{
  //! \brief Constructs a Provider that returns null pointers.
  Provider() :
    ptr_(Teuchos::rcp(new ProviderImpl<T, NullProviderFunctor<T> >))
  {}

  //! \brief Constructs an uninitialized Provider.
  //! \details A valid Provider \e must be assigned to <tt>*this</tt> before Provider::operator() is called.
  //!          This constructor exists mostly for debugging purposes and its use is discouraged.
  //!          The \link Provider::Provider() default constructor \endlink should be preferred in most cases.
  /* implicit */ Provider(Teuchos::ENull) :
    ptr_(Teuchos::null)
  {}

  //! \brief Constructs a Provider returning an already existing instance.
  //! \details The template parameter \c U must refer to a type that can be converted to \c T.
  template <typename U>
  /* implicit */ Provider(const Teuchos::RCP<U> &instance) :
    ptr_(Teuchos::rcp(new ProviderImpl<T, SharingProviderFunctor<T> >(makeSharingProviderFunctor<T>(instance))))
  {}

  //! \brief Constructs a Provider from a function object.
  //! \details The template parameter \c P must refer to a copy-constructible function object type
  //!          such that the following code is valid
  //!          \code Teuchos::RCP<T> instance = p(params); \endcode
  //!          where \c params denotes a variable of type <tt>Teuchos::RCP<Teuchos::ParameterList><tt/>.
  template <typename P>
  /* implicit */ Provider(const P &p) :
    ptr_(Teuchos::rcp(new ProviderImpl<T, P>(p)))
  {}

  //! \brief Constructs a Provider handle that wraps the provided implementation.
  /* implicit */ Provider(const Teuchos::RCP<ProviderBase<T> > &ptr_in) :
    ptr_(ptr_in)
  {}
  //@}


  //! \name Core functionality
  //@{
  //! \brief Returns an instance of a subclass of the type T.
  //! \details \pre <tt>this->nonnull()</tt>
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) {
    return ptr_->getInstance(params);
  }
  //@}

  //! \name Validity check
  //@{
  //! \brief Checks that the Provider has been initialized
  bool nonnull() const { return ptr_.nonnull(); }

  //! \brief Checks whether the Provider has been left uninitialized
  bool is_null() const { return ptr_.is_null(); }
  //@}

  //! \name Access to internals
  //@{
  //! \brief Returns a <tt>const</tt> pointer to the internal implementation.
  Teuchos::RCP<const ProviderBase<T> > ptr() const { return ptr_; }

  //! \brief Returns a non-<tt>const</tt> pointer to the internal implementation.
  Teuchos::RCP<ProviderBase<T> > ptr() { return ptr_; }
  //@}

private:
  Teuchos::RCP<ProviderBase<T> > ptr_;
};

//! \name Validity check
//@{

//! \brief Returns \c true if \c handle is initialized.
//! \relates Provider
template <typename T>
inline
bool nonnull(const Provider<T> &handle)
{
  return handle.nonnull();
}

//! \brief Returns \c true if \c handle is uninitialized.
//! \relates Provider
template <typename T>
inline
bool is_null(const Provider<T> &handle)
{
  return handle.is_null();
}

//@}

} // namespace Piro

#endif /*PIRO_PROVIDER_H*/
