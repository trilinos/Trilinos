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

#ifndef PIRO_EXTENSIBLEFACTORY_H
#define PIRO_EXTENSIBLEFACTORY_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Teuchos_Ptr.hpp"

#include <string>
#include <map>
#include <functional>

namespace Piro {

template <typename T>
class ProviderBase {
public:
  virtual Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) = 0;

  ProviderBase() {}
  virtual ~ProviderBase() {}

private:
  // Disallow copy & assignment
  ProviderBase(const ProviderBase &);
  ProviderBase &operator=(const ProviderBase &);
};


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

  virtual Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) {
    return functor_(params);
  }

private:
  Functor functor_;
};


template <typename T>
class SharingProviderFunctor :
  public std::unary_function<Teuchos::RCP<Teuchos::ParameterList>, Teuchos::RCP<T> > {
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


template <typename T, typename F>
class CachingProviderFunctor :
  public std::unary_function<const Teuchos::RCP<Teuchos::ParameterList> &, Teuchos::RCP<T> > {
public:
  CachingProviderFunctor() :
    otherFunctor_(),
    instance_(Teuchos::null)
  {}

  explicit CachingProviderFunctor(F otherFunctor) :
    otherFunctor_(otherFunctor),
    instance_(Teuchos::null)
  {}

  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) {
    if (Teuchos::is_null(instance_)) {
      instance_ = otherFunctor_(params);
    }
    return instance_;
  }

private:
  F otherFunctor_;
  Teuchos::RCP<T> instance_;
};

template <typename T, typename F>
CachingProviderFunctor<T, F>
makeCachingProviderFunctor(F otherFunctor)
{
  return CachingProviderFunctor<T, F>(otherFunctor);
}


template <typename T>
struct ConstructorProviderFunctor {
public:
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return Teuchos::rcp(new T(params));
  }
};

template <typename T>
struct DefaultConstructorProviderFunctor {
public:
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &/*params*/) const {
    return Teuchos::rcp(new T);
  }
};

template <typename T>
struct ReferenceAcceptingConstructorProviderFunctor {
public:
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return Teuchos::rcp(new T(*params));
  }
};


template <typename T>
class Provider {
public:
  Provider() :
    ptr_(Teuchos::null)
  {}

  /* implicit */ Provider(Teuchos::ENull) :
    ptr_(Teuchos::null)
  {}

  /* implicit */ Provider(const Teuchos::RCP<ProviderBase<T> > &ptr) :
    ptr_(ptr)
  {}

  template <typename U>
  /* implicit */ Provider(const Teuchos::RCP<U> &instance) :
    ptr_(Teuchos::rcp(new ProviderImpl<T, SharingProviderFunctor<T> >(makeSharingProviderFunctor<T>(instance))))
  {}

  template <typename P>
  /* implicit */ Provider(const P &p) :
    ptr_(Teuchos::rcp(new ProviderImpl<T, P>(p)))
  {}

  bool is_null() const { return ptr_.is_null(); }
  Teuchos::RCP<const ProviderBase<T> > ptr() const { return ptr_; }
  Teuchos::RCP<ProviderBase<T> > ptr() { return ptr_; }

  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) {
    return (*ptr_)(params);
  }

private:
  Teuchos::RCP<ProviderBase<T> > ptr_;
};

template <typename T>
inline
bool is_null(const Provider<T> &h)
{
  return h.is_null();
}

template <typename T>
inline
bool nonnull(const Provider<T> &h)
{
  return !is_null(h);
}


template <typename T>
inline
Provider<T> providerFromInstance(const Teuchos::RCP<T> &instance)
{
  return SharingProviderFunctor<T>(instance);
}


template <typename T>
inline
Provider<T> providerFromConstructor()
{
  return ConstructorProviderFunctor<T>();
}

template <typename T>
inline
Provider<T> providerFromDefaultConstructor()
{
  return DefaultConstructorProviderFunctor<T>();
}

template <typename T>
inline
Provider<T> providerFromReferenceAcceptingConstructor()
{
  return ReferenceAcceptingConstructorProviderFunctor<T>();
}


template <typename T>
inline
Provider<T> cachingProvider(const Provider<T> &p)
{
  return CachingProviderFunctor<T, Provider<T> >(p);
}


template <typename T>
struct NullProviderFunctor :
  public std::unary_function<Teuchos::RCP<Teuchos::ParameterList>, Teuchos::RCP<T> > {
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &/*params*/) const { return Teuchos::null; }
};

template <typename T>
class ExtensibleFactory {
public:
  ExtensibleFactory() :
    defaultProvider_(NullProviderFunctor<T>()),
    selectorToken_("Type")
  {}

  const Provider<T> &defaultProvider() const { return defaultProvider_; }
  void setDefaultProvider(const Provider<T> &p) { defaultProvider_ = p; }

  const std::string &selectorToken() const { return selectorToken_; }
  void setSelectorToken(const std::string &t) { selectorToken_ = t; }

  Teuchos::RCP<const ProviderBase<T> > provider(const std::string &key) const;
  void setProvider(const std::string &key, const Provider<T> &p);

  Teuchos::RCP<T> create(const Teuchos::RCP<Teuchos::ParameterList> &params);

private:
  Provider<T> defaultProvider_;

  std::string selectorToken_;

  Provider<T> getProvider(const std::string &key) const;
  typedef std::map<std::string, Provider<T> > ProviderMap;
  ProviderMap providers_;
};

template <typename T>
Teuchos::RCP<T>
ExtensibleFactory<T>::create(const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  const Teuchos::Ptr<const std::string> type(Teuchos::getParameterPtr<std::string>(*params, selectorToken_));
  if (Teuchos::nonnull(type)) {
    const std::string &key = *type;
    Provider<T> typeProvider = this->getProvider(key);
    if (nonnull(typeProvider)) {
      const Teuchos::RCP<Teuchos::ParameterList> &providerParams = Teuchos::sublist(params, key);
      return typeProvider(providerParams);
    }
  }
  return defaultProvider_(params);
}

template <typename T>
Provider<T>
ExtensibleFactory<T>::getProvider(const std::string &key) const
{
  const typename ProviderMap::const_iterator it = providers_.find(key);
  return (it != providers_.end()) ? it->second : Provider<T>();
}

template <typename T>
Teuchos::RCP<const ProviderBase<T> >
ExtensibleFactory<T>::provider(const std::string &key) const
{
  return this->getProvider(key).ptr();
}

template <typename T>
void
ExtensibleFactory<T>::setProvider(const std::string &key, const Provider<T> &h)
{
  if (nonnull(h)) {
    providers_[key] = h;
  } else {
    providers_.erase(key);
  }
}

} // namespace Piro

#endif /*PIRO_EXTENSIBLEFACTORY_H*/
