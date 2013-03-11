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
#include "Teuchos_TypeTraits.hpp"

#include <string>
#include <map>
#include <functional>

namespace Piro {

template <typename T>
class ProviderBase {
public:
  virtual Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const = 0;
  ~ProviderBase() {}
};

template <typename T, typename Functor>
class ProviderImpl : public ProviderBase<T> {
public:
  explicit ProviderImpl(const Functor &functor) :
    functor_(functor)
  {}

  virtual Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
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


template <bool b> struct TruthType {};
typedef TruthType<true> TrueType;
typedef TruthType<false> FalseType;

template <typename T, typename ArgumentHandlingStrategy>
class CreatingProviderFunctor :
  public std::unary_function<Teuchos::RCP<Teuchos::ParameterList>, Teuchos::RCP<T> > {
public:
  CreatingProviderFunctor() :
    strategy_()
  {}

  explicit CreatingProviderFunctor(const ArgumentHandlingStrategy &strategy) :
    strategy_(strategy)
  {}

  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return Teuchos::rcp(
        this->create(
          params,
          TruthType<Teuchos::TypeTraits::is_same<typename ArgumentHandlingStrategy::result_type, void>::value>()));
  }

private:
  T* create(
      const Teuchos::RCP<Teuchos::ParameterList> &params,
      FalseType /*strategy_returns_void*/) const {
    return new T(strategy_(params));
  }

  T* create(
      const Teuchos::RCP<Teuchos::ParameterList> &/*params*/,
      TrueType /*strategy_returns_void*/) const {
    return new T;
  }

  ArgumentHandlingStrategy strategy_;
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

  /* implicit */ Provider(const Teuchos::RCP<const ProviderBase<T> > &ptr) :
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

  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return (*ptr_)(params);
  }

private:
  Teuchos::RCP<const ProviderBase<T> > ptr_;
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


struct Forward :
  public std::unary_function<Teuchos::RCP<Teuchos::ParameterList>, Teuchos::RCP<Teuchos::ParameterList> > {
  const Teuchos::RCP<Teuchos::ParameterList> &operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return params;
  }
};

struct Dereference :
  public std::unary_function<Teuchos::RCP<Teuchos::ParameterList>, Teuchos::ParameterList> {
  Teuchos::ParameterList &operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return *params;
  }
};

struct Ignore :
  public std::unary_function<Teuchos::RCP<Teuchos::ParameterList>, void> {
  void operator()(const Teuchos::RCP<Teuchos::ParameterList> &/*params*/) const {
    return;
  }
};


template <typename T>
inline
Provider<T> makeProvider(const Teuchos::RCP<T> &instance)
{
  return SharingProviderFunctor<T>(instance);
}

template <typename T, typename ArgumentHandlingStrategy>
inline
Provider<T> makeProvider()
{
  return CreatingProviderFunctor<T, ArgumentHandlingStrategy>();
}

template <typename T>
inline
Provider<T> makeProvider()
{
  return CreatingProviderFunctor<T, Forward>();
}


template <typename T>
class ExtensibleFactory {
public:
  explicit ExtensibleFactory(const std::string &selectorToken = "Type") :
    selectorToken_(selectorToken)
  {}

  const std::string &selectorToken() const { return selectorToken_; }

  Teuchos::RCP<const ProviderBase<T> > provider(const std::string &key) const;

  void setProvider(const std::string &key, const Provider<T> &h);

  Teuchos::RCP<T> create(const Teuchos::RCP<Teuchos::ParameterList> &params) const;

private:
  Provider<T> getProvider(const std::string &key) const;

  std::string selectorToken_;

  typedef std::map<std::string, Provider<T> > ProviderMap;
  ProviderMap providers_;
};

template <typename T>
Teuchos::RCP<T>
ExtensibleFactory<T>::create(const Teuchos::RCP<Teuchos::ParameterList> &params) const
{
  Teuchos::RCP<T> result;
  const Teuchos::Ptr<const std::string> type(Teuchos::getParameterPtr<std::string>(*params, selectorToken_));
  if (Teuchos::nonnull(type)) {
    const std::string &key = *type;
    const Provider<T> typeProvider = this->getProvider(key);
    if (nonnull(typeProvider)) {
      const Teuchos::RCP<Teuchos::ParameterList> &providerParams = Teuchos::sublist(params, key);
      result = typeProvider(providerParams);
    }
  }
  return result;
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
