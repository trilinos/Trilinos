// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_PROVIDERHELPERS_H
#define PIRO_PROVIDERHELPERS_H

#include "Piro_Provider.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include <functional>

namespace Piro {

//! \cond DETAILS

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


template <typename T, typename F>
class CachingProviderFunctor {
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


template <typename T, typename NullaryFunctor>
class NullaryProviderFunctorAdapter {
public:
  NullaryProviderFunctorAdapter() :
    nullaryFunctor_()
  {}

  explicit NullaryProviderFunctorAdapter(NullaryFunctor nf) :
    nullaryFunctor_(nf)
  {}

  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &/*params*/)
  {
    return nullaryFunctor_();
  }

public:
  NullaryFunctor nullaryFunctor_;
};

template <typename T, typename NullaryFunctor>
NullaryProviderFunctorAdapter<T, NullaryFunctor>
makeNullaryProviderFunctorAdapter(NullaryFunctor nf)
{
  return NullaryProviderFunctorAdapter<T, NullaryFunctor>(nf);
}


template <typename T, typename ReferenceAcceptingFunctor>
class ReferenceAcceptingProviderFunctorAdapter {
public:
  ReferenceAcceptingProviderFunctorAdapter() :
    referenceAcceptingFunctor_()
  {}

  explicit ReferenceAcceptingProviderFunctorAdapter(ReferenceAcceptingFunctor raf) :
    referenceAcceptingFunctor_(raf)
  {}

  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
    return referenceAcceptingFunctor_(*params);
  }

public:
  ReferenceAcceptingFunctor referenceAcceptingFunctor_;
};

template <typename T, typename ReferenceAcceptingFunctor>
ReferenceAcceptingProviderFunctorAdapter<T, ReferenceAcceptingFunctor>
makeReferenceAcceptingProviderFunctorAdapter(ReferenceAcceptingFunctor raf)
{
  return ReferenceAcceptingProviderFunctorAdapter<T, ReferenceAcceptingFunctor>(raf);
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
Provider<T> providerFromInstance(const Teuchos::RCP<T> &instance)
{
  return SharingProviderFunctor<T>(instance);
}


template <typename T>
inline
Provider<T> cachingProvider(const Provider<T> &p)
{
  return CachingProviderFunctor<T, Provider<T> >(p);
}


template <typename T, typename NullaryFunctor>
inline
Provider<T> providerFromNullary(NullaryFunctor nf)
{
  return NullaryProviderFunctorAdapter<T, NullaryFunctor>(nf);
}

template <typename T, typename NullaryFunctor>
inline
Provider<T> providerFromNullary()
{
  return NullaryProviderFunctorAdapter<T, NullaryFunctor>();
}


template <typename T, typename ReferenceAcceptingFunctor>
inline
Provider<T> providerFromReferenceAccepting(ReferenceAcceptingFunctor raf)
{
  return NullaryProviderFunctorAdapter<T, ReferenceAcceptingFunctor>(raf);
}

template <typename T, typename ReferenceAcceptingFunctor>
inline
Provider<T> providerFromReferenceAccepting()
{
  return NullaryProviderFunctorAdapter<T, ReferenceAcceptingFunctor>();
}

//! \endcond

} // namespace Piro

#endif /*PIRO_PROVIDERHELPERS_H*/
