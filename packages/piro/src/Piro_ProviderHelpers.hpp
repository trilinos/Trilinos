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

#ifndef PIRO_PROVIDERHELPERS_H
#define PIRO_PROVIDERHELPERS_H

#include "Piro_Provider.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include <functional>

namespace Piro {

template <typename T>
struct ConstructorProviderFunctor : public ProviderFunctorBase<T> {
public:
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return Teuchos::rcp(new T(params));
  }
};

template <typename T>
struct DefaultConstructorProviderFunctor : public ProviderFunctorBase<T> {
public:
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &/*params*/) const {
    return Teuchos::rcp(new T);
  }
};

template <typename T>
struct ReferenceAcceptingConstructorProviderFunctor : public ProviderFunctorBase<T> {
public:
  Teuchos::RCP<T> operator()(const Teuchos::RCP<Teuchos::ParameterList> &params) const {
    return Teuchos::rcp(new T(*params));
  }
};


template <typename T, typename F>
class CachingProviderFunctor : public ProviderFunctorBase<T> {
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
class NullaryProviderFunctorAdapter : public ProviderFunctorBase<T> {
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
class ReferenceAcceptingProviderFunctorAdapter : public ProviderFunctorBase<T> {
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

} // namespace Piro

#endif /*PIRO_PROVIDERHELPERS_H*/
