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

#include "Piro_Provider.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"

#include <string>
#include <map>

namespace Piro {

template <typename T>
class ExtensibleFactory {
public:
  ExtensibleFactory() :
    defaultProvider_(),
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
