// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_FACTORY_HPP
#define MUELU_FACTORY_HPP

#include <string>
#include <map>

#include "Teuchos_RCP.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryAcceptor.hpp"
#include "MueLu_ParameterListAcceptor.hpp"

#include "MueLu_Level.hpp"

namespace MueLu {

  class Factory : public FactoryBase, public FactoryAcceptor, public ParameterListAcceptorImpl {

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    Factory() { }

    //! Destructor.
    virtual ~Factory() { }
    //@}

    //@{
    //! Configuration

    //! SetFactory is for expert users only. To change configuration of the preconditioner, use a factory manager.
    virtual void SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory) {
      RCP<const FactoryBase> f = factory;
      SetParameter(varName, ParameterEntry(f)); // parameter validation done in ParameterListAcceptorImpl
    }

    const RCP<const FactoryBase> GetFactory(const std::string & varName) const {
      return GetParameterList().get< RCP<const FactoryBase> >(varName);
    }

    // SetParameterList(...);

    // GetParameterList(...);

    //@}

    virtual RCP<const ParameterList> GetValidParameterList(const ParameterList& paramList = ParameterList()) const {
      return rcp(new ParameterList()); // no parameter by default - TMP
    }

  protected:

    void Input(Level & level, const std::string & varName) const {
      level.DeclareInput(varName, GetFactory(varName).get(), this);
    }

    template <class T>
    T Get(Level & level, const std::string & varName) const {
      return level.Get<T>(varName, GetFactory(varName).get());
    }

    template <class T>
    void Set(Level & level, const std::string & varName, const T & data) const {
      return level.Set<T>(varName, data, this);
    }

    bool IsAvailable(Level & level, const std::string & varName) const {
      return level.IsAvailable(varName, GetFactory(varName).get());
    }

  }; //class Factory

} //namespace MueLu

#define MUELU_FACTORY_SHORT
#endif //ifndef MUELU_FACTORY_HPP
