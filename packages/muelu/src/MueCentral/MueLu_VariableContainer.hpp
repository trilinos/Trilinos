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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_VARIABLECONTAINER_HPP
#define MUELU_VARIABLECONTAINER_HPP

#include <map>
#include <Teuchos_ParameterEntry.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_KeepType.hpp"

namespace MueLu {

  /*!
    @class VariableContainer
    @brief Class that stores all relevant data for a variable

    Maintains all data for a variable, that is, the data itself, a boolean flag for the "Keep" status,
    a boolean flag for the "Available" status, a reference counter for all requests and a list with
    all requesting factories.
  */
  class VariableContainer : public BaseClass {
  public:
    typedef std::map<const FactoryBase*,int> request_container;

  private:
    Teuchos::ParameterEntry           data_;        ///< the data itself
    bool                              available_;   ///< is data available?
    KeepType                          keep_;        ///< keep flag
    int                               count_;       ///< number of requests by all factories

    request_container                 requests_;    ///< requesting factories

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Default constructor.
    VariableContainer() : available_(false), keep_(false), count_(0) { }
    ~VariableContainer() { }

    //@}

    //! @name Data access
    //@{

    //! Store data in container class and set the "Available" status true.
    void SetData(const Teuchos::ParameterEntry& entry) {
      data_      = entry;
      available_ = true;
    }

    //! Return const reference to data stored in container
    //! NOTE: we do not check if data is available
    const Teuchos::ParameterEntry& GetData() const                              { return data_; }

    //! Return reference to data stored in container
    //! NOTE: we do not check if data is available
    Teuchos::ParameterEntry& GetData()                                          { return data_; }

    //! Returns true if data is available, i.e. SetData has been called before
    bool IsAvailable() const                                                    { return available_; }

    //@}

    //! @name Request/Release
    //@{

    //! Request data
    void Request(const FactoryBase* reqFactory) {
      request_container::iterator it = requests_.find(reqFactory);
      if (it == requests_.end())
        requests_[reqFactory] = 1;
      else
        (it->second)++;
      count_++;   // increment request counter
    }

    //! Release data
    void Release(const FactoryBase* reqFactory) {
      request_container::iterator it = requests_.find(reqFactory);
      TEUCHOS_TEST_FOR_EXCEPTION(it == requests_.end(), Exceptions::RuntimeError, "MueLu::VariableContainer::Release():"
                                 "cannot call Release if factory has not been requested before by factory " << reqFactory);
      if (--(it->second) == 0)
        requests_.erase(it);
      count_--;
    }

    //! Return the number of times the data has been requested by a specific factory
    int NumRequests(const FactoryBase* reqFactory) const {
      request_container::const_iterator it = requests_.find(reqFactory);
      return (it != requests_.end()) ? it->second : 0;
    }

    //! Returns the number of times the data has been requested
    int NumAllRequests() const                                                  { return count_; }

    //! Returns true, if data is requested by reqFactory
    bool IsRequested(const FactoryBase* reqFactory) const                       { return (NumRequests(reqFactory) > 0); }

    //! Returns true, if data is requested by at least one factory
    bool IsRequested() const                                                    { return (count_ > 0); }

    const request_container& Requests() const                                   { return requests_; }
    //@}

    //! @name Keep status
    //@{

    //! Returns true if at least one keep flag is set
    bool IsKept(KeepType keep) const                                            { return keep_ & keep; }

    //! Adds a keep flag to the flag combination
    void AddKeepFlag(KeepType keep = UserData)                                  { keep_ |= keep; }

    //! Removes a keep flag to the flag combination
    void RemoveKeepFlag(KeepType keep = UserData)                               { keep_ = keep_ & (keep_ ^ keep); }

    //! Returns the keep flag combination
    KeepType GetKeepFlag() const                                                { return keep_; }

    //@}
  };

}

#endif /* MUELU_VARIABLECONTAINER_HPP */

//TODO: move implementation to .cpp file + fwd decl of this class
