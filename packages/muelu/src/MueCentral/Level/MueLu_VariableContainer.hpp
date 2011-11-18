/*
 * MueLu_VariableContainer.hpp
 *
 *  Created on: 17.11.2011
 *      Author: tobias
 */

#ifndef MUELU_VARIABLECONTAINER_HPP_
#define MUELU_VARIABLECONTAINER_HPP_

#include <Teuchos_ParameterEntry.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {

  /*!
    @class VariableContainer
    @brief Class that stores all relevant data for a variable

    Maintains all data for a variable, that is, the data itself, a boolean flag for the "Keep" status,
    a factory pointer of the generating factory, a reference counter for all requests and a list with
    all requesting factories.
  */
  class VariableContainer : public BaseClass {
  public:
    //! @name Constructors/Destructors.
    //@{

    //! Default constructor.
    VariableContainer() :
      available_(false), keep_(false), count_(0)
    {};

    virtual ~VariableContainer(){};

    //@}

    //! Store need label and its associated data. This does not increment the storage counter.
    void SetData(const Teuchos::ParameterEntry & entry) {
      data_ = entry;
      available_ = true;
    } //Set

    const Teuchos::ParameterEntry & GetData() const {
      return data_;
    }

    Teuchos::ParameterEntry & GetData() {
      return data_;
    }

    void Request(const FactoryBase* reqFactory) {
      count_++;   // increment request counter
      if(requests_.count(reqFactory)==0) {
        requests_[reqFactory] = 1;
      } else
      {
        int cnt = requests_[reqFactory];
        requests_[reqFactory] = ++cnt;
      }
    }

    void Release(const FactoryBase* reqFactory) {
      if(requests_.count(reqFactory) > 0) {
        int cnt = requests_[reqFactory];
        requests_[reqFactory] = --cnt;
        if(cnt == 0 && IsKept() == false) {
          requests_.erase(reqFactory);
        }
      }
      count_--; // decrement request counter
    }

    int NumRequests(const FactoryBase* reqFactory) const {
      if(requests_.count(reqFactory)>0) {
        return requests_.find(reqFactory)->second;
      }
      return 0;
    }

    int NumAllRequests() const {
      return count_;
    }

    bool IsRequested(const FactoryBase* reqFactory) const {
      if (NumRequests(reqFactory) > 0) return true;
      return false;
    }

    bool IsRequested() const {
      if (count_ > 0) return true;
      return false;
    }

    bool IsKept() const {
      return keep_;
    }

    void Keep(bool bKeep = true) {
      keep_ = bKeep;
    }

    bool IsAvailable() const {
      return available_;
    }
  private:
    Teuchos::ParameterEntry           data_;        ///< the data itself
    bool                              available_;   ///< is data available?
    bool                              keep_;        ///< keep flag
    int                               count_;       ///< number of requests by all factories
    std::map<const FactoryBase*,int>  requests_;    ///< requesting factories
  };

}

#endif /* MUELU_VARIABLECONTAINER_HPP_ */
