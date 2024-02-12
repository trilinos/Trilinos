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

#include <Teuchos_TypeNameTraits.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Operator.hpp>

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
 private:
  // Motivated by Teuchos_any.hpp
  class DataBase {
   public:
    virtual ~DataBase() {}
    virtual const std::type_info& type() const = 0;
    virtual std::string typeName() const       = 0;
  };

  template <typename T>
  class Data : public DataBase {
   public:
    Data(const T& data)
      : data_(data) {}
    const std::type_info& type() const { return typeid(T); }
    std::string typeName() const { return Teuchos::TypeNameTraits<T>::name(); }
    T data_;
  };

  template <typename T>
  struct Getter {
    static T& get(DataBase* data_, DataBase*& /* datah_ */) {
      if ((data_ == NULL) || (data_->type() != typeid(T)))  // NVR added guard to avoid determining typeName unless we will use it
      {
        const std::string typeName = Teuchos::TypeNameTraits<T>::name();  // calls Teuchos::demangleName(), which can be expensive
        TEUCHOS_TEST_FOR_EXCEPTION(data_ == NULL, Teuchos::bad_any_cast,
                                   "Error, cast to type Data<" << typeName << "> failed since the content is NULL");

        TEUCHOS_TEST_FOR_EXCEPTION(data_->type() != typeid(T), Teuchos::bad_any_cast,
                                   "Error, cast to type Data<" << typeName << "> failed since the actual underlying type is "
                                                                              "\'"
                                                               << data_->typeName() << "!");
      }

      Data<T>* data = dynamic_cast<Data<T>*>(data_);
      if (!data)  // NVR added guard to avoid determining typeName unless we will use it
      {
        const std::string typeName = Teuchos::TypeNameTraits<T>::name();  // calls Teuchos::demangleName(), which can be expensive
        TEUCHOS_TEST_FOR_EXCEPTION(!data, std::logic_error,
                                   "Error, cast to type Data<" << typeName << "> failed but should not have and the actual underlying type is "
                                                                              "\'"
                                                               << data_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!");
      }
      return data->data_;
    }
  };

  template <typename T>
  struct Checker {
    static bool check(DataBase* data_, DataBase*& /* datah_ */) {
      if ((data_ == NULL) || (data_->type() != typeid(T)))  // NVR added guard to avoid determining typeName unless we will use it
      {
        return false;
      }

      Data<T>* data = dynamic_cast<Data<T>*>(data_);
      if (!data)  // NVR added guard to avoid determining typeName unless we will use it
      {
        return false;
      }
      return true;
    }
  };

 public:
  typedef std::map<const FactoryBase*, int> request_container;

 private:
  DataBase* data_;           ///< the data itself
  mutable DataBase* datah_;  ///< temporary data storage (need to get a reference
                             ///< to RCP to a base class (like Operator)
  bool available_;           ///< is data available?
  KeepType keep_;            ///< keep flag
  int count_;                ///< number of requests by all factories

  request_container requests_;  ///< requesting factories

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Default constructor.
  VariableContainer()
    : data_(NULL)
    , datah_(NULL)
    , available_(false)
    , keep_(false)
    , count_(0) {}
  ~VariableContainer() {
    delete data_;
    data_ = NULL;
    delete datah_;
    datah_ = NULL;
  }

  //@}

  //! @name Data access
  //@{

  //! Store data in container class and set the "Available" status true.
  template <typename T>
  void SetData(const T& entry) {
    delete data_;
    delete datah_;
    data_      = new Data<T>(entry);
    datah_     = NULL;
    available_ = true;
  }

  //! Return const reference to data stored in container
  //! NOTE: we do not check if data is available
  template <typename T>
  const T& GetData() const {
    return Getter<T>::get(data_, datah_);
  }

  //! Return reference to data stored in container
  //! NOTE: we do not check if data is available
  template <typename T>
  T& GetData() {
    return Getter<T>::get(data_, datah_);
  }

  //! Return reference to data stored in container
  //! NOTE: we do not check if data is available
  template <typename T>
  bool CheckType() {
    return Checker<T>::check(data_, datah_);
  }

  std::string GetTypeName() {
    if (data_ == NULL)
      return std::string("");
    return data_->typeName();
  }

  //! Returns true if data is available, i.e.
  //  if SetData has been called before
  bool IsAvailable() const { return available_; }

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
    count_++;  // increment request counter
  }

  //! Release data
  void Release(const FactoryBase* reqFactory) {
    request_container::iterator it = requests_.find(reqFactory);
    TEUCHOS_TEST_FOR_EXCEPTION(it == requests_.end(), Exceptions::RuntimeError,
                               "MueLu::VariableContainer::Release(): "
                               "cannot call Release if factory has not been requested before by factory "
                                   << reqFactory);
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
  int NumAllRequests() const { return count_; }

  //! Returns true, if data is requested by reqFactory
  bool IsRequested(const FactoryBase* reqFactory) const { return (NumRequests(reqFactory) > 0); }

  //! Returns true, if data is requested by at least one factory
  bool IsRequested() const { return (count_ > 0); }

  const request_container& Requests() const { return requests_; }
  //@}

  //! @name Keep status
  //@{

  //! Returns true if at least one keep flag is set
  bool IsKept(KeepType keep) const { return keep_ & keep; }

  //! Adds a keep flag to the flag combination
  void AddKeepFlag(KeepType keep = UserData) { keep_ |= keep; }

  //! Removes a keep flag to the flag combination
  void RemoveKeepFlag(KeepType keep = UserData) { keep_ = keep_ & (keep_ ^ keep); }

  //! Returns the keep flag combination
  KeepType GetKeepFlag() const { return keep_; }

  //@}
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct VariableContainer::Getter<Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > {
  typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> Operator;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;

  static Teuchos::RCP<Operator>& get(DataBase* data_, DataBase*& datah_) {
    typedef Teuchos::RCP<Operator> TO;
    typedef Teuchos::RCP<Matrix> TM;

    if (data_ == NULL)  // NVR added guard to avoid unnecessary call to TypeNameTraits<TO>::name(), which calls Teuchos::demangleName(), which can be expensive
    {
      TEUCHOS_TEST_FOR_EXCEPTION(data_ == NULL, Teuchos::bad_any_cast,
                                 "Error, cast to type Data<" << Teuchos::TypeNameTraits<TO>::name() << "> failed since the content is NULL");
    }
    if (data_->type() == typeid(TO)) {
      Data<TO>* data = dynamic_cast<Data<TO>*>(data_);
      if (!data)  // NVR added guard to avoid unnecessary call to TypeNameTraits<TO>::name(), which calls Teuchos::demangleName(), which can be expensive
      {
        TEUCHOS_TEST_FOR_EXCEPTION(!data, std::logic_error,
                                   "Error, cast to type Data<" << Teuchos::TypeNameTraits<TO>::name() << "> failed but should not have and the actual underlying type is "
                                                                                                         "\'"
                                                               << data_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!");
      }
      return data->data_;
    }

    if (data_->type() != typeid(TM))  // NVR added guard to avoid unnecessary call to TypeNameTraits<TO>::name(), which calls Teuchos::demangleName(), which can be expensive
    {
      TEUCHOS_TEST_FOR_EXCEPTION(data_->type() != typeid(TM), Teuchos::bad_any_cast,
                                 "Error, cast to type Data<" << Teuchos::TypeNameTraits<TM>::name() << "> failed since the actual underlying type is "
                                                                                                       "\'"
                                                             << data_->typeName() << "!");
    }
    Data<TM>* data = dynamic_cast<Data<TM>*>(data_);
    if (!data)  // NVR added guard to avoid unnecessary call to TypeNameTraits<TO>::name(), which calls Teuchos::demangleName(), which can be expensive
    {
      TEUCHOS_TEST_FOR_EXCEPTION(!data, std::logic_error,
                                 "Error, cast to type Data<" << Teuchos::TypeNameTraits<TM>::name() << "> failed but should not have and the actual underlying type is "
                                                                                                       "\'"
                                                             << data_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!");
    }
    if (datah_ == NULL)
      datah_ = new Data<TO>(Teuchos::rcp_dynamic_cast<Operator>(data->data_));
    Data<TO>* datah = dynamic_cast<Data<TO>*>(datah_);

    if (!datah)  // NVR added guard to avoid unnecessary call to TypeNameTraits<TO>::name(), which calls Teuchos::demangleName(), which can be expensive
    {
      TEUCHOS_TEST_FOR_EXCEPTION(!datah, std::logic_error,
                                 "Error, cast to type Data<" << Teuchos::TypeNameTraits<TO>::name() << "> failed but should not have and the actual underlying type is "
                                                                                                       "\'"
                                                             << datah_->typeName() << "! The problem might be related to incompatible RTTI systems in static and shared libraries!");
    }
    return datah->data_;
  }
};

}  // namespace MueLu

#endif /* MUELU_VARIABLECONTAINER_HPP */
