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
#ifndef MUELU_NEEDS_HPP
#define MUELU_NEEDS_HPP

#include <string>
#include <Teuchos_ParameterEntry.hpp>
#include "MueLu_VariableContainer.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_Exceptions.hpp"

#include "MueLu_TwoKeyMap.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_NoFactory.hpp"

#ifdef HAVE_MUELU_BOOST
#include "boost/graph/graphviz.hpp"
#endif

namespace MueLu {

  /*!
    @class Needs
    @brief Class that allows cross-factory communication of data needs.

    Maintains a list of 'Needs' for a given Level. For example, a restriction factory that
    transposes the tentative prolongator 'Needs' the prolongator factory to save this.

    The data is stored using a variable name and a pointer to the generating factory. We use
    a TwoKeyMap object with the pointer to the generating factory as primary key and the variable
    name as secondary key. The pointer to the generating factory is only used as primary "key".
    For the Needs class it doesn't matter if the given factory pointer is valid or not.
    So the NULL pointer can also be used.

    The data itself is stored within a VariableContainer object.
    The 'Needs' class only provides basic routines for handling the TwoKeyMap object with all
    VariableContainer objects.
    A reference counter keeps track of the storage and automatically frees the memory if
    the data is not needed any more. In the standard mode, the data first has to be requested
    by calling the Request function. Then the data can be set by calling Set.
    With Get the user can fetch data when needed. Release decrements the reference counter
    for the given variable.
  */
  class Needs : public BaseClass {

  private:

    UTILS::TwoKeyMap<const FactoryBase*, std::string, RCP<MueLu::VariableContainer> > dataTable_;

    //TODO: Key1 = const std::string?

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Default constructor.
    Needs();

    virtual ~Needs();

    //@}

    //! @name Set
    //! @brief functions for setting data in data storage
    //@{

    //      void Set(const Key1 & key1, const Key2 & key2, const Value & entry)

    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void Set(const std::string & ename, const T & entry, const FactoryBase* factory) {
      // Store entry only if data have been requested (or any keep flag)
      if (IsRequested(ename, factory) ||
          GetKeepFlag(ename, factory) != 0) {
        TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory, ename), Exceptions::RuntimeError, "MueLu::Needs::Set(): " + ename + " not found in dataTable_");
        const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
        var->SetData( Teuchos::ParameterEntry(entry));
      }
    } //Set

    //@}

    //! @name Request/Release data

    //@{

    //! Indicate that an object is needed. This increments the storage counter.
    //void Request(const std::string & ename, const FactoryBase* factory); //Request
    void Request(const std::string & ename, const FactoryBase* factory, const FactoryBase* requestedBy = MueLu::NoFactory::get()); //Request

    //! Decrement the storage counter.
    //void Release(const std::string & ename, const FactoryBase* factory); //Release
    void Release(const std::string & ename, const FactoryBase* factory, const FactoryBase* requestedBy = MueLu::NoFactory::get()); //Release

    //@}

    //! @name Get data
    //@{

    //! @brief Get data without decrementing associated storage counter (i.e., read-only access)
    // Usage: Level->Get< RCP<Matrix> >("A", factoryPtr)
    template <class T>
    const T & Get(const std::string & ename, const FactoryBase* factory) const {
      TEUCHOS_TEST_FOR_EXCEPTION(!dataTable_.IsKey(factory, ename), Exceptions::RuntimeError, "MueLu::Needs::Get(): " + ename + " not found in dataTable_");
      const Teuchos::RCP<MueLu::VariableContainer> & var = dataTable_.Get(factory,ename);
      return Teuchos::getValue<T>(var->GetData());
    }

    //! @brief Get data without decrementing associated storage counter (i.e., read-only access)
    // Usage: Level->Get< RCP<Matrix> >("A", factoryPtr)
    template <class T>
    T & Get(const std::string & ename, const FactoryBase* factory) {
      return const_cast<T &>(const_cast<const Needs &>(*this).Get<T>(ename, factory)); // Valid cast. See Effective C++, Item 3.
    }

    //@}

    //! @name Permanent storage
    //@{

    //! See documentation of Level::IsKept
    bool IsKept(KeepType keep, const std::string& ename, const FactoryBase* factory = NoFactory::get()) const { return GetKeepFlag(ename, factory) & keep; }

    //! See documentation of Level::AddKeepFlag
    void AddKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep = MueLu::Keep);

    //! See documentation of Level::RemoveKeepFlag
    void RemoveKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep = MueLu::All);

    //! See documentation of Level::GetKeepFlag
    KeepType GetKeepFlag(const std::string& ename, const FactoryBase* factory) const;

    //@}

    //! @name Utilities.
    //@{

    //! Test whether some information about (ename, factory) are stored
    bool IsKey(const std::string & ename, const FactoryBase* factory) const;

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string & ename, const FactoryBase* factory) const;

    //! Test whether some data has been requested.
    // Note1: IsRequested() can be used inside of a factory to know if it worth it to build some data.
    // (if a factory is called, at least one ename of the factory is requested but maybe not all of them)
    // Note2: this tells nothing about whether the data actually exists.
    bool IsRequested(const std::string & ename, const FactoryBase* factory) const;

    bool IsRequestedBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const;

    //! Test whether a factory is generating factory of a requested variable in Needs
    bool IsRequestedFactory(const FactoryBase* factory);

    //! returns how often factory is generating factory of requested variables in Needs
    //! used to decide whether Level.Release(factory) for releasing factory dependencies can
    //! be called safely
    int CountRequestedFactory(const FactoryBase* factory);

    //! Test whether a factory is among the generating factories of data that is already available
    bool IsAvailableFactory(const FactoryBase* factory);

    //! @brief Return the number of outstanding requests for a need.
    //!  Throws a <tt>Exceptions::RuntimeError</tt> exception if the need either hasn't been requested or hasn't been saved.
    int NumRequests(const std::string & ename, const FactoryBase* factory) const;

    int NumRequestsBy(const FactoryBase* factory, const std::string & ename, const FactoryBase* requestedBy) const;

    //! @name Helper functions
    //@{

    //! Returns a vector of strings containing all key names of requested variables

    std::vector<std::string> RequestedKeys(const FactoryBase* factory) const;

    std::vector<const FactoryBase*> RequestedFactories() const;

    std::string GetType(const std::string & ename, const FactoryBase* factory) const;

    //@}

    //! @name I/O Functions
    //@{

    //! Printing method
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

#if defined(HAVE_MUELU_BOOST) && defined(BOOST_VERSION) && (BOOST_VERSION >= 104400)
    void UpdateGraph(std::map<const FactoryBase*, long unsigned>&                           vindices,
                     std::map<std::pair<long unsigned,long unsigned>, std::string>&         edges,
                     boost::dynamic_properties&                                             dp,
                     boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                     boost::property<boost::vertex_name_t, std::string,
                     boost::property<boost::vertex_color_t, std::string,
                     boost::property<boost::vertex_index_t, std::string> > >,
                     boost::property<boost::edge_name_t, std::string,
                     boost::property<boost::edge_color_t, std::string> > >&                 graph) const;
#endif

    //@}

  private:

    //! Copy constructor
    Needs(const Needs & source);

  }; //class Needs

} //namespace MueLu

#define MUELU_NEEDS_SHORT
#endif // MUELU_NEEDS_HPP
