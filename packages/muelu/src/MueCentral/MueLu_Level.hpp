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
#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include "MueLu_BoostGraphviz.hpp"

#include "MueLu_BaseClass.hpp"

#include "MueLu_Needs.hpp"
#include "MueLu_NoFactory.hpp"

#include "MueLu_FactoryManagerBase_fwd.hpp"

namespace MueLu {

  /*!
    @class Level
    @brief Class that holds all level-specific information.

    All data is stored in an associative list. See the Needs class for more information.

    The Level class uses the functionality of the Needs class with the extended hashtables and
    adds the handling of default factories.
    All data that is stored in the <tt>Level</tt> class need a variable name (e.g. "A", "P", ...) and
    a pointer to the generating factory. Only with both the variable name and the generating factory
    the data can be accessed.

    If no pointer to the generating factory is provided (or it is NULL) then the Level class
    uses the information from a factory manager, which stores default factories for different
    variable names.
  */
  class Level : public BaseClass {

  public:
    //@{

    //! @name Constructors / Destructors
    Level();

    //! Constructor
    Level(RCP<FactoryManagerBase> & factoryManager);

    //@}

    //@{
    //! @name Build methods
    //! Builds a new Level object.
    RCP<Level> Build();

    //@}

    //! Destructor
    virtual ~Level();

    //@{
    //! @name Level handling

    //! @brief Return level number.
    int GetLevelID() const;

    //! @brief Set level number.
    void SetLevelID(int levelID);

    //! Previous level
    RCP<Level> & GetPreviousLevel();

    //! Set previous level object
    //! @\param[in] const RCP<Level>& previousLevel
    void SetPreviousLevel(const RCP<Level> & previousLevel);

    //! @name Set/Get factory manager
    //@{
    //! Set default factories (used internally by Hierarchy::SetLevel()).
    // Users should not use this method.
    void SetFactoryManager(const RCP<const FactoryManagerBase> & factoryManager);

    //! returns the current factory manager
    // Users should not use this method
    const RCP<const FactoryManagerBase> GetFactoryManager();
    //@}

    //@{
    //! @name Get functions

    //! Store need label and its associated data. This does not increment the storage counter.
    //! - If factory is not specified or factory == NoFactory::get(), mark data as user-defined by default
    //! - If factory is == NULL, use defaultFactory (FactoryManager have to be available).
    template <class T>
    void Set(const std::string & ename, const T &entry, const FactoryBase* factory = NoFactory::get()) {
      const FactoryBase* fac = GetFactory(ename, factory);

      if (fac == NoFactory::get()) {
        // user defined data
        // keep data
        AddKeepFlag(ename, NoFactory::get(), MueLu::UserData);
      }

      needs_.Set<T>(ename, entry, fac);

    } // Set

    //@}

    //! @name Get functions
    //! @brief Get functions for accessing stored data

    //@{
    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).
     *   Usage: Level->Get< RCP<Matrix> >("A", factory)
     *   if factory == NULL => use default factory
     *
     *  @param[in] const std::string& ename
     *  @param[in] const FactoryBase* factory
     *  @return data (templated)
     * */
    template <class T>
    T & Get(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
      const FactoryBase* fac = GetFactory(ename, factory);

      if (!IsAvailable(ename, fac)) {

        TEUCHOS_TEST_FOR_EXCEPTION(needs_.NumRequests(ename, fac) < 1 && needs_.GetKeepFlag(ename, fac) == 0, Exceptions::RuntimeError,
                                   "MueLu::Level::Get(): " << ename << " has not been requested (counter = " << needs_.NumRequests(ename, fac) << ", KeepFlag = " << needs_.GetKeepFlag(ename, fac) << "). " << std::endl << "Generating factory:" << *fac << " NoFactory="<<NoFactory::get());

        fac->CallBuild(*this);
        Release(*fac);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(!IsAvailable(ename, fac), Exceptions::RuntimeError, "MueLu::Level::Get(): factory did not produce expected output. " << ename << " has not been generated by " << *fac);

      return needs_.Get<T>(ename, fac);
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).*/
    template <class T>
    void Get(const std::string& ename, T& Value, const FactoryBase* factory = NoFactory::get()) {
      Value = Get<T>(ename, factory);
    }

    //TODO: add a const version of Get()?

    //@}

    //! @name Permanent storage
    //@{

    //! Request to keep variable 'ename' generated by 'factory' after the setup phase.
    //  This method is intented to be used by user drivers for printing, debugging or to keep some computed data for a next run of the setup phase.
    //
    //  This method is an alias for: AddKeepFlag(ename, factory, MueLu::Keep)
    //  See also the description of KeepEnum for more information.
    //
    //  To undo a keep request, one can use:
    //  - Delete(ename, factory) to delete the data and remove the "Keep" flag
    //  - or RemoveKeepFlag(ename, factory, MueLu::Keep) to go back to previous condition (data are kept only if internal MueLu logic need so).
    //
    // Note: Level variables tagged using this methods are also keep on the levels built using this.Build().
    //       This means that is you request to keep a specific variable on a fine level, all the coarser level that are created automatically during the setup phase will also retain the same variable.
    void Keep(const std::string & ename, const FactoryBase* factory) { AddKeepFlag(ename, factory, MueLu::Keep); } // Note: do not add default value for input parameter 'factory'

    //! Delete data that have been retained after the setup phase (using Keep(), AddKeepFlag(), or internal MueLu logic).
    // Special cases:
    // - If entry (ename, factory) does not exist, nothing is done.
    // - If entry exists but counter !=, entry cannot be desallocated before counter set to 0 (using Release()) so an exeption is thrown.
    void Delete(const std::string& ename, const FactoryBase* factory) { // Note: do not add default value for input parameter 'factory'
      if (!IsKey(ename, factory)) { return; } // if entry (ename, factory) does not exist, nothing is done.

      // Precondition:
      // Delete() should only be called if counter == 0
      // Note: It better to throw an exception rather than deleting the data if counter != 0 because users are not supposed to manipulate data with counter != 0
      TEUCHOS_TEST_FOR_EXCEPTION(IsRequested(ename, factory) == true, Exceptions::RuntimeError, "MueLu::Level::Delete(): IsRequested() == true. Ref counter != 0. You are not allowed to delete data that are still in use.");
      // If counter == 0 and entry exists, this means that a keep flag is set. Or there is an internal logic problem.
      TEUCHOS_TEST_FOR_EXCEPTION(GetKeepFlag(ename, factory) == 0, Exceptions::RuntimeError, "MueLu::Level::Delete(), Keep flag == 0?");

      RemoveKeepFlag(ename, factory, MueLu::All); // will delete the data if counter == 0

      // Post condition: data must have been deleted
      TEUCHOS_TEST_FOR_EXCEPTION(IsAvailable(ename, factory) == true, Exceptions::RuntimeError, "MueLu::Level::Delete(): Internal error (Post condition). Data have not been deleted.");
    }

    //! Test if a keep flag is set for variable 'ename' generated by 'factory'
    //! The input parameter keep can be a combination of flags. IsKept() will then return true if at least one of the flag is set.
    //! Note: There is no default parameter for IsKept() because it might be confusing (user generally wants to test IsKept with keep=MueLu::Keep but classes Level and Needs generally use keep = All)
    bool IsKept(const std::string& ename, const FactoryBase* factory, KeepType keep) const { return GetKeepFlag(ename, factory) & keep; }

    //! Alias for IsKept(ename, factory, keep) with factory = NoFactory::get()
    bool IsKept(const std::string& ename, KeepType keep) const { return GetKeepFlag(ename, NoFactory::get()) & keep; }

    //! Add a keep flag for variable 'ename' generated by 'factory'
    //! A variable can cumulate several keep flags (UserData+Final for example). This function just add a flag to the current flag combination.
    //! By default, the flag 'Keep' is added. So Keep(ename, factory) is the same as AddKeepFlag(ename, factory).
    //! See also the description of KeepEnum for more information.
    void AddKeepFlag(const std::string & ename, const FactoryBase* factory = NoFactory::get(), KeepType keep = MueLu::Keep); // TODO: remove default value for input parameter 'factory'?

    //! Remove a keep flag for variable 'ename' generated by 'factory'
    //! A variable can cumulate several keep flags (UserData+Final for example). This function just add a flag to the current flag combination.
    //! By default, all the flags are removed.
    void RemoveKeepFlag(const std::string & ename, const FactoryBase* factory, KeepType keep = MueLu::All);

    //! Get the flag combination set for variable 'ename' generated by 'factory'
    KeepType GetKeepFlag(const std::string& ename, const FactoryBase* factory = NoFactory::get()) const; // TODO: remove default value for input parameter 'factory'?

    //@}

    //! @name Request/Release functions
    //! @brief Request and Release for incrementing/decrementing the reference count pointer for a specific variable.
    //@{

    //! Increment the storage counter for all the inputs of a factory
    void Request(const FactoryBase& factory);

    //! Decrement the storage counter for all the inputs of a factory
    void Release(const FactoryBase& factory);

    //! Callback from FactoryBase::CallDeclareInput() and FactoryBase::DeclareInput()
    void DeclareInput(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy = NoFactory::get() );

    //! Callback from FactoryBase::CallDeclareInput() and FactoryBase::DeclareInput() to declare factory dependencies
    void DeclareDependencies(const FactoryBase* factory, bool bRequestOnly = false, bool bReleaseOnly = false);

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string& ename, const FactoryBase* factory = NoFactory::get(), const FactoryBase* requestedBy = NoFactory::get());

    //! Decrement the storage counter.
    void Release(const std::string& ename, const FactoryBase* factory = NoFactory::get(), const FactoryBase* requestedBy = NoFactory::get());

    //@}

    //! @name Utility functions
    //@{

    //! Test whether some information about (ename, factory) are stored
    bool IsKey(const std::string & ename, const FactoryBase* factory = NoFactory::get()) const;

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string & ename, const FactoryBase* factory = NoFactory::get()) const;

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string & ename, const FactoryBase* factory = NoFactory::get()) const;

    //@}

    //! @name I/O Functions
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Printing method
    // TODO: print only shows requested variables. check if we also list kept factories with ref counter=0?
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

#if defined(HAVE_MUELU_BOOST) && defined(BOOST_VERSION) && (BOOST_VERSION >= 104400)
    void UpdateGraph(std::map<const FactoryBase*, BoostVertex>&                   vindices,
                     std::map<std::pair<BoostVertex, BoostVertex>, std::string>&  edges,
                     BoostProperties&                                             dp,
                     BoostGraph&                                                  graph) const {
      needs_.UpdateGraph(vindices, edges, dp, graph);
    }
#endif

    //@}

  private:

    //! Copy constructor.
    Level(const Level& source);
    //
    // explicit Level(const Level& source);

    //! If input factory == NULL, returns the default factory. Else, return input factory.
    //
    //  If factory == NULL, the default factory is defined as follow:
    // - If user data is available, it is considered as the default and the factory manager is ignored.
    //   => The default factory is then NoFactory.
    // - Else, the factory manager is used to get the default factory.
    //
    // This strategy allows to use the same factory manager on the fine and coarse level without any trouble.
    //     Example :
    //
    //     FineLevel:
    //     ----------
    //     A          -> User provided
    //     Nullspace  -> User provided
    //
    //     CoarseLevel:
    //     ------------
    //     A         -> RAPFactory
    //     NullSpace -> NullspaceFactory
    //
    const FactoryBase* GetFactory(const std::string& varname, const FactoryBase* factory) const;

    enum   RequestMode { REQUEST, RELEASE, UNDEF }; //EI TODO
    static RequestMode requestMode_;                //EI TODO

    //
    //
    //

    int levelID_; // id number associated with level
    RCP<const FactoryManagerBase> factoryManager_;
    RCP<Level> previousLevel_;  // linked list of Level

    Needs needs_;

  }; //class Level

} //namespace MueLu

//TODO: Caps should not matter

#endif // MUELU_LEVEL_HPP
