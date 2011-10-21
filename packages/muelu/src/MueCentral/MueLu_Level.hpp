#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <iostream>
#include <sstream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Needs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#undef HEAVY_DEBUG_OUTPUT

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
    Level() : levelID_(-1) { }

    //! Constructor
    Level(RCP<FactoryManagerBase> & factoryManager) : levelID_(-1), factoryManager_(factoryManager) { }

    //@}

    //@{
    //! @name Build methods
    //! Builds a new Level object.
    RCP<Level> Build() {
      RCP<Level> newLevel = rcp( new Level() );

      // Copy 'keep' status of variables
      // TODO: this only concerns needs_. so a function in Needs class should be provided to do that!
      typedef std::vector<std::string> KeyList;
      
      KeyList ekeys = needs_.RequestedKeys();
      for (KeyList::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
        std::vector<const MueLu::FactoryBase*> ehandles = needs_.RequestedFactories(*it);
        for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
          const std::string & ename = *it;
          const MueLu::FactoryBase* fac = *kt;
          if (IsKept(ename, fac)) {
            if (fac == NULL) 
              newLevel->Keep(ename);
            else 
              newLevel->Keep(ename, fac);
          }
        }
      }
      
      return newLevel;
    }

    //@}

    //! Destructor
    virtual ~Level() {}

    //@{
    //! @name Level handling

    //! @brief Return level number.
    int GetLevelID() const { return levelID_; }

    //! @brief Set level number.
    void SetLevelID(int levelID) {
      if (levelID_ != -1 && levelID_ != levelID)
        GetOStream(Warnings1, 0) << "Warning: Level::SetLevelID(): Changing an already defined LevelID (previousID=" << levelID_ << "newID=" << levelID << std::endl;

      levelID_ = levelID;
    }

    //! Previous level
    RCP<Level> & GetPreviousLevel() { return previousLevel_; }

    //! Set previous level object
    //! @\param[in] const RCP<Level>& previousLevel
    void SetPreviousLevel(const RCP<Level> & previousLevel) {
      if (previousLevel_ != Teuchos::null && previousLevel_ != previousLevel)
        GetOStream(Warnings1, 0) << "Warning: Level::SetPreviousLevel(): PreviousLevel was already defined" << std::endl;

      previousLevel_ = previousLevel;
    }

    //! @name Set factory manager
    //@{
    //! Set default factories (used internally by Hierarchy::SetLevel()).
    // Users should not use this method.
    void SetFactoryManager(const RCP<const FactoryManagerBase> & factoryManager) {
      factoryManager_ = factoryManager;
    }

    //@}

    //@{
    //! @name Get functions

    //! Store need label and its associated data. This does not increment the storage counter.
    //! - If factory is not specified, mark data as user-defined
    //! - If factory == NULL, use defaultFactory (if available).
    template <class T>
    void Set(const std::string & ename, const T &entry, const FactoryBase* factory = NoFactory::get()) {
      const FactoryBase* fac = GetFactory(ename, factory);
      
      if (fac == NoFactory::get()) {
        // user defined data
        // keep data
        Keep(ename, NoFactory::get());
      }
      
      needs_.Set<T>(ename, entry, fac);

    } // Set

    //@}

    //! @name Get functions
    //! @brief Get functions for accessing stored data

    //@{
    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).
     *   Usage: Level->Get< RCP<Operator> >("A", factory)
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
        TEST_FOR_EXCEPTION(needs_.NumRequests(ename, fac) < 1 && !needs_.IsKept(ename, fac), Exceptions::RuntimeError, 
                           "MueLu::Level::Get(): " << ename << " has not been requested (counter = " << needs_.NumRequests(ename, fac) << ", IsKept() = " << needs_.IsKept(ename, fac) << "). " << std::endl << "Generating factory:" << *fac);
        
        fac->CallBuild(*this);
        Release(*fac);
      }
      
      TEST_FOR_EXCEPTION(!IsAvailable(ename, fac), Exceptions::RuntimeError, "MueLu::Level::Get(): factory did not produce expected output. " << ename << " has not been generated by " << *fac);
      
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

    //! keep variable 'ename' generated by 'factory'
    void Keep(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
      /*TEST_FOR_EXCEPTION(factory == NULL && factoryManager_ == null, Exceptions::RuntimeError, "MueLu::Level::Keep(): no FactoryManager");
      const FactoryBase* fac = GetFactory(ename, factory);

      // call Request for factory fac only if
      // the factory has not been requested before and no data has been generated by fact (independent of ename) OR
      if (needs_.IsAvailableFactory(fac) == false && needs_.IsRequestedFactory(fac) == false) {
        Request(*fac);
      }*/
      needs_.Keep(ename, GetFactory(ename, factory));
    }

    //! returns true, if 'ename' generated by 'factory' is marked to be kept
    bool IsKept(const std::string& ename, const FactoryBase* factory = NoFactory::get()) const {
      return needs_.IsKept(ename, GetFactory(ename, factory));
    }

    /*! @brief remove the permanently stored variable 'ename' generated by 'factory' */
    void Delete(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
      needs_.Delete(ename, GetFactory(ename, factory));
    }

    //@}

    //! @name Request/Release functions
    //! @brief Request and Release for incrementing/decrementing the reference count pointer for a specific variable.
    //@{

    //! Increment the storage counter for all the inputs of a factory
    void Request(const FactoryBase& factory) {
      RequestMode prev = requestMode_;
      requestMode_ = REQUEST;
      factory.CallDeclareInput(*this);
      requestMode_ = prev;
    }

    //! Decrement the storage counter for all the inputs of a factory
    void Release(const FactoryBase& factory) {
      RequestMode prev = requestMode_;
      requestMode_ = RELEASE;
      factory.CallDeclareInput(*this);
      requestMode_ = prev;
    }

    //! Callback from FactoryBase::CallDeclareInput() and FactoryBase::DeclareInput()
    void DeclareInput(const std::string& ename, const FactoryBase* factory) {
      if (requestMode_ == REQUEST)
        Request(ename, factory);
      else if (requestMode_ == RELEASE)
        Release(ename, factory);
      else
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareInput(): requestMode_ undefined.");
    }

    //! Callback from FactoryBase::CallDeclareInput() and FactoryBase::DeclareInput() to declare factory dependencies
    void DeclareDependencies(const FactoryBase* factory, bool bRequestOnly = false, bool bReleaseOnly = false) { //TODO: replace bReleaseOnly, bReleaseOnly by one RequestMode enum
      if (bRequestOnly && bReleaseOnly)
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareDependencies(): Both bRequestOnly and bReleaseOnly set to true makes no sense.");

      if (requestMode_ == REQUEST) {

        if (bReleaseOnly == false) Request(*factory);

      } else if (requestMode_ == RELEASE) {

        if (bRequestOnly == false) Release(*factory);

      } else TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareDependencies(): requestMode_ undefined.");
    }

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
#if defined(HEAVY_DEBUG_OUTPUT)
      std::cout << "call Request(" << ename << ", " << factory << ")" << std::endl;
#endif

      const FactoryBase* fac = GetFactory(ename, factory);
      
      // call Request for factory fac only if
      // 1) the factory has not been requested before and no data has been generated by fact (independent of ename) OR
      // 2) the factory has not been used for generating data AND a variable name (ename) is requested to be kept!
      // TODO: is needs_.IsAvailableFactory(fac) == false enough? do we need to check the IsRequested status at all?
      if ((needs_.IsAvailableFactory(fac) == false && needs_.IsRequestedFactory(fac) == false)) { // ||
          //(needs_.IsAvailableFactory(fac) == false && needs_.IsRequestedFactory(fac) == true && needs_.IsKept(ename,fac) == true)) {
#if defined(HEAVY_DEBUG_OUTPUT)
        std::cout << "call Request(" << fac << ") [for declareInput]" << std::endl;
#endif
        Request(*fac);
      }

      // 2) request data 'ename' generated by 'fac'
      needs_.Request(ename, fac);
    }

    //! Decrement the storage counter.
    void Release(const std::string& ename, const FactoryBase* factory = NoFactory::get()) {
#if defined(HEAVY_DEBUG_OUTPUT)
      std::cout << "call Release(" << ename << ", " << factory << ")" << std::endl;
#endif

      const FactoryBase* fac = GetFactory(ename, factory);

      if(IsAvailable(ename,fac) == false && IsRequested(ename,fac) == true) {
        // data has been requested but never built
        // can we release the dependencies of fac safely?
        int cnt = needs_.CountRequestedFactory(fac);
        if(cnt == 1) {
          // factory is only generating factory of current variable (ename,factory)
          // Release(fac) can be called safely
          // TODO: not sure if this works correctly for kept data. Should not be a big problem now, but have to check this
          Release(*fac);
        }
        needs_.Release(ename, fac);
      }
      else if(IsAvailable(ename,fac) == true && IsRequested(ename,fac) == true) {
        needs_.Release(ename, fac);  // standard
      }
      else {
        // error!
        std::cout << "Try to release " << ename << " generated by " << fac << std::endl;
        TEST_FOR_EXCEPTION(IsAvailable(ename,fac) == false && IsRequested(ename,fac) == false && IsKept(ename,fac) == true, Exceptions::RuntimeError,
            "MueLu::Level::Release(): Data not requested and not available, but kept. Call of Release does not make sense. You should check this!");
        TEST_FOR_EXCEPTION(IsAvailable(ename,fac) == true && IsRequested(ename,fac) == false && IsKept(ename,fac) == false, Exceptions::RuntimeError,
            "MueLu::Level::Release(): Strange! Data is available but not requested and not kept??");
      }
    }

    //@}

    //! @name Utility functions
    //@{

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string & ename, const FactoryBase* factory = NoFactory::get()) {
      return needs_.IsAvailable(ename, GetFactory(ename, factory));
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string & ename, const FactoryBase* factory = NoFactory::get()) {
      return needs_.IsRequested(ename, GetFactory(ename, factory));
    }

    //@}

    //! @name I/O Functions
    //@{

    //! Printing method
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
      Teuchos::TabularOutputter outputter(out);
      outputter.pushFieldSpec("name",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 32);
      outputter.pushFieldSpec("gen. factory addr.", Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
      outputter.pushFieldSpec("req",                Teuchos::TabularOutputter::INT,    Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 3);
      outputter.pushFieldSpec("keep",               Teuchos::TabularOutputter::STRING,    Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 5);
      outputter.pushFieldSpec("type",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 10);
      outputter.pushFieldSpec("data",               Teuchos::TabularOutputter::STRING, Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 20);
      outputter.outputHeader();

      std::vector<std::string> ekeys = needs_.RequestedKeys();
      for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++) {
        std::vector<const MueLu::FactoryBase*> ehandles = needs_.RequestedFactories(*it);
        for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++) {
          outputter.outputField(*it);   // variable name
          outputter.outputField(*kt);   // factory ptr

          int reqcount = needs_.NumRequests(*it, *kt); // request counter
          outputter.outputField(reqcount);
          if (needs_.IsKept(*it, *kt)) {
            outputter.outputField("true");
          }
          else
          {
            outputter.outputField("false");
          }

          if (needs_.IsAvailable(*it, *kt)) {
            std::string strType = needs_.GetType(*it, *kt); // Variable type
            if (strType.find("Xpetra::Operator") != std::string::npos) {
              outputter.outputField("Operator" );
              outputter.outputField(" ");
            } else if (strType.find("Xpetra::MultiVector") != std::string::npos) {
              outputter.outputField("Vector");
              outputter.outputField("");
            } else if (strType.find("MueLu::SmootherBase") != std::string::npos) {
              outputter.outputField("SmootherBase");
              outputter.outputField("");
            } else if (strType == "int") {
              outputter.outputField(strType);
              int data = needs_.Get<int>(*it, *kt);
              outputter.outputField(data);
            } else if (strType == "double") {
              outputter.outputField(strType);
              double data = needs_.Get<double>(*it, *kt);
              outputter.outputField(data);
            } else if (strType == "string") {
              outputter.outputField(strType);
              std::string data = needs_.Get<std::string>(*it, *kt);
              outputter.outputField(data);
            } else {
              outputter.outputField(strType);
              outputter.outputField("");
            }
          } else {
            outputter.outputField("unknown");
            outputter.outputField("not available");
          }


          outputter.nextRow();
        }
      }

    }

    //@}

  private:
    
    //! Copy constructor.
    Level(const Level& source) { }
    //
    // explicit Level(const Level& source) {
    //   levelID_ = source.levelID_;
    //   factoryManager_ = source.factoryManager_;
    //   previousLevel_ = source.previousLevel_;
    //   needs_ = source.needs_; // TODO: deep copy
    //   // TODO factorize with Build()
    // }

    //! If input factory == NULL, returns the default factory. Else, return input factory.
    const FactoryBase* GetFactory(const std::string& varname, const FactoryBase* factory) const {
      if (factory == NULL) {
        TEST_FOR_EXCEPTION(factoryManager_ == null, Exceptions::RuntimeError, "MueLu::Level::GetFactory(): no FactoryManager");
        const FactoryBase* fac = factoryManager_->GetFactory(varname).get();
        TEST_FOR_EXCEPTION(fac == NULL, Exceptions::RuntimeError, "MueLu::Level::GetFactory(): Default factory cannot be NULL");
        return fac;
      } else {
        return factory;
      }
    }

    enum   RequestMode { REQUEST, RELEASE, UNDEF };
    static RequestMode requestMode_;

    //
    //
    //
  
    int levelID_; // id number associated with level
    RCP<const FactoryManagerBase> factoryManager_;
    RCP<Level> previousLevel_;  // linked list of Level

    Needs needs_;

  }; //class Level

} //namespace MueLu

#define MUELU_LEVEL_SHORT
#endif //ifndef MUELU_LEVEL_HPP

//TODO: Caps should not matter
