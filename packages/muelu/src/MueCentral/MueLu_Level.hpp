#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <iostream>
#include <sstream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Needs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_DefaultFactoryHandlerBase.hpp"

namespace MueLu {

/*!
    @class Level
    @brief Class that holds all level-specific information.

    All data is stored in an associative list. See the Needs class for more information.

    The Level class uses the functionality of the Needs class with the extended hashtables and
    adds the handling of default factories.
    All data that is stored in the <tt>Level</tt> class need a variable name (e.g. "A", "P",...) and
    a pointer to the generating factory. Only with both the variable name and the generating factory
    the data can be accessed.

    If no pointer to the generating factory is provided (or it is NULL) then the Level class
    uses the information from a default factory handler, which stores default factories for different
    variable names.
 */
class Level : public Needs {

private:
    mutable int levelID_; // id number associated with level
    RCP<DefaultFactoryHandlerBase> defaultFactoryHandler_;

    // linked list of Level
    RCP<Level> previousLevel_;

protected:
    RCP<Teuchos::FancyOStream> out_;

public:

    //@{
    //! @name Constructors / Destructors
    Level() : levelID_(-1), out_(this->getOStream()) {
        //Teuchos::OSTab tab(out_); MueLu_cout(Teuchos::VERB_HIGH) << "Instantiating new uninitialized Level" << std::endl;
    }

    //! Constructor
    Level(RCP<DefaultFactoryHandlerBase>& defaultFactoryHandler) : levelID_(-1), defaultFactoryHandler_(defaultFactoryHandler), out_(this->getOStream()) {
    }

    //! Copy constructor.
    explicit Level(const Level& source) {
        levelID_ = source.levelID_;
        defaultFactoryHandler_ = source.defaultFactoryHandler_;
    }

    //@}

    //@{
    //! @name Build methods //TODO: merge with copy constructor?
    //! Builds a new Level object.
    RCP<Level> Build(std::ostream &os) { //TODO: why ostream in argument?
      // todo copy keep status of variables.
      RCP<Level> newLevel = rcp( new Level(defaultFactoryHandler_));

      // copy keep status of variables
      copyKeepStatus(newLevel);

      return newLevel;
    }
    //@}

    virtual ~Level() {}

    void Print(std::ostream &os) {
        os << this << std::endl;
    }

    //@{
    //! @name Level handling

    //! @brief Set level number.
    void SetLevelID(int i) const {
        levelID_ = i;
    }

    //! @brief Return level number.
    int GetLevelID() const { return levelID_; }

    void SetPreviousLevel(const RCP<Level> & previousLevel) {
        previousLevel_ = previousLevel;
    }

    //! Previous level
    RCP<Level> & GetPreviousLevel() { return previousLevel_; }

    //@}

    //@{
    //! @name Set methods.

    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void Set(const std::string ename, const T &entry, const FactoryBase* factory) {
      //      TEST_FOR_EXCEPTION(not null, or reference instead);
        const FactoryBase* fac = factory;
        if (factory == NULL)
        {
            fac = GetDefaultFactoryPtr(ename);
        }
        Needs::SetData<T>(ename, entry, fac);
    } //Set

    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void Set(const std::string& ename, const T &entry, Teuchos::RCP<const FactoryBase> factory) {
        Set<T>(ename, entry, factory.get());
    }

    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void Set(const std::string& ename, const T &entry) {
        Set<T>(ename, entry, NULL);
    }

    //@}

    //! @name Get functions
    //! @brief Get functions for accessing stored data

    //@{
    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A", factory)
    // factory == NULL => use default factory
    template <class T>
    T & Get(const std::string& ename, const FactoryBase* factory)
    {
        // data not available
        // if no generating factory given, use DefaultFactoryHandler
        if (factory == NULL)
        {
            const FactoryBase* defaultFactory = GetDefaultFactoryPtr(ename);
            if( defaultFactory == NULL)
            {
              return Needs::GetData<T>(ename,defaultFactory);
            }
            // check if data for default factory has already been generated
            if(!Needs::IsAvailable(ename,defaultFactory))
            {
              //defaultFactory->callDeclareInput(*this);
              defaultFactory->NewBuild(*this);
            }

            TEST_FOR_EXCEPTION(! Needs::IsAvailable(ename,defaultFactory), Exceptions::RuntimeError, "MueLu::Level::Get(): factory did not produce expected output. no " << ename << " generated by " << defaultFactory);
            return Needs::GetData<T>(ename,defaultFactory);
        }
        else
        {
            // variable 'ename' generated by 'factory' available in Level
            if (  !IsAvailable(ename, factory) )
            {
              //factory->callDeclareInput(*this);
              factory->NewBuild(*this);
            }

            TEST_FOR_EXCEPTION(! IsAvailable(ename,factory), Exceptions::RuntimeError, "MueLu::Level::Get(): factory did not produce expected output.");
            return Needs::GetData<T>(ename,factory);
        }
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A", factory)
    // factory == Teuchos::null => use default factory
    template <class T>
    T & Get(const std::string& ename, Teuchos::RCP<const FactoryBase> factory)
    {
        return Get<T>(ename,factory.get());
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    template <class T>
    T & Get(const std::string& ename)
    {
        return Get<T>(ename, NULL);
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).*/
    template <class T>
    void Get(const std::string& ename, T& Value, Teuchos::RCP<const FactoryBase> factory)
    {
        Value = Get<T>(ename,factory.get());
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).*/
    template <class T>
    void Get(const std::string& ename, T& Value, const FactoryBase* factory)
    {
        Value = Get<T>(ename,factory);
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    template <class T>
    void Get(const std::string& ename, T& Value)
    {
        Value = Get<T>(ename,NULL); // todo fix me (call Needs::GetData directly)
    }

    //@}

    //! @name Permanent storage
    //@{

    ///! keep variable 'ename' generated by 'factory'
    void Keep(const std::string& ename, RCP<const FactoryBase> factory)
    {
      Keep(ename,factory.get());
    }

    ///! keep variable 'ename' generated by 'factory'
    virtual void Keep(const std::string& ename, const FactoryBase* factory)
    {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      Needs::Keep(ename,fac);
    }

    ///! keep variable 'ename' generated by no factory
    virtual void Keep(const std::string& ename)
    {
      //Needs::Keep(ename);
        Keep(ename,NULL);
    }

    ///! returns true, if 'ename' generated by 'factory' is marked to be kept
    bool isKept(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null)
    {
      return isKept(ename,factory.get());
    }

    ///! returns true, if 'ename' generated by 'factory' is marked to be kept
    virtual bool isKept(const std::string& ename, const FactoryBase* factory)
    {
      return Needs::isKept(ename,factory);
    }

    //@}

    //! @name Request/Release functions
    //! @brief Request and Release for incrementing/decrementing the reference count pointer for a specific variable.
    //@{
    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string& ename) {
      //Needs::Request(ename,NULL);
        Request(ename,NULL);
    } //Request

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string& ename, RCP<const FactoryBase> factory) {
      Request(ename,factory.get());
    } //Request

    //!
    void Request(const std::string& ename, const FactoryBase* factory) {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      Needs::Request(ename,fac);
    }

    //! Decrement the storage counter.
    void Release(const std::string& ename, RCP<const FactoryBase> factory)
    {
      Release(ename,factory.get());
    } //Release

    //! Decrement the storage counter.
    void Release(const std::string& ename)
    {
      //Needs::Release(ename,NULL);
        Release(ename,NULL);
    } //Release

    //! Decrement the storage counter.
    void Release(const std::string& ename, const FactoryBase* factory)
    {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      Needs::Release(ename,fac);
    }

    //@}

    //! @name Utility functions
    //@{
    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string ename, RCP<const FactoryBase> factory) {
      return IsAvailable(ename,factory.get());
    }

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string ename) {
      //return Needs::IsAvailable(ename,NULL);
        return IsAvailable(ename,NULL);
    }

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string ename, const FactoryBase* factory) {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      return Needs::IsAvailable(ename,fac);
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename, RCP<const FactoryBase> factory) {
      return IsRequested(ename,factory.get());
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename) {
      //return Needs::IsRequested(ename,NULL);
        return IsRequested(ename,NULL);
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename, const FactoryBase* factory) {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      return Needs::IsRequested(ename,fac);
    }
    //@}

    //! @name Default factory handler
    //@{
    //! Set default factories (used internally by Hierarchy::SetLevel()).
    // Users should not use this method.
    void SetDefaultFactoryHandler(RCP<DefaultFactoryHandlerBase>& defaultFactoryHandler) {
        defaultFactoryHandler_ = defaultFactoryHandler;
    }

    //! Get ptr to default factory. // TODO: make me private again
    const FactoryBase* GetDefaultFactoryPtr(const std::string& varname) {
        TEST_FOR_EXCEPTION(defaultFactoryHandler_ == null, Exceptions::RuntimeError, "MueLu::Level::GetDefaultFactory(): no DefaultFactoryHandler.");
        return defaultFactoryHandler_->GetDefaultFactoryRCP(varname).get();
    }
    //@}
private:
    //! Get RCP to default factory for given variable 'varname'.
    RCP<const FactoryBase> GetDefaultFactoryRCP(const std::string& varname) {
        TEST_FOR_EXCEPTION(defaultFactoryHandler_ == null, Exceptions::RuntimeError, "MueLu::Level::GetDefaultFactory(): no DefaultFactoryHandler.");
        return defaultFactoryHandler_->GetDefaultFactoryRCP(varname);
    }

    //! Get default factory.
    const FactoryBase & GetDefaultFactory(const std::string& varname) {
        TEST_FOR_EXCEPTION(defaultFactoryHandler_ == null, Exceptions::RuntimeError, "MueLu::Level::GetDefaultFactory(): no DefaultFactoryHandler.");
        return defaultFactoryHandler_->GetDefaultFactory(varname);
    }

}; //class Level

std::ostream& operator<<(std::ostream& os, Level const &level);

} //namespace MueLu

#define MUELU_LEVEL_SHORT
#endif //ifndef MUELU_LEVEL_HPP
