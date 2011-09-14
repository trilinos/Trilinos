#ifndef MUELU_LEVEL_HPP
#define MUELU_LEVEL_HPP

#include <iostream>
#include <sstream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Needs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_NoFactory.hpp"
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
class Level : public BaseClass { //public Needs {

private:
    mutable int levelID_; // id number associated with level
    RCP<DefaultFactoryHandlerBase> defaultFactoryHandler_;
    RCP<Level> previousLevel_;  // linked list of Level

    RCP<Needs> needs_;
public:

    //@{

    //! @name Constructors / Destructors
    Level() : levelID_(-1) {
    	needs_ = rcp(new Needs());
    }

    //! Constructor
    Level(RCP<DefaultFactoryHandlerBase> & defaultFactoryHandler) : levelID_(-1), defaultFactoryHandler_(defaultFactoryHandler) {
    	needs_ = rcp(new Needs());
    }

    //! Copy constructor.
    explicit Level(const Level& source) {
        levelID_ = source.levelID_;
        defaultFactoryHandler_ = source.defaultFactoryHandler_;
        previousLevel_ = source.previousLevel_;
        needs_ = source.needs_; // TODO: deep copy
        // TODO factorize with Build()
    }

    //@}

    //@{
    //! @name Build methods
    //! Builds a new Level object.
    RCP<Level> Build(std::ostream &os) { //TODO: why ostream in argument?
      // todo copy keep status of variables.
      RCP<Level> newLevel = rcp( new Level(defaultFactoryHandler_));

      // copy keep status of variables
      std::vector<std::string> ekeys = needs_->RequestedKeys();
      for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++)
      {
          std::vector<const MueLu::FactoryBase*> ehandles = needs_->RequestedHandles(*it);
          for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++)
          {
            const std::string ename = *it;
            const MueLu::FactoryBase* fac = *kt;
            if(isKept(ename,fac))
            {
              if(fac == NULL) newLevel->Keep(ename);
              else newLevel->Keep(ename,fac);
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

    //! @brief Set level number.
    void SetLevelID(int i) const { levelID_ = i; }

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
    //! If factory == NULL, use defaultFactory (if available).
    template <class T>
    void Set(const std::string ename, const T &entry, const FactoryBase* factory) {
      //      TEST_FOR_EXCEPTION(not null, or reference instead);
        const FactoryBase* fac = factory;
        if (factory == NULL)
        {
            fac = GetDefaultFactoryPtr(ename);
        }
        needs_->SetData<T>(ename, entry, fac);
    } //Set

    //! Store need label and its associated data. This does not increment the storage counter.
    //! no generating factory available. mark it as user-generated data
    template <class T>
    void Set(const std::string& ename, const T &entry) {
        Set<T>(ename, entry, MueLu::NoFactory::get());
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
              return needs_->GetData<T>(ename,defaultFactory);
            }
            // check if data for default factory has already been generated
            if(!needs_->IsAvailable(ename,defaultFactory))
            {
              //defaultFactory->callDeclareInput(*this);
              defaultFactory->NewBuild(*this);
            }

            TEST_FOR_EXCEPTION(! needs_->IsAvailable(ename,defaultFactory), Exceptions::RuntimeError, "MueLu::Level::Get(): factory did not produce expected output. no " << ename << " generated by " << defaultFactory);
            return needs_->GetData<T>(ename,defaultFactory);
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
            return needs_->GetData<T>(ename,factory);
        }
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    template <class T>
    T & Get(const std::string& ename)
    {
        return Get<T>(ename, MueLu::NoFactory::get());
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
        Value = Get<T>(ename,MueLu::NoFactory::get()); // todo fix me (call Needs::GetData directly)
    }

    //@}

    //! @name Permanent storage
    //@{

    ///! keep variable 'ename' generated by 'factory'
    virtual void Keep(const std::string& ename, const FactoryBase* factory)
    {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      needs_->Keep(ename,fac);
    }

    ///! keep variable 'ename' generated by no factory
    virtual void Keep(const std::string& ename)
    {
      //Needs::Keep(ename);
        Keep(ename,MueLu::NoFactory::get());
    }

    ///! returns true, if 'ename' generated by 'factory' is marked to be kept
    bool isKept(const std::string& ename) const
    {
      return isKept(ename,MueLu::NoFactory::get());
    }

    ///! returns true, if 'ename' generated by 'factory' is marked to be kept
    virtual bool isKept(const std::string& ename, const FactoryBase* factory) const
    {
      return needs_->isKept(ename,factory);
    }

    //@}

    //! @name Request/Release functions
    //! @brief Request and Release for incrementing/decrementing the reference count pointer for a specific variable.
    //@{
    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string& ename) {
    	needs_->Request(ename,MueLu::NoFactory::get());
    } //Request

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string& ename, const FactoryBase* factory, bool bCallDeclareInput = true) {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }

      TEST_FOR_EXCEPTION(fac == NULL, Exceptions::RuntimeError, "MueLu::Level::Request(): ptr to generating factory must not be null! ERROR.");

      // 1) check if 'ename' has already been requested or is already available
      //    if not: call DeclareInput of generating factory 'fac'
      if( bCallDeclareInput == true && //ename != "A" && // hack for RAPFactory test
    	 !IsRequested(ename,fac) &&
    	 !IsAvailable(ename,fac))
      {
    	  fac->callDeclareInput(*this);
      }

      // 2) request data 'ename' generated by 'fac'
      needs_->Request(ename,fac);
    }

    //! Decrement the storage counter.
    void Release(const std::string& ename)
    {
      //Needs::Release(ename,NULL);
        Release(ename,MueLu::NoFactory::get());
    } //Release

    //! Decrement the storage counter.
    void Release(const std::string& ename, const FactoryBase* factory)
    {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      needs_->Release(ename,fac);
    }

    //@}

    //! @name Utility functions
    //@{

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string ename) {
      //return Needs::IsAvailable(ename,NULL);
        return IsAvailable(ename,MueLu::NoFactory::get());
    }

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string ename, const FactoryBase* factory) {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      return needs_->IsAvailable(ename,fac);
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename) {
      //return Needs::IsRequested(ename,NULL);
        return IsRequested(ename,MueLu::NoFactory::get());
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename, const FactoryBase* factory) {
      const FactoryBase* fac = factory;
      if (factory == NULL)
      {
          fac = GetDefaultFactoryPtr(ename);
      }
      return needs_->IsRequested(ename,fac);
    }
    //@}

    //! @name Default factory handler
    //@{
    //! Set default factories (used internally by Hierarchy::SetLevel()).
    // Users should not use this method.
    void SetDefaultFactoryHandler(RCP<DefaultFactoryHandlerBase>& defaultFactoryHandler) {
        defaultFactoryHandler_ = defaultFactoryHandler;
    }

    //@}

    //! @name I/O Functions
    //@{

    /*! \brief Printing method for Needs class.*/
    std::ostream& print(std::ostream& os) const
    {
        Teuchos::TabularOutputter outputter(os);
        outputter.pushFieldSpec("name", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT,Teuchos::TabularOutputter::GENERAL,32);
        outputter.pushFieldSpec("gen. factory addr.", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
        outputter.pushFieldSpec("gen by", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 6);
        outputter.pushFieldSpec("req", Teuchos::TabularOutputter::INT,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 3);
        outputter.pushFieldSpec("type", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 10);
        outputter.pushFieldSpec("data", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 20);
        outputter.outputHeader();

        std::vector<std::string> ekeys = needs_->RequestedKeys();
        for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++)
        {
            std::vector<const MueLu::FactoryBase*> ehandles = needs_->RequestedHandles(*it);
            for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++)
            {
                outputter.outputField(*it);   // variable name
                outputter.outputField(*kt);   // factory ptr

                if(defaultFactoryHandler_->IsAvailable(*it)) // check if default factory is available for variable name
                {
                    if(GetDefaultFactoryPtr(*it)==*kt)
                       	outputter.outputField("def"); // factory ptr (deault factory)
                }
                else
                {
                    if (*kt == MueLu::NoFactory::get())
                    	outputter.outputField("user"); // factory ptr (user generated)
                    else
                    	outputter.outputField(" ");
                }

                int reqcount = 0;             // request counter
                reqcount = needs_->NumRequests(*it, *kt);
                outputter.outputField(reqcount);
                                              // variable type
                std::string strType = needs_->GetDataType(*it,*kt);
                if(strType.find("Xpetra::Operator")!=std::string::npos)
                {
                    outputter.outputField("Operator" );
                    outputter.outputField(" ");
                }
                else if(strType.find("Xpetra::MultiVector")!=std::string::npos)
                {
                  outputter.outputField("Vector");
                  outputter.outputField("");
                }
                else if(strType == "int")
                {
                    outputter.outputField(strType);
                    int data = 0; needs_->GetData<int>(*it,data,*kt);
                    outputter.outputField(data);
                }
                else if(strType == "double")
                {
                    outputter.outputField(strType);
                    double data = 0.0; needs_->GetData<double>(*it,data,*kt);
                    outputter.outputField(data);
                }
                else if(strType == "string")
                {
                    outputter.outputField(strType);
                    std::string data = ""; needs_->GetData<std::string>(*it,data,*kt);
                    outputter.outputField(data);
                }
                else
                {
                    outputter.outputField(strType);
                    outputter.outputField("unknown");
                }

                outputter.nextRow();
            }
        }
        return os;
    }

    //@}

private:

    //! Get ptr to default factory.
    const FactoryBase* GetDefaultFactoryPtr(const std::string& varname) const
    {
  		TEST_FOR_EXCEPTION(defaultFactoryHandler_ == null, Exceptions::RuntimeError, "MueLu::Level::GetDefaultFactory(): no DefaultFactoryHandler.");
   		return defaultFactoryHandler_->GetDefaultFactoryRCP(varname).get();
    }

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
