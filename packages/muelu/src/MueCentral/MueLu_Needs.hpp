#ifndef MUELU_NEEDS_HPP
#define MUELU_NEEDS_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Exceptions.hpp"

#include <MueLu_ExtendedHashtable.hpp>

#define MueLu_cout(minimumVerbLevel)                                    \
  if (this->getVerbLevel() >= minimumVerbLevel) *(this->getOStream())

namespace MueLu {
  
  /*!
    @class Needs
    @brief Class that allows cross-factory communication of data needs.

    Maintains a list of 'Needs' for a given Level. For example, a restriction factory that
    transposes the tentative prolongator 'Needs' the prolongator factory to save this.

  */
  class Needs : public BaseClass {

    /*
    //! Friendly print.  FIXME pretty print doesn't work ;(
    friend std::ostream& operator<<(std::ostream &os, Needs const &foo);
    */

  private:

    //! Prior output level
    Teuchos::EVerbosityLevel priorOutputLevel_;
    //! Stores number of outstanding requests for a need.
    MueLu::UTILS::ExtendedHashtable countTable_;
    //! Stores data associated with a need.
    MueLu::UTILS::ExtendedHashtable dataTable_;

    //! flag: setup phase
    /*!
     * if setupPhase=true, store always data, even if not requested
     * Note: handling slightly different than in MueMat!
     */
    bool setupPhase_;
  protected:
    RCP<Teuchos::FancyOStream> out_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Default constructor.
    Needs() : out_(this->getOStream()) {
      setupPhase_ = true; //TODO fix me.
    }

    virtual ~Needs() {}
    //@}

    //! @name Access methods
    //@{
    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void SetData(const std::string& ename, const T &entry, RCP<const FactoryBase> factory = Teuchos::null) {
        SetData(ename,entry,factory.get());
    } //Set


    template <class T>
    void SetData(const std::string& ename, const T &entry, const FactoryBase* factory) {
        // check if in setupPhase (=store all data) or data is requested
        if(setupPhase_ ||
            (countTable_.isKey(ename, factory) && countTable_.Get<int>(ename, factory) != 0))
        {
          if(!countTable_.isKey(ename, factory))
            countTable_.Set(ename,0,factory);    // make sure that 'ename' is counted
          dataTable_.Set(ename,entry,factory);
        }
    }

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string ename, RCP<const FactoryBase> factory = Teuchos::null) {
      Request(ename,factory.get());
    } //Request

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string ename, const FactoryBase* factory) {
      // if it's the first request for 'ename', create a new key in the hashtable
      if(!countTable_.isKey(ename, factory))
        countTable_.Set(ename,0,factory);

      // increment counter
      IncrementCounter(ename,factory);
    } //Request

    //! Decrement the storage counter.
    void Release(const std::string ename, RCP<const FactoryBase> factory = Teuchos::null)
    {
      Release(ename,factory.get());
    } //Release

    //! Decrement the storage counter.
    void Release(const std::string ename, const FactoryBase* factory)
    {
      // test: we can only release data if the key 'ename' exists in countTable (and dataTable)
      if (!countTable_.isKey(ename,factory) ||
          !dataTable_.isKey(ename,factory)) // why releasing data that has never been stored?
      {
        std::string msg =  "Release: " + ename + " not found. You must first request " + ename;
        throw(Exceptions::RuntimeError(msg));
      }

      // decrement reference counter
      DecrementCounter(ename,factory);

      // desallocation if counter gets zero
      int numReq = -1;
      countTable_.Get<int>(ename,numReq,factory);
      if(numReq==-1) throw(Exceptions::RuntimeError("Release: error reading countTable_(" + ename + ")"));
      if (numReq == 0) // todo: handle keepAll
      {
        countTable_.Remove(ename,factory);
        dataTable_.Remove(ename,factory);
      }
    } //Release

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A")
    template <class T>
    T & GetData(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null) {
        return GetData<T>(ename,factory.get());
    }

    template <class T>
    void GetData(const std::string& ename, T &value, RCP<const FactoryBase> factory = Teuchos::null) {
        GetData<T>(ename,value,factory.get());
    }

    template <class T>
    void GetData(const std::string& ename, T &value, const FactoryBase* factory) {
        dataTable_.Get<T>(ename,value,factory);
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A")
    template <class T>
    T & GetData(const std::string& ename, const FactoryBase* factory)
    {
      if (dataTable_.isKey(ename,factory))
      {
        return dataTable_.Get<T>(ename,factory);
      }
      else
        throw(Exceptions::RuntimeError("Get: " + ename + " not found in dataTable_"));
    }

    //@}

    //! @name Access methods (without reference counter mechanism)
    //@{

    ///! keep variable 'ename' generated by 'factory'
    void Keep(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null)
    {
      countTable_.Set(ename,-1,factory.get());
    }

    ///! keep variable 'ename' generated by 'factory'
    void Keep(const std::string& ename, const FactoryBase* factory)
    {
      countTable_.Set(ename,-1,factory);
    }

    ///! returns true, if 'ename' generated by 'factory' is marked to be kept
    bool isKept(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null)
    {
      if (!countTable_.isKey(ename,factory.get())) return false;
      return countTable_.Get<int>(ename,factory.get()) == -1;
    }

    ///! returns true, if 'ename' generated by 'factory' is marked to be kept
    bool isKept(const std::string& ename, const FactoryBase* factory)
    {
      if (!countTable_.isKey(ename,factory)) return false;
      return countTable_.Get<int>(ename,factory) == -1;
    }

    void Delete(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null)
    {
      if (!countTable_.isKey(ename,factory.get())) return; // data not available?
      if (countTable_.Get<int>(ename,factory.get()) != -1)
      {
        std::cout << "Needs::Delete(): This method is intended to be used when the automatic garbage collector is disabled. Use Release instead to decrement the reference counter!" << std::endl;
      }

      countTable_.Remove(ename,factory.get());
      if (dataTable_.isKey(ename,factory.get()))
          dataTable_.Remove(ename,factory.get());
    }

    void Delete(const std::string& ename, const FactoryBase* factory)
    {
      if (!countTable_.isKey(ename,factory)) return; // data not available?
      if (countTable_.Get<int>(ename,factory) != -1)
      {
        std::cout << "Needs::Delete(): This method is intended to be used when the automatic garbage collector is disabled. Use Release instead to decrement the reference counter!" << std::endl;
      }

      countTable_.Remove(ename,factory);
      if (dataTable_.isKey(ename,factory))
          dataTable_.Remove(ename,factory);
    }

    //@}

    //! @name Utilities.
    //@{
    //! NeedsTable raw print
    //FIXME not covered by unit test right now
    void RawPrint(std::ostream &os, Needs const &foo) {
      //std::cout << "name | #requests" << std::endl;
      //std::cout << "=================" << std::endl;
      //std::cout << foo.countTable_ << std::endl << std::endl;
      //std::cout << "name | value" << std::endl;
      //std::cout << "============" << std::endl;
      //std::cout << foo.dataTable_ << std::endl;
    } //RawPrint

    void Print(std::ostream &out)
    {
      Teuchos::TabularOutputter outputter(out);
      outputter.pushFieldSpec("name", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT,Teuchos::TabularOutputter::GENERAL,12);
      outputter.pushFieldSpec("gen. factory addr.", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
      outputter.pushFieldSpec("req", Teuchos::TabularOutputter::INT,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
      outputter.outputHeader();

      std::vector<std::string> ekeys = countTable_.keys();
      for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++)
      {
        std::vector<const MueLu::FactoryBase*> ehandles = countTable_.handles(*it);
        for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++)
        {
          outputter.outputField(*it);
          outputter.outputField(*kt);
          int reqcount = countTable_.Get<int>(*it,*kt);
          outputter.outputField(reqcount);  // TODO: fix me
          outputter.nextRow();
        }
      }
    }

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string ename, RCP<const FactoryBase> factory = Teuchos::null) {
      return IsAvailable(ename,factory.get());
    }

    //! Test whether a need's value has been saved.
    bool IsAvailable(const std::string ename, const FactoryBase* factory) {
      return dataTable_.isKey(ename,factory);
    }


    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename, RCP<const FactoryBase> factory = Teuchos::null) {
      return IsRequested(ename,factory.get());
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename, const FactoryBase* factory) {
      return countTable_.isKey(ename,factory);
    }

    /*! @brief Return the number of outstanding requests for a need.

    Throws a <tt>Exceptions::RuntimeError</tt> exception if the need either hasn't been requested or
    hasn't been saved.
    */
    int NumRequests(const std::string ename, RCP<const FactoryBase> factory = Teuchos::null) {
    	return NumRequests(ename,factory.get());
    } //NumRequests

    /*! @brief Return the number of outstanding requests for a need.

    Throws a <tt>Exceptions::RuntimeError</tt> exception if the need either hasn't been requested or
    hasn't been saved.
    */
    int NumRequests(const std::string ename, const FactoryBase* factory) {
      //FIXME should we return 0 instead of throwing an exception?

      if (!countTable_.isKey(ename,factory))
      {
            std::string msg =  "NumRequests: " + ename + " not found in countTable_";
            throw(Exceptions::RuntimeError(msg));
      }

      int numReq = -1;
      countTable_.Get<int>(ename,numReq,factory);
      if(numReq==-1) throw(Exceptions::RuntimeError("NumRequests: number of requests for " + ename + " cannot be determined. Error."));
      return numReq;
    } //

    bool SetupPhase(bool bSetup)
    {
      setupPhase_ = bSetup;

      if(!setupPhase_)
      {
        // desallocation of data, that has not been requested
        std::vector<string> dataKeys = dataTable_.keys();
        for(std::vector<string>::iterator it = dataKeys.begin(); it!=dataKeys.end(); ++it)
        {
          for(std::vector<const MueLu::FactoryBase*>::iterator ith = dataTable_.handles(*it).begin(); ith!=dataTable_.handles(*it).end(); ++ith)
          {
            RCP<const MueLu::FactoryBase> rcpFacBase = Teuchos::rcp(*ith);
            if(!countTable_.isKey(*it,/*rcpFacBase*/*ith) || countTable_.Get<int>(*it,/*rcpFacBase*/*ith) == 0)
            {
              std::cout << "Warning: SetupPhaseDesalloc: " << *it << " is present at the Level structure but has not been requested (-> removed)" << std::endl;
              dataTable_.Remove(*it,/*rcpFacBase*/ *ith);
              countTable_.Remove(*it,/*rcpFacBase*/ *ith);
            }
          }
        }
      }

      return setupPhase_;
    }

    //@}

  private:

    //! @name Helper functions
    //@{

    /*! @brief increments counter for variable <tt>ename</tt> (generated by the given factory <tt>factory</tt>)
     * */
    void IncrementCounter(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null)
    {
      IncrementCounter(ename,factory.get());
    }

    /*! @brief increments counter for variable <tt>ename</tt> (generated by the given factory <tt>factory</tt>)
     * */
    void IncrementCounter(const std::string& ename, const FactoryBase* factory)
    {
      int currentCount = -2;
      countTable_.Get<int>(ename,currentCount,factory);
      if(currentCount == -2) std::cout << "ERROR in IncrementCounter " << std::endl;
      if(currentCount != -1) // counter not disabled (no debug mode)
      {
        countTable_.Set(ename,++currentCount,factory);
      }
    }

    /*! @brief decrements counter for variable <tt>ename</tt> (generated by the given factory <tt>factory</tt>)
     * */
    void DecrementCounter(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null)
    {
      DecrementCounter(ename,factory.get());
    }

    /*! @brief decrements counter for variable <tt>ename</tt> (generated by the given factory <tt>factory</tt>)
     * */
    void DecrementCounter(const std::string& ename, const FactoryBase* factory)
    {
        int currentCount = -2;
        countTable_.Get<int>(ename,currentCount,factory);
        if(currentCount == -2) std::cout << "ERROR in DecrementCounter " << std::endl;
        if(currentCount != -1) // counter not disabled (no debug mode)
        {
          countTable_.Set(ename,--currentCount,factory);
        }
    }

    //@}

  }; //class Needs

  /*
  //! NeedsTable pretty print
  std::ostream& operator<<(std::ostream &os, Needs const &foo) {
  std::cout << "name  |  #requests  |  value" << std::endl;
  std::cout << "============================" << std::endl;
  for (Teuchos::ParameterList::ConstIterator param=foo.countTable_.begin(); param!=foo.countTable_.end() ; param++) {
  const std::string pname = foo.countTable_.name(param);
  const int& numRequests = foo.countTable_.get<int>(pname);
  const Teuchos::ParameterEntry &someEntry = foo.dataTable_.getEntry(pname);
  std::cout << pname << " | " << numRequests << " | " << someEntry << std::endl;
  }
  return os;
  }
  */

} //namespace MueLu

#define MUELU_NEEDS_SHORT

#endif //ifndef MUELU_NEEDS_HPP
