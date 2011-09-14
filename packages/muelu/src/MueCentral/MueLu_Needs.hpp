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

    The data is stored using a variable name and a pointer to the generating factory. The
    pointer to the generating factory is only used as "key". For the Needs class it doesn't
    matter if the given factory pointer is valid or not. So the NULL pointer can also be used.

    A reference counter keeps track of the storage and automatically frees the memory if
    the data is not needed any more. In the standard mode, the data first has to be requested
    by calling the Request function. Then the data can be set by calling SetData.
    With GetData the user can fetch data when needed. Release decrements the reference counter
    for the given variable.
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

  protected:
    RCP<Teuchos::FancyOStream> out_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Default constructor.
    Needs() : out_(this->getOStream()) {
    }

    virtual ~Needs() {}
    //@}

    //! @name SetData
    //! @brief functions for setting data in data storage
    //@{

    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void SetData(const std::string& ename, const T &entry, RCP<const FactoryBase> factory = Teuchos::null) {
        SetData(ename,entry,factory.get());
    } //Set

    //! Store need label and its associated data. This does not increment the storage counter.
    template <class T>
    void SetData(const std::string& ename, const T &entry, const FactoryBase* factory) {
        // check if data is requested
        if(countTable_.isKey(ename, factory) && countTable_.Get<int>(ename, factory) != 0)
        {
          if(!countTable_.isKey(ename, factory))
            countTable_.Set(ename,0,factory);    // make sure that 'ename' is counted
          dataTable_.Set(ename,entry,factory);
        }
    } // Set

    //@}

    //! @name Request/Release data

    //@{

    //! Indicate that an object is needed. This increments the storage counter.
    virtual void Request(const std::string& ename, const FactoryBase* factory) {
      // if it's the first request for 'ename', create a new key in the hashtable
      if(!countTable_.isKey(ename, factory))
        countTable_.Set(ename,0,factory);

      // increment counter
      IncrementCounter(ename,factory);
    } //Request

    //! Decrement the storage counter.
    virtual void Release(const std::string& ename, const FactoryBase* factory)
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
      int numReq = -2;
      countTable_.Get<int>(ename,numReq,factory);
      if(numReq==-2) throw(Exceptions::RuntimeError("Release: error reading countTable_(" + ename + ")"));
      if (numReq == 0) // todo: handle keepAll
      {
        countTable_.Remove(ename,factory);
        dataTable_.Remove(ename,factory);
      }
    } //Release

    //@}

    //! @name Get data
    //@{

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A", factory)
    template <class T>
    T & GetData(const std::string& ename, RCP<const FactoryBase> factory = Teuchos::null) {
        return GetData<T>(ename,factory.get());
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A", A, factory)
    template <class T>
    void GetData(const std::string& ename, T &value, RCP<const FactoryBase> factory = Teuchos::null) {
        GetData<T>(ename,value,factory.get());
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A", A, factoryPtr)
    template <class T>
    void GetData(const std::string& ename, T &value, const FactoryBase* factory) {
        dataTable_.Get<T>(ename,value,factory);
    }

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access). */
    // Usage: Level->Get< RCP<Operator> >("A", factoryPtr)
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

    //! @name Permanent storage
    //@{

    /*! @brief Keep variable 'ename' generated by 'factory'. This variable is not handled by the internal reference counter. */
    virtual void Keep(const std::string& ename, const FactoryBase* factory)
    {
      countTable_.Set(ename,-1,factory);
    }

    /*! @brief Keep variable 'ename' generated by 'factory'. This variable is not handled by the internal reference counter. */
    virtual void Keep(const std::string& ename)
    {
      countTable_.Set(ename,-1,NULL);
    }

    /*! @brief returns true, if variable 'ename' generated by 'factory' is permanently stored */
    virtual bool isKept(const std::string& ename, const FactoryBase* factory) const
    {
      //if (!countTable_.isKey(ename,factory)) return false;
      int nCnt = 0;
      countTable_.Get<int>(ename,nCnt,factory);
      if (nCnt == -1) return true;
      return false;
    }

    /*! @brief returns true, if variable 'ename' generated by 'factory' is permanently stored */
    void isKept(const std::string& ename, bool& bKept, const FactoryBase* factory) const
    {
      //if (!countTable_.isKey(ename,factory)) bKept = false;
      int nCnt = 0;
      countTable_.Get<int>(ename,nCnt,factory);
      if(nCnt == -1) bKept = true;
      else bKept = false;
    }

    /*! @brief remove the permanently stored variable 'ename' generated by 'factory' */
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

    /*! @brief remove the permanently stored variable 'ename' generated by 'factory' */
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

    //! Test whether a need's value has been saved.
    virtual bool IsAvailable(const std::string ename, const FactoryBase* factory) {
      return dataTable_.isKey(ename,factory);
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    virtual bool IsRequested(const std::string ename, const FactoryBase* factory) {
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

    //! @name I/O Functions
    //@{

    /*! \brief Printing method for Needs class.*/
    std::ostream& print(std::ostream& os) const
    {
        Teuchos::TabularOutputter outputter(os);
        outputter.pushFieldSpec("name", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT,Teuchos::TabularOutputter::GENERAL,32);
        outputter.pushFieldSpec("gen. factory addr.", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 18);
        outputter.pushFieldSpec("req", Teuchos::TabularOutputter::INT,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 3);
        outputter.pushFieldSpec("type", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 10);
        outputter.pushFieldSpec("data", Teuchos::TabularOutputter::STRING,Teuchos::TabularOutputter::LEFT, Teuchos::TabularOutputter::GENERAL, 20);
        outputter.outputHeader();

        std::vector<std::string> ekeys = countTable_.keys();
        for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++)
        {
            std::vector<const MueLu::FactoryBase*> ehandles = countTable_.handles(*it);
            for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++)
            {
                outputter.outputField(*it);   // variable name
                outputter.outputField(*kt);   // factory ptr
                int reqcount = 0;             // request counter
                countTable_.Get<int>(*it,reqcount,*kt);
                outputter.outputField(reqcount);
                                              // variable type
                std::string strType = dataTable_.GetType(*it,*kt);
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
                    int data = 0; dataTable_.Get<int>(*it,data,*kt);
                    outputter.outputField(data);
                }
                else if(strType == "double")
                {
                    outputter.outputField(strType);
                    double data = 0.0; dataTable_.Get<double>(*it,data,*kt);
                    outputter.outputField(data);
                }
                else if(strType == "string")
                {
                    outputter.outputField(strType);
                    std::string data = ""; dataTable_.Get<std::string>(*it,data,*kt);
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



    //! @name Helper functions
    //@{

    /*bool SetupPhase(bool bSetup)
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
            if(!countTable_.isKey(*it,*ith) || countTable_.Get<int>(*it,*ith) == 0)
            {
              std::cout << "Warning: SetupPhaseDesalloc: " << *it << " is present at the Level structure but has not been requested (-> removed)" << std::endl;
              dataTable_.Remove(*it, *ith);
              countTable_.Remove(*it, *ith);
            }
          }
        }

      }

      return setupPhase_;
    }*/

    void copyKeepStatus(const RCP<Needs>& newNeeds) const
    {
      std::vector<std::string> ekeys = countTable_.keys();
      for (std::vector<std::string>::iterator it = ekeys.begin(); it != ekeys.end(); it++)
      {
          std::vector<const MueLu::FactoryBase*> ehandles = countTable_.handles(*it);
          for (std::vector<const MueLu::FactoryBase*>::iterator kt = ehandles.begin(); kt != ehandles.end(); kt++)
          {
            const std::string ename = *it;
            const MueLu::FactoryBase* fac = *kt;
            if(isKept(ename,fac))
            {
              if(fac == NULL) newNeeds->Keep(ename);
              else newNeeds->Keep(ename,fac);
            }
          }
      }
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
        if(currentCount != -1 && currentCount != 0) // counter not disabled (no debug mode)
        {
          countTable_.Set(ename,--currentCount,factory);
        }
    }

    //@}

  }; //class Needs


  //! NeedsTable pretty print
  /*std::ostream& operator<<(std::ostream &os, Needs const &foo) {
  std::cout << "name  |  #requests  |  value" << std::endl;
  std::cout << "============================" << std::endl;
  for (Teuchos::ParameterList::ConstIterator param=foo.countTable_.begin(); param!=foo.countTable_.end() ; param++) {
  const std::string pname = foo.countTable_.name(param);
  const int& numRequests = foo.countTable_.get<int>(pname);
  const Teuchos::ParameterEntry &someEntry = foo.dataTable_.getEntry(pname);
  std::cout << pname << " | " << numRequests << " | " << someEntry << std::endl;
  }*/
  //foo.print(os);
  //return os;
  //}


} //namespace MueLu

#define MUELU_NEEDS_SHORT

#endif //ifndef MUELU_NEEDS_HPP
