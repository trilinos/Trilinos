#ifndef MUELU_FACTORYMANAGERBASE_HPP
#define MUELU_FACTORYMANAGERBASE_HPP

#include <string>
#include <Teuchos_Hashtable.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
//#include "MueLu_Exceptions.hpp"

namespace MueLu {
  class FactoryBase;

  //! Class that provides default factories within Needs class.
  class FactoryManagerBase : public BaseClass {

  public:
    //@{ Constructors/Destructors.

    //! Destructor.
    virtual ~FactoryManagerBase() { }

    //@}

    //@{ Get/Set functions.

    //! Get
    // Return ref because user also give ref to the Hierarchy.
    // Factory freed at the end of FillHierarchy() //->TODO
    virtual const FactoryBase & GetFactory(const std::string & varName) const {
      // TODO: try/catch + better exception msg if not found
      return *factoryTable_.get(varName);
    }    

    void SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory) const { //TODO: remove const, remame SetFactory()
      // TODO: if (varName already exist) ...
      factoryTable_.put(varName, factory);
    }

    bool IsAvailable(const std::string & varName) const {
      return factoryTable_.containsKey(varName);
    }

    //@}

  protected:
    mutable Teuchos::Hashtable<std::string, RCP<const FactoryBase> > factoryTable_; //TODO: use std lib hashtable instead (Teuchos::Hashtable is deprecated).
        
  }; // class FactoryManagerBase

} // namespace MueLu

#define MUELU_FACTORYMANAGERBASE_SHORT
#endif //ifndef MUELU_FACTORYMANAGERBASE_HPP

//TODO: factoryTable_ must be cleaned at the end of hierarchy Populate() (because Hierarchy is not holding any factories after construction)

