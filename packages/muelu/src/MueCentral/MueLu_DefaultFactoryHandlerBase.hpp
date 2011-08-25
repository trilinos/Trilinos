#ifndef MUELU_DEFAULTFACTORYHANDLERBASE_HPP
#define MUELU_DEFAULTFACTORYHANDLERBASE_HPP

#include <string>
#include <Teuchos_Hashtable.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
//#include "MueLu_Exceptions.hpp"

namespace MueLu {
  class FactoryBase;

  //! Class that provides default factories within Needs class.
  class DefaultFactoryHandlerBase : public BaseClass {

  public:
    //@{ Constructors/Destructors.

    //! Destructor.
    virtual ~DefaultFactoryHandlerBase() { }

    //@}

    //@{ Get/Set functions.

    //! Get
    // Return ref because user also give ref to the Hierarchy.
    // Factory freed at the end of FillHierarchy() //->TODO
    virtual const FactoryBase & GetDefaultFactory(const std::string & varName) {
      // TODO: try/catch + better exception msg if not found
      return *factoryTable_.get(varName);
    }    

    void SetDefaultFactory(const std::string & varName, const RCP<FactoryBase> & factory) {
      // TODO: if (varName already exist) ...
      factoryTable_.put(varName, factory);
    }

    bool IsAvailable(const std::string & varName) {
      return factoryTable_.containsKey(varName);
    }

    //@}

  protected:
    Teuchos::Hashtable<std::string, RCP<FactoryBase> > factoryTable_;
        
  }; // class DefaultFactoryHandlerBase

} // namespace MueLu

#define MUELU_DEFAULTFACTORYHANDLERBASE_SHORT
#endif //ifndef MUELU_DEFAULTFACTORYHANDLERBASE_HPP

//TODO: factoryTable_ must be cleaned at the end of hierarchy Populate() (because Hierarchy is not holding any factories after construction)

