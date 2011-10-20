#ifndef MUELU_NOFACTORY_HPP
#define MUELU_NOFACTORY_HPP

#include <Teuchos_TestForException.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryBase.hpp"

namespace MueLu {

  class Level;
  
  /*!
    @class NoFactory class.
    @brief NoFactory that is used for data stored in level class for that no generating factory is available/necessary.

    Use Singleton pattern.
  */
  class NoFactory : public FactoryBase {

    //! Constructor.
    NoFactory() { }

  public:

    //! Destructor.
    virtual ~NoFactory() { }

    //! Implementation of FactoryBase interface
    //{@
    
    //!
    void CallBuild(Level & requestedLevel) const {  
      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::NoFactory::Build(): this method cannot be called.");
    }

    //!
    void CallDeclareInput(Level & requestedLevel) const {  }

    //@}
    
    //! Static Get() functions
    //{@
    
    //! 
    static const RCP<const NoFactory> getRCP() {
      if(noFactory_ == Teuchos::null) {
        noFactory_ = rcp(new NoFactory());
      }

      return noFactory_;
    }

    //! 
    static const NoFactory* get() {
      return getRCP().get();
    }

    //@}

  private:
    static RCP<const NoFactory> noFactory_; // static NoFactory instance for user defined "factories"

  }; // class NoFactory

} // namespace MueLu

#endif
