#ifndef MUELU_FACTORYMANAGERBASE_HPP
#define MUELU_FACTORYMANAGERBASE_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

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
    const virtual RCP<const FactoryBase> & GetFactory(const std::string & varName) const = 0;
    //@}

    // Free temporarily hold data at the end of Hierarchy::Setup()
    // This method is const because the clean concerns only mutable data.
    virtual void Clean() const { }

  }; // class FactoryManagerBase

} // namespace MueLu

#define MUELU_FACTORYMANAGERBASE_SHORT
#endif //ifndef MUELU_FACTORYMANAGERBASE_HPP
