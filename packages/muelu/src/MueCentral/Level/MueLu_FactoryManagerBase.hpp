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
    // Factory freed at the end of FillHierarchy() //->TODO
    virtual const FactoryBase & GetFactory(const std::string & varName) const = 0;

    virtual bool IsAvailable(const std::string & varName) const = 0;

    //@}
        
  }; // class FactoryManagerBase

} // namespace MueLu

#define MUELU_FACTORYMANAGERBASE_SHORT
#endif //ifndef MUELU_FACTORYMANAGERBASE_HPP
