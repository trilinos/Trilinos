#ifndef MUELU_HIERARCHYFACTORYBASE_HPP
#define MUELU_HIERARCHYFACTORYBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {

  //! Class storing a strategy for building a Hierachy
  class HierarchyFactoryBase : public BaseClass {

  public:
    //@{ Constructors/Destructors.

    //! Destructor.
    virtual ~HierarchyFactoryBase() { }

    //@}

    //@{ Get/Set functions.

    //! Get

    //   virtual RCP<const FactoryManagerBase> GetFactoryManager(int level = 0) const = 0;

    //@}

    //  protected:


  }; // class HierarchyFactoryBase

} // namespace MueLu

#define MUELU_HIERARCHYFACTORYBASE_SHORT
#endif //ifndef MUELU_HIERARCHYFACTORYBASE_HPP
