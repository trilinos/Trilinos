#ifndef MUELU_FACTORYBASE_HPP
#define MUELU_FACTORYBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {
  class Level;

  //! Base class for factories (e.g., R, P, and A_coarse).
  class FactoryBase : public BaseClass {

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    FactoryBase() {}

    //! Destructor.
    virtual ~FactoryBase() {}
    //@}

    //@{
    //! @name Build methods.

    virtual bool NewBuild(Level & requestedLevel) const = 0;

    //@}

  }; //class FactoryBase

} //namespace MueLu

#define MUELU_FACTORYBASE_SHORT
#endif //ifndef MUELU_FACTORYBASE_HPP
