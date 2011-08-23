#ifndef MUELU_FACTORY_HPP
#define MUELU_FACTORY_HPP

#include <Teuchos_VerboseObject.hpp>
#include "MueLu_ConfigDefs.hpp"

namespace MueLu {
  class Level;

  //! Base class for factories (e.g., R, P, and A_coarse).
  class FactoryBase : public Teuchos::VerboseObject<FactoryBase> {

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    FactoryBase() {}

    //! Destructor.
    virtual ~FactoryBase() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    //virtual bool Build(Level & fineLevel, Level & coarseLevel) const = 0;

    //@}

  }; //class FactoryBase

} //namespace MueLu

#define MUELU_FACTORY_SHORT
#endif //ifndef MUELU_FACTORY_HPP
