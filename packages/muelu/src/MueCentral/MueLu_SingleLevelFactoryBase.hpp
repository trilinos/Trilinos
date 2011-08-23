#ifndef MUELU_SINGLELEVELFACTORY_HPP
#define MUELU_SINGLELEVELFACTORY_HPP

#include "Teuchos_VerboseObject.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  /*!
    @class Base class for factories (e.g., Aggregation, Smoothers).
    @brief Base class for factories that use one levels (currentLevel).
  */

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class SingleLevelFactoryBase : public Teuchos::VerboseObject<SingleLevelFactoryBase<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps> > {

#include "MueLu_UseShortNames.hpp"
    
  private:

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    SingleLevelFactoryBase() {}

    //! Destructor.
    virtual ~SingleLevelFactoryBase() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    virtual bool Build(Level & currentLevel) const = 0;

    //@}

  }; //class SingleLevelFactoryBase

} //namespace MueLu

#define MUELU_SINGLELEVELFACTORY_SHORT
#endif //ifndef MUELU_SINGLELEVELFACTORY_HPP
