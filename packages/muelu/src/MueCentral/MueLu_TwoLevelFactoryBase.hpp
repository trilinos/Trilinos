#ifndef MUELU_TWOLEVELFACTORY_HPP
#define MUELU_TWOLEVELFACTORY_HPP

#include "Teuchos_VerboseObject.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  /*!
    @class Base class for factories (e.g., R, P, and A_coarse).
    @brief Base class for factories that use two levels (fineLevel and coarseLevel).
  */

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class TwoLevelFactoryBase : public Teuchos::VerboseObject<TwoLevelFactoryBase<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps> > {

#include "MueLu_UseShortNames.hpp"
    
  private:

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    TwoLevelFactoryBase() {}

    //! Destructor.
    virtual ~TwoLevelFactoryBase() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    virtual bool Build(Level & fineLevel, Level & coarseLevel) const = 0;

    //@}

  }; //class TwoLevelFactoryBase

} //namespace MueLu

#define MUELU_TWOLEVELFACTORY_SHORT
#endif //ifndef MUELU_TWOLEVELFACTORY_HPP
