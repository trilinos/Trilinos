#ifndef MUELU_RFACTORY_HPP
#define MUELU_RFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class RFactory
  @brief Factory that provides an interface for a concrete implementation of a restriction operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

class RFactory : public TwoLevelFactoryBase {

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    RFactory() { }

    //! Destructor.
    virtual ~RFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    virtual void BuildR(Level &fineLevel, Level &coarseLevel) const = 0;
    //@}

}; //class RFactory

} //namespace MueLu

#define MUELU_RFACTORY_SHORT

#endif //ifndef MUELU_RFACTORY_HPP
