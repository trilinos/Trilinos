#ifndef MUELU_PFACTORY_HPP
#define MUELU_PFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {

  class Level;

/*!
  @class PFactory
  @brief Factory that provides an interface for a concrete implementation of a prolongation operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

class PFactory : public TwoLevelFactoryBase {

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PFactory() { }

    //! Destructor.
    virtual ~PFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    virtual bool BuildP(Level &fineLevel, Level &coarseLevel) const = 0;
    //@}

}; //class PFactory

} //namespace MueLu

#define MUELU_PFACTORY_SHORT

#endif //ifndef MUELU_PFACTORY_HPP
