#ifndef MUELU_RFACTORY_HPP
#define MUELU_RFACTORY_HPP

#include "MueLu_OperatorFactory.hpp"
#include "MueLu_Exceptions.hpp"

#include <iostream>

namespace MueLu {

/*!
  @class RFactory
  @brief Factory that provides an interface for a concrete implementation of a restriction operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class RFactory : public OperatorFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
#include "MueLu_UseShortNames.hpp"

  private:

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    RFactory() : reUseGraph_(false), reUseAggregates_(false)
    {
      Teuchos::OSTab tab(this->out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "RFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~RFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    bool BuildR(Level &fineLevel, Level &coarseLevel) = 0;
    //@}

}; //class RFactory

} //namespace MueLu

#define MUELU_RFACTORY_SHORT

#endif //ifndef MUELU_RFACTORY_HPP
