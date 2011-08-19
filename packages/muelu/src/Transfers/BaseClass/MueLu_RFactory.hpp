#ifndef MUELU_RFACTORY_HPP
#define MUELU_RFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp" //TODO: inheritence of RFactory
#include "MueLu_Exceptions.hpp"

#define MueLu_cout(minimumVerbLevel) \
    if (this->getVerbLevel() >= minimumVerbLevel) *(this->getOStream())

#include <iostream>

namespace MueLu {

/*!
  @class RFactory
  @brief Factory that provides an interface for a concrete implementation of a restriction operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class RFactory : public Teuchos::VerboseObject<RFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > {
#include "MueLu_UseShortNames.hpp"

  private:

  protected:
     RCP<Teuchos::FancyOStream> out_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    RFactory() : out_(this->getOStream())
    {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "RFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~RFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    virtual bool BuildR(Level &fineLevel, Level &coarseLevel) = 0;
    //@}

}; //class RFactory

} //namespace MueLu

#define MUELU_RFACTORY_SHORT

#endif //ifndef MUELU_RFACTORY_HPP
