#ifndef MUELU_SMOOTHERFACTORY_HPP
#define MUELU_SMOOTHERFACTORY_HPP

#include <iostream>

#include "MueLu_Needs.hpp"
#include "MueLu_Smoother.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

/*!
  @class Smoother factory base class.
  @brief Base class for smoother factories.
*/

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class SmootherFactory : public Needs {

#include "MueLu_UseShortNames.hpp"

  public:
    //@{ Constructors/Destructors.
    SmootherFactory() {Teuchos::OSTab tab(this->out_); *(this->out_) << "Instantiating a new SmootherFactory" << std::endl;}

    virtual ~SmootherFactory() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build pre-smoother and/or post-smoother
    bool Build(Level &level /*,Teuchos::ParameterList Specs*/) {
      Teuchos::OSTab tab(this->out_);
      *(this->out_) << "Building pre-smoother and/or post-smoother" << std::endl;
      return true;
    }

    //@}

}; //class SmootherFactory

} //namespace MueLu

#define MUELU_SMOOTHERFACTORY_SHORT

#endif //ifndef MUELU_SMOOTHERFACTORY_HPP
