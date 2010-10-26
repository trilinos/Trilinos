#ifndef MUELU_SMOOTHERFACTORY_HPP
#define MUELU_SMOOTHERFACTORY_HPP

#include <iostream>

#include "MueLu_BaseFactory.hpp"
#include "MueLu_Smoother.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

/*!
  @class Smoother factory base class.
  @brief Base class for smoother factories.
*/

template <class Scalar,class LO,class GO,class Node>
class SmootherFactory : public BaseFactory {

  typedef MueLu::Smoother<Scalar,LO,GO,Node> Smoother;

  public:
    //@{ Constructors/Destructors.
    SmootherFactory() {Teuchos::OSTab tab(this->out_); *(this->out_) << "Instantiating a new SmootherFactory" << std::endl;}

    virtual ~SmootherFactory() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build pre-smoother and/or post-smoother
    bool Build(Level<Scalar,LO,GO,Node> &level /*,Teuchos::ParameterList Specs*/) {
      Teuchos::OSTab tab(this->out_);
      *(this->out_) << "Building pre-smoother and/or post-smoother" << std::endl;
      return true;
    }

    //@}

}; //class SmootherFactory

} //namespace MueLu

#endif //ifndef MUELU_SMOOTHERFACTORY_HPP
