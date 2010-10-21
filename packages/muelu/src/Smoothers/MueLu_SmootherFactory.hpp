#ifndef MUELU_SMOOTHERFACTORY_HPP
#define MUELU_SMOOTHERFACTORY_HPP

#include <iostream>

#include "MueLu_BaseFactory.hpp"
#include "MueLu_Smoother.hpp"
#include "MueLu_Level.hpp"

/*!
  @class Smoother factory base class.
  @brief Base class for smoother factories.
*/

namespace MueLu {

template <class Scalar,class LO,class GO,class Node>
class SmootherFactory : public BaseFactory {

  typedef MueLu::Smoother<Scalar,LO,GO,Node> Smoother;

  private:
    int outputLevel_;
    int priorOutputLevel_;

  public:
    //@{ Constructors/Destructors.
    SmootherFactory() {std::cout << "Instantiating a new SmootherFactory" << std::endl;}

    virtual ~SmootherFactory() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build a smoother
    bool Build(Teuchos::RCP<Smoother> &preSm,
               Teuchos::RCP<Smoother> &postSm,
               Level<Scalar,LO,GO,Node> &level) {
      std::cout << "Building ";
      if (preSm != Teuchos::null) {
        std::cout << "a pre-";
        if (postSm != Teuchos::null)
          std::cout << " and post-";
        std::cout << "smoother" << std::endl;
      }
      else if (postSm != Teuchos::null) std::cout << "a post-smoother." << std::endl;
      else std::cout << "no smoothers." << std::endl;
      return true;
    }

    //@}

/*
//FIXME The needs mechanism hasn't been decided upon.
    void AddNeeds(Teuchos::ParameterList const &newNeeds) {
      CrossFactory.MergeNeeds(newNeeds,Needs_);
    }

    Teuchos::ParameterList& GetNeeds(Teuchos::ParameterList &newNeeds) {
      return Needs_;
    }
*/

    void SetOutputLevel(int outputLevel) {
      outputLevel_ = outputLevel;
    }

    int GetOutputLevel() {
      return outputLevel_;
    }

}; //class SmootherFactory

} //namespace MueLu

#endif //ifndef MUELU_SMOOTHERFACTORY_HPP
