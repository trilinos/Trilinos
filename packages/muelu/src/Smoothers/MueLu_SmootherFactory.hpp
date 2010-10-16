#ifndef MUELU_SMOOTHERFACTORY_HPP
#define MUELU_SMOOTHERFACTORY_HPP

#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "MueLu_Level.hpp"

/*!
  @class Smoother factory base class.
  @brief Base class for smoother factories.
*/

namespace MueLu {

template <class Scalar,class LO,class GO,class Node>
class SmootherFactory : public BaseFactory {
  private:
    Teuchos::ParameterList Needs_;
    int outputLevel_;
    int priorOutputLevel_;

  public:
    //@{ Constructors/Destructors.
    SmootherFactory() {}

    virtual ~SmootherFactory() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build a smoother
    bool Build(Teuchos::RCP<Smoother> &preSm,
               Teuchos::RCP<Smoother> &postSm,
               Level<Scalar,LO,GO,Node> &level) {
      std::cout << "Building a pre- and post-smoother." << std::endl;
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

/*
//FIXME The needs mechanism hasn't been decided upon.

    void TempOutputLevel(int outputLevel) {
    }

    void RestoreOutputLevel(int outputLevel) {
    }
*/

}; //class SmootherFactory

} //namespace MueLu

#endif //ifndef MUELU_SMOOTHERFACTORY_HPP
