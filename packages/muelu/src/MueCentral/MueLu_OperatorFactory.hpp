#ifndef MUELU_OPERATORFACTORY_HPP
#define MUELU_OPERATORFACTORY_HPP

#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "MueLu_BaseFactory.hpp"
#include "MueLu_Level.hpp"

/*!
  @class Base class for operator factories (e.g., R, P, and A_coarse).
  @brief Base class for factories that build operators.
  Very similar to BaseFactory, but with an additional virtual Build method.
*/

namespace MueLu {

template<class Scalar, class LO, class GO, class Node>
class OperatorFactory : public BaseFactory {

  private:
    Teuchos::ParameterList Needs_;
    int outputLevel_;
    int priorOutputLevel_;

  public:
    //@{ Constructors/Destructors.
    OperatorFactory() {}

    virtual ~OperatorFactory() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    virtual bool Build(Level<Scalar, LO, GO, Node> &i,Level<Scalar, LO, GO, Node> &j) = 0;

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

}; //class OperatorFactory

} //namespace MueLu

#endif //ifndef MUELU_OPERATORFACTORY_HPP
