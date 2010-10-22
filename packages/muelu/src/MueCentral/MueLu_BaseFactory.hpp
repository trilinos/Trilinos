#ifndef MUELU_BASEFACTORY_HPP
#define MUELU_BASEFACTORY_HPP

//#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

/*!
  @class Factory base class.
  @brief Base class for factories.
  Maintains just 2 things:
   - Ouput level status
   - A list of 'Needs' for the factory. For example, a restriction factory that transposes the tentative
     prolongator 'Needs' the prolongator factory to save this.
*/

namespace MueLu {

class BaseFactory : public Teuchos::VerboseObject<BaseFactory> {
  private:
    Teuchos::ParameterList Needs_;
    int outputLevel_;
    int priorOutputLevel_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:
    //@{ Constructors/Destructors.
    BaseFactory() : out_(this->getOStream()) {}

    virtual ~BaseFactory() {}
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

}; //class BaseFactory

} //namespace MueLu

#endif //ifndef MUELU_BASEFACTORY_HPP
