#ifndef MUELU_PFACTORY_HPP
#define MUELU_PFACTORY_HPP

#include "Teuchos_VerboseObject.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"

#define MueLu_cout(minimumVerbLevel) \
    if (this->getVerbLevel() >= minimumVerbLevel) *(this->getOStream())

#include <iostream>

namespace MueLu {

/*!
  @class PFactory
  @brief Factory that provides an interface for a concrete implementation of a prolongation operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

class PFactory : public TwoLevelFactoryBase {

  protected:

     bool reUseGraph_;
     bool reUseAggregates_;
     RCP<Teuchos::FancyOStream> out_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PFactory() : reUseGraph_(false), reUseAggregates_(false), out_(this->getOStream())
    {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "PFactory: Instantiating a new factory" << std::endl;
    }

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

    //! @name Set methods.
    //@{

    void ReUseAggregates(bool const &value) {
      reUseAggregates_ = value;
    }

    void ReUseGraph(bool const &value) {
      reUseGraph_ = value;
    }
    //@}

    //! @name Get methods.
    //@{

    bool ReUseAggregates() const {
      return reUseAggregates_;
    }

    bool ReUseGraph() const {
      return reUseGraph_;
    }

    //@}

}; //class PFactory

} //namespace MueLu

#define MUELU_PFACTORY_SHORT

#endif //ifndef MUELU_PFACTORY_HPP
