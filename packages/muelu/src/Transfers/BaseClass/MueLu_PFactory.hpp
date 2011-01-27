#ifndef MUELU_PFACTORY_HPP
#define MUELU_PFACTORY_HPP

#include "MueLu_OperatorFactory.hpp"
#include "MueLu_Exceptions.hpp"

#include <iostream>

namespace MueLu {

/*!
  @class PFactory
  @brief Factory that provides an interface for a concrete implementation of a prolongation operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class PFactory : public OperatorFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
#include "MueLu_UseShortNames.hpp"

  private:

     bool reUseGraph_;
     bool reUseAggregates_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PFactory() : reUseGraph_(false), reUseAggregates_(false)
    {
      Teuchos::OSTab tab(this->out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "PFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~PFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    bool BuildP(Level &fineLevel, Level &coarseLevel) = 0;
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
