#ifndef MUELU_PRFACTORY_HPP
#define MUELU_PRFACTORY_HPP

#include <iostream>

#include "Xpetra_ConfigDefs.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
  @class PRFactory
  @brief Factory that provides an interface for multigrid transfer operators (both restriction and prolongation).

  The user has to implement the Build function. The default implementation is GenericPRFactory.
*/

class PRFactory : public TwoLevelFactoryBase {

  protected:

     Xpetra::global_size_t maxCoarseSize_;
     bool reUseGraph_;
     bool reUseAggregates_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PRFactory() : maxCoarseSize_(50), reUseGraph_(false), reUseAggregates_(false)
    {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "PRFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~PRFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    bool Build(Level &fineLevel, Level &coarseLevel) const = 0;
    //@}

    //! @name Set methods.
    //@{
    void SetMaxCoarseSize(Xpetra::global_size_t const &maxCoarseSize) {
      maxCoarseSize_ = maxCoarseSize;
    }

    void ReUseAggregates(bool const &value) {
      reUseAggregates_ = value;
    }

    void ReUseGraph(bool const &value) {
      reUseGraph_ = value;
    }
    //@}

    //! @name Get methods.
    //@{

  Xpetra::global_size_t GetMaxCoarseSize() const {
      return maxCoarseSize_;
    }

    bool ReUseAggregates() const {
      return reUseAggregates_;
    }

    bool ReUseGraph() const {
      return reUseGraph_;
    }

    //@}

}; //class PRFactory

} //namespace MueLu

#define MUELU_PRFACTORY_SHORT

#endif //ifndef MUELU_PRFACTORY_HPP
