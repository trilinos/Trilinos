#ifndef MUELU_PRFACTORY_HPP
#define MUELU_PRFACTORY_HPP

#include <Xpetra_ConfigDefs.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {

  /*!
    @class PRFactory
    @brief Factory that provides an interface for multigrid transfer operators (both restriction and prolongation).

    The user has to implement the Build function. The default implementation is GenericPRFactory.
  */

  class PRFactory : public TwoLevelFactoryBase {

  protected:

    Xpetra::global_size_t maxCoarseSize_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PRFactory() : maxCoarseSize_(50)
    { }

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

    //@}

    //! @name Get methods.
    //@{
  
    Xpetra::global_size_t GetMaxCoarseSize() const {
      return maxCoarseSize_;
    }

    //@}

  }; //class PRFactory

} //namespace MueLu

#define MUELU_PRFACTORY_SHORT

#endif //ifndef MUELU_PRFACTORY_HPP
