#ifndef MUELU_SMOOTHERBASE_HPP
#define MUELU_SMOOTHERBASE_HPP

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Needs.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  /*!
    @class SmootherBase
    @brief Base class for smoothers

    This has the signature for the required Apply function and contains data that is generic across all
    smoothers.
  */

  template <class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
  class SmootherBase : public Needs {

#include "MueLu_UseShortNames.hpp"

  private:
    std::string Type_;

  public:
    //@{ Constructors/Destructors.
    SmootherBase() {}

    virtual ~SmootherBase() {}
    //@}

    //! @name Apply methods.
    //@{

    //! Apply smoother.
    //virtual void Apply(RCP<MultiVector> x, RCP<MultiVector> const rhs, bool InitialGuessIsZero) = 0;
    virtual void Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero) = 0;

    //@}

    //! @name Set/Get methods.
    //@{

    //! Prototype of the method to setup the number of iteration of the smoother.
    virtual void SetNIts(LO const &Nits) { //TODO: should be removed from the interface
       TEST_FOR_EXCEPTION(1, MueLu::Exceptions::RuntimeError, "SetNIts() not implemented for this smoother");
    };

    //! Get the smoother type.
    std::string GetType() const {
      return Type_;
    }

    /*! @brief Set the smoother type.

    This method must be called by constructors of derived classes.
    */
    void SetType(std::string type) {
      Type_ = type;
    }

    //! @name Utilities.
    //@{

    //! @brief Print information about the smoother
    virtual void Print(std::string prefix) const = 0;

    //@}

  }; //class SmootherBase

} //namespace MueLu

#define MUELU_SMOOTHERBASE_SHORT

#endif //ifndef MUELU_SMOOTHERBASE_HPP
