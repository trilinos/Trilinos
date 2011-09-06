#ifndef MUELU_SMOOTHERPROTOTYPE_HPP
#define MUELU_SMOOTHERPROTOTYPE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

  class Level;

  /*!
    @class AdvSmootherPrototype

    'Advanced Smoother prototypes' can be fully copied using the Copy() method.
    They can also copy the parameters of
    another smoother object of the same type (CopyParameters()). Both
    capabilities are used by the SmootherFactory.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class AdvSmootherPrototype : public SmootherPrototypex<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  public:
    //!@nameConstructors/Destructors.
    //@{
    AdvSmootherPrototype() 
      : type_("undefined")
    {}

    virtual ~AdvSmootherPrototype() {}
    //@}

    //! @name Build methods.
    //@{

    virtual void CopyParameters(const AdvSmootherPrototype& smootherPrototype) = 0;

    //@}

    //! @name Get/Set methods.
    //@{

    //! Get the smoother type.
    std::string GetType() const { return type_; }

    /*! @brief Set the smoother type.
      This method must be called by constructors of derived classes.
    */
    //TODO: remove, type_ should be const
    void SetType(std::string & type) { type_ = type; }
  
    //@}

  private:
    std::string type_;

  }; //class AdvSmootherPrototype

} //namespace MueLu

#define MUELU_SMOOTHERPROTOTYPE_SHORT

#endif //ifndef MUELU_SMOOTHERPROTOTYPE_HPP
