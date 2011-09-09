#ifndef MUELU_SMOOTHERPROTOTYPE_HPP
#define MUELU_SMOOTHERPROTOTYPE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherBase.hpp"

namespace MueLu {

  class Level;

  /*!
    @class SmootherPrototype
    @brief Base class for smoother prototypes

    A smoother prototype is a smoother which can be in two states:
    - ready to be duplicated (parameters defined)
    - ready to be used (setup phase completed)

    'Smoother prototypes' can be fully copied using the Copy() method.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SmootherPrototype : public SmootherBase<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  public:
    //!@nameConstructors/Destructors.
    //@{

    SmootherPrototype() : isSetup_(false) {}

    virtual ~SmootherPrototype() {}

    //@}

    //! @name Build methods.
    //@{

    virtual void Setup(Level &) = 0;

    virtual RCP<SmootherPrototype> Copy() const = 0;

    //@}

    //! @name Get/Set methods.
    //@{

    //! Get the state of a smoother prototype.
    bool IsSetup() const { return isSetup_; }
    
    //! Set the state of a smoother prototype.
    // Developpers: this method must be called by our Setup() method.
    void IsSetup(bool const &ToF) { isSetup_ = ToF; }
    
    //@}

  private:
    bool isSetup_;
  
  }; //class SmootherPrototype

} //namespace MueLu

#define MUELU_SMOOTHERPROTOTYPE_SHORT
#endif //ifndef MUELU_SMOOTHERPROTOTYPE_HPP

//TODO: private copy constructor
//TODO: update comments
