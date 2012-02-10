#ifndef MUELU_SMOOTHERPROTOTYPEBASE_HPP
#define MUELU_SMOOTHERPROTOTYPEBASE_HPP



#include "MueLu_ConfigDefs.hpp"
//#include "MueLu_Level.hpp"

namespace MueLu {
  class Level;
  
  /*!
    @class SmootherPrototypeBase
    @brief Base class for smoother prototypes

    This has the signature for the required DeclareInput.
  */

  class SmootherPrototypeBase  {
  public:
    //@{ Constructors/Destructors.
    SmootherPrototypeBase() {}

    virtual ~SmootherPrototypeBase() {}
    //@}

  public:
    
    //! @name DeclareInput methods.
    //@{

    //! Declare Input for smoother.
    virtual void DeclareInput(Level &currentLevel) const = 0;

    //@}

  }; //class SmootherPrototypeBase

} //namespace MueLu

#endif //ifndef MUELU_SMOOTHERPROTOTYPEBASE_HPP
