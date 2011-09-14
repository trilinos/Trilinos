#ifndef MUELU_BASECLASS_HPP
#define MUELU_BASECLASS_HPP

#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_Describable.hpp>
#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

  //! Base class for MueLu classes
  class BaseClass
    : public Teuchos::VerboseObject<BaseClass>, public Teuchos::Describable
  {
    
  public:
    
    //! Destructor.
    virtual ~BaseClass() {}

  }; // class BaseClass

} // namespace MueLu

#define MUELU_BASECLASS_SHORT
#endif // ifndef MUELU_BASECLASS_HPP

//TODO: remove default implementation to force developpers to write helpful description/describe functions
