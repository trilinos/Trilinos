#ifndef MUELU_BASECLASS_HPP
#define MUELU_BASECLASS_HPP

#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_Describable.hpp>
#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

  //! Base class for MueLu classes
  class BaseClass
    : public Teuchos::VerboseObject<BaseClass>, Teuchos::Describable
  {
    
  public:
    
    //! Destructor.
    virtual ~BaseClass() {}

    //! Return a simple one-line description of this object.
    virtual std::string description() const {
      return "BaseClass::description() not implemented for this object.";
    }

    //! Print the object with some verbosity level to an FancyOStream object. 
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {
      if (verbLevel != Teuchos::VERB_NONE)
        out << "BaseClass::describe() not implemented for this object." << std::endl;
    }

  }; // class BaseClass

} // namespace MueLu

#define MUELU_BASECLASS_SHORT
#endif // ifndef MUELU_BASECLASS_HPP

//TODO: remove default implementation to force developpers to write helpful description/describe functions
