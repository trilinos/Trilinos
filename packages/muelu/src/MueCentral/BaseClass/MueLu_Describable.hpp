#ifndef MUELU_DESCRIBABLE_HPP
#define MUELU_DESCRIBABLE_HPP

#include <Teuchos_Describable.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerbosityLevel.hpp"

namespace MueLu {

  //! Base class for MueLu classes
  class Describable
    : public Teuchos::Describable
  {
    
  public:
    
    //! Destructor.
    virtual ~Describable() {}

    //! @name MueLu Describe
    //@{

    virtual void describe(Teuchos::FancyOStream &out_arg, const VerbLevel verbLevel = Default) const {
      RCP<Teuchos::FancyOStream> out = rcp(&out_arg,false); //JG: no idea why we have to do that, but it's how Teuchos::Describable::describe() is implemented
      Teuchos::OSTab tab(out);
      *out << this->description() << std::endl;
    }

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    virtual std::string description() const {
      std::string str = Teuchos::Describable::description();

      // remove template parameters
      size_t found = str.find_first_of("<");
      if (found != std::string::npos)
        return str.substr(0, found);
      
      return str;
    }
    
    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const { describe(out, toMueLuVerbLevel(verbLevel)); }

    //@}

  }; // class Describable

} // namespace MueLu

#define MUELU_DESCRIBABLE_SHORT
#endif // ifndef MUELU_DESCRIBABLE_HPP
