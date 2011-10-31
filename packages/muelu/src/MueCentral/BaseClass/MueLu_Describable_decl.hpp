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
    virtual ~Describable() ;

    //! @name MueLu Describe
    //@{

    virtual void describe(Teuchos::FancyOStream &out_arg, const VerbLevel verbLevel = Default) const ;

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    virtual std::string description() const ;
    
    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const ;

    //@}

  }; // class Describable

} // namespace MueLu

#define MUELU_DESCRIBABLE_SHORT
#endif // ifndef MUELU_DESCRIBABLE_HPP
