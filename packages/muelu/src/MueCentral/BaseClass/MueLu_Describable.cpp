#ifndef MUELU_DESCRIBABLE_HPP
#define MUELU_DESCRIBABLE_HPP

#include "MueLu_Describable.hpp"

namespace MueLu {

    Describable::~Describable() { }

    void Describable::describe(Teuchos::FancyOStream &out_arg, const VerbLevel verbLevel) const {
      RCP<Teuchos::FancyOStream> out = rcp(&out_arg,false); //JG: no idea why we have to do that, but it's how Teuchos::Describable::describe() is implemented
      Teuchos::OSTab tab(out);
      *out << this->description() << std::endl;
    }

    std::string Describable::description() const {
      std::string str = Teuchos::Describable::description();

      // remove template parameters
      size_t found = str.find_first_of("<");
      if (found != std::string::npos)
        return str.substr(0, found);
      
      return str;
    }
    
    void Describable::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { describe(out, toMueLuVerbLevel(verbLevel)); }

    std::string Describable::ShortClassName() const {
      std::string str = Teuchos::Describable::description();
      
      // remove template parameters
      {
        size_t found = str.find_first_of("<");
        if (found != std::string::npos)
          str = str.substr(0, found);
      }

      // remove namespace
      {
        size_t found = str.find_last_of(":");
        if (found != std::string::npos)
          str = str.substr(found+1);
      }

      return str;
    }

} // namespace MueLu

#define MUELU_DESCRIBABLE_SHORT
#endif // MUELU_DESCRIBABLE_HPP
