#ifndef XPETRA_PARAMETERS_HPP
#define XPETRA_PARAMETERS_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Xpetra_Map.hpp> // for UnderlyingLib definition
#include <Xpetra_Utils.hpp> // for toString(lib_)

namespace Xpetra {
  
  class Parameters
    : public Teuchos::VerboseObject<Parameters>, public Teuchos::Describable
  {

  public:
    
    Parameters(Teuchos::CommandLineProcessor& clp) {
      setCLP(clp);
    }
      
    void setCLP(Teuchos::CommandLineProcessor& clp) {
      int nOptions=0;                                  // Gives the number of possible option values to select 
      const int maxOptions=2;                          // No ore than 2 libraries are supported right now
      Xpetra::UnderlyingLib optionValues[maxOptions]; // Array that gives the numeric values for each option.
      const char*            optionNames [maxOptions]; // Array that gives the name used in the commandline for each option.

      std::stringstream documentation; // documentation for the option
      documentation << "linear algebra library (0=Epetra, 1=Tpetra)";

      // Default is Tpetra if available. If not, default is Epetra
#if defined(HAVE_XPETRA_EPETRA)
      lib_ = Xpetra::UseEpetra; // set default (if Tpetra support is missing)
      optionValues[nOptions] = Xpetra::UseEpetra;
      //optionValues[nOptions] = "epetra"; //TODO: do not break compatibility right now
      optionNames[nOptions] = "0";            
      nOptions++;
#endif
#if defined(HAVE_XPETRA_TPETRA)
      lib_ = Xpetra::UseTpetra; // set default
      optionValues[nOptions] = Xpetra::UseTpetra;
      //optionsValues[nOptions] = "tpetra"; //TODO: do not break compatibility right now
      optionNames[nOptions] = "1";
      nOptions++;
#endif
        
      clp.setOption<Xpetra::UnderlyingLib>("linAlgebra", &lib_, nOptions, optionValues, optionNames, documentation.str().c_str());

    }

    void check() const {
      //TODO with ifdef...
    }
      
    Xpetra::UnderlyingLib GetLib() const {
      check();
      return lib_;
    }
   
    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const {
      std::ostringstream out;
      out << Teuchos::Describable::description();
      out << "{lib = "  << toString(lib_) << "} ";
      return out.str();
    }
    
    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const {
      using std::endl;
      int vl = (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;
      if (vl == Teuchos::VERB_NONE) return;
      
      if (vl == Teuchos::VERB_LOW) { out << description() << endl; } else { out << Teuchos::Describable::description() << endl; }
      
      if (vl == Teuchos::VERB_MEDIUM || vl == Teuchos::VERB_HIGH || vl == Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab1(out);
        out << "Linear algebra library: " << toString(lib_) << endl;
      }
    }

    //@}
    
  private:
    Xpetra::UnderlyingLib lib_;
  };
  
}

#endif
