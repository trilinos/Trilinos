#ifndef CTHULHU_PARAMETERS_HPP
#define CTHULHU_PARAMETERS_HPP

#include <Teuchos_CommandLineProcessor.hpp>
#include <Cthulhu_Map.hpp> // for UnderlyingLib definition

namespace Cthulhu {
  
  class Parameters
  {

  public:
    
    Parameters(Teuchos::CommandLineProcessor& clp) {
      setCLP(clp);
    }
      
    void setCLP(Teuchos::CommandLineProcessor& clp) {
      int nOptions=0;                                  // Gives the number of possible option values to select 
      const int maxOptions=2;                          // No ore than 2 libraries are supported right now
      Cthulhu::UnderlyingLib optionValues[maxOptions]; // Array that gives the numeric values for each option.
      const char*            optionNames [maxOptions]; // Array that gives the name used in the commandline for each option.

      std::stringstream documentation; // documentation for the option
      documentation << "linear algebra library (0=Epetra, 1=Tpetra)";

      // Default is Tpetra if available. If not, default is Epetra
#if defined(HAVE_CTHULHU_EPETRA)
      lib_ = Cthulhu::UseEpetra; // set default (if Tpetra support is missing)
      optionValues[nOptions] = Cthulhu::UseEpetra;
      //optionValues[nOptions] = "epetra"; //TODO: do not break compatibility right now
      optionNames[nOptions] = "0";            
      nOptions++;
#endif
#if defined(HAVE_CTHULHU_TPETRA)
      lib_ = Cthulhu::UseTpetra; // set default
      optionValues[nOptions] = Cthulhu::UseTpetra;
      //optionsValues[nOptions] = "tpetra"; //TODO: do not break compatibility right now
      optionNames[nOptions] = "1";
      nOptions++;
#endif
        
      clp.setOption<Cthulhu::UnderlyingLib>("linAlgebra", &lib_, nOptions, optionValues, optionNames, documentation.str().c_str());

    }

    void check() {
      //TODO with ifdef...
    }
      
    Cthulhu::UnderlyingLib GetLib() {
      check();
      return lib_;
    }
     
    void print() { // TODO: Teuchos::Describable 
      check();
     
      std::cout << "Cthulhu parameters: " << std::endl;

      if (lib_ == Cthulhu::UseTpetra)
        std::cout << " * linAlgebra = Tpetra" << std::endl;
      else if (lib_ == Cthulhu::UseEpetra)
        std::cout << " * linAlgebra = Epetra" << std::endl;

      std::cout << std::endl;
    }

  private:
    Cthulhu::UnderlyingLib lib_;
  };
  
}

#endif
