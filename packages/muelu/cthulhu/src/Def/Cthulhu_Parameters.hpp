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
      // default
#if defined(HAVE_CTHULHU_TPETRA)
      lib_ = Cthulhu::UseTpetra;
#elif defined(HAVE_CTHULHU_EPETRA)
      lib_ = Cthulhu::UseEpetra;
#else
      //throw error ?
#endif

      // add option --linAlgebra==
#if defined(HAVE_CTHULHU_EPETRA) && defined (HAVE_CTHULHU_TPETRA)
      std::stringstream description;
      description << "use Tpetra (==" << Cthulhu::UseTpetra << ") or Epetra (==" << Cthulhu::UseEpetra << ")";
      clp.setOption("linAlgebra",reinterpret_cast<int*>(&lib_),description.str().c_str()); //TODO: use templated method setOption<>
#endif
    }

    void check() {
      //TODO with ifdef...
    }
      
    Cthulhu::UnderlyingLib GetLib() {
      check();
      return lib_;
    }
     
    void print() { // TODO: Teuchos::Describale 
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
