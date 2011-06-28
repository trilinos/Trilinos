#ifndef MUELU_GALLERYPARAMETERS_HPP
#define MUELU_GALLERYPARAMETERS_HPP

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>

#include "MueLu_Exceptions.hpp"

namespace MueLu {
  
  namespace Gallery {
    
    // TODO nx/ny/nz == GO or global_size_t ? But what is the best to do?

    template<typename GO>
    class Parameters
    {

    public:
    
      Parameters(Teuchos::CommandLineProcessor& clp, GO nx=16, GO ny=-1, GO nz=-1, const std::string & matrixType="Laplace1D"): nx_(nx), ny_(ny), nz_(nz), matrixType_(matrixType) {
        clp.setOption("nx", &nx_, "mesh points in x-direction.");
        clp.setOption("ny", &ny_, "mesh points in y-direction.");
        clp.setOption("nz", &nz_, "mesh points in z-direction.");
        clp.setOption("matrixType", &matrixType_, "matrix type: Laplace1D, Laplace2D, Laplace3D"); //TODO: Star2D, numGlobalElements=...
      }
      
      void check() {
        
        //if (nx < 0) ...

      }
      
      GO GetNumGlobalElements() {
        check();

        GO numGlobalElements=-1;
        if (matrixType_ == "Laplace1D")
          numGlobalElements = static_cast<GO>(nx_);
        else if (matrixType_ == "Laplace2D")
          numGlobalElements = static_cast<GO>(nx_*ny_);
        else if (matrixType_ == "Laplace3D")
          numGlobalElements = static_cast<GO>(nx_*ny_*nz_);
        //TODO else throw

        if (numGlobalElements < 0) throw Exceptions::Overflow("Gallery: numGlobalElements < 0 (did you forget --nx for 2D problems?)");

        return numGlobalElements;
      }

      const std::string & GetMatrixType() {
        check();
        return matrixType_;
      }

      Teuchos::ParameterList & GetParameterList() {
        check();
        //TODO: check if already on it...
        paramList.set("nx", static_cast<GO>(nx_));
        paramList.set("ny", static_cast<GO>(ny_));
        paramList.set("nz", static_cast<GO>(nz_));

        return paramList;
      }

      void print() { // TODO: Teuchos::Describale 
        check();
        std::cout << "Gallery parameters: " << std::endl
                  << " * matrix type  = " << matrixType_ << std::endl
                  << " * problem size = " << GetNumGlobalElements() << std::endl
                  << std::endl;
      }

    private:
      Teuchos::ParameterList paramList;

      // See Teuchos BUG 5249: https://software.sandia.gov/bugzilla/show_bug.cgi?id=5249
      double nx_;
      double ny_;
      double nz_;
      // GO nx_;
      // GO ny_;
      // GO nz_;

      std::string matrixType_;
    };
  
  }
}

#endif


//TODO: add capability to read from file ??

//GOAL: link between InputReader and ParameterList + Hide everything from examples
