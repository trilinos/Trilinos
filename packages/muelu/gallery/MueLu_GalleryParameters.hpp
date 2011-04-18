#ifndef MUELU_GALLERYPARAMETERS_HPP
#define MUELU_GALLERYPARAMETERS_HPP

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>

namespace MueLu {
  
  namespace Gallery {
    
    class Parameters
    {

    public:
    
      Parameters(Teuchos::CommandLineProcessor& cmdp, int nx=16, int ny=-1, int nz=-1, const std::string & matrixType="Laplace1D"): nx_(nx), ny_(ny), nz_(nz), matrixType_(matrixType) {
        cmdp.setOption("nx", &nx_, "mesh points in x-direction.");
        cmdp.setOption("ny", &ny_, "mesh points in y-direction.");
        cmdp.setOption("nz", &nz_, "mesh points in z-direction.");
        cmdp.setOption("matrixType", &matrixType_, "matrix type: Laplace1D, Laplace2D, Laplace3D"); //TODO: Star2D, numGlobalElements=...
      }
      
      void check() {
        
        // if (nx < 0) ...

      }
      
      int GetNumGlobalElements() {
        check();

        int numGlobalElements;
        if (matrixType_ == "Laplace1D")
          numGlobalElements = nx_;
        else if (matrixType_ == "Laplace2D")
          numGlobalElements = nx_*ny_;
        else if (matrixType_ == "Laplace3D")
          numGlobalElements = nx_*ny_*nz_;
        //TODO else throw

        return numGlobalElements;
      }

      const std::string & GetMatrixType() {
        check();
        return matrixType_;
      }

      Teuchos::ParameterList & GetParameterList() {
        check();
        //TODO: check if already on it...
        paramList.set("nx", nx_);
        paramList.set("ny", ny_);
        paramList.set("nz", nz_);

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

      int nx_;
      int ny_;
      int nz_;
      std::string matrixType_;
    };
  
  }
}

#endif


//TODO: add capability to read from file ??

//GOAL: link between InputReader and ParameterList + Hide everything from examples
