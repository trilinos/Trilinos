#ifndef MUELU_TEST_HELPERS_H
#define MUELU_TEST_HELPERS_H

#include <iostream>

// Teuchos
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

// Xpetra
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Parameters.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_CrsOperator.hpp"

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Hierarchy.hpp"

// Gallery
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

namespace MueLu {

  namespace TestHelpers {
    
    using Xpetra::global_size_t;

    class Parameters {

    private:
      Parameters() {} // static class
      
    public:
      
      static Xpetra::Parameters xpetraParameters;
      
      inline static RCP<const Teuchos::Comm<int> > getDefaultComm() {
        return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
      }
      
      inline static Xpetra::UnderlyingLib getLib() {
        return MueLu::TestHelpers::Parameters::xpetraParameters.GetLib();
      }
    };

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    class Factory {
#include "MueLu_UseShortNames.hpp"

    private: 
      Factory() {} // static class

    public:

      //
      // Method that creates a map containing a specified number of local elements per process.
      //
      static const RCP<const Map> BuildMap(LO numElementsPerProc) { 

        RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  
        const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  
        return MapFactory::Build(TestHelpers::Parameters::getLib(), INVALID, numElementsPerProc, 0, comm);
  
      } // BuildMap()

      // Create a matrix as specified by parameter list options
      static RCP<Operator> BuildMatrix(Teuchos::ParameterList &matrixList) {
        RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

        int nx,ny,nz; //global_size_t
        nx = ny = nz = 5;
        nx = matrixList.get("nx",nx);
        ny = matrixList.get("ny",ny);
        nz = matrixList.get("nz",nz);

        std::string matrixType = matrixList.get("matrixType","Laplace1D");
        GO numGlobalElements; //global_size_t
        if (matrixType == "Laplace1D")
          numGlobalElements = nx;
        else if (matrixType == "Laplace2D")
          numGlobalElements = nx*ny;
        else if (matrixType == "Laplace3D")
          numGlobalElements = nx*ny*nz;
        else {
          std::string msg = matrixType + " is unsupported (in unit testing)";
          throw(MueLu::Exceptions::RuntimeError(msg));
        }

        RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, 0, comm);

        RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>(matrixType,map,matrixList);
        return Op;
      } // BuildMatrix()

      // Create a 1D Poisson matrix with the specified number of rows
      // nx: global number of rows
      static RCP<Operator> Build1DPoisson(int nx) { //global_size_t
        Teuchos::ParameterList matrixList;
        matrixList.set("nx", nx);
        matrixList.set("matrixType","Laplace1D");
        RCP<Operator> A = BuildMatrix(matrixList);
        return A;
      } // Build1DPoisson()

      // Create a 2D Poisson matrix with the specified number of rows
      // nx: global number of rows
      // ny: global number of rows
      static RCP<Operator> Build2DPoisson(int nx, int ny=-1) { //global_size_t
        Teuchos::ParameterList matrixList;
        if (ny==-1) ny=nx;
        matrixList.set("nx", nx);
        matrixList.set("ny", ny);
        matrixList.set("matrixType","Laplace2D");
        RCP<Operator> A = BuildMatrix(matrixList);
        return A;
      } // Build2DPoisson()
 
      // Needed to initialize correctly a level used for testing SingleLevel factory Build() methods.
      // This method initializes LevelID and DefaultFactoryHandlers
      static void createSingleLevelHierarchy(Level& currentLevel) {
        Hierarchy H;
        H.SetLevel(rcpFromRef(currentLevel)); // rcpFromRef safe here
      }
      
      // Needed to initialize correctly levels used for testing TwoLevel factory Build() methods.
      // This method initializes LevelIDs, DefaultFactoryHandlers and references to previous level
      static void createTwoLevelHierarchy(Level& coarseLevel, Level& fineLevel) {
        Hierarchy H;
        H.SetLevel(rcpFromRef(fineLevel));   // rcpFromRef safe here
        H.SetLevel(rcpFromRef(coarseLevel)); // rcpFromRef safe here
      }
      
    }; // class Factory

  } // namespace TestHelpers

} // namespace MueLu


// Macro to skip a test when UnderlyingLib==Epetra or Tpetra 
#define MUELU_TEST_ONLY_FOR(UnderlyingLib) \
  if (MueLu::TestHelpers::Parameters::getLib() == UnderlyingLib)

// Macro to skip a test when Epetra is used with Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal) \
  if (!(MueLu::TestHelpers::Parameters::getLib() == Xpetra::UseEpetra && (Teuchos::OrdinalTraits<LocalOrdinal>::name() != string("int") || Teuchos::OrdinalTraits<GlobalOrdinal>::name() != string("int"))))

// Macro to skip a test when Epetra is used with Scalar != double or Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) \
  if (!(MueLu::TestHelpers::Parameters::getLib() == Xpetra::UseEpetra && Teuchos::ScalarTraits<Scalar>::name() != string("double"))) \
    MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal)






//

//TODO: add directly to Teuchos ?
#include "../xpetra/test/Xpetra_UnitTestHelpers.hpp" // declaration of TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL


//


//! Namespace for MueLu test classes
namespace MueLuTests {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arcp;
  using Teuchos::arcpFromArrayView;
  using Teuchos::rcpFromRef;
  using Teuchos::null;
  using Teuchos::arcp_reinterpret_cast;
  using Teuchos::Array;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcpFromRef;

  using namespace MueLu::TestHelpers;
}

#endif // ifndef MUELU_TEST_HELPERS_H
