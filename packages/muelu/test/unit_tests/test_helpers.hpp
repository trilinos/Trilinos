#ifndef MUELU_UNITTEST_HELPERS_H
#define MUELU_UNITTEST_HELPERS_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Cthulhu_ConfigDefs.hpp" //TODO: use Cthulhu
#include "Cthulhu_DefaultPlatform.hpp"

#include "Cthulhu_Map.hpp"
#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#endif
#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_CrsOperator.hpp"

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"
#include <iostream>

#include "MueLu_Exceptions.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLu {
namespace UnitTest {

  using Cthulhu::global_size_t;
  using Teuchos::RCP;

  typedef Cthulhu::DefaultPlatform::DefaultPlatformType::NodeType Node;

  inline
  RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    return Cthulhu::DefaultPlatform::getDefaultPlatform().getComm(); //TODO: use Cthulhu here
  }

#ifdef HAVE_CTHULHU_TPETRA
  //
  // Function that creates a map containing a specified number of local elements per process.
  //
  template<class LocalOrdinal,class GlobalOrdinal,class Node>
  const RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> >
  create_map(LocalOrdinal num_elements_per_proc)
  { 
    RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    const LocalOrdinal indexBase = 0;
  
    //TODO: use CthulhuMapFactory here
    return Teuchos::rcp(new Cthulhu::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_elements_per_proc, indexBase, comm));
  
  } // create_map()
#endif

  //create a matrix as specified by parameter list options
  template<class SC, class LO,class GO>
  RCP<CrsOperator> create_test_matrix(Teuchos::ParameterList &matrixList)
  {
    RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

    GO nx,ny,nz;
    nx = ny = nz = 5;
    nx = matrixList.get("nx",nx);
    ny = matrixList.get("ny",ny);
    nz = matrixList.get("nz",nz);

    std::string matrixType = matrixList.get("matrixType","Laplace1D");
    GO numGlobalElements;
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
    LO indexBase = 0;
    RCP<const Map > map;
    map = rcp( new Cthulhu::EpetraMap(numGlobalElements, indexBase, comm) );

    RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map,
                                                    CrsOperator>(matrixType,map,matrixList);
    return Op;
  } //create_test_matrix

  //create a 1D Poisson matrix with the specified number of rows
  template<class SC, class LO, class GO>
  RCP<CrsOperator> create_1d_poisson_matrix(GO numRows)
  {
    Teuchos::ParameterList matrixList;
    matrixList.set("nx",numRows);
    matrixList.set("matrixType","Laplace1D");
    RCP<CrsOperator> A = create_test_matrix<SC,LO,GO>(matrixList);
    return A;
  }

} // namespace UnitTest
} // namespace MueLu

#endif // ifndef MUELU_UNITTEST_HELPERS_H
