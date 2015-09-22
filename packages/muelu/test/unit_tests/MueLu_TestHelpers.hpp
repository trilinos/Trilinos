// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_TEST_HELPERS_H
#define MUELU_TEST_HELPERS_H

#include <string>
#ifndef _MSC_VER
#include <dirent.h>
#endif

// Teuchos
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

// Xpetra
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Parameters.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_CrsGraph.hpp"

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Level.hpp"

// Galeri
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"

#include "MueLu_NoFactory.hpp"

// Conditional Tpetra stuff
#ifdef HAVE_MUELU_TPETRA
#include "Xpetra_TpetraCrsGraph.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"
#include "Xpetra_TpetraBlockCrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Experimental_BlockCrsMatrix.hpp"
#endif

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
        return TestHelpers::Parameters::xpetraParameters.GetLib();
      }
    };

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    class TestFactory {
#include "MueLu_UseShortNames.hpp"

    private:
      TestFactory() {} // static class

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
      static RCP<Matrix> BuildMatrix(Teuchos::ParameterList &matrixList, Xpetra::UnderlyingLib lib) {
        RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

        if (lib == Xpetra::NotSpecified)
          lib = TestHelpers::Parameters::getLib();

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

        RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, 0, comm);
        RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
            Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixType, map, matrixList);
        RCP<Matrix> Op = Pr->BuildMatrix();

        return Op;
      } // BuildMatrix()

      // Create a tridiagonal matrix (stencil = [b,a,c]) with the specified number of rows
      // dofMap: row map of matrix
      static RCP<Matrix> BuildTridiag(RCP<const Map> dofMap, Scalar a, Scalar b, Scalar c, Xpetra::UnderlyingLib lib=Xpetra::NotSpecified) { //global_size_t

        if (lib == Xpetra::NotSpecified)
          lib = TestHelpers::Parameters::getLib();

        RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

        Teuchos::RCP<Matrix> mtx = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(dofMap, 3);

        LocalOrdinal NumMyElements = dofMap->getNodeNumElements();
        Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = dofMap->getNodeElementList();
        GlobalOrdinal indexBase = dofMap->getIndexBase();

        GlobalOrdinal NumEntries;
        LocalOrdinal  nnz = 3;
        std::vector<Scalar>        Values(nnz);
        std::vector<GlobalOrdinal> Indices(nnz);

        // access information from strided block map
        std::vector<size_t> strInfo = std::vector<size_t>();
        size_t blockSize = 1;
        LocalOrdinal blockId = -1;
        GlobalOrdinal offset = 0;

        Teuchos::RCP<const StridedMap> strdofMap = Teuchos::rcp_dynamic_cast<const StridedMap>(dofMap);
        if(strdofMap != Teuchos::null) {
          strInfo = strdofMap->getStridingData();
          blockSize = strdofMap->getFixedBlockSize();
          blockId = strdofMap->getStridedBlockId();
          offset = strdofMap->getOffset();
          TEUCHOS_TEST_FOR_EXCEPTION(blockId > -1 && strInfo[blockId]==1, MueLu::Exceptions::RuntimeError,
                                             "MueLu::TestHelpers::BuildTridiag: strInfo block size must be > 1.");
          // todo: write one more special case for a row being first and last bock row
        } else {
          // no striding information. emulate block matrix
          blockId = 0;
          strInfo.push_back(blockSize); // default block size = 1
        }

        GlobalOrdinal rrmax = (dofMap->getMaxAllGlobalIndex()-offset-indexBase) / blockSize;

        // loop over all rows
        for (LocalOrdinal i = 0; i < NumMyElements; i++) {

          GlobalOrdinal rr = (MyGlobalElements[i]-offset-indexBase) / blockSize + indexBase;  // node index

          // distinguish 5 different cases

          GlobalOrdinal blockOffset = 0;
          for (LocalOrdinal k=0; k<blockId; k++)
            blockOffset += Teuchos::as<GlobalOrdinal>(strInfo[k]);

          if (MyGlobalElements[i] == blockOffset + offset + indexBase) {
            // very first row
            Indices[0] = MyGlobalElements[i];
            Values [0] = a;
            Indices[1] = MyGlobalElements[i] + 1;
            Values [1] = c;
            NumEntries = 2;
          } else if (MyGlobalElements[i] == rrmax * Teuchos::as<GO>(blockSize) + blockOffset + Teuchos::as<GO>(strInfo[blockId]) - 1 + offset + indexBase) {
            // very last row
            Indices[0] = MyGlobalElements[i] - 1;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            NumEntries = 2;
          } else if (MyGlobalElements[i] == rr * Teuchos::as<GO>(blockSize) + blockOffset + Teuchos::as<GO>(strInfo[blockId]) - 1 + offset + indexBase) {
            // last row in current node block
            Indices[0] = MyGlobalElements[i] - 1;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            Indices[2] = (rr+1)*blockSize + blockOffset + offset + indexBase;
            Values [2] = c;
            NumEntries = 3;
          } else if (MyGlobalElements[i] == rr * Teuchos::as<GO>(blockSize) + blockOffset + offset + indexBase) {
            // first row in current node block
            Indices[0] = (rr-1)*blockSize + blockOffset + strInfo[blockId] - 1 + offset + indexBase;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            Indices[2] = MyGlobalElements[i] + 1;
            Values [2] = c;
            NumEntries = 3;
          } else {
            // usual row entries in block rows
            Indices[0] = MyGlobalElements[i] - 1;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            Indices[2] = MyGlobalElements[i] + 1;
            Values [2] = c;
            NumEntries = 3;
          }

          // debug output
          /*std::cout << "--------> Proc " << comm->getRank() << " lrow " << i << " grow " << MyGlobalElements[i] << " node " << rr << " values = " ;
          for (size_t k = 0; k < NumEntries; k++) {
            std::cout << " " << Indices[k] << "->" << Values[k] << "   ";
          }
          std::cout << std::endl;*/

          // put the off-diagonal entries
          // Xpetra wants ArrayViews (sigh)
          Teuchos::ArrayView<Scalar>        av(&Values [0], NumEntries);
          Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
          mtx->insertGlobalValues(MyGlobalElements[i], iv, av);
        }

        mtx->fillComplete();

        return mtx;
      } // BuildTridiag()

      // Create a 1D Poisson matrix with the specified number of rows
      // nx: global number of rows
      static RCP<Matrix> Build1DPoisson(int nx, Xpetra::UnderlyingLib lib=Xpetra::NotSpecified) { //global_size_t
        Teuchos::ParameterList matrixList;
        matrixList.set("nx", nx);
        matrixList.set("matrixType","Laplace1D");
        RCP<Matrix> A = BuildMatrix(matrixList,lib);
        return A;
      } // Build1DPoisson()

      // Create a 2D Poisson matrix with the specified number of rows
      // nx: global number of rows
      // ny: global number of rows
      static RCP<Matrix> Build2DPoisson(int nx, int ny=-1, Xpetra::UnderlyingLib lib=Xpetra::NotSpecified) { //global_size_t
        Teuchos::ParameterList matrixList;
        if (ny==-1) ny=nx;
        matrixList.set("nx", nx);
        matrixList.set("ny", ny);
        matrixList.set("matrixType","Laplace2D");
        RCP<Matrix> A = BuildMatrix(matrixList,lib);
        return A;
      } // Build2DPoisson()


     // Create a matrix as specified by parameter list options
     static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList &matrixList, Xpetra::UnderlyingLib lib=Xpetra::NotSpecified) {
       RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
       RCP<Matrix> Op;

        if (lib == Xpetra::NotSpecified)
          lib = TestHelpers::Parameters::getLib();

        // This only works for Tpetra
        if(lib!=Xpetra::UseTpetra) return Op;

        // Make the graph
        RCP<Matrix> FirstMatrix = BuildMatrix(matrixList,lib);
        RCP<const Xpetra::CrsGraph<LO,GO,NO> > Graph = FirstMatrix->getCrsGraph();

#if defined(HAVE_MUELU_TPETRA)
        // Thanks for the code, Travis!
        int blocksize = 3;
        RCP<const Xpetra::TpetraCrsGraph<LO,GO,NO> > TGraph = rcp_dynamic_cast<const Xpetra::TpetraCrsGraph<LO,GO,NO> >(Graph);
        RCP<const Tpetra::CrsGraph<LO,GO,NO> > TTGraph = TGraph->getTpetra_CrsGraph();

        RCP<Tpetra::Experimental::BlockCrsMatrix<SC,LO,GO,NO> > bcrsmatrix = rcp(new Tpetra::Experimental::BlockCrsMatrix<SC,LO,GO,NO> (*TTGraph, blocksize));

        const Tpetra::Map<LO,GO,NO>& meshRowMap = *bcrsmatrix->getRowMap();
        const Scalar zero   = Teuchos::ScalarTraits<SC>::zero();
        const Scalar one   = Teuchos::ScalarTraits<SC>::one();
        const Scalar two   = one+one;
        const Scalar three = two+one;

        Teuchos::Array<SC> basematrix(blocksize*blocksize, zero);
        basematrix[0] = two;
        basematrix[2] = three;
        basematrix[3] = three;
        basematrix[4] = two;
        basematrix[7] = three;
        basematrix[8] = two;
        Teuchos::Array<LO> lclColInds(1);
        for (LocalOrdinal lclRowInd = meshRowMap.getMinLocalIndex (); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
          lclColInds[0] = lclRowInd;
          bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &basematrix[0], 1);
        }
        bcrsmatrix->computeDiagonalGraph(); // Needs to get done to smooth for some reason

        RCP<Xpetra::CrsMatrix<SC,LO,GO,NO> > temp = rcp(new Xpetra::TpetraBlockCrsMatrix<SC,LO,GO,NO>(bcrsmatrix));
        Op = rcp(new Xpetra::CrsMatrixWrap<SC,LO,GO,NO>(temp));
#endif
        return Op;
     } // BuildMatrix()


      // Needed to initialize correctly a level used for testing SingleLevel factory Build() methods.
      // This method initializes LevelID and linked list of level
      static void createSingleLevelHierarchy(Level& currentLevel) {
        RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
        currentLevel.SetFactoryManager(factoryHandler);

        currentLevel.SetLevelID(0);
#ifdef HAVE_MUELU_TIMER_SYNCHRONIZATION
        currentLevel.SetComm(TestHelpers::Parameters::getDefaultComm());
#endif
      }

      // Needed to initialize correctly levels used for testing TwoLevel factory Build() methods.
      // This method initializes LevelID and linked list of level
      static void createTwoLevelHierarchy(Level& fineLevel, Level& coarseLevel) {
        RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
        fineLevel.SetFactoryManager(factoryHandler);
        coarseLevel.SetFactoryManager(factoryHandler);

        coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));

        fineLevel.SetLevelID(0);
        coarseLevel.SetLevelID(1);
#ifdef HAVE_MUELU_TIMER_SYNCHRONIZATION
        fineLevel.SetComm(TestHelpers::Parameters::getDefaultComm());
        coarseLevel.SetComm(TestHelpers::Parameters::getDefaultComm());
#endif
      }

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
      static RCP<SmootherPrototype> createSmootherPrototype(const std::string& type="Gauss-Seidel", LO sweeps=1) {
        Teuchos::ParameterList  ifpackList;
        ifpackList.set("relaxation: type", type);
        ifpackList.set("relaxation: sweeps", (LO) sweeps);
        ifpackList.set("relaxation: damping factor", (SC) 1.0);
        return rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
      }
#endif

    }; // class Factory



    //! Return the list of files in the directory. Only files that are matching '*filter*' are returned.
    ArrayRCP<std::string> GetFileList(const std::string & dirPath, const std::string & filter);

  } // namespace TestHelpers

} // namespace MueLu


// Macro to skip a test when UnderlyingLib==Epetra or Tpetra
#define MUELU_TEST_ONLY_FOR(UnderlyingLib) \
  if (TestHelpers::Parameters::getLib() == UnderlyingLib)

// Macro to skip a test when Epetra is used with Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal) \
  if (!(TestHelpers::Parameters::getLib() == Xpetra::UseEpetra && (Teuchos::OrdinalTraits<LocalOrdinal>::name() != string("int") || Teuchos::OrdinalTraits<GlobalOrdinal>::name() != string("int"))))

// Macro to skip a test when Epetra is used with Scalar != double or Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) \
  if (!(TestHelpers::Parameters::getLib() == Xpetra::UseEpetra && Teuchos::ScalarTraits<Scalar>::name() != string("double"))) \
    MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal)






//

//TODO: add directly to Teuchos ?
//#include "../xpetra/test/Xpetra_UnitTestHelpers.hpp" // declaration of TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL


//


//! Namespace for MueLu test classes
namespace MueLuTests {

  using namespace TestHelpers;
}

#endif // ifndef MUELU_TEST_HELPERS_H
