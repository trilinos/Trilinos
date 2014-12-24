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
#include "MueLu_Utilities_def.hpp"

#include <string>

#ifdef HAVE_MUELU_EPETRAEXT
#include "EpetraExt_Transpose_RowMatrix.h"
#endif

namespace MueLu {

  RCP<Xpetra::Matrix<double, int, int> > Utils2<double, int, int>::Transpose(Matrix& Op, bool optimizeTranspose) {
   typedef double                                           Scalar;
   typedef int                                              LocalOrdinal;
   typedef int                                              GlobalOrdinal;
   typedef KokkosClassic::DefaultNode::DefaultNodeType      Node;

   Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ZZ Entire Transpose"));

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    try {
      Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(Op);
    }
    catch (...) {
      TorE = "tpetra";
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    if (TorE == "tpetra") {
      try {
        const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(Op);

        // Compute the transpose A of the Tpetra matrix tpetraOp.
        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
        {
          Teuchos::TimeMonitor tmm (*Teuchos::TimeMonitor::getNewCounter ("ZZ Tpetra Transpose Only"));
          Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp));
          A = transposer.createTranspose();
        }
        RCP<Xpetra::TpetraCrsMatrix<SC> > AA   = rcp(new Xpetra::TpetraCrsMatrix<SC>(A));
        RCP<Xpetra::CrsMatrix<SC> >       AAA  = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
        RCP<Xpetra::CrsMatrixWrap<SC> >   AAAA = rcp( new Xpetra::CrsMatrixWrap<SC> (AAA));

        return AAAA;
      }
      catch (std::exception& e) {
        std::cout << "threw exception '" << e.what() << "'" << std::endl;
        throw Exceptions::RuntimeError("Utils::Transpose failed, perhaps because matrix is not a Crs matrix");
      }
    } //if
#endif

    if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
#else
      throw Exceptions::RuntimeError("Tpetra");
#endif // HAVE_MUELU_TPETRA

    } else {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      // Epetra case
      Epetra_CrsMatrix& epetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Op);
      EpetraExt::RowMatrix_Transpose transposer;
      Epetra_CrsMatrix * A = dynamic_cast<Epetra_CrsMatrix*>(&transposer(epetraOp));
      transposer.ReleaseTranspose(); // So we can keep A in Muelu...

      RCP<Epetra_CrsMatrix> rcpA(A);
      RCP<EpetraCrsMatrix>            AA   = rcp(new EpetraCrsMatrix(rcpA));
      RCP<Xpetra::CrsMatrix<SC> >     AAA  = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
      RCP<Xpetra::CrsMatrixWrap<SC> > AAAA = rcp( new Xpetra::CrsMatrixWrap<SC>(AAA));
      AAAA->fillComplete(Op.getRangeMap(), Op.getDomainMap());

      return AAAA;
#else
      throw Exceptions::RuntimeError("Epetra (Err. 2)");
#endif
    }
    return Teuchos::null;
  } //Transpose

  // -- ------------------------------------------------------- --

  void Utils2<double,int,int>::MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage) {
#ifdef HAVE_MUELU_EPETRA
    try {
      const Epetra_CrsMatrix& epOp = Utils<double,int,int>::Op2NonConstEpetraCrs(Op);

      Epetra_Map const &rowMap = epOp.RowMap();
      int nnz;
      double *vals;
      int *cols;

      for (int i = 0; i < rowMap.NumMyElements(); ++i) {
        epOp.ExtractMyRowView(i, nnz, vals, cols);
        for (int j = 0; j < nnz; ++j)
          vals[j] *= scalingVector[i];
      }

    } catch (...){
      throw Exceptions::RuntimeError("Only Epetra_CrsMatrix types can be scaled");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Epetra has not been enabled.");
#endif // HAVE_MUELU_EPETRA
  } //Utils2::MyOldScaleMatrix_Epetra()

  // -- ------------------------------------------------------- --

  void Utils2<double, int, int>::TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta) {
   typedef double                                           Scalar;
   typedef int                                              LocalOrdinal;
   typedef int                                              GlobalOrdinal;
   typedef KokkosClassic::DefaultNode::DefaultNodeType      Node;

    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      const Epetra_CrsMatrix& epA = Utils<double,int,int>::Op2EpetraCrs(A);
      Epetra_CrsMatrix&       epB = Utils<double,int,int>::Op2NonConstEpetraCrs(B);

      //FIXME is there a bug if beta=0?
      int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, beta);

      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value of " + toString(rv));
        std::ostringstream buf;
#else
      throw Exceptions::RuntimeError("MueLu must be compiled with EpetraExt.");
#endif

    } else if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = Utils<Scalar,LocalOrdinal,GlobalOrdinal>::Op2TpetraCrs(A);
      Tpetra::CrsMatrix<SC, LO, GO, NO>&       tpB = Utils<Scalar,LocalOrdinal,GlobalOrdinal>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, beta);
#else
      throw Exceptions::RuntimeError("MueLu must be compiled with Tpetra.");
#endif
    }
  } //Utils2::TwoMatrixAdd() (specialization)

  // -- ------------------------------------------------------- --

  void Utils2<double,int,int>::TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha,
                                            const Matrix& B, bool transposeB, SC beta,
                                            RCP<Matrix>& C,  Teuchos::FancyOStream &fos, bool AHasFixedNnzPerRow) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (C == Teuchos::null) {
      if (!A.isFillComplete() || !B.isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Global statistics are not available for estimates.");

      size_t maxNzInA     = A.getGlobalMaxNumRowEntries();
      size_t maxNzInB     = B.getGlobalMaxNumRowEntries();
      size_t numLocalRows = A.getNodeNumRows();

      if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
        // first check if either A or B has at most 1 nonzero per row
        // the case of both having at most 1 nz per row is handled by the ``else''
        Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);

        if ((maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = B.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;

        } else {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = A.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
        }

        fos << "Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
             << ", using static profiling" << std::endl;
        C = rcp(new Xpetra::CrsMatrixWrap<double,int,int,NO>(A.getRowMap(), exactNnzPerRow, Xpetra::StaticProfile));

      } else {
        // general case
        double nnzPerRowInA = Teuchos::as<double>(A.getGlobalNumEntries()) / A.getGlobalNumRows();
        double nnzPerRowInB = Teuchos::as<double>(B.getGlobalNumEntries()) / B.getGlobalNumRows();
        LO    nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

        LO maxPossible = A.getGlobalMaxNumRowEntries() + B.getGlobalMaxNumRowEntries();
        //Use static profiling (more efficient) if the estimate is at least as big as the max
        //possible nnz's in any single row of the result.
        Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

        fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
        fos << "Utils::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
             << ", max possible nnz per row in sum = " << maxPossible
             << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
             << std::endl;

        C = rcp(new Xpetra::CrsMatrixWrap<double,int,int,NO>(A.getRowMap(), nnzToAllocate, pft));
      }
      if (transposeB)
        fos << "Utils::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }

    if (C == Teuchos::null) {
      if (!A.isFillComplete() || !B.isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"Global statistics are not available for estimates.");

      double nnzPerRowInA = Teuchos::as<double>(A.getGlobalNumEntries()) / A.getGlobalNumRows();
      double nnzPerRowInB = Teuchos::as<double>(B.getGlobalNumEntries()) / B.getGlobalNumRows();
      LO    nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

      LO maxPossible = A.getGlobalMaxNumRowEntries() + B.getGlobalMaxNumRowEntries();
      //Use static profiling (more efficient) if the estimate is at least as big as the max possible nnz's in any single row of the result.
      Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

      fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
      fos << "Utils::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
           << ", max possible nnz per row in sum = " << maxPossible
           << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
           << std::endl;

      C = rcp(new Xpetra::CrsMatrixWrap<double,int,int,NO>(A.getRowMap(), nnzToAllocate, pft));

      if (transposeB)
        fos << "Utils::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      const Epetra_CrsMatrix& epA = Utils<double,int,int>::Op2EpetraCrs(A);
      const Epetra_CrsMatrix& epB = Utils<double,int,int>::Op2EpetraCrs(B);
      RCP<Epetra_CrsMatrix>   epC = Utils<double,int,int>::Op2NonConstEpetraCrs(C);
      Epetra_CrsMatrix* ref2epC = &*epC; //to avoid a compiler error...

      //FIXME is there a bug if beta=0?
      int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, transposeB, beta, ref2epC);

      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value of " + toString(rv));
#else
      throw Exceptions::RuntimeError("MueLu must be compile with EpetraExt.");
#endif

    } else if (C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      const Tpetra::CrsMatrix<SC, LO, GO, NO>&  tpA = Utils<double,int,int>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<SC, LO, GO, NO>&  tpB = Utils<double,int,int>::Op2TpetraCrs(B);
      RCP<  Tpetra::CrsMatrix<SC, LO, GO, NO> > tpC = Utils<double,int,int>::Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, transposeB, beta, tpC);
#else
      throw Exceptions::RuntimeError("MueLu must be compile with Tpetra.");
#endif
    }

    ///////////////////////// EXPERIMENTAL
    if (A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(A));
    if (B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(B));
    ///////////////////////// EXPERIMENTAL

  } //TwoMatrixAdd()

  // -- ------------------------------------------------------- --

  RCP<Xpetra::MultiVector<double,int,int> > Utils2<double,int,int>::ReadMultiVector(const std::string& fileName, const RCP<const Map>& map) {
    Xpetra::UnderlyingLib lib = map->lib();

    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      Epetra_MultiVector * MV;
      EpetraExt::MatrixMarketFileToMultiVector(fileName.c_str(), toEpetra(map), MV);
      return Xpetra::toXpetra<int>(rcp(MV));
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<SC,LO,GO,NO>                    sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>  reader_type;
      typedef Tpetra::Map<LO,GO,NO>                             map_type;
      typedef Tpetra::MultiVector<SC,LO,GO,NO>                  multivector_type;

      RCP<const map_type>   temp = toTpetra(map);
      RCP<multivector_type> TMV  = reader_type::readDenseFile(fileName,map->getComm(),map->getNode(),temp);
      RCP<MultiVector>      rmv  = Xpetra::toXpetra(TMV);
      return rmv;
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    return Teuchos::null;
  }

  RCP<const Xpetra::Map<int,int> > Utils2<double,int,int>::ReadMap(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm) {
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
        Epetra_Map *eMap;
        int rv = EpetraExt::MatrixMarketFileToMap(fileName.c_str(), *(Xpetra::toEpetra(comm)), eMap);
        if (rv != 0)
          throw Exceptions::RuntimeError("Error reading matrix with EpetraExt::MatrixMarketToMap (returned " + toString(rv) + ")");

        RCP<Epetra_Map> eMap1 = rcp(new Epetra_Map(*eMap));
        return Xpetra::toXpetra<int>(*eMap1);
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<double,int,int,NO> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>                          reader_type;

      RCP<NO> node = rcp(new NO());

      RCP<const Tpetra::Map<int,int,NO> > tMap = reader_type::readMapFile(fileName, comm, node);
      if (tMap.is_null())
        throw Exceptions::RuntimeError("The Tpetra::Map returned from readSparseFile() is null.");

      return Xpetra::toXpetra(tMap);
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }
  }


  /* Removes the following non-serializable data (A,P,R,Nullspace,Coordinates)
     from level-specific sublists from inList
     and moves it to nonSerialList.  Everything else is copied to serialList.
     This function returns the level number of the highest level for which
     non-serializable data was provided.
  */
  long ExtractNonSerializableData(const Teuchos::ParameterList& inList, Teuchos::ParameterList& serialList, Teuchos::ParameterList& nonSerialList) {
    using Teuchos::ParameterList;

    ParameterList dummy;
    long maxLevel = 0;

    for (ParameterList::ConstIterator it = inList.begin(); it != inList.end(); it++) {
      const std::string& levelName = it->first;

      // Check for mach of the form "level X" where X is a positive integer
      if (inList.isSublist(levelName) && levelName.find("level ") == 0 && levelName.size() > 6) {
        int levelID = strtol(levelName.substr(6).c_str(), 0, 0);
        if (maxLevel < levelID)
          maxLevel = levelID;

        // Split the sublist
        const ParameterList& levelList = inList.sublist(levelName);
        for (ParameterList::ConstIterator it2 = levelList.begin(); it2 != levelList.end(); it2++) {
          const std::string& name = it2->first;
          if (name == "A" || name == "P" || name == "R" || name == "Nullspace" || name == "Coordinates")
            nonSerialList.sublist(levelName).setEntry(name, it2->second);
          else
            serialList   .sublist(levelName).setEntry(name, it2->second);
        }

      } else {
        std::cout << levelName << std::endl;
        serialList.setEntry(it->first, it->second);
      }
    }

    return maxLevel;
  }


} // namespace MueLu
