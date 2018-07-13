// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_HPP_

#include "Xpetra_ConfigDefs.hpp"

// #include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_StridedMap.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include <TpetraExt_TripleMatrixMultiply.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
// #include <Xpetra_TpetraMultiVector.hpp>
// #include <Xpetra_TpetraVector.hpp>
#endif // HAVE_XPETRA_TPETRA

namespace Xpetra {

  template <class Scalar,
            class LocalOrdinal  /*= int*/,
            class GlobalOrdinal /*= LocalOrdinal*/,
            class Node          /*= KokkosClassic::DefaultNode::DefaultNodeType*/>
  class TripleMatrixMultiply {
#undef XPETRA_TRIPLEMATRIXMULTIPLY_SHORT
#include "Xpetra_UseShortNames.hpp"

  public:

    /** Given CrsMatrix objects A, B and C, form the product C = A*B.
      In a parallel setting, A and B need not have matching distributions,
      but C needs to have the same row-map as A (if transposeA is false).
      At this time C=AT*B and C=A*BT are known to not work. However,
      C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

      @param A Input, must already have had 'FillComplete()' called.
      @param transposeA Input, whether to use transpose of matrix A.
      @param B Input, must already have had 'FillComplete()' called.
      @param transposeB Input, whether to use transpose of matrix B.
      @param C Result. On entry to this method, it doesn't matter whether
      FillComplete() has already been called on C or not. If it has,
      then C's graph must already contain all nonzero locations that
      will be produced when forming the product A*B. On exit,
      C.FillComplete() will have been called, unless the last argument
      to this function is specified to be false.
      @param call_FillComplete_on_result Optional argument, defaults to true.
      Power users may specify this argument to be false if they *DON'T*
      want this function to call C.FillComplete. (It is often useful
      to allow this function to call C.FillComplete, in cases where
      one or both of the input matrices are rectangular and it is not
      trivial to know which maps to use for the domain- and range-maps.)

*/
    static void MultiplyRAP(const Matrix& R, bool transposeR,
                            const Matrix& A, bool transposeA,
                            const Matrix& P, bool transposeP,
                            Matrix& Ac,
                            bool call_FillComplete_on_result = true,
                            bool doOptimizeStorage           = true,
                            const std::string & label        = std::string(),
                            const RCP<ParameterList>& params = null) {

      TEUCHOS_TEST_FOR_EXCEPTION(transposeR == false && Ac.getRowMap()->isSameAs(*R.getRowMap()) == false,
        Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as row map of R");
      TEUCHOS_TEST_FOR_EXCEPTION(transposeR == true && Ac.getRowMap()->isSameAs(*R.getDomainMap()) == false,
        Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as domain map of R");

      TEUCHOS_TEST_FOR_EXCEPTION(!R.isFillComplete(), Exceptions::RuntimeError, "R is not fill-completed");
      TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
      TEUCHOS_TEST_FOR_EXCEPTION(!P.isFillComplete(), Exceptions::RuntimeError, "P is not fill-completed");

      bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

      if (Ac.getRowMap()->lib() == Xpetra::UseEpetra) {
        throw(Xpetra::Exceptions::RuntimeError("Xpetra::TripleMatrixMultiply::MultiplyRAP is only implemented for Tpetra"));
      } else if (Ac.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpR = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(R);
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpP = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(P);
        Tpetra::CrsMatrix<SC,LO,GO,NO> &       tpAc = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Ac);

        // 18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
        // Previously, Tpetra's matrix matrix multiply did not support fillComplete.
        Tpetra::TripleMatrixMultiply::MultiplyRAP(tpR, transposeR, tpA, transposeA, tpP, transposeP, tpAc, haveMultiplyDoFillComplete, label, params);
#else
        throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
      }

      if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
        RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
        fillParams->set("Optimize Storage", doOptimizeStorage);
        Ac.fillComplete((transposeP) ? P.getRangeMap() : P.getDomainMap(),
                        (transposeR) ? R.getDomainMap() : R.getRangeMap(),
                        fillParams);
      }

      // transfer striding information
      RCP<const Map> domainMap = Teuchos::null;
      RCP<const Map> rangeMap  = Teuchos::null;

      const std::string stridedViewLabel("stridedMaps");
      const size_t        blkSize = 1;
      std::vector<size_t> stridingInfo(1, blkSize);
      LocalOrdinal        stridedBlockId = -1;
      
      if (R.IsView(stridedViewLabel)) {
        rangeMap  = transposeR ? R.getColMap(stridedViewLabel) : R.getRowMap(stridedViewLabel);
      } else {
        rangeMap  = transposeR ? R.getDomainMap()       : R.getRangeMap();
        rangeMap  = StridedMapFactory::Build(rangeMap,  stridingInfo, stridedBlockId);
      }
      
      if (P.IsView(stridedViewLabel)) {
          domainMap = transposeP ? P.getRowMap(stridedViewLabel) : P.getColMap(stridedViewLabel);
      } else {
        domainMap = transposeP ? P.getRangeMap()        : P.getDomainMap();
        domainMap = StridedMapFactory::Build(domainMap, stridingInfo, stridedBlockId);
      }
      Ac.CreateView(stridedViewLabel, rangeMap, domainMap);

    } // end Multiply

  }; // class TripleMatrixMultiply

#ifdef HAVE_XPETRA_EPETRA
  // specialization TripleMatrixMultiply for SC=double, LO=GO=int
  template <>
  class TripleMatrixMultiply<double,int,int,EpetraNode> {
    typedef double          Scalar;
    typedef int             LocalOrdinal;
    typedef int             GlobalOrdinal;
    typedef EpetraNode      Node;
#include "Xpetra_UseShortNames.hpp"

  public:

    static void MultiplyRAP(const Matrix& R, bool transposeR,
                            const Matrix& A, bool transposeA,
                            const Matrix& P, bool transposeP,
                            Matrix& Ac,
                            bool call_FillComplete_on_result = true,
                            bool doOptimizeStorage           = true,
                            const std::string & label        = std::string(),
                            const RCP<ParameterList>& params = null) {

      TEUCHOS_TEST_FOR_EXCEPTION(transposeR == false && Ac.getRowMap()->isSameAs(*R.getRowMap()) == false,
        Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as row map of R");
      TEUCHOS_TEST_FOR_EXCEPTION(transposeR == true && Ac.getRowMap()->isSameAs(*R.getDomainMap()) == false,
        Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as domain map of R");

      TEUCHOS_TEST_FOR_EXCEPTION(!R.isFillComplete(), Exceptions::RuntimeError, "R is not fill-completed");
      TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
      TEUCHOS_TEST_FOR_EXCEPTION(!P.isFillComplete(), Exceptions::RuntimeError, "P is not fill-completed");

      bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

      if (Ac.getRowMap()->lib() == Xpetra::UseEpetra) {
        throw(Xpetra::Exceptions::RuntimeError("Xpetra::TripleMatrixMultiply::MultiplyRAP is only implemented for Tpetra"));
      } else if (Ac.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
# if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
      (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
        throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra <double,int,int> ETI enabled."));
# else
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpR = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(R);
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpP = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(P);
        Tpetra::CrsMatrix<SC,LO,GO,NO> &       tpAc = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Ac);

        // 18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
        // Previously, Tpetra's matrix matrix multiply did not support fillComplete.
        Tpetra::TripleMatrixMultiply::MultiplyRAP(tpR, transposeR, tpA, transposeA, tpP, transposeP, tpAc, haveMultiplyDoFillComplete, label, params);
# endif
#else
        throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
        if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
          RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
          fillParams->set("Optimize Storage", doOptimizeStorage);
          Ac.fillComplete((transposeP) ? P.getRangeMap() : P.getDomainMap(),
                          (transposeR) ? R.getDomainMap() : R.getRangeMap(),
                          fillParams);
        }

        // transfer striding information
        RCP<const Map> domainMap = Teuchos::null;
        RCP<const Map> rangeMap  = Teuchos::null;

        const std::string stridedViewLabel("stridedMaps");
        const size_t        blkSize = 1;
        std::vector<size_t> stridingInfo(1, blkSize);
        LocalOrdinal        stridedBlockId = -1;
      
        if (R.IsView(stridedViewLabel)) {
          rangeMap  = transposeR ? R.getColMap(stridedViewLabel) : R.getRowMap(stridedViewLabel);
        } else {
          rangeMap  = transposeR ? R.getDomainMap()       : R.getRangeMap();
          rangeMap  = StridedMapFactory::Build(rangeMap,  stridingInfo, stridedBlockId);
        }
      
        if (P.IsView(stridedViewLabel)) {
          domainMap = transposeP ? P.getRowMap(stridedViewLabel) : P.getColMap(stridedViewLabel);
        } else {
          domainMap = transposeP ? P.getRangeMap()        : P.getDomainMap();
          domainMap = StridedMapFactory::Build(domainMap, stridingInfo, stridedBlockId);
        }
        Ac.CreateView(stridedViewLabel, rangeMap, domainMap);
      }

    } // end Multiply

  }; // end specialization on SC=double, GO=int and NO=EpetraNode

  // specialization TripleMatrixMultiply for SC=double, GO=long long and NO=EpetraNode
  template <>
  class TripleMatrixMultiply<double,int,long long,EpetraNode> {
    typedef double          Scalar;
    typedef int             LocalOrdinal;
    typedef long long       GlobalOrdinal;
    typedef EpetraNode      Node;
#include "Xpetra_UseShortNames.hpp"

  public:

    static void MultiplyRAP(const Matrix& R, bool transposeR,
                            const Matrix& A, bool transposeA,
                            const Matrix& P, bool transposeP,
                            Matrix& Ac,
                            bool call_FillComplete_on_result = true,
                            bool doOptimizeStorage           = true,
                            const std::string & label        = std::string(),
                            const RCP<ParameterList>& params = null) {

      TEUCHOS_TEST_FOR_EXCEPTION(transposeR == false && Ac.getRowMap()->isSameAs(*R.getRowMap()) == false,
        Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as row map of R");
      TEUCHOS_TEST_FOR_EXCEPTION(transposeR == true && Ac.getRowMap()->isSameAs(*R.getDomainMap()) == false,
        Exceptions::RuntimeError, "XpetraExt::TripleMatrixMultiply::MultiplyRAP: row map of Ac is not same as domain map of R");

      TEUCHOS_TEST_FOR_EXCEPTION(!R.isFillComplete(), Exceptions::RuntimeError, "R is not fill-completed");
      TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
      TEUCHOS_TEST_FOR_EXCEPTION(!P.isFillComplete(), Exceptions::RuntimeError, "P is not fill-completed");

      bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

      if (Ac.getRowMap()->lib() == Xpetra::UseEpetra) {
        throw(Xpetra::Exceptions::RuntimeError("Xpetra::TripleMatrixMultiply::MultiplyRAP is only implemented for Tpetra"));
      } else if (Ac.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
# if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
      (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra <double,int,long long,EpetraNode> ETI enabled."));
# else
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpR = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(R);
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
        const Tpetra::CrsMatrix<SC,LO,GO,NO> & tpP = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(P);
        Tpetra::CrsMatrix<SC,LO,GO,NO> &       tpAc = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(Ac);

        // 18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
        // Previously, Tpetra's matrix matrix multiply did not support fillComplete.
        Tpetra::TripleMatrixMultiply::MultiplyRAP(tpR, transposeR, tpA, transposeA, tpP, transposeP, tpAc, haveMultiplyDoFillComplete, label, params);
# endif
#else
        throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
        if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
          RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
          fillParams->set("Optimize Storage", doOptimizeStorage);
          Ac.fillComplete((transposeP) ? P.getRangeMap() : P.getDomainMap(),
                          (transposeR) ? R.getDomainMap() : R.getRangeMap(),
                          fillParams);
        }

        // transfer striding information
        RCP<const Map> domainMap = Teuchos::null;
        RCP<const Map> rangeMap  = Teuchos::null;

        const std::string stridedViewLabel("stridedMaps");
        const size_t        blkSize = 1;
        std::vector<size_t> stridingInfo(1, blkSize);
        LocalOrdinal        stridedBlockId = -1;
      
        if (R.IsView(stridedViewLabel)) {
          rangeMap  = transposeR ? R.getColMap(stridedViewLabel) : R.getRowMap(stridedViewLabel);
        } else {
          rangeMap  = transposeR ? R.getDomainMap()       : R.getRangeMap();
          rangeMap  = StridedMapFactory::Build(rangeMap,  stridingInfo, stridedBlockId);
        }
      
        if (P.IsView(stridedViewLabel)) {
          domainMap = transposeP ? P.getRowMap(stridedViewLabel) : P.getColMap(stridedViewLabel);
        } else {
          domainMap = transposeP ? P.getRangeMap()        : P.getDomainMap();
          domainMap = StridedMapFactory::Build(domainMap, stridingInfo, stridedBlockId);
        }
        Ac.CreateView(stridedViewLabel, rangeMap, domainMap);
      }

    } // end Multiply

  }; // end specialization on GO=long long and NO=EpetraNode
#endif

} // end namespace Xpetra

#define XPETRA_TRIPLEMATRIXMULTIPLY_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_HPP_ */
