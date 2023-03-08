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
#ifndef MUELU_UTILITIES_KOKKOS_DECL_HPP
#define MUELU_UTILITIES_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <string>

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Xpetra_BlockedCrsMatrix_fwd.hpp"
#include "Xpetra_CrsMatrix_fwd.hpp"
#include "Xpetra_CrsMatrixWrap_fwd.hpp"
#include "Xpetra_ExportFactory.hpp"
#include "Xpetra_ImportFactory_fwd.hpp"
#include "Xpetra_MapFactory_fwd.hpp"
#include "Xpetra_Map_fwd.hpp"
#include "Xpetra_MatrixFactory_fwd.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_MultiVector_fwd.hpp"
#include "Xpetra_Operator_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_Vector_fwd.hpp"

#include "Xpetra_IO.hpp"

#include "Kokkos_ArithTraits.hpp"

#ifdef HAVE_MUELU_EPETRA
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Xpetra_EpetraCrsMatrix_fwd.hpp"
#include "Xpetra_EpetraMultiVector_fwd.hpp"
#endif

#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_UtilitiesBase.hpp"

#ifdef HAVE_MUELU_TPETRA
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Xpetra_TpetraCrsMatrix_fwd.hpp"
#include "Xpetra_TpetraMultiVector_fwd.hpp"
#endif


namespace MueLu {

  /*!
    @class Utilities
    @brief MueLu utility class.

    This class provides a number of static helper methods. Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
  */
  template <class Scalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class Utilities_kokkos : public MueLu::UtilitiesBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
#undef MUELU_UTILITIES_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    using TST                   = Teuchos::ScalarTraits<SC>;
    using Magnitude             = typename TST::magnitudeType;
    using CoordinateType        = typename TST::coordinateType;
    using RealValuedMultiVector = Xpetra::MultiVector<CoordinateType,LO,GO,NO>;

     /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in RCP<Vector>.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixDiagonal(const Matrix& A); // FIXME

    /*! @brief Extract Matrix Diagonal

    Returns inverse of the Matrix diagonal in RCP<Vector>.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = TST::eps()*100, const bool doLumped = false); // FIXME



    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix& A); // FIXME

  }; // class Utils


  /*!
        @class Utilities
        @brief MueLu utility class (specialization SC=double and LO=GO=int).

        This class provides a number of static helper methods. Some are temporary and will eventually
        go away, while others should be moved to Xpetra.

        Note: this is the implementation for Epetra. Tpetra throws if TPETRA_INST_INT_INT is disabled!
   */
  template <class Node>
  class Utilities_kokkos<double,int,int,Node> : public UtilitiesBase<double,int,int,Node> {
  public:
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    using CoordinateType        = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
    using RealValuedMultiVector = Xpetra::MultiVector<CoordinateType,LocalOrdinal,GlobalOrdinal,Node>;

  private:
#undef MUELU_UTILITIES_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    static RCP<Vector> GetMatrixDiagonal(const Matrix& A) {
      const auto rowMap = A.getRowMap();
      auto diag = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowMap,true);

      A.getLocalDiagCopy(*diag);

      return diag;
    }

    static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<SC>::eps()*100, const bool doLumped=false);


  }; // class Utilities (specialization SC=double LO=GO=int)




} //namespace MueLu

#define MUELU_UTILITIES_KOKKOS_SHORT

#endif // MUELU_UTILITIES_KOKKOS_DECL_HPP
