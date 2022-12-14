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
//                    Jonathan Hu        (jhu@sandia.gov)
//                    Ray Tuminaro       (rstumin@sandia.gov)
//                    Chris Siefert      (csiefer@sandia.gov)
//                    Luc Berger-Vergoat (lberge@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_XPETRA_CRSMATRIX_UTILS_HPP_
#define PACKAGES_XPETRA_CRSMATRIX_UTILS_HPP_

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"
#include "Xpetra_Map.hpp"          // definition of UnderlyingLib

#ifdef HAVE_XPETRA_EPETRA
#include "Epetra_Util.h"
#endif

#include "Tpetra_Import_Util2.hpp"

namespace Xpetra {

  /*!
    @class CrsMatrixUtils
    @brief Xpetra utility class for CrsMatrix-related routines

    The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.

  */
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class CrsMatrixUtils {
#undef XPETRA_CRSMATRIXUTILS_SHORT

  public:

/// \brief Sort the entries of the (raw CSR) matrix by column index
///   within each row.
    static void sortCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                               const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                               const Teuchos::ArrayView<Scalar>&CRS_vals,
                               const UnderlyingLib lib) {

      if(lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
        throw(Xpetra::Exceptions::RuntimeError("Xpetra::CrsMatrixUtils::sortCrsEntries only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)"));
#endif // HAVE_XPETRA_EPETRA
      } else if(lib == Xpetra::UseTpetra) {
        Tpetra::Import_Util::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
      }

      return;
    }

/// \brief Sort and merge the entries of the (raw CSR) matrix by column index
///   within each row.
    static void sortAndMergeCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                                       const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                                       const Teuchos::ArrayView<Scalar>&CRS_vals,
                                       const UnderlyingLib lib) {

      if(lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
        throw(Xpetra::Exceptions::RuntimeError("Xpetra::CrsMatrixUtils::sortAndMergeCrsEntries only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)"));
#endif // HAVE_XPETRA_EPETRA
      } else if(lib == Xpetra::UseTpetra) {
        Tpetra::Import_Util::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
      }

      return;
    }

  }; // end class CrsMatrixUtils

#ifdef HAVE_XPETRA_EPETRA
// Specialization for double, int, int, EpetraNode
  template <>
  class CrsMatrixUtils<double,int,int,EpetraNode> {
    typedef double          Scalar;
    typedef int             LocalOrdinal;
    typedef int             GlobalOrdinal;
    typedef EpetraNode      Node;
#undef XPETRA_CRSMATRIXUTILS_SHORT

  public:

/// \brief Sort the entries of the (raw CSR) matrix by column index
///   within each row.
    static void sortCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                               const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                               const Teuchos::ArrayView<Scalar>&CRS_vals,
                               const UnderlyingLib lib) {

      if(lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
        int rv = Epetra_Util::SortCrsEntries(Teuchos::as<int>(CRS_rowptr.size()-1),
                                             CRS_rowptr.getRawPtr(),
                                             CRS_colind.getRawPtr(),
                                             CRS_vals.getRawPtr());
        if (rv != 0) {
          throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of "
                                         + Teuchos::toString(rv));
        }
#endif // HAVE_XPETRA_EPETRA
      } else if(lib == Xpetra::UseTpetra) {
        Tpetra::Import_Util::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
      }

      return;
    }

/// \brief Sort and merge the entries of the (raw CSR) matrix by column index
///   within each row.
    static void sortAndMergeCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                                       const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                                       const Teuchos::ArrayView<Scalar>&CRS_vals,
                                       const UnderlyingLib lib) {

      if(lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
        int rv = Epetra_Util::SortAndMergeCrsEntries(Teuchos::as<int>(CRS_rowptr.size()-1),
                                                     CRS_rowptr.getRawPtr(),
                                                     CRS_colind.getRawPtr(),
                                                     CRS_vals.getRawPtr());
        if (rv != 0) {
          throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of "
                                         + Teuchos::toString(rv));
        }
#endif // HAVE_XPETRA_EPETRA
      } else if(lib == Xpetra::UseTpetra) {
        Tpetra::Import_Util::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
      }

      return;
    }

  }; // end class CrsMatrixUtils


// Specialization for double, int, long long, EpetraNode
  template <>
  class CrsMatrixUtils<double,int,long long,EpetraNode> {
    typedef double          Scalar;
    typedef int             LocalOrdinal;
    typedef long long       GlobalOrdinal;
    typedef EpetraNode      Node;
#undef XPETRA_CRSMATRIXUTILS_SHORT

  public:

/// \brief Sort the entries of the (raw CSR) matrix by column index
///   within each row.
    static void sortCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                               const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                               const Teuchos::ArrayView<Scalar>&CRS_vals,
                               const UnderlyingLib lib) {

      if(lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
        int rv = Epetra_Util::SortCrsEntries(Teuchos::as<int>(CRS_rowptr.size()-1),
                                             CRS_rowptr.getRawPtr(),
                                             CRS_colind.getRawPtr(),
                                             CRS_vals.getRawPtr());
        if (rv != 0) {
          throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of "
                                         + Teuchos::toString(rv));
        }
#endif // HAVE_XPETRA_EPETRA
      } else if(lib == Xpetra::UseTpetra) {
        Tpetra::Import_Util::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
      }

      return;
    }

/// \brief Sort and merge the entries of the (raw CSR) matrix by column index
///   within each row.
    static void sortAndMergeCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                                       const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                                       const Teuchos::ArrayView<Scalar>&CRS_vals,
                                       const UnderlyingLib lib) {

      if(lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
        int rv = Epetra_Util::SortAndMergeCrsEntries(Teuchos::as<int>(CRS_rowptr.size()-1),
                                                     CRS_rowptr.getRawPtr(),
                                                     CRS_colind.getRawPtr(),
                                                     CRS_vals.getRawPtr());
        if (rv != 0) {
          throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of "
                                         + Teuchos::toString(rv));
        }
#endif // HAVE_XPETRA_EPETRA
      } else if(lib == Xpetra::UseTpetra) {
        Tpetra::Import_Util::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
      }

      return;
    }

  }; // end class CrsMatrixUtils
#endif // HAVE_XPETRA_EPETRA for Epetra scpecialization

} // end namespace Xpetra

#define XPETRA_CRSMATRIXUTILS_SHORT

#endif // PACKAGES_XPETRA_CRSMATRIX_UTILS_HPP_
