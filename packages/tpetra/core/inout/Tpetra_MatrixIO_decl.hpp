// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIX_IO_DECL
#define TPETRA_MATRIX_IO_DECL

/// \file Tpetra_MatrixIO_decl.hpp
/// \brief Declarations of functions to read a Tpetra::CrsMatrix
///   sparse matrix from a Harwell-Boeing file.
/// \warning DON'T TRUST ANYTHING IN THIS HEADER FILE.  It's probably
///   broken for anything other than real unsymmetric Harwell-Boeing
///   files with a very simple format.

#include <Tpetra_CrsMatrix.hpp>

namespace Tpetra {
  namespace Utils {

    bool parseIfmt (Teuchos::ArrayRCP<char> fmt, int &perline, int &width);
    bool parseRfmt (Teuchos::ArrayRCP<char> fmt, int &perline, int &width, int &prec, char &flag);
    void readHBInfo (const std::string &filename, int &M, int &N, int &nz, Teuchos::ArrayRCP<char> &Type, int &Nrhs);

    void
    readHBHeader (std::ifstream &in_file, Teuchos::ArrayRCP<char> &Title,
                  Teuchos::ArrayRCP<char> &Key, Teuchos::ArrayRCP<char> &Type,
                  int &Nrow, int &Ncol, int &Nnzero, int &Nrhs,
                  Teuchos::ArrayRCP<char> &Ptrfmt,
                  Teuchos::ArrayRCP<char> &Indfmt,
                  Teuchos::ArrayRCP<char> &Valfmt,
                  Teuchos::ArrayRCP<char> &Rhsfmt,
                  int &Ptrcrd, int &Indcrd, int &Valcrd, int &Rhscrd,
                  Teuchos::ArrayRCP<char> &Rhstype);

    void
    readHBMatDouble (const std::string &filename, int &M, int &N, int &nonzeros,
                     std::string &Type, Teuchos::ArrayRCP<int> &colptr,
                     Teuchos::ArrayRCP<int> &rowind,
                     Teuchos::ArrayRCP<double> &val);



    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void
    readHBMatrix (const std::string &filename,
                  const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &A,
                  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = Teuchos::null,
                  const Teuchos::RCP<Teuchos::ParameterList> &params = Teuchos::null);

  } // end of Tpetra::Utils namespace
} // end of Tpetra namespace

#endif
