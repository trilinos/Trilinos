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
//
// ***********************************************************************
//
// @HEADER
// Get rid of template parameters

// New definition of types using the types Scalar, LocalOrdinal, GlobalOrdinal, Node of the current context.

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by Scalar, LocalOrdinal, GlobalOrdinal, Node).

#ifdef XPETRA_CRSMATRIX_SHORT
typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrix;
#endif

#ifdef XPETRA_VECTOR_SHORT
typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Vector;
#endif

#ifdef XPETRA_MULTIVECTOR_SHORT
typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
#endif

#ifdef XPETRA_MATRIX_SHORT
typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
#endif

#ifdef XPETRA_OPERATOR_SHORT
typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> Operator;
#endif

#ifdef XPETRA_TPETRAOPERATOR_SHORT
typedef Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraOperator;
#endif

#ifdef XPETRA_BLOCKEDCRSMATRIX_SHORT
typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedCrsMatrix;
#endif

#ifdef XPETRA_CRSMATRIXWRAP_SHORT
typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixWrap;
#endif

#ifdef XPETRA_VECTORFACTORY_SHORT
typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorFactory;
#endif

#ifdef XPETRA_CRSMATRIXFACTORY_SHORT
typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> CrsMatrixFactory;
#endif

#ifdef XPETRA_MULTIVECTORFACTORY_SHORT
typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVectorFactory;
#endif

#ifdef XPETRA_MATRIXFACTORY_SHORT
typedef Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixFactory;
#endif

#ifdef XPETRA_MATRIXFACTORY2_SHORT
typedef Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node> MatrixFactory2;
#endif

#ifdef XPETRA_TPETRACRSMATRIX_SHORT
typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMatrix;
#endif

#ifdef XPETRA_EPETRACRSMATRIX_SHORT
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
typedef Xpetra::EpetraCrsMatrix EpetraCrsMatrix;
#endif
#endif

#ifdef XPETRA_TPETRAMULTIVECTOR_SHORT
typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraMultiVector;
#endif

#ifdef XPETRA_TPETRAVECTOR_SHORT
typedef Xpetra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraVector;
#endif

#ifdef XPETRA_MAPEXTRACTOR_SHORT
typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;
#endif

#ifdef XPETRA_MAPEXTRACTORFACTORY_SHORT
typedef Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorFactory;
#endif

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
typedef Scalar    SC;
// TODO: do the same for Epetra object (problem of namespace)
