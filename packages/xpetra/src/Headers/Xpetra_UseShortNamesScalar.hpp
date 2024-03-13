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
using CrsMatrix [[maybe_unused]] = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_IO_SHORT
using IO [[maybe_unused]] = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_ITERATOROPS_SHORT
using IteratorOps [[maybe_unused]] = Xpetra::IteratorOps<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_VECTOR_SHORT
using Vector [[maybe_unused]] = Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_BLOCKEDVECTOR_SHORT
using BlockedVector [[maybe_unused]] = Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MULTIVECTOR_SHORT
using MultiVector [[maybe_unused]] = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MATRIX_SHORT
using Matrix [[maybe_unused]] = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MATRIXMATRIX_SHORT
using MatrixMatrix [[maybe_unused]] = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TRIPLEMATRIXMULTIPLY_SHORT
using TripleMatrixMultiply [[maybe_unused]] = Xpetra::TripleMatrixMultiply<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MATRIXUTILS_SHORT
using MatrixUtils [[maybe_unused]] = Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_OPERATOR_SHORT
using Operator [[maybe_unused]] = Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TPETRAOPERATOR_SHORT
using TpetraOperator [[maybe_unused]] = Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TPETRAHALFPRECISIONOPERATOR_SHORT
using TpetraHalfPrecisionOperator [[maybe_unused]] = Xpetra::TpetraHalfPrecisionOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_BLOCKEDCRSMATRIX_SHORT
using BlockedCrsMatrix [[maybe_unused]] = Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_BLOCKEDMULTIVECTOR_SHORT
using BlockedMultiVector [[maybe_unused]] = Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_REORDEREDBLOCKEDMULTIVECTOR_SHORT
using ReorderedBlockedMultiVector [[maybe_unused]] = Xpetra::ReorderedBlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_REORDEREDBLOCKEDCRSMATRIX_SHORT
using ReorderedBlockedCrsMatrix [[maybe_unused]] = Xpetra::ReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef HAVE_XPETRA_THYRA
#ifdef XPETRA_THYRAUTILS_SHORT
using ThyraUtils [[maybe_unused]] = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif
#endif

#ifdef XPETRA_CRSMATRIXWRAP_SHORT
using CrsMatrixWrap [[maybe_unused]] = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_VECTORFACTORY_SHORT
using VectorFactory [[maybe_unused]] = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_CRSMATRIXFACTORY_SHORT
using CrsMatrixFactory [[maybe_unused]] = Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MULTIVECTORFACTORY_SHORT
using MultiVectorFactory [[maybe_unused]] = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MATRIXFACTORY_SHORT
using MatrixFactory [[maybe_unused]] = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MATRIXFACTORY2_SHORT
using MatrixFactory2 [[maybe_unused]] = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TPETRACRSMATRIX_SHORT
using TpetraCrsMatrix [[maybe_unused]] = Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TPETRABLOCKCRSMATRIX_SHORT
using TpetraBlockCrsMatrix [[maybe_unused]] = Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

// TODO remove this
#ifdef XPETRA_EPETRACRSMATRIX_SHORT
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
using EpetraCrsMatrix64 [[maybe_unused]] = Xpetra::EpetraCrsMatrixT<long long, Xpetra::EpetraNode>;
#endif
using EpetraCrsMatrix [[maybe_unused]] = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;  // do we need this???
#endif
// TODO remove above entries

#ifdef XPETRA_TPETRAMULTIVECTOR_SHORT
using TpetraMultiVector [[maybe_unused]] = Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_TPETRAVECTOR_SHORT
using TpetraVector [[maybe_unused]] = Xpetra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MAPEXTRACTOR_SHORT
using MapExtractor [[maybe_unused]] = Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

#ifdef XPETRA_MAPEXTRACTORFACTORY_SHORT
using MapExtractorFactory [[maybe_unused]] = Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#endif

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
using SC [[maybe_unused]] = Scalar;
// TODO: do the same for Epetra object (problem of namespace)
