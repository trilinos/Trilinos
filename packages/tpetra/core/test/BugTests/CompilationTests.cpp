/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>

// include all the Tpetra headers twice to check if headers are protected by #ifndef FILE_HPP #define ... #endif
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrixSolveOp.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Directory.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <TpetraExt_MMHelpers.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_MatrixIO.hpp>

#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrixSolveOp.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Directory.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <TpetraExt_MMHelpers.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_MatrixIO.hpp>

// This test does not compile if Tpetra objects are declared outside
// of the Tpetra namespace.
typedef int BinaryFunctorAdapter;
typedef int BinaryFunctorAdapterWithAlphaBeta;
typedef int BinaryOp;
typedef int CrsGraph;
typedef int CrsMatrix;
typedef int CrsMatrixMultiplyOp;
typedef int CrsMatrixSolveOp;
typedef int DefaultPlatform;
typedef int Directory;
typedef int DistObject;
typedef int Distributor;
typedef int Export;
typedef int FDStencil;
typedef int Import;
typedef int KernelOp;
typedef int Map;
typedef int MpiPlatform;
typedef int mprec_mult;
typedef int MultiVector;
typedef int OneOp;
typedef int Operator;
typedef int ReductionGlob;
typedef int RowGraph;
typedef int RowMatrix;
typedef int RowMatrixTransposer;
typedef int RTIPreTransformReductionAdapter;
typedef int RTIReductionAdapter;
typedef int ScaleKernel;
typedef int SerialPlatform;
typedef int StdOpKernel;
typedef int TeuchosValueTypeReductionOpAdapter;
typedef int TransformReductionGlob;
typedef int UnaryFunctorAdapter;
typedef int Vector;
typedef int ZeroOp;

TEUCHOS_UNIT_TEST( Compilation, Bug_ClassDeclarationOutsideOfNamespace )
{

}



