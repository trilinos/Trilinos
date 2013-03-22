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

#include "Tpetra_ConfigDefs.hpp"

#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)

// definitions
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Map_def.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_MultiVector_def.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_Vector_def.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>
#include <Tpetra_CrsGraph_def.hpp>
#include <Tpetra_CrsMatrixMultiplyOp_decl.hpp>
#include <Tpetra_CrsMatrixMultiplyOp_def.hpp>
#include <Tpetra_MatrixIO_decl.hpp>
#include <Tpetra_MatrixIO_def.hpp>
#include <Ifpack2_BorderedOperator_decl.hpp>
#include <Ifpack2_BorderedOperator_def.hpp>
#include <Ifpack2_Relaxation_decl.hpp>
#include <Ifpack2_Relaxation_def.hpp>
#include <Ifpack2_ILUT_decl.hpp>
#include <Ifpack2_ILUT_def.hpp>
#include <Ifpack2_RILUK_decl.hpp>
#include <Ifpack2_RILUK_def.hpp>
#include <Ifpack2_Diagonal_decl.hpp>
#include <Ifpack2_Diagonal_def.hpp>
#include <Ifpack2_Chebyshev_decl.hpp>
#include <Ifpack2_Chebyshev_def.hpp>

/* these examples require some explicit instantiations that is not enabled in the build of the library 
*/

#if defined(HAVE_TPETRA_QD)
#include <qd/dd_real.h>

#define INSTANT_ALL( NODE ) \
  namespace Tpetra { \
  TPETRA_MULTIVECTOR_INSTANT(dd_real,int,int,NODE) \
  TPETRA_VECTOR_INSTANT(     dd_real,int,int,NODE) \
  TPETRA_CRSMATRIX_INSTANT(dd_real,int,int,NODE) \
  TPETRA_CRSMATRIX_MULTIPLYOP_INSTANT(dd_real,dd_real,int,int,NODE) \
  namespace Utils { \
    TPETRA_MATRIXIO_INSTANT(dd_real,int,int,NODE) \
  }  \
  }  \
  namespace Ifpack2 { \
    template class BorderedOperator<dd_real,int,int,NODE>; \
    template class Relaxation<Tpetra::CrsMatrix<dd_real,int,int,NODE> >; \
    template class ILUT<Tpetra::CrsMatrix<dd_real,int,int,NODE> >; \
    template class RILUK<Tpetra::CrsMatrix<dd_real,int,int,NODE> >; \
    template class Diagonal<Tpetra::CrsMatrix<dd_real,int,int,NODE> >; \
    template class Chebyshev<Tpetra::CrsMatrix<dd_real,int,int,NODE> >; \
  }

INSTANT_ALL(Kokkos::DefaultNode::DefaultNodeType)
#endif

namespace Ifpack2 {
  template class Diagonal<Tpetra::CrsMatrix<double,int,int,Kokkos::SerialNode> >;
}

#endif
