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


#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#include <Kokkos_View.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>

typedef int LocalOrdinal;
typedef long long GlobalOrdinal;
typedef double Scalar;

typedef Tpetra::Details::DefaultTypes::node_type NT;

typedef Kokkos::DefaultExecutionSpace ExecutionSpace;

typedef Kokkos::View<GlobalOrdinal*,ExecutionSpace> global_ordinal_view_type;

// NOTE: Arrays are hardwired for QUAD4
typedef Kokkos::View<LocalOrdinal*[4],ExecutionSpace> local_ordinal_2d_array_type;
typedef Kokkos::View<GlobalOrdinal*[4],ExecutionSpace> global_ordinal_2d_array_type;
typedef Kokkos::View<Scalar*[4],ExecutionSpace> scalar_2d_array_type;

typedef typename Tpetra::Map<LocalOrdinal, GlobalOrdinal, NT> MapType;
typedef typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, NT> GraphType;
typedef typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, NT> MatrixType;

typedef typename Tpetra::Export<LocalOrdinal, GlobalOrdinal, NT> ExportType;

#endif


