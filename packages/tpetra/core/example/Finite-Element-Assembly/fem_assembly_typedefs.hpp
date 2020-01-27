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
#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_TYPEDEFS_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_TYPEDEFS_HPP

#include "Kokkos_View.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_FECrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_FECrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_FEMultiVector.hpp"

namespace TpetraExamples {

using local_ordinal_type = Tpetra::Map<>::local_ordinal_type;
using global_ordinal_type = Tpetra::Map<>::global_ordinal_type;
using execution_space = Tpetra::Map<>::device_type::execution_space;

using map_type = Tpetra::Map<>;
using crs_graph_type = Tpetra::CrsGraph<>;
using fe_graph_type = Tpetra::FECrsGraph<>;
using Scalar = Tpetra::CrsMatrix<>::scalar_type;
using crs_matrix_type = Tpetra::CrsMatrix<Scalar>;
using fe_matrix_type = Tpetra::FECrsMatrix<Scalar>;

using import_type = Tpetra::Import<>;
using export_type = Tpetra::Export<>;
using multivector_type = Tpetra::MultiVector<Scalar>;
using fe_multivector_type = Tpetra::FEMultiVector<Scalar>;


using global_ordinal_view_type =
  Kokkos::View<global_ordinal_type*, execution_space>;
using local_ordinal_view_type =
  Kokkos::View<local_ordinal_type*, execution_space>;
using scalar_1d_array_type = Kokkos::View<Scalar*, execution_space>;
using bool_1d_array_type = Kokkos::View<bool*, execution_space>;

// NOTE: Arrays are hardwired for QUAD4
using local_ordinal_2d_array_type =
  Kokkos::View<local_ordinal_type*[4], execution_space>;
using global_ordinal_2d_array_type =
  Kokkos::View<global_ordinal_type*[4], execution_space>;
using scalar_2d_array_type = Kokkos::View<Scalar*[4], execution_space>;


}

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_TYPEDEFS_HPP

