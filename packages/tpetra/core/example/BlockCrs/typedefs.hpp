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


#ifndef __TYPEDEFS_HPP__
#define __TYPEDEFS_HPP__

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockMultiVector.hpp>

namespace BlockCrsTest {
  
  typedef double value_type;
  typedef double magnitude_type;
  
  typedef Tpetra::Map<> map_type;
  typedef typename map_type::local_ordinal_type local_ordinal_type;
  typedef typename map_type::global_ordinal_type global_ordinal_type;
  typedef typename map_type::node_type node_type;
  
  typedef Tpetra::Import<> tpetra_import_type;
  typedef Tpetra::Export<> tpetra_export_type;
  typedef Tpetra::MultiVector<value_type> tpetra_multivector_type;
  typedef Tpetra::MultiVector<magnitude_type> tpetra_multivector_magnitude_type;

  typedef Tpetra::CrsGraph<> tpetra_crs_graph_type;
  typedef Tpetra::RowMatrix<value_type> tpetra_rowmatrix_type;
  typedef Tpetra::BlockCrsMatrix<value_type> tpetra_blockcrs_type;

  typedef Kokkos::DefaultExecutionSpace exec_space;
  typedef Kokkos::DefaultHostExecutionSpace host_space;

  typedef Kokkos::pair<local_ordinal_type,local_ordinal_type> local_ordinal_range_type;
  typedef Kokkos::pair<global_ordinal_type,global_ordinal_type> global_ordinal_range_type;

  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;
}

#endif
  

