/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#include "Kokkos_Sparse_impl_MKL.hpp"

namespace KokkosSparse {
namespace Impl {
namespace Mkl {

//! Human-readable string representation of MKL's index base enum.
std::string
sparseIndexBaseToString (const sparse_index_base_t base)
{
#ifdef HAVE_TPETRAKERNELS_MKL
  if (base == SPARSE_INDEX_BASE_ZERO) {
    return "SPARSE_INDEX_BASE_ZERO";
  }
  else if (base == SPARSE_INDEX_BASE_ONE) {
    return "SPARSE_INDEX_BASE_ONE";
  }
  else {
    throw std::invalid_argument ("Invalid sparse_index_base_t value!");
  }
#else
  return "DID NOT BUILD WITH MKL";
#endif // HAVE_TPETRAKERNELS_MKL
}

std::string
sparseOperationToString (const sparse_operation_t op)
{
#ifdef HAVE_TPETRAKERNELS_MKL
  if (op == SPARSE_OPERATION_NON_TRANSPOSE) {
    return "SPARSE_OPERATION_NON_TRANSPOSE";
  }
  else if (op == SPARSE_OPERATION_TRANSPOSE) {
    return "SPARSE_OPERATION_TRANSPOSE";
  }
  else if (op == SPARSE_OPERATION_CONJUGATE_TRANSPOSE) {
    return "SPARSE_OPERATION_CONJUGATE_TRANSPOSE";
  }
  else {
    throw std::invalid_argument ("Invalid sparse_operation_t value!");
  }
#else
  return "DID NOT BUILD WITH MKL";
#endif // HAVE_TPETRAKERNELS_MKL
}

std::string
sparseMatrixTypeToString (const sparse_matrix_type_t type)
{
#ifdef HAVE_TPETRAKERNELS_MKL
  if (type == SPARSE_MATRIX_TYPE_GENERAL) {
    return "SPARSE_MATRIX_TYPE_GENERAL";
  }
  else if (type == SPARSE_MATRIX_TYPE_SYMMETRIC) {
    return "SPARSE_MATRIX_TYPE_SYMMETRIC";
  }
  else if (type == SPARSE_MATRIX_TYPE_HERMITIAN) {
    return "SPARSE_MATRIX_TYPE_HERMITIAN";
  }
  else if (type == SPARSE_MATRIX_TYPE_TRIANGULAR) {
    return "SPARSE_MATRIX_TYPE_TRIANGULAR";
  }
  else if (type == SPARSE_MATRIX_TYPE_DIAGONAL) {
    return "SPARSE_MATRIX_TYPE_DIAGONAL";
  }
  else if (type == SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR) {
    return "SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR";
  }
  else if (type == SPARSE_MATRIX_TYPE_BLOCK_DIAGONAL) {
    return "SPARSE_MATRIX_TYPE_BLOCK_DIAGONAL";
  }
  else {
    throw std::invalid_argument ("Invalid sparse_matrix_type_t value!");
  }
#else
  return "DID NOT BUILD WITH MKL";
#endif // HAVE_TPETRAKERNELS_MKL
}

std::string
sparseFillModeToString (const sparse_fill_mode_t mode)
{
#ifdef HAVE_TPETRAKERNELS_MKL
  if (mode == SPARSE_FILL_MODE_LOWER) {
    return "SPARSE_FILL_MODE_LOWER";
  }
  else if (mode == SPARSE_FILL_MODE_UPPER) {
    return "SPARSE_FILL_MODE_UPPER";
  }
  else {
    throw std::invalid_argument ("Invalid sparse_fill_mode_t value!");
  }
#else
  return "DID NOT BUILD WITH MKL";
#endif // HAVE_TPETRAKERNELS_MKL
}

std::string
sparseDiagTypeToString (const sparse_diag_type_t diag)
{
#ifdef HAVE_TPETRAKERNELS_MKL
  if (diag == SPARSE_DIAG_NON_UNIT) {
    return "SPARSE_DIAG_NON_UNIT";
  }
  else if (diag == SPARSE_DIAG_UNIT) {
    return "SPARSE_DIAG_UNIT";
  }
  else {
    throw std::invalid_argument ("Invalid sparse_diag_type_t value!");
  }
#else
  return "DID NOT BUILD WITH MKL";
#endif // HAVE_TPETRAKERNELS_MKL
}

std::string
matrixDescriptorToString (const matrix_descr& descr,
                          const std::string inputIndent,
                          const bool oneLine,
                          const bool printHeader)
{
  using std::endl;
  std::ostringstream out;

  if (printHeader) {
    out << inputIndent << "MKL matrix descriptor:";
  }
  std::string indent = inputIndent;
  if (oneLine) {
    if (printHeader) {
      out << " {";
    }
    else {
      out << inputIndent << "{";
    }
  }
  else {
    if (printHeader) {
      out << endl;
    }
    indent = inputIndent + " ";
  }

  out << indent << "Built with MKL: ";
#ifdef HAVE_TPETRAKERNELS_MKL
  out << "YES";
#else
  out << "NO";
#endif // HAVE_TPETRAKERNELS_MKL
  if (oneLine) {
    out << ", ";
  }
  else {
    out << endl;
  }

  out << indent << "Type: " << sparseMatrixTypeToString (descr.type);
  if (oneLine) {
    out << ", ";
  }
  else {
    out << endl;
  }

  out << indent << "Mode: " << sparseFillModeToString (descr.mode);
  if (oneLine) {
    out << ", ";
  }
  else {
    out << endl;
  }

  out << indent << "Diag: " << sparseDiagTypeToString (descr.diag);
  if (oneLine) {
    out << "}";
  }
  else {
    out << endl;
  }

  return out.str ();
}

std::string
sparseLayoutToString (const sparse_layout_t layout)
{
#ifdef HAVE_TPETRAKERNELS_MKL
  if (layout == SPARSE_LAYOUT_COLUMN_MAJOR) {
    return "SPARSE_LAYOUT_COLUMN_MAJOR";
  }
  if (layout == SPARSE_LAYOUT_ROW_MAJOR) {
    return "SPARSE_LAYOUT_ROW_MAJOR";
  }
  else {
    throw std::invalid_argument ("Invalid sparse_layout_t value!");
  }
#else
  return "DID NOT BUILD WITH MKL";
#endif // HAVE_TPETRAKERNELS_MKL
}

} // namespace Mkl
} // namespace Impl
} // namespace KokkosSparse



