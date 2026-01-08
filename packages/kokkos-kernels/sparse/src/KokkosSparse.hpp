// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
/// \file KokkosSparse.hpp
/// \brief Public interface to local computational kernels on sparse
///   matrices.
///
/// KokkosSparse::spmv implements local sparse matrix-vector multiply.
/// It computes y = beta*y + alpha*Op(A)*x, where x and y are either
/// both rank 1 (single vectors) or rank 2 (multivectors) Kokkos::View
/// instances, A is a KokkosSparse::CrsMatrix, and Op(A) is determined
/// by the \c mode input (either no transpose, transpose, or conjugate
/// transpose).  If beta == 0, ignore and overwrite the initial
/// entries of y; if alpha == 0, ignore the entries of A and x.
///
/// KokkosSparse::trsv implements local sparse triangular solve.
/// It solves Ax=b, where A is either upper or lower triangular.
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_trsv.hpp"
#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_gauss_seidel.hpp"
#include "KokkosSparse_par_ilut.hpp"
#include "KokkosSparse_gmres.hpp"
