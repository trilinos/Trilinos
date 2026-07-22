// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MATRIX_MARKET_HPP
#define STOKHOS_MATRIX_MARKET_HPP

namespace Stokhos {

template < typename MatrixType > class MatrixMarketWriter;

template < typename MatrixType >
void write_matrix_market(const MatrixType& A ,
                         const std::string& filename)
{
  MatrixMarketWriter<MatrixType>::write(A, filename);
}

} // namespace Stokhos

#endif
