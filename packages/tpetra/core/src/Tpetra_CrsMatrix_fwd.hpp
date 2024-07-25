// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIX_FWD_HPP
#define TPETRA_CRSMATRIX_FWD_HPP

#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_CrsMatrix_fwd.hpp
/// \brief Forward declaration of Tpetra::CrsMatrix

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
template<class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
         class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
         class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
         class Node = ::Tpetra::Details::DefaultTypes::node_type>
class CrsMatrix;

/// \class is_crs_matrix
/// \brief is_crs_matrix<T>::value is true if T is a CrsMatrix<...>, false
/// otherwise
template <typename>
struct is_crs_matrix : public std::false_type {};
template <typename... P>
struct is_crs_matrix<CrsMatrix<P...>> : public std::true_type {};
template <typename... P>
struct is_crs_matrix<const CrsMatrix<P...>> : public std::true_type {};

/// \brief Equivalent to is_crs_matrix<T>::value.
template <typename T>
inline constexpr bool is_crs_matrix_v = is_crs_matrix<T>::value;

} // namespace Tpetra
#endif // DOXYGEN_SHOULD_SKIP_THIS



#endif // TPETRA_CRSMATRIX_FWD_HPP
