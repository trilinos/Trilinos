// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_AlgebraicPermutationStrategy_decl.hpp
 *
 *  Created on: Feb 20, 2013
 *      Author: tobias
 */

#ifndef MUELU_ALGEBRAICPERMUTATIONSTRATEGY_DECL_HPP_
#define MUELU_ALGEBRAICPERMUTATIONSTRATEGY_DECL_HPP_

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_Export_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {

// template struct for comparing pairs
template <class Scalar = DefaultScalar, class LocalOrdinal = DefaultGlobalOrdinal>
struct CompPairs {
  CompPairs(const std::vector<Scalar>& v)
    : vinternal_(v) {}
  std::vector<Scalar> vinternal_;
  bool operator()(LocalOrdinal a, LocalOrdinal b) {
    // return vinternal_[a] < vinternal_[b];
    return Teuchos::ScalarTraits<Scalar>::magnitude(vinternal_[a]) > Teuchos::ScalarTraits<Scalar>::magnitude(vinternal_[b]);
  }
};

// template function for comparison
template <class Scalar, class LocalOrdinal>
CompPairs<Scalar, LocalOrdinal> CreateCmpPairs(const std::vector<Scalar>& v) {
  return CompPairs<Scalar, LocalOrdinal>(v);
}

// template function for sorting permutations
template <class Scalar, class LocalOrdinal>
void sortingPermutation(const std::vector<Scalar>& values, std::vector<LocalOrdinal>& v) {
  size_t size = values.size();
  v.clear();
  v.reserve(size);
  for (size_t i = 0; i < size; ++i)
    v.push_back(i);

  std::sort(v.begin(), v.end(), MueLu::CreateCmpPairs<Scalar, LocalOrdinal>(values));
}

//! @brief Algebraic permutation strategy
/*!
   This class permutes columns of a input matrix A trying to make A
   a diagonal dominant matrix.

  */

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AlgebraicPermutationStrategy : public BaseClass {
#undef MUELU_ALGEBRAICPERMUTATIONSTRATEGY_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  /*!
    @class AlgebraicPermutationStrategy class.
    @brief Class which defines local permutations of matrix columns.
    */

  //! @name build permutation methods.
  //@{

  //! build permutation operators
  /*!
   *  The following variables produced
   *  "A"     :      permuted and scaled A
   *  "permA" :      permuted A without scaling
   *  "permP" :      permutation opertor (should be identity)
   *  "permQT":      transpose permutation operators
   *  "permScaling": scaling operator
   *
   * \param A: input matrix (input)
   * \param permRowMap: Dof row map permutation shall be restricted on (input)
   * \param currentLevel: only for output of variables
   * \param genFactory: const pointer to generating (calling) PermutationFactory
   // TODO avoid this, not very elegant.
   // Decide which variables have to be generated, give them back per reference to the PermutationFactory.
   */
  void BuildPermutation(const Teuchos::RCP<Matrix>& A, const Teuchos::RCP<const Map>& permRowMap,
                        Level& currentLevel, const FactoryBase* genFactory) const;

  //@}

 private:
};

}  // namespace MueLu

#define MUELU_ALGEBRAICPERMUTATIONSTRATEGY_SHORT

#endif /* MUELU_ALGEBRAICPERMUTATIONSTRATEGY_DECL_HPP_ */
