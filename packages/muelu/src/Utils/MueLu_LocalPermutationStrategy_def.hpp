/*
 * MueLu_LocalPermutationStrategy_def.hpp
 *
 *  Created on: Feb 19, 2013
 *      Author: tobias
 */

#ifndef MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_
#define MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsGraph.hpp>

#include "MueLu_LocalPermutationStrategy_decl.hpp"

namespace MueLu {

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildPermutation(Level & currentLevel) const {
  }

} // namespace MueLu


#endif /* MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_ */
