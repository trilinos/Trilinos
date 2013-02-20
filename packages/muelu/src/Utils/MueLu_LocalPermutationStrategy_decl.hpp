/*
 * MueLu_LocalPermutationStrategy_decl.hpp
 *
 *  Created on: Feb 19, 2013
 *      Author: tobias
 */

#ifndef MUELU_LOCALPERMUTATIONSTRATEGY_DECL_HPP_
#define MUELU_LOCALPERMUTATIONSTRATEGY_DECL_HPP_

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {

  //! @brief Local permutation strategy
  /*!
     This class permutes columns of a input matrix A. No inter-node permutations are allowed,
     only permutations of columns that correspond to DOFs of the same node.
    */

  template<class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::sparseOps>
  class LocalPermutationStrategy : public BaseClass {
#undef MUELU_LOCALPERMUTATIONSTRATEGY_SHORT
#include "MueLu_UseShortNames.hpp"
  public:

  /*!
    @class LocalPermutationStrategy class.
    @brief Class which defines local permutations of matrix columns which correspond to DOFs of the same node.
    */

    //! @name Apply methods.
    //@{

    //! Apply constraint.
    void BuildPermutation(Level & currentLevel) const;

    //@}


  private:
   };

} // namespace MueLu

#define MUELU_LOCALPERMUTATIONSTRATEGY_SHORT

#endif /* MUELU_LOCALPERMUTATIONSTRATEGY_DECL_HPP_ */
