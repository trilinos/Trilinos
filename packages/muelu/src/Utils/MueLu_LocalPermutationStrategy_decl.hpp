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
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>

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
     * \param genFactory: const pointer to generating (calling) PermutationFactory // TODO avoid this, not very elegant. Decide which variables have to be generated, give them back per reference to the PermutationFactory.
     */
    void BuildPermutation(const Teuchos::RCP<Matrix> & A, const Teuchos::RCP<const Map> permRowMap, Level & currentLevel, const FactoryBase* genFactory) const;



  private:

    GlobalOrdinal getGlobalDofId(const Teuchos::RCP<Matrix> & A, LocalOrdinal localNodeId, LocalOrdinal localDof) const;
    GlobalOrdinal globalDofId2globalNodeId( const Teuchos::RCP<Matrix> & A, GlobalOrdinal grid ) const;
   };

} // namespace MueLu

#define MUELU_LOCALPERMUTATIONSTRATEGY_SHORT

#endif /* MUELU_LOCALPERMUTATIONSTRATEGY_DECL_HPP_ */
