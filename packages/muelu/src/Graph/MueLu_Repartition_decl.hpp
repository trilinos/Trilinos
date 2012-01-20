#ifndef MUELU_REPARTITION_DECL_HPP
#define MUELU_REPARTITION_DECL_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_OperatorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Repartition_fwd.hpp"

namespace MueLu {

  /*!
    @class Repartition class.
    @brief Factory for building permutation matrix that can be be used to shuffle data (matrices, vectors) among processes

    This factory acts on both the number of partitions and a vector (usually created by Zoltan) that indicates to which partitions
    the current level's system matrix's DOFS belong.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class Repartition : public SingleLevelFactoryBase {
#undef MUELU_REPARTITION_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    Repartition();

    //! Destructor.
    virtual ~Repartition();

    //@}

    //! @name Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Build methods.
    //@{

    //! Build an object with this factory.
    void Build(Level & currentLevel) const;

    //@}


    //! @name Helper methods.
    //@{

    /*! Determine which process should own each partition.

      For right now, we assign the partitions 0..N to pids 0..N, respectively.

      TODO: This should be changed in order to minimize data movement.  A good choice for partition owner is to
      choose the pid that already has the greatest number of nonzeros for a particular partition.

      @param myPartitionNumber On output, either the partition number this PID owns, otherwise -1.
      @param partitionOwner    On output, an Array (fancy std::vector) such that PID partitionOwner[i] is the owner of partition i.
    */
    void DeterminePartitionPlacement(Level & currentLevel, GlobalOrdinal &myPartitionNumber, Array<int> &partitionOwner) const;

  private:

  }; // class Repartition

} // namespace MueLu

#define MUELU_REPARTITION_SHORT
#endif // MUELU_REPARTITION_DECL_HPP
