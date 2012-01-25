#ifndef MUELU_REPARTITION_DECL_HPP
#define MUELU_REPARTITION_DECL_HPP

// Some classes are only used in the definition (_def.hpp) of this class 
// but forward declarations are needed here to enable the UseShortNames mechanism.
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Export_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_OperatorFactory_fwd.hpp>

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

    /*!
      @brief Determine which process should own each partition.

      Partitions are assigned to processes in order to minimize data movement.  The basic idea is that a good choice for partition
      owner is to choose the pid that already has the greatest number of nonzeros for a particular partition.

      @param[in]  currentLevel      The current multigrid level's Level object.
      @param[out] myPartitionNumber The partition number this PID owns, otherwise -1.
      @param[out] partitionOwner    An Array (fancy std::vector) such that the PID of the process that owns partition i is given by partitionOwner[i].
    */
    void DeterminePartitionPlacement(Level & currentLevel, GlobalOrdinal &myPartitionNumber, Array<int> &partitionOwner) const;

  private:

  }; // class Repartition

} // namespace MueLu

#define MUELU_REPARTITION_SHORT
#endif // MUELU_REPARTITION_DECL_HPP
