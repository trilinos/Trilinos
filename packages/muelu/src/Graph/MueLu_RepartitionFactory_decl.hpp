#ifndef MUELU_REPARTITIONFACTORY_DECL_HPP
#define MUELU_REPARTITIONFACTORY_DECL_HPP

// Some classes are only used in the definition (_def.hpp) of this class 
// but forward declarations are needed here to enable the UseShortNames mechanism.
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Export_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_OperatorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
 
namespace MueLu {

  /*!
    @class RepartitionFactory class.
    @brief Factory for building permutation matrix that can be be used to shuffle data (matrices, vectors) among processes

    This factory acts on both the number of partitions and a vector (usually created by Zoltan) that indicates to which partitions
    the current level's system matrix's DOFS belong.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class RepartitionFactory : public SingleLevelFactoryBase {
#undef MUELU_REPARTITIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    RepartitionFactory(RCP<const FactoryBase> loadBalancer=Teuchos::null,
                       RCP<const FactoryBase> AFact=Teuchos::null,
                       LO minRowsPerProcessor=1000, SC nnzMaxMinRatio=1.2, GO startLevel=1, LO useDiffusiveHeuristic=0, GO minNnzPerProcessor=-1);

    //! Destructor.
    virtual ~RepartitionFactory();

    //@}

    //! @name Input
    //@{

    /*! @brief Determines the data that RepartitionFactory needs, and the factories that generate that data.

        If this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */
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
  //! @name Get/Set methods.
  //@{

  /*! @brief Suppress any possible repartitioning until specified level.

      Setting this to a very large number will prevent repartitioning from ever happening.
  */
  void SetStartLevel(int startLevel);

  /*! @brief Set imbalance threshold, below which repartitioning is initiatied.
  
  Imbalance is measured by \f$\max_k{N_k} / min_k{N_k}\f$, where \f$N_k\f$ is the number of nonzeros in the local matrix on process \f$k\f$.
  */
  void SetImbalanceThreshold(double threshold);

  /*! @brief Set minimum allowable number of rows on any single process, below which repartitioning is initiated.
      
      This option takes precedence over SetMinNnzPerProcessor.
  */
  void SetMinRowsPerProcessor(GO threshold);

  /*! @brief Set minimum allowable number of nonzeros on any single process, below which repartitioning is initiated.

      This option is ignored if SetMinRowPerProcessor is set.
  */

  /*! @brief @todo Currently does nothing.
  */
  void SetMinNnzPerProcessor(GO threshold);

  //@}

  private:
    //! Load-balancing factory.
    RCP<const FactoryBase> loadBalancer_;
    RCP<const FactoryBase> AFact_;
    //! Minimum number of rows over all processes.  If any process falls below this, repartitioning is initiated.
    LO     minRowsPerProcessor_;
    //! Imbalance threshold, below which repartitioning is initiated.  Imbalance is measured by ratio of maximum nonzeros over all processes to minimum number of nonzeros over all processes.
    double nnzMaxMinRatio_;
    //! First level at which repartitioning can possibly occur.  Repartitioning at finer levels is suppressed.
    int    startLevel_;

    mutable LO useDiffusiveHeuristic_; //FIXME HACK!!!
    //! Minimum number of nonzeros over all processes.  If any process falls below this, repartitioning is initiated.
    GO     minNnzPerProcessor_;

  }; // class RepartitionFactory

} // namespace MueLu

#endif //ifdef HAVE_MPI

#define MUELU_REPARTITIONFACTORY_SHORT
#endif // MUELU_REPARTITIONFACTORY_DECL_HPP
