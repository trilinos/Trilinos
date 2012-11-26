// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_LOCALAGGREGATIONALGORITHM_DECL_HPP
#define MUELU_LOCALAGGREGATIONALGORITHM_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_LocalAggregationAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"

#include "MueLu_AggOptions.hpp" // includes Ordering enum

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {

  using namespace AggOptions; // necessary

  /* ************************************************************************* */
  /* definition of the structure from ML for holding aggregate information     */
  /* ------------------------------------------------------------------------- */
  typedef struct MueLu_SuperNode_Struct
  {
    int    length;
    int    maxLength;
    int    index;
    Teuchos::ArrayRCP<int> list;
    struct MueLu_SuperNode_Struct *next;
  } MueLu_SuperNode;

  /* In the algorithm, aggStat[]=READY/NOTSEL/SELECTED indicates whether a node has been aggregated. */
  enum NodeState {
    READY   = -11,   /* indicates that a node is available to be */
    /* selected as a root node of an aggregate  */

    NOTSEL  = -12,   /* indicates that a node has been rejected  */
    /* as a root node. This could perhaps be    */
    /* because if this node had been selected a */
    /* small aggregate would have resulted.     */

    SELECTED = -13   /* indicates that a node has been assigned  */
    /* to an aggregate.                         */
  };

  /*!
    @class LocalAggregationAlgorithm class.
    @brief Algorithm for coarsening a graph with uncoupled aggregation.

    This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
    that can include unknowns from more than one process.

  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class LocalAggregationAlgorithm : public BaseClass {
#undef MUELU_LOCALAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t; //TODO

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    LocalAggregationAlgorithm();

    //! Destructor.
    virtual ~LocalAggregationAlgorithm() { }

    //@}

    //! @name Set/get methods.
    //@{

    void SetOrdering(Ordering ordering)                          { ordering_                = ordering;                }
    void SetMinNodesPerAggregate(int minNodesPerAggregate)       { minNodesPerAggregate_    = minNodesPerAggregate;    }
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { maxNeighAlreadySelected_ = maxNeighAlreadySelected; }

    Ordering GetOrdering()                const { return ordering_;                }
    int      GetMinNodesPerAggregate()    const { return minNodesPerAggregate_;    }
    int      GetMaxNeighAlreadySelected() const { return maxNeighAlreadySelected_; }

    //@}

    //! @name Aggregation methods.
    //@{

    /*! @brief Local aggregation. */
    void CoarsenUncoupled(Graph const & graph, Aggregates & aggregates) const; // CoarsenUncoupled

  private:
    //! Aggregation options (TODO: Teuchos::ParameterList?)
    Ordering ordering_;                /**<  natural, random, graph           */
    int      minNodesPerAggregate_;    /**<  aggregate size control           */
    int      maxNeighAlreadySelected_; /**<  complexity control               */

    //! @name Utilities
    //@{

    /*! @brief Utility to take a list of integers and reorder them randomly (by using a local permutation).
      @param list On input, a bunch of integers. On output, the same integers in a different order
      that is determined randomly.
    */
    void RandomReorder(Teuchos::ArrayRCP<LO> list) const;

    /*! @brief Generate a random number in the range [min, max] */
    int RandomOrdinal(int min, int max) const;

    //@}

  }; //class LocalAggregationFactory

} //namespace MueLu

#define MUELU_LOCALAGGREGATIONALGORITHM_SHORT
#endif // MUELU_LOCALAGGREGATIONALGORITHM_DECL_HPP
