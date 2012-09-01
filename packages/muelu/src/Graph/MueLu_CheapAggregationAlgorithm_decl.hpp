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
/*
 * MueLu_CheapAggregationAlgorithm_decl.hpp
 *
 *  Created on: Jul 25, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CHEAPAGGREGATIONALGORITHM_DECL_HPP_
#define MUELU_CHEAPAGGREGATIONALGORITHM_DECL_HPP_

#include <bitset>
#include <vector>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_CheapAggregationAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {

  namespace AggOptions {
    /* Options defining how to pick-up the next root node in the local aggregation procedure */
    enum Ordering {
      NATURAL = 0, /* node ordering   */
      RANDOM  = 1, /* random ordering */
      GRAPH   = 2  /* graph ordering  */
    };
  } // namespace AggOptions

  using namespace AggOptions;

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

  class Aggregate {
  public:
    int length;                   // current size of aggregate
    int maxLength;                // max size of aggregate
    int index;                    // local aggregate id
    std::vector<int> list;  // list of node ids in aggregate
  };

  /* In the algorithm, aggStat[]=READY/NOTSEL/SELECTED indicates whether a node has been aggregated. */
  enum NodeState {
    READY   = -11,   /* indicates that a node is available to be */
    /* selected as a root node of an aggregate  */

    NOTSEL  = -12,   /* indicates that a node has been rejected  */
    /* as a root node. This could perhaps be    */
    /* because if this node had been selected a */
    /* small aggregate would have resulted.     */

    SELECTED = -13,  /* indicates that a node has been assigned  */
    /* to an aggregate.                         */

    BDRY = -15, /* indicates that a node is a Dirichlet bdry node */

    READY_1PT = -16, /* indicates that a node is ready to be aggregates (as a single point aggregate only) */
    SELECTED_1PT = -17 /* indicates that a 1pt aggregate node is already aggregated */
  };

  enum NodeState2 {
	NODEAGGREGATED = 0x01,
	NODENOTSEL     = 0x02,
	NODEONEPT      = 0x04,
	NODESELECTED   = (NODENOTSEL) | (NODEAGGREGATED)
  };

  /*!
    @class CheapAggregationAlgorithm class.
    @brief Algorithm for coarsening a graph with uncoupled aggregation.
  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class CheapAggregationAlgorithm : public BaseClass {
#undef MUELU_CHEAPAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t; //TODO

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    CheapAggregationAlgorithm(RCP<const FactoryBase> const &graphFact = Teuchos::null);

    //! Destructor.
    virtual ~CheapAggregationAlgorithm() { }

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

    LocalOrdinal PhaseOnePt(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat, Teuchos::ArrayRCP<unsigned int> & coarse_aggStat) const;
    LocalOrdinal Phase1(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat, Teuchos::ArrayRCP<unsigned int> & coarse_aggStat) const; // local uncoupled coarsening (Phase 1b)
    LocalOrdinal Phase2_maxlink(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat, Teuchos::ArrayRCP<unsigned int> & coarse_aggStat) const;
    LocalOrdinal Phase3(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat, Teuchos::ArrayRCP<unsigned int> & coarse_aggStat) const; // local uncoupled coarsening (Phase 3)


  private:

    void PrintAggregationInformation(const std::string phase, Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat) const;

    //! Aggregation options (TODO: Teuchos::ParameterList?)
    Ordering ordering_;                /**<  natural, random, graph           */
    int      minNodesPerAggregate_;    /**<  aggregate size control           */
    int      maxNeighAlreadySelected_; /**<  complexity control               */

    /*! @brief Utility to take a list of integers and reorder them randomly (by using a local permutation).
      @param list On input, a bunch of integers. On output, the same integers in a different order
      that is determined randomly.
    */
    void RandomReorder(Teuchos::ArrayRCP<LO> list) const; 

    /*! @brief Generate a random number in the range [min, max] */
    int RandomOrdinal(int min, int max) const;

    //@}

  }; //class CheapAggregationAlgorithm

} //namespace MueLu

#endif /* MUELU_CHEAPAGGREGATIONALGORITHM_DECL_HPP_ */
