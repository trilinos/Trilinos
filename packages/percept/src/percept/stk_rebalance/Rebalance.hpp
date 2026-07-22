// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef stk_rebalance_Rebalance_hpp
#define stk_rebalance_Rebalance_hpp

#include <string>
#include <percept/PerceptMesh.hpp>


#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <percept/stk_rebalance/Partition.hpp>

/** \file Rebalance.hpp
 *   \ingroup stk_rebalance_module
 *
 * \brief Static functions for dynamic load balancing.
 *
 *  The rebalance namespace is intended to provide an application
 *  the top level functions to perform a mesh redistribution.
 *  The action is controlled by instances of the Partition
 *  class that is passed into the rebalance function.
 */

namespace stk {
  namespace rebalance {

    class Rebalance {
    public:
      virtual ~Rebalance() = default;

      /** \brief Rebalance with a Partition object.
       *
       * \param bulk_data      BulkData must be in a parallel consistent state.
       *
       * \param selector       Used to select a subset of mesh entities to compute measure.
       *
       * \param coord_ref      The field containing the nodal coordinates. For the default
       *                       ZoltanPartition class in stk::reblance, this should be non-NULL.
       *
       * \param elem_weight_ref This field will be used by the \a Partition class and
       *                        can be NULL.
       *
       * \param Partition       The base class of a derived class that is used to
       *                        determine the new partition.  See the \a ZoltanPartition
       *                        class for an example.
       *
       * \param rank            Rank of the entities \a elem_weight_ref is defined on.
       *
       * \param entities_list   If non-zero, use this list instead of BulkData::buckets(rank)
       *
       * \param do_rebal        You can call rebalance first with do_rebal=true, then again
       *                          with do_rebal=false and debug=true to print the imbalance
       *                          factors from within Zoltan.
       *
       * This \a rebalance function will use the \a Partition object passed
       * to perform the rebalancing.  It will be necessary to use one of the
       * pre-defined derived classes in stk::rebalance, like \a ZoltanPartition,
       * or to define your own.
       */
      virtual bool rebalance(mesh::BulkData & bulk_data ,
                             const mesh::Selector & selector ,
                             const stk::mesh::FieldBase * coord_ref ,
                             const stk::mesh::FieldBase * elem_weight_ref,
                             Partition & partition,
                             const stk::mesh::EntityRank rank = stk::mesh::InvalidEntityRank,
                             stk::mesh::EntityVector *entities_list = 0,
                             bool do_rebal = true,
                             bool debug = false);

      /** \} */

    protected:
      virtual bool balance_comm_spec_domain( Partition * partition,
                                     mesh::EntityProcVec & rebal_spec );

      virtual void rebalance_dependent_entities( const mesh::BulkData    & bulk_data ,
                                         const Partition         * partition,
                                         const mesh::EntityRank  & dep_rank,
                                         mesh::EntityProcVec     & entity_procs,
                                         const stk::mesh::EntityRank rank);

      virtual bool full_rebalance(mesh::BulkData  & bulk_data ,
                                  Partition       * partition,
                                  const stk::mesh::EntityRank rank);

    };

  }
} // namespace stk

#endif
