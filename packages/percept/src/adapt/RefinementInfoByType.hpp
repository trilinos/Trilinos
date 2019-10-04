// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinementInfoByType_hpp
#define adapt_RefinementInfoByType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <stdint.h>
#include <limits>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>



  namespace percept {

    class Refiner;
    typedef uint64_t RefinementInfoCount ;

    struct RefinementInfoByType
    {
      RefinementInfoByType() :
        m_numOrigElems(0),
        m_numNewElems(0),
        m_numOrigNodes(0),
        m_numNewNodes(0),

        m_numOrigElemsLast(0),
        m_numNewElemsLast(0),
        m_rank(0)
      {}

      RefinementInfoCount m_numOrigElems;
      RefinementInfoCount m_numNewElems;
      shards::CellTopology m_topology;
      RefinementInfoCount m_numOrigNodes;
      RefinementInfoCount m_numNewNodes;

      RefinementInfoCount m_numOrigElemsLast;
      RefinementInfoCount m_numNewElemsLast;
      unsigned m_rank;

    };

    struct RefinementInfo
    {
      Refiner *m_refiner;
      std::vector<RefinementInfoByType> m_refinementInfoByType;

      RefinementInfo(Refiner *ref);

      // following is for full stats
      RefinementInfoCount m_totalNumElementsBeforeRefine;
      RefinementInfoCount m_totalNumSidesBeforeRefine;
      RefinementInfoCount m_totalNumElementsAfterRefine;
      RefinementInfoCount m_totalNumSidesAfterRefine;
      RefinementInfoCount m_numMarkedSides;
      RefinementInfoCount m_numMarkedElements;

      RefinementInfoCount m_numberOfMarkedSubDimEntities[4]; // edge=1, face=2, elem=3 (by topo rank)
      RefinementInfoCount m_numberOfSubDimEntities[4];

      RefinementInfoCount m_totalNumEntitiesBeforeFilter[4];
      RefinementInfoCount m_totalNumEntitiesAfterFilter[4];

      RefinementInfoCount m_totalNumEntitiesBeforeEstimate[4];
      RefinementInfoCount m_totalNumEntitiesAfterEstimate[4];

      bool m_full_stats;

      void estimateNew(int iRefinePass);
      void printTable(std::ostream& os, int iRefinePass, bool printAll = false, const std::string& extra_title="");
      void countCurrentNodes(percept::PerceptMesh& eMesh);

      void full_stats_before_mark();
      void full_stats_after_mark();
      void full_stats_before_refine();
      void full_stats_after_refine();
      void full_stats_before_filter(stk::mesh::EntityRank rank, RefinementInfoCount count);
      void full_stats_after_filter(stk::mesh::EntityRank rank, RefinementInfoCount count);
      void full_stats_before_estimate(stk::mesh::EntityRank rank, RefinementInfoCount count);
      void full_stats_after_estimate(stk::mesh::EntityRank rank, RefinementInfoCount count);
      void countRefinedEntities(stk::mesh::EntityRank rank, RefinementInfoCount& num_elem, RefinementInfoCount * num_elem_marked);

    };

  }

#endif
