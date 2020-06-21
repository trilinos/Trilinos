// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef adapt_NodeRegistryTypes_hpp
#define adapt_NodeRegistryTypes_hpp


#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <utility>
#include <math.h>
#include <unordered_map>
#include <set>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_mesh/base/EntityKey.hpp>

#include <percept/stk_mesh.hpp>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <percept/NoMallocArray.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>

#include <percept/PerceptBoostArray.hpp>

#include <tuple>

#include <adapt/SubDimCell.hpp>
#include <adapt/SDCEntityType.hpp>
#include <adapt/NodeIdsOnSubDimEntityType.hpp>

namespace percept {
// pair of rank and number of entities of that rank needed on a SubDimCell
   struct NeededEntityType
   {
     stk::mesh::EntityRank first;  // e.g. EDGE_RANK if edges are needed to be marked
     unsigned              second; // number of new nodes needed on this rank (e.g. a quadratic edge needs 2)
     std::vector<int>      third;  // special case: for non-homogeneous topos, like wedge/pyramid, say which sub-dim entities get marked
     NeededEntityType(stk::mesh::EntityRank f = stk::topology::NODE_RANK, unsigned s = 0)
     : first(f), second(s), third(0) {}
   };

   inline std::ostream &operator<<(std::ostream& out, const std::array<stk::mesh::EntityId, 1>& arr)
   {
     out << arr[0];
     return out;
   }

#define DOUBLE2LEN 2
   typedef std::array<double, DOUBLE2LEN> Double2;

   inline std::ostream &operator<<(std::ostream& out, const Double2 arr)
   {
       for(unsigned i=0;i<DOUBLE2LEN;i++)
           out << arr[i] << " ";
       return out;
   }

   // tuple storage: SDC_DATA_GLOBAL_NODE_IDS, SDC_DATA_OWNING_ELEMENT_KEY,  SDC_DATA_OWNING_SUBDIM_RANK, SDC_DATA_OWNING_SUBDIM_ORDINAL, SDC_DATA_SPACING
   typedef std::tuple<NodeIdsOnSubDimEntityType, stk::mesh::EntityKey, unsigned char, unsigned char, Double2> SubDimCellData;

   enum { MAX_NODES_ON_A_FACE = 4 };
   typedef MySubDimCell<SDCEntityType, MAX_NODES_ON_A_FACE, CompareSDCEntityType> SubDimCell_SDCEntityType;

   template<>
   inline
   SubDimCell_SDCEntityType& SubDimCell_SDCEntityType::operator=(const SubDimCell_SDCEntityType& from)
   {
       m_eMesh =  from.m_eMesh;
       base_type::m_hash = from.getHash();
       base_type::m_HashCode = from.m_HashCode;
       base_type::m_CompareClass = from.m_CompareClass;
       m_size = from.m_size;
       for(unsigned i=0;i<=m_size;i++)
           m_data[i] = from.m_data[i];

       return *this;
   }

   inline std::ostream& operator<<(std::ostream& out,  SubDimCellData& val)
   {
     std::ostringstream ostr;
     ostr << "SDC:: node ids= [";
     for (unsigned ii=0; ii < std::get<SDC_DATA_GLOBAL_NODE_IDS>(val).size(); ++ii)
       {
         ostr << " " << std::get<SDC_DATA_GLOBAL_NODE_IDS>(val).m_entity_id_vector[ii];
       }
     ostr << "] owning element rank= " << std::get<SDC_DATA_OWNING_ELEMENT_KEY>(val).rank()
         << " id= " << std::get<SDC_DATA_OWNING_ELEMENT_KEY>(val).id()
         << " subDim-ord= " << (int)std::get<SDC_DATA_OWNING_SUBDIM_ORDINAL>(val)
         << " subDim-rank= " << (int)std::get<SDC_DATA_OWNING_SUBDIM_RANK>(val)
         << " spacing info= " << std::get<SDC_DATA_SPACING>(val)[0] << " " << std::get<SDC_DATA_SPACING>(val)[1];
     out << ostr.str() << std::endl;
     return out;
   }

   /// map of the node ids on a sub-dim entity to the data on the sub-dim entity
   typedef std::unordered_map<SubDimCell_SDCEntityType, SubDimCellData, my_fast_hash<SDCEntityType, 4>, my_fast_equal_to<SDCEntityType, 4> > SubDimCellToDataMap;

   // Size and rank of sub-dim cells needing new nodes, actual nodes' EntityKeys stored in a SubDimCell<EntityKey>
   enum CommDataTypeEnum {
     CDT_SUBDIM_ENTITY_SIZE,
     CDT_SUBDIM_ENTITY_RANK,
     CDT_SUBDIM_ENTITY
   };

   // decode: size of SubDimCell, SubDimCell, sub-dim entity rank, non-owning element RANK
   typedef std::tuple<unsigned, stk::mesh::EntityRank, SubDimCell<stk::mesh::EntityKey> > CommDataType;

   enum NodeRegistryState {
     NRS_NONE,
     NRS_START_REGISTER_NODE,
     NRS_END_REGISTER_NODE,
     NRS_START_CHECK_FOR_REMOTE,
     NRS_END_CHECK_FOR_REMOTE,
     NRS_START_PARALLEL_OPS,
     NRS_START_GET_FROM_REMOTE,
     NRS_END_GET_FROM_REMOTE,
     NRS_END_PARALLEL_OPS
   };
} //percept
#endif

