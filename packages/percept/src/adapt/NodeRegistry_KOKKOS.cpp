// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <adapt/NodeRegistry_KOKKOS.hpp>
#include <adapt/FindValidCentroid.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>
#include <adapt/RefinerUtil.hpp>

#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <stk_mesh/base/DataTraits.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

#include <set>
#include <typeinfo>
#include <Kokkos_Macros.hpp>
  namespace percept {

    void NodeRegistry_KOKKOS::initialize()
    {
      m_cell_2_data_map.clear();
    }

    void NodeRegistry_KOKKOS::
    beginRegistration(int ireg, int nreg)
    {
//      if (CHECK_DB) checkDB("beginRegistration");

      m_nodes_to_ghost.resize(0);
      m_pseudo_entities.clear();
      m_state = NRS_START_REGISTER_NODE;
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry_KOKKOS::beginRegistration" << std::endl;
    }

    void NodeRegistry_KOKKOS::resizeMap()
    {   //attempts to resize the map to fit only it's current key/val pairs. The kokkos rehash doesn't allow for downsizing the map which is problematic for our use case of overestimating the number of nodes needed to avoid resizing the map
            //in the node registration process
        //madbrew : Note that this calls the size() function which is NOT a device function and cannot be called from a parallel kernel. Further, the Kokkos documentation says that this function has undefined behavior when erasable() is true
        //It appears to be able to reliably tell me how many entries are in the map in a unit test, so I'm guessing the undefined behavior arises when the size function is called simultaneously with an insert/removal
        unsigned current_cap = m_cell_2_data_map.capacity();
        unsigned current_size = m_cell_2_data_map.size();
        if(current_size>current_cap)
            throw std::runtime_error("Something went horribly wrong, map size is greater than capacity");
        if( ( (current_cap - current_size) < m_waste_tolerance ) || current_size == 128) //since the map always sizes itself to the closest multiple of 128, 128 is the minimum size. No need to resize
            return; //if the map is already full or an acceptable amount of space is being wasted, then do not resize

        m_cell_2_data_map.rehash(current_size);
    }

    void NodeRegistry_KOKKOS::
    endRegistration(int ireg, int nreg)
    {
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry_KOKKOS::endRegistration start" << std::endl;


      communicate_marks();

      removeUnmarkedSubDimEntities();

      std::cout << "map's size = " << m_cell_2_data_map.size() << std::endl;
      resizeMap();
      std::cout << "map's capacity = " << m_cell_2_data_map.capacity() << std::endl;

      mod_begin();
      this->createNewNodesInParallel();

      {
        m_nodes_to_ghost.resize(0);

        if (m_debug)
          std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry_KOKKOS::endRegistration end" << std::endl;

//        if (CHECK_DB) checkDB("endRegistration - after rehash");
      }
      m_state = NRS_END_REGISTER_NODE;
    }

    /// when a sub-dim entity is visited during node registration but is flagged as not being marked, and thus not requiring
    ///   any new nodes, we flag it with NR_MARK_NONE, then remove it here
    void NodeRegistry_KOKKOS::removeUnmarkedSubDimEntities()
    {
      EXCEPTWATCH;
//      SubDimCellToDataMap_KOKKOS::iterator iter;
      bool debug = false;

//      for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
//        {
//          SubDimCellData& nodeId_elementOwnerId = (*iter).second;
      for(unsigned iMapIndx = 0; iMapIndx<m_cell_2_data_map.capacity();iMapIndx++)
      {
          const SubDimCell_SDCEntityType& subDimEntity = m_cell_2_data_map.key_at(iMapIndx);
          unsigned iData = m_cell_2_data_map.find(subDimEntity);
          if(iData == Kokkos::UnorderedMapInvalidIndex)
              continue;
          SubDimCellData& nodeId_elementOwnerId = m_cell_2_data_map.value_at(iData);
          NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
          if (debug) std::cout << "tmp SRK nodeIds_onSE.size= " << nodeIds_onSE.size() << std::endl;
          if (nodeIds_onSE.size())
            {
              unsigned mark = nodeIds_onSE.m_mark;
              unsigned is_marked = mark & NR_MARK;
              unsigned is_not_marked = mark & NR_MARK_NONE;

              if (debug)
                std::cout << "tmp SRK is_marked= " << is_marked << " is_not_marked= " << is_not_marked << std::endl;
              if (!is_marked && is_not_marked)
                {
                  // check if the node is a "hanging node" in which case it has relations
                  bool found = false;
                  for (unsigned in=0; in < nodeIds_onSE.size(); ++in)
                    {
                      if(!m_eMesh.is_valid(nodeIds_onSE[in])) continue;
                      if (debug) std::cout << "is valid = " << m_eMesh.is_valid(nodeIds_onSE[in]) << std::endl;
                      size_t num_rels = m_eMesh.get_bulk_data()->count_relations(nodeIds_onSE[in]);
                      if (num_rels)
                        {
                          if (debug) std::cout << "tmp SRK num_rels is non-zero in removeUnmarkedSubDimEntities, id= " << m_eMesh.identifier(nodeIds_onSE[in]) <<  std::endl;
                          found = true;
                          break;
                        }
                    }

                  if (debug)
                    std::cout << "P[" << m_eMesh.get_rank() << " removeUnmarkedSubDimEntities:: tmp SRK for node= " << m_eMesh.identifier(nodeIds_onSE[0])  << " FOUND mark = " << mark
                              << " NR_MARK =" << NR_MARK << " NR_MARK_NONE= " << NR_MARK_NONE
                              << " resize to 0 (delete unmarked entity) = " << (!found)
                              << std::endl;
                  if (!found)
                    {
                      nodeIds_onSE.resize(0);
                    }
                }
            }
        }
    }

    /// Register the need for a new node on the sub-dimensional entity @param subDimEntity on element @param element.
    /// If the element is a ghost element, the entity is still registered: the locality/ownership of the new entity
    /// can be determined by the locality of the element (ghost or not).
    bool NodeRegistry_KOKKOS::registerNeedNewNode(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes,const CellTopologyData * const bucket_topo_data)
    {
      bool ret_val = false;
      SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
      bool foundGhostNode = getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd, bucket_topo_data);
      if (foundGhostNode)
        return ret_val;

      if (s_use_new_ownership_check && m_eMesh.isGhostElement(element))
        return false;

      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnerId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnerId = (nodeId_elementOwnerId_ptr ? *nodeId_elementOwnerId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnerId_ptr == 0;
      bool is_not_empty_but_data_cleared = (!is_empty && std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId).size() == 0);

      // if empty or if my id is the smallest, make this element the owner
      stk::mesh::EntityId db_id = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId).id();
      stk::mesh::EntityRank db_rank = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId).rank();

      stk::mesh::EntityId element_id = m_eMesh.identifier(element);
      stk::mesh::EntityId element_rank = m_eMesh.entity_rank(element);
      bool should_put_in_id = (element_id < db_id);
      bool should_put_in_rank_gt = (element_rank > db_rank);
      bool should_put_in_rank_gte = (element_rank >= db_rank);
      bool should_put_in = should_put_in_rank_gt || (should_put_in_id && should_put_in_rank_gte);

      unsigned smark=0;
      if (!is_empty)
        {
          unsigned& mark = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId).m_mark;
          if (needNodes)
            mark |= NR_MARK;
          else
            mark |= NR_MARK_NONE;
          smark = mark;
        }

      /// once it's in, the assertion should be:
      ///   owning_elementId < non_owning_elementId && owning_elementRank >= non_owning_elementRank
      ///
      if (is_empty || is_not_empty_but_data_cleared || should_put_in)
        {
          // create SubDimCellData for the map rhs
          // add one to iSubDimOrd for error checks later
          VERIFY_OP_ON(m_eMesh.identifier(element), !=, 0, "hmmm registerNeedNewNode #1");
          VERIFY_OP_ON(m_eMesh.is_valid(element), ==, true, "hmmm registerNeedNewNode #2");
          if (is_empty || is_not_empty_but_data_cleared)
            {
              unsigned numNewNodes = needed_entity_rank.second;

              SubDimCellData data = std::forward_as_tuple( NodeIdsOnSubDimEntityType(numNewNodes, stk::mesh::Entity(), smark), stk::mesh::EntityKey(m_eMesh.entity_rank(element), m_eMesh.identifier(element)), (unsigned char)needed_entity_rank.first, (unsigned char)(iSubDimOrd+1), Double2() );
              NodeIdsOnSubDimEntityType& nid_new = std::get<SDC_DATA_GLOBAL_NODE_IDS>(data);
              if (needNodes)
                nid_new.m_mark |= NR_MARK;
              else
                nid_new.m_mark |= NR_MARK_NONE;

              smark = nid_new.m_mark;
              putInMap(subDimEntity,  data);
            }
          else
            {
              stk::mesh::Entity owning_element = m_eMesh.get_entity(db_rank, db_id);
              VERIFY_OP_ON(m_eMesh.is_valid(owning_element), ==, true, "hmmm");
              NodeIdsOnSubDimEntityType& nid = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
              SubDimCellData data = std::forward_as_tuple(nid, stk::mesh::EntityKey(m_eMesh.entity_rank(element), m_eMesh.identifier(element)), (unsigned char)needed_entity_rank.first, (unsigned char)(iSubDimOrd+1), Double2() );
              putInMap(subDimEntity,data);
            }
          ret_val = true;
        }

      // all debug code from here on

      bool debug = false;
      if (debug && !is_empty)
        {
      unsigned originalMark = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId).m_mark;
          std::cout << m_eMesh.rank() << " registerNeedNewNode element= " << m_eMesh.print_entity_compact(element) << "\n needNodes= " << needNodes << " iSubDimOrd= " << iSubDimOrd
                    << " is_empty= " << is_empty << " is_not_empty_but_data_cleared= " << is_not_empty_but_data_cleared << " should_put_in= " << should_put_in
                    << " originalMark= " << originalMark << " currentMark= " << smark << " NR_MARK= " << NR_MARK
                    << m_eMesh.demangled_stacktrace(20)
                    << std::endl;
        }
      if (debug)
        {
          SubDimCellData* a_nodeId_elementOwnerId_ptr = getFromMapPtr(subDimEntity);
          SubDimCellData& a_nodeId_elementOwnerId = (a_nodeId_elementOwnerId_ptr ? *a_nodeId_elementOwnerId_ptr : empty_SubDimCellData);
          bool a_is_empty = a_nodeId_elementOwnerId_ptr == 0;
          unsigned& gotMark = std::get<SDC_DATA_GLOBAL_NODE_IDS>(a_nodeId_elementOwnerId).m_mark;
          unsigned nidSize = std::get<SDC_DATA_GLOBAL_NODE_IDS>(a_nodeId_elementOwnerId).size();

          std::ostringstream sout;
          sout << "P[" << m_eMesh.get_rank() << "] registerNeedNewNode:: element= " << m_eMesh.identifier(element) << " nidSize= " << nidSize
               << " nid= " << (nidSize ? (int)m_eMesh.identifier(std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId)[0]) : -1);

          sout << " smark= " << smark << " gotMark= " << gotMark << " needNodes= " << needNodes << " isG= " << m_eMesh.isGhostElement(element)
               << " is_empty= " << is_empty
               << " a_is_empty= " << a_is_empty
               << " should_put_in= " << should_put_in
               << " should_put_in_id= " << should_put_in_id
               << " should_put_in_rank_gt= " << should_put_in_rank_gt
               << " should_put_in_rank_gte= " << should_put_in_rank_gte
               << " needed_entity_rank= "
               << needed_entity_rank.first << " subDimEntity= ";

          for (unsigned k=0; k < subDimEntity.size(); k++)
            {
              sout << " " << m_eMesh.identifier(subDimEntity[k]) << " ";
            }
          std::cout << sout.str() << std::endl;
        }

      //deletethis
/*      SubDimCellData * tester = getFromMapPtr(subDimEntity);
      if(!tester)
          std::cout << "    Cannot access the entity just registered in the map\n\n";
      else
          std::cout << "    Can access the entity just registered in the map\n\n";*/
      //deletethis
      return ret_val;
    }


    /**
     *  For any given subDimEntity (edge, tri or quad face), we can define an "owner" for it as the processor
     *    that owns the node with the minimum ID.  This is implemented below through use of the ordered std::map.
     *  However, a proc can "own" a subDimEntity, but not have any elements use that subDimEntity.
     *
     *  Thus, we use the more global rule that the for all subDimEntity's that are used by an element, the
     *    proc with minimum rank is the one that actually "owns" the subDimEntity - this is implemented in the
     *    routines that unpack the node id's that get sent to all sharing procs.
     *
     *  So, the main function of this routine is to tell the caller all the sharing procs of the subDimEntity
     *    and thus all procs that have an element that use this subDimEntity will send their info to all other
     *    sharing procs, but only the one with the min proc rank will not unpack and override its node ID.
     *    Other procs will unpack and override if the proc rank is smaller than their rank.
     *
     *  Note: for consistency, and for the unpack to work correctly and assign ownership of the data, we
     *    return the current proc in the list of sharing procs @param other_procs
     *
     */

    int NodeRegistry_KOKKOS::proc_owning_subdim_entity(const SubDimCell_SDCEntityType& subDimEntity, std::vector<int>& other_procs, bool& all_shared)
    {
      other_procs.resize(0);
      all_shared = false;

      VERIFY_OP_ON(subDimEntity.size(), >=, 1, "bad size");

      int my_proc = m_eMesh.get_rank();

      bool debug = false;

      std::vector<int> sharing_procs;
      std::map<int, unsigned> count_per_proc;
      all_shared = true;
      for (unsigned i = 0; i < subDimEntity.size(); ++i)
        {
          if (m_eMesh.shared(subDimEntity[i]))
            {
              m_eMesh.get_bulk_data()->comm_shared_procs(m_eMesh.entity_key(subDimEntity[i]), sharing_procs);

              sharing_procs.push_back(my_proc);

              for (unsigned j=0; j < sharing_procs.size(); ++j)
                {
                  ++count_per_proc[sharing_procs[j]];
                }
            }
          else
            {
              all_shared = false;
            }
        }
      if (debug)
        {
          std::cerr << m_eMesh.rank() << "tmp srk all_shared= " << all_shared
              << "\n n0= " << m_eMesh.print_entity_compact(subDimEntity[0])
              << "\n n1= " << m_eMesh.print_entity_compact(subDimEntity[1])
              << std::endl;
        }
      if (!all_shared)
        return -1;

      VERIFY_OP_ON(count_per_proc.size(), >=, 1, "bad size");
      int proc_owner = -1;
      for (auto& counts : count_per_proc)
        {
          if (counts.second == subDimEntity.size())
            {
              if (proc_owner < 0)
                proc_owner = counts.first;
              other_procs.push_back(counts.first);
            }
        }
      if (debug) {
        std::cerr << m_eMesh.rank() << "tmp srk proc_owner= " << proc_owner << std::endl;
      }
      VERIFY_OP_ON(proc_owner, >=, 0, "bad proc_owner");
      return proc_owner;
    }


    /// check the newly registered node from the registry, which does one of three things, depending on what mode we are in:
    ///   1. counts buffer in prep for sending (just does a pack)
    ///   2. packs the buffer (after buffers are alloc'd)
    ///   3. returns the new node after all communications are done


    void NodeRegistry_KOKKOS::
    doForAllSubEntities(ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks,const CellTopologyData * const bucket_topo_data)
    {
      const CellTopologyData * const cell_topo_data = (bucket_topo_data ? bucket_topo_data : m_eMesh.get_cell_topology(element));

      shards::CellTopology cell_topo(cell_topo_data);

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_rank == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_rank == stk::topology::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          curr_funct = function;
          doForAllSubEntitiesFunctor for_each(this,element,needed_entity_ranks,bucket_topo_data,numSubDimNeededEntities,ineed_ent);
          for_each.run();

/*
          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              (this ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true, bucket_topo_data);

            } // iSubDimOrd
*/
        } // ineed_ent
    }

    unsigned NodeRegistry_KOKKOS::local_size()
    {
      unsigned sz=0;
//      for (SubDimCellToDataMap_KOKKOS::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
//        {
//          SubDimCellData& data = (*cell_iter).second;
      for(unsigned iMapIndx = 0; iMapIndx<m_cell_2_data_map.capacity();iMapIndx++)
      {
          const SubDimCell_SDCEntityType& subDimEntity = m_cell_2_data_map.key_at(iMapIndx);
          unsigned iData = m_cell_2_data_map.find(subDimEntity);
          if(iData == Kokkos::UnorderedMapInvalidIndex)
              continue;
          SubDimCellData& nodeId_elementOwnerId = m_cell_2_data_map.value_at(iData);
          stk::mesh::EntityId owning_elementId = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId).id();

          NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
          if (nodeIds_onSE.size())
            {
              stk::mesh::EntityRank erank = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId).rank();
              stk::mesh::Entity owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);

              if (!m_eMesh.is_valid(owning_element))
                {
                  if (!s_use_new_ownership_check)
                    {
                      std::cout << "tmp owning_element = null, owning_elementId= " << owning_elementId
                                  << std::endl;
                      throw std::logic_error("logic: hmmm #5.2");
                    }
                  else
                    continue;
                }
              if (!m_eMesh.isGhostElement(owning_element))
                {
                  sz += nodeIds_onSE.size();
                }
            }
        }
      return sz;
    }

    /// after registering all needed nodes, this method is used to request new nodes on this processor
    void NodeRegistry_KOKKOS::createNewNodesInParallel()
    {
      stk::mesh::Part* new_nodes_part = m_eMesh.get_non_const_part("refine_new_nodes_part");
      VERIFY_OP_ON(new_nodes_part, !=, 0, "new_nodes_part");
      unsigned num_nodes_needed = local_size();

      // assert( bulk data is in modifiable mode)
      // create new entities on this proc
      vector<stk::mesh::Entity> new_nodes;

      if (m_useAddNodeSharing)
        {
          stk::diag::Timer *timerCE_ = 0;
          if (m_refiner)
            {
              timerCE_ = new stk::diag::Timer("NR_CreateEnt", m_refiner->rootTimer());
              timerCE_->start();
            }

#if USE_CREATE_ENTITIES
          m_eMesh.createEntities( stk::topology::NODE_RANK, num_nodes_needed, new_nodes);
#else
          m_eMesh.initializeIdServer();
          stk::mesh::Part& nodePart = m_eMesh.get_fem_meta_data()->get_topology_root_part(stk::topology::NODE);
          stk::mesh::PartVector nodeParts(1, &nodePart);
          m_eMesh.getEntitiesUsingIdServer( m_eMesh.node_rank(), num_nodes_needed, new_nodes, nodeParts);
#endif
          if (timerCE_)
            {
              timerCE_->stop();
              delete timerCE_;
            }
        }
      std::vector<stk::mesh::EntityId> ids(num_nodes_needed);

      if (new_nodes_part)
        {
          stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->locally_owned_part() );
          std::vector<stk::mesh::Part*> add_parts(1, new_nodes_part);
          std::vector<stk::mesh::Part*> remove_parts;
          for (unsigned ind = 0; ind < new_nodes.size(); ind++)
            {
              if (m_eMesh.m_new_nodes_field)
                {
                  NewNodesType_type *ndata = stk::mesh::field_data(*m_eMesh.m_new_nodes_field, new_nodes[ind]);
                  if (ndata)
                    {
                      ndata[0] = static_cast<NewNodesType_type>(1);
                    }
                }
              m_eMesh.get_bulk_data()->change_entity_parts( new_nodes[ind], add_parts, remove_parts );
            }
        }
      // set map values to new node id's
//      unsigned inode=0;

      Kokkos::View<unsigned*,DataLayout,MemSpace>threadSafeiNode("threadSafeiNode",1);
      Kokkos::View<unsigned*,DataLayout,MemSpace>::HostMirror threadSafeiNode_mir("threadSafeiNode_mir",1);
      threadSafeiNode_mir(0) = 0;
      Kokkos::deep_copy(threadSafeiNode,threadSafeiNode_mir);

      // Malachi Phillips Jun 21, 2018
      //
      // This needs a capture by *this in order to work on __device__
      //
      // Capture by *this cannot work on __host__, unless using a complaint c++17 compiler
      //
      Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,m_cell_2_data_map.capacity()), [&] (const unsigned int iMapIndx)
      {
          const SubDimCell_SDCEntityType& subDimEntity = m_cell_2_data_map.key_at(iMapIndx);
          unsigned iData = m_cell_2_data_map.find(subDimEntity);
          if(iData == Kokkos::UnorderedMapInvalidIndex)
              //continue; //continue doesn't work the same way in pf
              (void)iData;
          else{
              SubDimCellData& data = m_cell_2_data_map.value_at(iData);
              NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(data);
              if (!nodeIds_onSE.size())
                  //            continue; //continue doesn't work the same way in pf
                  (void)iData;
              else{
                  stk::mesh::EntityId owning_elementId = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(data).id();

#ifndef __CUDACC__
                  if (!owning_elementId)
                  {
                    throw std::logic_error("logic: hmmm #5.4.0");
                  }
#endif

                  stk::mesh::EntityRank erank = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(data).rank();
                  stk::mesh::Entity owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);

                  if (!m_eMesh.is_valid(owning_element))
                  {
#ifndef __CUDACC__
                      if (!s_use_new_ownership_check)
                          throw std::logic_error("logic: hmmm #5.4");
#endif
                  }

                  else if (s_use_new_ownership_check && m_eMesh.isGhostElement(owning_element))
                  {
                      //std::cerr << m_eMesh.rank() << " owning_element= " << m_eMesh.print_entity_compact(owning_element) << std::endl;
                      // FIXME
                      //VERIFY_MSG("found ghost");
                      //              continue; //continue doesn't work the same way in pf
                  }

                  else if (!m_eMesh.isGhostElement(owning_element))
                  {
                      if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
                      {
#ifndef __CUDACC__
                          throw std::logic_error("NodeRegistry_KOKKOS:: createNewNodesInParallel logic err #0.0");
#endif
                      }

                      for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
                      {
                          unsigned localiNode = Kokkos::atomic_fetch_add<unsigned>(&threadSafeiNode(0),1);  //fetches old value then increments the index in a thread safe manner

                          VERIFY_OP(localiNode, < , num_nodes_needed, "UniformRefiner::doBreak() too many nodes");
#ifndef __CUDACC__
                          if ( DEBUG_NR_UNREF)
                          {
                              std::cout << "tmp createNewNodesInParallel: old node id= " << (m_eMesh.is_valid(nodeIds_onSE[ii]) ? toString(m_eMesh.identifier(nodeIds_onSE[ii])) : std::string("null")) << std::endl;
                              std::cout << "tmp createNewNodesInParallel: new node=";
                              m_eMesh.print_entity(std::cout, new_nodes[localiNode]);
                          }
#endif

                          // if already exists from a previous iteration/call to doBreak, don't reset it and just use the old node
                          if (m_eMesh.is_valid(nodeIds_onSE[ii]))
                          {
#ifndef __CUDACC__
                              if (DEBUG_NR_UNREF)
                              {
                                  std::cout << "tmp createNewNodesInParallel: old node id is no-null, re-using it= " << (m_eMesh.is_valid(nodeIds_onSE[ii]) ? toString(m_eMesh.identifier(nodeIds_onSE[ii])) : std::string("null")) << std::endl;
                                  std::cout << "tmp createNewNodesInParallel: new node=";
                                  m_eMesh.print_entity(std::cout, new_nodes[localiNode]);
                              }
#endif
                          }
                          else
                          {
                              nodeIds_onSE[ii] = new_nodes[localiNode];
                              nodeIds_onSE.m_entity_id_vector[ii] = m_eMesh.identifier(new_nodes[localiNode]);
                          }

                          //                  inode++;
                      }
                  }
              }
          }
        });//iMapIndx
    }

    void NodeRegistry_KOKKOS::mod_begin()
    {
      if (m_refiner)
        {
          m_refiner->mod_begin();
        }
      else
        {
          m_eMesh.get_bulk_data()->modification_begin();
        }
    }

    void NodeRegistry_KOKKOS::mod_end(const std::string& msg)
    {
      if (m_refiner)
        {
          m_refiner->mod_end(0,"NReg:"+msg);
        }
      else
        {
          m_eMesh.get_bulk_data()->modification_end();
        }
    }

    void NodeRegistry_KOKKOS::communicate_marks()
    {

      // only one stage now - each proc sends to other sharing procs and each can |= and accumulate locally

      for (int stage = 0; stage < 1; ++stage)
      {
        stk::CommSparse commAll (m_eMesh.parallel());

        communicate_marks_pack(commAll, stage);

        commAll.allocate_buffers();

        communicate_marks_pack(commAll, stage);
        commAll.communicate();
        communicate_marks_unpack(commAll);
      }
    }

    void NodeRegistry_KOKKOS::communicate_marks_pack(stk::CommSparse& commAll, int stage)
    {
      CommDataType buffer_entry;

      SubDimCellToDataMap_KOKKOS& map = m_cell_2_data_map;

//      for (SubDimCellToDataMap_KOKKOS::iterator iter = map.begin(); iter != map.end(); ++iter)
//        {
//          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
//          SubDimCellData& nodeId_elementOwnerId = (*iter).second;
      for(unsigned iMapIndx = 0; iMapIndx<map.capacity();iMapIndx++)
      {
          const SubDimCell_SDCEntityType& subDimEntity = map.key_at(iMapIndx);
          unsigned iData = map.find(subDimEntity);
          if(iData == Kokkos::UnorderedMapInvalidIndex)
              continue;
          SubDimCellData& nodeId_elementOwnerId = map.value_at(iData);
          stk::mesh::EntityRank owning_element_subDim_rank = static_cast<stk::mesh::EntityRank>(std::get<SDC_DATA_OWNING_SUBDIM_RANK>(nodeId_elementOwnerId));

          static std::vector<int> procs_to_send_to;

          bool all_shared = false;
          int new_owning_proc = proc_owning_subdim_entity(subDimEntity, procs_to_send_to, all_shared);
          (void)new_owning_proc;

          bool need_to_send = false;
          if (stage == 0)
            {
              need_to_send = all_shared && (procs_to_send_to.size() != 0);
            }

          if (!need_to_send)
            continue;

          NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);

          if (nodeIds_onSE.size())
            {
              unsigned mark = nodeIds_onSE.m_mark;
              unsigned is_marked = mark & NR_MARK;
              unsigned is_not_marked = mark & NR_MARK_NONE;
              // only need to send if it's got non-zero mark info
              if (is_marked || is_not_marked)
                {
                  std::get<CDT_SUBDIM_ENTITY_SIZE>(buffer_entry) = subDimEntity.size();
                  std::get<CDT_SUBDIM_ENTITY_RANK>(buffer_entry) = owning_element_subDim_rank;

                  for (unsigned inode=0; inode < subDimEntity.size(); ++inode)
                    std::get<CDT_SUBDIM_ENTITY>(buffer_entry)[inode] = m_eMesh.entity_key(subDimEntity[inode]);

                  for (unsigned jprocs = 0; jprocs < procs_to_send_to.size(); ++jprocs)
                    {
                      if (procs_to_send_to[jprocs] != m_eMesh.get_rank())
                        {
                          commAll.send_buffer( procs_to_send_to[jprocs] ).pack< CommDataType > (buffer_entry);
                          commAll.send_buffer( procs_to_send_to[jprocs] ).pack<unsigned>(mark);
                        }
                    }
                }
            }
        }
    }

    void NodeRegistry_KOKKOS::communicate_marks_unpack(stk::CommSparse& commAll)
    {
      unsigned proc_size = m_eMesh.get_parallel_size();

      CommDataType buffer_entry;

      //try
      {
        for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
          {
            stk::CommBuffer & recv_buffer = commAll.recv_buffer( from_proc );

            while ( recv_buffer.remaining() )
              {
                recv_buffer.unpack< CommDataType >( buffer_entry );
                unsigned mark=0;
                recv_buffer.unpack< unsigned > (mark);

                {
                  unsigned              subDimEntitySize   = std::get<CDT_SUBDIM_ENTITY_SIZE>(buffer_entry);

                  SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
                  //getSubDimEntity(subDimEntity, owning_element, needed_entity_rank, iSubDimOrd);
                  for (unsigned inode = 0; inode < subDimEntitySize; ++inode)
                    {
                      stk::mesh::Entity node = m_eMesh.get_entity(std::get<CDT_SUBDIM_ENTITY>(buffer_entry)[inode]);
                      VERIFY_OP_ON(m_eMesh.is_valid(node), ==, true, "bad node");
                      subDimEntity.insert(node);
                    }

                  static SubDimCellData empty_SubDimCellData;

                  SubDimCellData* nodeId_elementOwnerId_ptr = getFromMapPtr(subDimEntity);
                  SubDimCellData& nodeId_elementOwnerId = (nodeId_elementOwnerId_ptr ? *nodeId_elementOwnerId_ptr : empty_SubDimCellData);
                  bool is_empty = nodeId_elementOwnerId_ptr == 0;

                  //VERIFY_OP_ON(is_empty, !=, true, "hmmm");
                  if (!is_empty)
                    {
                      NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);

                      // accumulation from all other procs
                      nodeIds_onSE.m_mark |= mark;
                    }
                }
              }
          }
      }
    }

    bool NodeRegistry_KOKKOS::verifyAllKeysInBoostNR(NodeRegistry * nr, SetOfEntities& nodesMappedTo, unsigned& noKeysNotInCommon)
    {   //checks to if all the nodes within the kokkos NR have the same mapping and nodes within the boost NR
        //madbrew: not sure if this is really a fair comparison since the node id's in NodeIdsOnSubDimEntityType can differ due to the fact that map loops differ between the registry's when adding id's
            //this causes the registry's to have the same number of id's and even the same sets of id's but a different mapping between parent (i.e. edge/face) and children node(s)
        nodesMappedTo.clear();
        noKeysNotInCommon = 0;
        if(!nr)
        {
            std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  Invalid NodeRegistry_KOKKOS pointer, aborting NR diff\n";
            return false;
        }

        bool bUOmap_contains_kUOmap = true;

        SubDimCellToDataMap_KOKKOS * map = &m_cell_2_data_map;
        for (unsigned iMapIndx = 0; iMapIndx < map->capacity(); iMapIndx++)
        {
            const SubDimCell_SDCEntityType& subDimEntity_KOKKOS = map->key_at(iMapIndx);
            unsigned iVal = map->find(subDimEntity_KOKKOS);
            if(iVal==Kokkos::UnorderedMapInvalidIndex)
                continue;
            SubDimCellData& data_KOKKOS = map->value_at(iVal);

            SubDimCellData * data_BOOST = nr->getFromMapPtr(subDimEntity_KOKKOS);

            if(!data_BOOST){
                bUOmap_contains_kUOmap = false;
                if(bUOmap_contains_kUOmap)
                    std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  key with m_HashCode = " << subDimEntity_KOKKOS.getHash() << " NOT FOUND kokkos map\n";
                noKeysNotInCommon++;
            }
            else
            {
              if ( !( std::get<0>(data_KOKKOS) == std::get<0>(*data_BOOST) ))
                {
//                    bUOmap_contains_kUOmap = false;
                    std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  NodeIdsOnSubDimEntityType on tuple's do not match between NodeRegistry and KokkosNodeRegistry\n";
                    std::cout << "  kOUmap NodeIdsOnSubDimEntityType  = "<< std::get<0>(data_KOKKOS) << " bOUmap NodeIdsOnSubDimEntityType  = " << std::get<0>(*data_BOOST) << "\n";
                }
              if(!( std::get<1>(data_KOKKOS) == std::get<1>(*data_BOOST) ))
                {
//                    bUOmap_contains_kUOmap = false;
                    std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  EntityKey on tuple's do not match between NodeRegistry and KokkosNodeRegistry\n";
                    std::cout << "  kOUmap EntityKey = "<< std::get<1>(data_KOKKOS) << " bOUmap EntityKey = " << std::get<1>(*data_BOOST) << "\n";
                }
                if(!( std::get<2>(data_KOKKOS) == std::get<2>(*data_BOOST) ) )
                {
//                    bUOmap_contains_kUOmap = false;
                    std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  unsigned char on tuple's do not match between NodeRegistry and KokkosNodeRegistry\n";
                    std::cout << "  kOUmap unsigned char = "<< std::get<2>(data_KOKKOS) << " bOUmap unsigned char = " << std::get<2>(*data_BOOST) << "\n";
                }
                if(!( std::get<3>(data_KOKKOS) == std::get<3>(*data_BOOST) ) )
                {
//                    bUOmap_contains_kUOmap = false;
                    std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  unsigned char on tuple's do not match between NodeRegistry and KokkosNodeRegistry\n";
                    std::cout << "  kOUmap unsigned char = "<< std::get<3>(data_KOKKOS) << " bOUmap unsigned char = " << std::get<3>(*data_BOOST) << "\n";
                }
                if(!( std::get<4>(data_KOKKOS) == std::get<4>(*data_BOOST) ) )
                {
//                    bUOmap_contains_kUOmap = false;
                    std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  Double2 on tuple's do not match between NodeRegistry and KokkosNodeRegistry\n";
                    std::cout << "  kOUmap Double2 = "<< std::get<4>(data_KOKKOS) << " bOUmap Double2 = " << std::get<4>(*data_BOOST) << "\n";
                }
            }

            if(!bUOmap_contains_kUOmap)
                std::cout << "NodeRegistry_KOKKOS::compare_to_boost_NR  :  key with m_HashCode = " << subDimEntity_KOKKOS.getHash() << " doesn't match anything in kokkos map\n\n";
            for(unsigned iEnt=0;iEnt<std::get<0>(data_KOKKOS).m_entity_id_vector.size();iEnt++){
                unsigned ent_val = std::get<0>(data_KOKKOS).m_entity_id_vector[iEnt];
                nodesMappedTo.insert(stk::mesh::Entity(ent_val));
            }
        }
        return noKeysNotInCommon==0;
    }

    bool NodeRegistry_KOKKOS::
    getSubDimEntity(SubDimCell_SDCEntityType& subDimEntity, const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd,
                    const CellTopologyData * const bucket_topo_data
                    )
    {
      subDimEntity.clear();
      // in the case of elements, we don't share any nodes so we just make a map of element id to node
      if (needed_entity_rank == stk::topology::ELEMENT_RANK)
        {
          subDimEntity.resize(1);
          subDimEntity[0] =  element ;
          subDimEntity.updateHashCode();
          return false;
        }

      const CellTopologyData * const cell_topo_data = (bucket_topo_data ? bucket_topo_data : m_eMesh.get_cell_topology(element) );

      stk::mesh::Entity const * const elem_nodes = m_eMesh.get_bulk_data()->begin_nodes(element);

      const unsigned *  inodes = 0;
      unsigned nSubDimNodes = 0;
      static const unsigned edge_nodes_2[2] = {0,1};
      static const unsigned face_nodes_3[3] = {0,1,2};
      static const unsigned face_nodes_4[4] = {0,1,2,3};

      // special case for faces in 3D
      if (needed_entity_rank == m_eMesh.face_rank() && needed_entity_rank == m_eMesh.entity_rank(element))
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          if (nSubDimNodes ==3 )
            inodes = face_nodes_3;
          else
            inodes = face_nodes_4;

        }
      // special case for edges in 2D
      else if (needed_entity_rank == m_eMesh.edge_rank() && needed_entity_rank == m_eMesh.entity_rank(element))
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          if (nSubDimNodes == 2 )
            {
              inodes = edge_nodes_2;
            }
          else
            {
              throw std::runtime_error("NodeRegistry_KOKKOS bad for edges");
            }
        }
      else if (needed_entity_rank == m_eMesh.edge_rank())
        {
          inodes = cell_topo_data->edge[iSubDimOrd].node;
          nSubDimNodes = 2;
        }
      else if (needed_entity_rank == m_eMesh.face_rank())
        {
          nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = cell_topo_data->side[iSubDimOrd].node;
        }

      subDimEntity.resize(nSubDimNodes);
      for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
        {
          subDimEntity[jnode] =  elem_nodes[inodes[jnode]] ;
        }
      bool foundGhostNode = false;
      if (m_checkForGhostedNodes)
        {
          for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
            {
              stk::mesh::Bucket& bucket = m_eMesh.bucket(subDimEntity[jnode]);
              if (!bucket.owned() && !bucket.shared())
                foundGhostNode = true;
            }
        }

      subDimEntity.sort();
      subDimEntity.updateHashCode();
      return foundGhostNode;
    }

  }//percept
