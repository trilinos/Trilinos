
/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <stdint.h>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/RunEnvironment.hpp>
#include <stk_percept/ProgressMeter.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/util/MallocUsed.h>

#include <stk_mesh/base/MemoryUsage.hpp>

#include <stk_adapt/RefinerUtil.hpp>

#include <stk_adapt/UniformRefiner.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>



namespace stk { 

  namespace adapt {


#if STK_ADAPT_HAVE_YAML_CPP

    /**  in the following, we are looking at one partition, iM, of M-partitioned mesh
     *
     */

    class SerializeNodeRegistry {
      PerceptMesh& m_eMesh;
      NodeRegistry* m_nodeRegistry;
      std::string m_input_mesh_name;
      std::string m_output_mesh_name;
      std::string m_filePrefix;
      int M, iM, P, iP;
      const std::vector<std::string> & m_entity_rank_names;
      std::vector<stk::mesh::EntityId> m_id_max;
      std::string m_globalIdFile;
      std::string m_nodeRegistryFile;
      std::string m_globalNodeRegistryFile;
      const unsigned FAMILY_TREE_RANK;
      static const bool m_debug = false;
      int m_spatialDim;

    public:
      enum { MaxPass = 3 };

      SerializeNodeRegistry(PerceptMesh& eMesh, NodeRegistry* nodeRegistry, std::string input_mesh_name, std::string output_mesh_name, int M, int iM, int P=1, int iP=0) : 
        m_eMesh(eMesh), m_nodeRegistry(nodeRegistry), m_input_mesh_name(input_mesh_name), m_output_mesh_name(output_mesh_name), m_filePrefix(input_mesh_name), 
        M(M), iM(iM), P(P), iP(iP),
        m_entity_rank_names(eMesh.getFEM_meta_data()->entity_rank_names()),
        m_id_max(m_entity_rank_names.size(), 0u),  FAMILY_TREE_RANK(eMesh.element_rank() + 1u)
      {
        size_t pos = m_filePrefix.find(".");
        if (pos != std::string::npos)
          m_filePrefix = m_filePrefix.substr(0, pos);
        m_globalIdFile = "streaming-refine-global-data."+m_filePrefix+".yaml";
        m_nodeRegistryFile = std::string("streaming-refine-nodeRegistry."+m_filePrefix+".yaml.")+boost::lexical_cast<std::string>(M)+"."+boost::lexical_cast<std::string>(iM);
        m_globalNodeRegistryFile = std::string("streaming-refine-nodeRegistry."+m_filePrefix+".yaml");
        m_spatialDim = eMesh.getSpatialDim();
      }

      void pass(int streaming_pass)
      {
        std::cout << "\n\n ------ SerializeNodeRegistry ----- pass number " << streaming_pass << "\n\n" << std::endl;
        if ( streaming_pass== 0)
          {
            pass0();
          }
        else if (streaming_pass == 1)
          {
            pass1();
          }
        else if (streaming_pass == 2)
          {
            pass2();
          }
        else if (streaming_pass == 3)
          {
            pass3();
          }
      }

      /**
       *   pass0: open unrefined mesh, refine it, find max id
       *   (iM = 0...M)
       *   1. if iM==0, setCurrentGlobalMaxId to [0,0,0,0]
       *   2. getCurrentGlobalMaxId()
       *   3. find new max id from current mesh
       *   4. setCurrentGlobalMaxId()
       */
      void pass0()
      {
        if (0 == iM)
          {
            // create initial file, write 0's in it  for m_id_max (it should be initialized to 0, but just to be sure, we reset it here)
            for (unsigned irank=0; irank < m_id_max.size(); irank++)
              m_id_max[irank]=0u;
            setCurrentGlobalMaxId();
          }

        getCurrentGlobalMaxId();
        findCurrentMaxId(m_eMesh);
        setCurrentGlobalMaxId();
        printCurrentGlobalMaxId("pass0");
      }

      void printCurrentGlobalMaxId(std::string msg="")
      {
        if (m_debug) std::cout << "SerializeNodeRegistry["<<M<<", "<<iM<<"]::printCurrentGlobalMaxId: = " << m_id_max << " for: " << msg << std::endl;
      }

      /**
       *   pass1: refine mesh, write local NodeRegistry, set new max id from refined mesh
       *   (iM = 0...M)
       *   1. open file.e.M.iM, refine mesh
       *   2. if iM==0, setCurrentGlobalMaxId to [0,0,0,0]
       *   3. getCurrentGlobalMaxId()
       *   4. resetNewElementIds() (resets new element ids by looking at child elements only,
       *   5. write NodeRegistry in name.yaml.M.iM
       *   6. setCurrentGlobalMaxId()
       *   7. save refined mesh
       */
      void pass1()
      {
        NodeRegistry& nodeRegistry = *m_nodeRegistry;
        writeNodeRegistry(nodeRegistry, m_nodeRegistryFile);

        /*
        if (0 == iM)
          {
            // create initial file, write 0's in it  for m_id_max (it should be initialized to 0, but just to be sure, we reset it here)
            for (unsigned irank=0; irank < m_id_max.size(); irank++)
              m_id_max[irank]=0u;
            setCurrentGlobalMaxId();
          }
        */

        getCurrentGlobalMaxId();
        resetNewElementIds(m_eMesh);
        setCurrentGlobalMaxId();
        printCurrentGlobalMaxId("pass1");
      }

      /**
       *   pass2 - create global NodeRegistry from each local one by "last one wins"
       *   (single call, no loop over iM)
       *   1. create new (global) NodeRegistry
       *   2. getCurrentGlobalMaxId
       *   3. loop iM
       *      a. read NodeRegistry from name.yaml.M.iM -> input values into new NodeRegistry
       *   4. for each key/value pair, increment idserver, save new id in value
       *   5. write new global NodeRegistry
       *
       */
      void pass2()
      {
        if (iM != 0) throw std::logic_error("SerializeNodeRegistry::pass2 logic error");

        getCurrentGlobalMaxId();
        //setCurrentGlobalMaxId();
        printCurrentGlobalMaxId("pass2 initial");

        PerceptMesh eMeshNew(m_spatialDim);
        eMeshNew.openEmpty();
        NodeRegistry globalNR(eMeshNew);

        PerceptMesh eMeshLocal(m_spatialDim);
        eMeshLocal.openEmpty();

        for (iM = 0; iM < M; iM++)
          {
            NodeRegistry newLocalNR(eMeshLocal);
            m_nodeRegistryFile = std::string("streaming-refine-nodeRegistry."+m_filePrefix+".yaml.")+boost::lexical_cast<std::string>(M)+"."+boost::lexical_cast<std::string>(iM);
            readNodeRegistry(newLocalNR, m_nodeRegistryFile);
            processNodeRegistry(newLocalNR, globalNR);
          }
        resetIds(globalNR);  // on the in-memory global NodeRegistry
        writeNodeRegistry(globalNR, m_globalNodeRegistryFile);
        setCurrentGlobalMaxId();
        printCurrentGlobalMaxId("pass2 done");
      }

      /**
       *   pass 3: 
       *   1. read global NodeRegistry
       *   2. read each refined file-ref.M.iM
       *   3. lookup edge/face/elem in global NodeRegistry, reset id to that found in NR
       *   4. write refined file-ref-id.M.iM
       */
      void pass3()
      {
        PerceptMesh eMeshGlobal(m_spatialDim);
        eMeshGlobal.openEmpty();
        // here we use the same eMesh to use the same BulkData - this could be changed if we didn't use pointers
        // but always used id's in the maps' keys
        //PerceptMesh eMeshLocal(m_spatialDim);
        //eMeshLocal.openEmpty();
        NodeRegistry globalNR(eMeshGlobal);
        //NodeRegistry localNR(eMeshLocal);
        NodeRegistry localNR(eMeshGlobal);
        readNodeRegistry(localNR, m_nodeRegistryFile);
        readNodeRegistry(globalNR, m_globalNodeRegistryFile);
        lookupAndSetNewNodeIds(localNR, globalNR);
      }

      void lookupAndSetNewNodeIds(NodeRegistry& localNR, NodeRegistry& globalNR)
      {
        m_eMesh.getBulkData()->modification_begin();
        SubDimCellToDataMap& map = localNR.getMap();
        //SubDimCellToDataMap& globalMap = globalNR.getMap();
        std::cout << " tmp SerializeNodeRegistry::lookupAndSetNewNodeIds map size: " << map.size() << std::endl;

        SubDimCellToDataMap::iterator iter;
        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            const SubDimCell_SDSEntityType& subDimEntity = iter->first;
            SubDimCellData& nodeId_elementOwnderId = iter->second;
            SubDimCellData* global_nodeId_elementOwnderId_ptr = globalNR.getFromMapPtr(subDimEntity);
            if (!global_nodeId_elementOwnderId_ptr)
              {
                std::cout << "SerializeNodeRegistry::lookupAndSetNewNodeIds couldn't find subDimEntity= " << subDimEntity;
                for (unsigned kk=0; kk < subDimEntity.size(); kk++)
                  {
                    std::cout << " [" << subDimEntity[kk]->identifier() << "] ";
                  }
                std::cout << std::endl;
              }
            else
              {
//                 if (m_debug)
//                   std::cout << "found it" << std::endl;
              }
            VERIFY_OP_ON(global_nodeId_elementOwnderId_ptr, !=, 0, "SerializeNodeRegistry::lookupAndSetNewNodeIds couldn't find subDimEntity");

            NodeIdsOnSubDimEntityType& global_nodeIds_onSE = global_nodeId_elementOwnderId_ptr->get<SDC_DATA_GLOBAL_NODE_IDS>();
            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            unsigned global_nnodes = global_nodeIds_onSE.size();
            unsigned nnodes = nodeIds_onSE.size();
            VERIFY_OP_ON(global_nnodes, ==, nnodes, "SerializeNodeRegistry::lookupAndSetNewNodeIds: mismatch in nnodes");

            for (unsigned inode=0; inode < nnodes; inode++)
              {
                stk::mesh::EntityId id_new = global_nodeIds_onSE.m_entity_id_vector[inode];
                stk::mesh::EntityId id_old = nodeIds_onSE.m_entity_id_vector[inode];
                stk::mesh::Entity* tmp_global_node = globalNR.getMesh().getBulkData()->get_entity(0, id_new);
                stk::mesh::Entity* local_node_to_change = m_eMesh.getBulkData()->get_entity(0, id_old);
                if (m_debug) std::cout << "iM= " << iM << " id_new= " << id_new << " id_old= " << id_old << std::endl;
                VERIFY_OP_ON(local_node_to_change, !=, 0, "SerializeNodeRegistry::lookupAndSetNewNodeIds null local_node_to_change");
                VERIFY_OP_ON(tmp_global_node, !=, 0, "SerializeNodeRegistry::lookupAndSetNewNodeIds null tmp_global_node");
                //VERIFY_OP_ON(local_node_to_change, ==, nodeIds_onSE[inode], "SerializeNodeRegistry::lookupAndSetNewNodeIds old node mistmatch");
                VERIFY_OP_ON(tmp_global_node, ==, global_nodeIds_onSE[inode], "SerializeNodeRegistry::lookupAndSetNewNodeIds new node mistmatch");
                //nodeIds_onSE[inode] = global_nodeIds_onSE[inode];
                if (id_new != id_old)
                  {
                    if (m_debug) std::cout << "DIFF iM= " << iM << " id_new= " << id_new << " id_old= " << id_old << std::endl;
                    //bool did_destroy = m_eMesh.getBulkData()->destroy_entity(new_node);
                    //VERIFY_OP_ON(did_destroy, !=, false, "SerializeNodeRegistry::lookupAndSetNewNodeIds couldn't destroy global node");
                    stk::mesh::Entity* local_node_new_id_check = m_eMesh.getBulkData()->get_entity(0, id_new);
                    //VERIFY_OP_ON(local_node_new_id_check, ==, 0, "SerializeNodeRegistry::lookupAndSetNewNodeIds couldn't change id of local node");
                    if (local_node_new_id_check)
                      {
                        //VERIFY_OP_ON(local_node_new_id_check, ==, 0, "SerializeNodeRegistry::lookupAndSetNewNodeIds couldn't change id of local node");
                        
                      }
                    if (!local_node_new_id_check)
                      {
                        m_eMesh.getBulkData()->change_entity_id(id_new, *local_node_to_change);
                      }
                  }
              }
          }
        m_eMesh.getBulkData()->modification_end();
      }

      void writeNodeRegistry(NodeRegistry& nodeRegistry, std::string filename)
      {
        YAML::Emitter yaml;
        std::cout << "\nnodeRegistry.serialize_write(yaml) to file= " << filename << std::endl;
        nodeRegistry.serialize_write(yaml);
        if (!yaml.good())
          {
            std::cout << "Emitter error: " << yaml.good() << " " <<yaml.GetLastError() << "\n";
            throw std::runtime_error("Emitter error");
          }
        std::ofstream file(filename.c_str());
        if (!file.is_open() || !file.good())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::writeNodeRegistry couldn't open file ")+filename);
          }
        file << yaml.c_str();
        file.close();
      }

      void readNodeRegistry(NodeRegistry& nodeRegistry, std::string filename)
      {
        std::cout << "\nnodeRegistry.serialize_read() from file= " << filename << std::endl;
        std::ifstream file(filename.c_str());
        if (!file.is_open() || !file.good())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::readNodeRegistry couldn't open file ")+filename);
          }
        
        nodeRegistry.serialize_read(file);
      }

      void processNodeRegistry(NodeRegistry& newLocalNR, NodeRegistry& globalNR)
      {
        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = newLocalNR.getMap();
        SubDimCellToDataMap& globalMap = globalNR.getMap();
        std::cout << " tmp SerializeNodeRegistry::processNodeRegistry map size: " << map.size() << std::endl;

        // key.serialized = { nodeid_0,... : set<EntityId> }
        // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            if (m_debug) 
              {
                std::cout << "SerializeNodeRegistry::processNodeRegistry inserting map entry = " << subDimEntity ;
                for (unsigned kk=0; kk < subDimEntity.size(); kk++)
                  {
                    std::cout << " [" << subDimEntity[kk]->identifier() << "] ";
                  }
                std::cout << " data= " << nodeId_elementOwnderId << " nid=" <<  nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().m_entity_id_vector[0] << std::endl;
              }
            /// clone subDimEntity...
            globalMap[subDimEntity] = nodeId_elementOwnderId;
          }        
        std::cout << "SerializeNodeRegistry::processNodeRegistry globalMap size= " << globalMap.size() << std::endl;
      }

      /// Using the current global max id, in a simple id-server manner, generate new id's and assign to 
      ///   the shared nodes.
      void resetIds(NodeRegistry& nodeRegistry) 
      {
        nodeRegistry.getMesh().getBulkData()->modification_begin();

        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = nodeRegistry.getMap();
        std::cout << " tmp SerializeNodeRegistry::resetIds map size: " << map.size() << std::endl;

        // key.serialized = { nodeid_0,... : set<EntityId> }
        // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            //const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;
            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            unsigned nnodes = nodeIds_onSE.size();
            for (unsigned inode=0; inode < nnodes; inode++)
              {
                stk::mesh::EntityId id_new = m_id_max[0]+1;
                m_id_max[0] = id_new;
                nodeIds_onSE.m_entity_id_vector[inode] = id_new;

                stk::mesh::Entity * node = nodeRegistry.getMesh().getBulkData()->get_entity(0, id_new);
                //key.insert(const_cast<stk::mesh::Entity*>(&element) );
                VERIFY_OP_ON(node, ==, 0, "SerializeNodeRegistry::resetIds found existing node with proposed new id");
                if (!node)
                  {
                    stk::mesh::PartVector parts(1, &nodeRegistry.getMesh().getFEM_meta_data()->universal_part());
                    node = &nodeRegistry.getMesh().getBulkData()->declare_entity(0, id_new, parts);
                  }
                nodeIds_onSE[inode] = node;
              }
          }
        nodeRegistry.getMesh().getBulkData()->modification_end();
      }

      ///Note: we read and write to a file instead of storing in memory to allow stk_adapt_exe to be called on a 
      //   single piece of the partition - useful for testing and for potentially parallelizing the global outer loop
      void getCurrentGlobalMaxId()
      {
        // get max node, element id, dump to global file, if global file exists, read from it, use those as start values
        std::fstream file;
        file.open(m_globalIdFile.c_str(), std::ios_base::in);  
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::getCurrentGlobalMaxId couldn't open file ")+m_globalIdFile);
          }
            
        YAML::Parser parser(file);
        YAML::Node doc;

        try {
          while(parser.GetNextDocument(doc)) {
            std::cout << "\n read doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (doc.Type() == YAML::NodeType::Map)
              {
                for (unsigned irank=0; irank < m_id_max.size(); irank++)
                  {
                    doc[m_entity_rank_names[irank]] >> m_id_max[irank];
                  }
              }
            else
              {
                throw std::runtime_error("bad streaming-refine-global-data file");
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << "\n";
          throw std::runtime_error( e.what());
        }
        file.close();
      }

      void setCurrentGlobalMaxId()
      {
        std::fstream file;
        file.open(m_globalIdFile.c_str(), std::ios_base::out | std::ios_base::trunc);
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::setCurrentGlobalMaxId couldn't open file ")+m_globalIdFile);
          }

        YAML::Emitter out;
        out << YAML::BeginMap;
        for (unsigned irank=0; irank < m_id_max.size(); irank++)
          {
            out << YAML::Key << m_entity_rank_names[irank] << YAML::Value << m_id_max[irank];
          }
        out << YAML::EndMap;
        file << out.c_str();
        file.close();
      }

      void findCurrentMaxId(PerceptMesh& eMesh)
      {
        for (unsigned irank=0; irank < m_id_max.size(); irank++)
          {
            // bucket loop
            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.getBulkData()->buckets( irank );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (selector(**k))
                {
                  stk::mesh::Bucket & bucket = **k ;
                  const unsigned num_entities_in_bucket = bucket.size();

                  for (unsigned iEntity = 0; iEntity < num_entities_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity& entity = bucket[iEntity];
                      stk::mesh::EntityId id = entity.identifier();
                      m_id_max[irank] = std::max(m_id_max[irank], id);
                    }
                }
              }
          }
      }

      // use id_max values as a simple id server to reset new element ids
      void resetNewElementIds(PerceptMesh& eMesh)
      {
        eMesh.getBulkData()->modification_begin();

        for (unsigned irank=1; irank < m_id_max.size(); irank++)
          {
            if (irank == FAMILY_TREE_RANK) continue;
            // bucket loop
            const std::vector<stk::mesh::Bucket*> & buckets = eMesh.getBulkData()->buckets( irank );

            for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (selector(**k))
                {
                  stk::mesh::Bucket & bucket = **k ;
                  const unsigned num_entities_in_bucket = bucket.size();

                  for (unsigned iEntity = 0; iEntity < num_entities_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity& entity = bucket[iEntity];
                      if (eMesh.isChildElement(entity))
                        {
                          stk::mesh::EntityId id = m_id_max[irank] + 1;
                          m_id_max[irank] = id;
                          if (m_debug) std::cout << "SerializeNodeRegistry::resetNewElementIds irank= " << irank 
                                                 << " old id= " << entity.identifier() << " new id= " << id << std::endl;
                          eMesh.getBulkData()->change_entity_id(id, entity);
                        }
                    }
                }
              }
          }
        eMesh.getBulkData()->modification_begin();
      }

    };

#else
    class SerializeNodeRegistry {
    public:
      SerializeNodeRegistry(PerceptMesh& eMesh, NodeRegistry* nodeRegistry, std::string input_mesh_name, std::string output_mesh_name, int M, int iM, int P=1, int iP=0) {}
    };
#endif

  }
}
