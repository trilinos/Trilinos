// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#if HAVE_YAML

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <stdint.h>
#include <tuple>
#include <unordered_map>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <percept/Percept.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/RunEnvironment.hpp>
#include <percept/ProgressMeter.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_topology/topology.hpp>

#include <adapt/RefinerUtil.hpp>

#include <adapt/UniformRefiner.hpp>
#include <adapt/UniformRefinerPattern.hpp>

#include <stk_io/IossBridge.hpp>
#include <Ioss_NullEntity.h>
#include <Ioss_Utils.h>

#define DEBUG_YAML 0

#include <yaml-cpp/yaml.h>

#define YAML_CHECK(emitter) do { if (1 && !emitter.good()) { std::cout << "Emitter error: " << __FILE__ << ":" << __LINE__ << " emitter.good()= " \
                                                                       << emitter.good() << " Error Message: " << emitter.GetLastError() << std::endl; return;} } while(0)

#define YAML_ERRCHECK YAML_CHECK(emitter)

#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept {

    inline stk::topology find_topology_by_name(const std::string & topologyName)
    {
      stk::topology t(stk::topology::BEGIN_TOPOLOGY);
      for ( ; t < stk::topology::END_TOPOLOGY; ++t) {
        if (topologyName == t.name()) {
          return t;
        }
      }
      return stk::topology(stk::topology::INVALID_TOPOLOGY);
    }

    /**  in the following, we are looking at one partition, m_iM, of M-partitioned mesh
     *
     */

    class SerializeNodeRegistry {

    public:
      typedef std::vector<stk::mesh::EntityId> IdVector;

      // PartMap maps a part name to data about the part: its primary entity rank, topology name, subset parts
      typedef std::string TopologyName;
      typedef std::string PartName;
      typedef std::vector<PartName> PartSubsets;
      typedef std::tuple<stk::mesh::EntityRank, TopologyName, PartSubsets> PartMapData;
      typedef std::map<PartName, PartMapData> PartMap;

      typedef stk::mesh::EntityId NodeMapKey;
      typedef int ProcRank;
      typedef std::vector<ProcRank> SharedProcs;
      typedef SharedProcs NodeMapValue;
      typedef std::unordered_map<NodeMapKey, NodeMapValue> NodeMap;

    private:
      PerceptMesh& m_eMesh;
      NodeRegistry* m_nodeRegistry;
      std::string m_input_mesh_name;
      std::string m_output_mesh_name;
      std::string m_filePrefix;
      // iM, M = pseudo processor (reads file.e.M.iM)
      // iP, P = worker, nworker
      int m_M, m_iM, m_W, m_iW;
      // for W==1, [start,end] = [0,M-1]
      //     W==2, iW==0, [M_0,M_1] = [0,M/2-1] ,
      //           iW==1, = [M/2,M-1]
      // etc..
      int m_M_0, m_M_1;
      const std::vector<std::string> & m_entity_rank_names;
      IdVector m_id_max;
      std::string m_globalIdFile;
      std::string m_localNodeRegistryFile;
      std::string m_globalNodeRegistryFile;
      std::string m_globalPartsFile;
      const stk::mesh::EntityRank FAMILY_TREE_RANK;
      static const bool m_debug = false;
      int m_spatialDim;
      PartMap *m_partMap;
      NodeMap *m_nodeMap;
      std::string m_geomFile;

    public:
      enum { MaxPass = 3 };

      SerializeNodeRegistry(PerceptMesh& eMesh, NodeRegistry* nodeRegistry, std::string input_mesh_name, std::string output_mesh_name,
                            int M, int iM, int W=1, int iW=0, int M_0=-1, int M_1=-1) :
        m_eMesh(eMesh), m_nodeRegistry(nodeRegistry), m_input_mesh_name(input_mesh_name), m_output_mesh_name(output_mesh_name), m_filePrefix(input_mesh_name),
        m_M(M), m_iM(iM), m_W(W), m_iW(iW), m_M_0(M_0 >= 0 ? M_0 : 0), m_M_1(M_1 >= 0 ? M_1 : M-1),
        m_entity_rank_names(eMesh.get_fem_meta_data()->entity_rank_names()),
        m_id_max(m_entity_rank_names.size(), 0u),  FAMILY_TREE_RANK(static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u)),
        m_partMap(0), m_nodeMap(0), m_geomFile("")
      {
        size_t pos = m_filePrefix.find(".");
        if (pos != std::string::npos)
          m_filePrefix = m_filePrefix.substr(0, pos);
        m_globalIdFile = "global_id_data.yaml.W."+toString(m_W)+"."+toString(m_iW);
        m_localNodeRegistryFile = std::string("local_nodeRegistry.yaml.")+toString(m_M)+"."+toString(m_iM)+"."+toString(m_W)+"."+toString(m_iW);
        m_globalNodeRegistryFile = std::string("global_nodeRegistry.yaml.W."+toString(m_W)+"."+toString(m_iW));
        m_globalPartsFile = "global_parts.yaml.W."+toString(m_W)+"."+toString(m_iW);
        m_spatialDim = eMesh.get_spatial_dim();
      }

      ~SerializeNodeRegistry() {
        if (m_partMap) delete m_partMap;
        if (m_nodeMap) delete m_nodeMap;
      }

      void set_geometry_file(std::string geomFile) { m_geomFile = geomFile; }

      void pass(int streaming_pass)
      {
        std::cout << "\n\nM[" << m_iM << ", " << m_M << "] W[" << m_iW << ", " << m_W << "] ------ SerializeNodeRegistry ----- pass number " << streaming_pass << "\n\n" << std::endl;

        switch(streaming_pass)
          {
          case -1: passM1(); return;
          case 0: pass0(); return;
          case 1: pass1(); return;
          case 2: pass2(); return;
          case 3: pass3(); return;
          default:
            throw std::logic_error("SerializeNodeRegistry::pass unknown pass");
          }
      }

      void pass_init(int streaming_pass)
      {
        std::cout << "\n\nM[" << m_iM << ", " << m_M << "] W[" << m_iW << ", " << m_W << "] ------ SerializeNodeRegistry ----- init pass number " << streaming_pass << "\n\n" << std::endl;

        switch(streaming_pass)
          {
          case -1: passM1_init(); return;
          case 0: pass0_init(); return;
            //           case 1: pass1_init(); return;
            //           case 2: pass2_init(); return;
            //           case 3: pass3_init(); return;
            //           default:
            //             throw std::logic_error("SerializeNodeRegistry::pass unknown pass");
          }
      }

      void pass_final(int streaming_pass)
      {
        std::cout << "\n\nM[" << m_iM << ", " << m_M << "] W[" << m_iW << ", " << m_W << "] ------ SerializeNodeRegistry ----- final pass number " << streaming_pass << "\n\n" << std::endl;

        switch(streaming_pass)
          {
          case -1: passM1_final(); return;
          case 0: pass0_final(); return;
//           case 1: pass1_final(); return;
          case 2: pass2_final(); return;
//           case 3: pass3_final(); return;
//           default:
//             throw std::logic_error("SerializeNodeRegistry::pass unknown pass");
          }
      }

      //CellTopology get_cell_topology( const Part & part) const;

      /**
       *   passM1: open unrefined mesh, get Part (exodus block) information, put to global yaml file
       *   (m_iM = 0...m_M)
       */
      // part primary_entity_rank, topology name
      void writeGlobalPartsFile()
      {
        m_globalPartsFile = "global_parts.yaml.W."+toString(m_W)+"."+toString(m_iW);
        writeGlobalPartsFile(m_globalPartsFile, *m_partMap);
      }

      void writeGlobalPartsFile(std::string fileName, PartMap& partMap)
      {
        std::fstream file;
        file.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc);
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::writeGlobalPartsFile couldn't open file ")+fileName);
          }
        YAML::Emitter out;
        out << YAML::BeginMap; YAML_CHECK(out);
        PartMap::iterator iter;
        for (iter = partMap.begin(); iter != partMap.end(); ++iter)
          {
            out << YAML::Key << iter->first;  YAML_CHECK(out);
            out << YAML::Value;               YAML_CHECK(out);
            out << YAML::Flow;                YAML_CHECK(out);
            out << YAML::BeginSeq;            YAML_CHECK(out);
            out << std::get<0>(iter->second);
            out << std::get<1>(iter->second);
            {
              out << YAML::Flow;                YAML_CHECK(out);
              out << YAML::BeginSeq;            YAML_CHECK(out);
              PartSubsets subsets = std::get<2>(iter->second);
              for (unsigned isub=0; isub < subsets.size(); isub++)
                {
                  out << subsets[isub];
                }
              out << YAML::EndSeq;              YAML_CHECK(out);
            }
            out << YAML::EndSeq;              YAML_CHECK(out);
          }
        out << YAML::EndMap; YAML_CHECK(out);
        file << out.c_str();
        file.close();
      }

      void declareGlobalParts()
      {
        if (!m_partMap) readGlobalPartsFile();
        PartMap::iterator iter;
        for (iter = m_partMap->begin(); iter != m_partMap->end(); ++iter)
          {
            PartName part_name = iter->first;
            stk::mesh::EntityRank part_rank = std::get<0>(iter->second);
            TopologyName topo_name = std::get<1>(iter->second);
            const stk::mesh::Part* c_part = m_eMesh.get_fem_meta_data()->get_part(part_name);
            stk::mesh::Part *part = const_cast<stk::mesh::Part *>(c_part);
            if (!part)
              {
                part = &m_eMesh.get_fem_meta_data()->declare_part(part_name, part_rank);
                stk::io::put_io_part_attribute(*part);
                stk::topology topo = find_topology_by_name(topo_name.c_str());
                if (!topo.is_valid())
                  {
                    std::cout << "bad cell topo SerializeNodeRegistry::declareGlobalParts topo_name= " << topo_name << std::endl;
                    throw std::runtime_error("bad cell topo SerializeNodeRegistry::declareGlobalParts");
                  }
                stk::mesh::set_topology(*part, topo);
              }
          }
        // once parts are declared, declare part subsets
        for (iter = m_partMap->begin(); iter != m_partMap->end(); ++iter)
          {
            PartName part_name = iter->first;
            const stk::mesh::Part* c_part = m_eMesh.get_fem_meta_data()->get_part(part_name);
            stk::mesh::Part *part = const_cast<stk::mesh::Part *>(c_part);
            if (!part)
              {
                throw std::runtime_error(std::string("no part found SerializeNodeRegistry::declareGlobalParts: part= ")+part_name);
              }
            PartSubsets subsets = std::get<2>(iter->second);
            for (unsigned isub=0; isub < subsets.size(); isub++)
              {
                const stk::mesh::Part* c_subset = m_eMesh.get_fem_meta_data()->get_part(subsets[isub]);
                stk::mesh::Part *subset = const_cast<stk::mesh::Part *>(c_subset);

                if (!subset)
                  {
                    throw std::runtime_error(std::string("no subset part found SerializeNodeRegistry::declareGlobalParts: subset part= ")+subsets[isub]);
                  }
                m_eMesh.get_fem_meta_data()->declare_part_subset(*part, *subset);
              }
          }
      }

      // part primary_entity_rank, topology name
      void readGlobalPartsFile()
      {
        if (m_partMap) delete m_partMap;
        m_partMap = new PartMap;

        readGlobalPartsFile(m_globalPartsFile, *m_partMap);
      }

      // merge parts if already there
      void readGlobalPartsFile(std::string fileName, PartMap& partMap)
      {
        std::fstream file;
        file.open(fileName.c_str(), std::ios_base::in);
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::readGlobalPartsFile couldn't open file ")+fileName);
          }

        //YAML::Parser parser(file);
        YAML::Node doc;

        try {
          //while(parser.GetNextDocument(doc)) {
          if (1) {
            doc = YAML::Load(file);
            //std::cout << "\n readGlobalPartsFile doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (doc.Type() == YAML::NodeType::Map)
              {
                for(YAML::const_iterator iter=doc.begin();iter!=doc.end();++iter)
                  {
                    const YAML::Node key = iter->first;
                    PartName part_name = key.as<PartName>();

                    const YAML::Node valSeq = iter->second;
                    UInt rank_input;
                    TopologyName topo_name;
                    YAML::const_iterator itv=valSeq.begin();
                    rank_input = itv->as<UInt>();
                    ++itv;
                    topo_name = itv->as<TopologyName>();
                    ++itv;
                    stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(rank_input);
                    const YAML::Node subsetSeq = *itv;
                    YAML::const_iterator iss;
                    PartSubsets subsets;
                    for (iss = subsetSeq.begin(); iss != subsetSeq.end(); ++iss)
                      {
                        PartName subset_name = iss->as<PartName>();
                        subsets.push_back(subset_name);
                      }
                    PartMapData pmd(rank, topo_name, subsets);
                    if (0)
                      {
                        partMap[part_name] = pmd;
                      }
                    else
                      {
                        PartMapData& pmd_existing = partMap[part_name];
                        std::get<0>(pmd_existing) = rank;
                        std::get<1>(pmd_existing) = topo_name;
                        PartSubsets& subsets_existing = std::get<2>(pmd_existing);
                        //std::cout << "pmd read= " << pmd << std::endl;
                        //std::cout << "pmd existing= " << pmd_existing << std::endl;
                        for (unsigned isub=0; isub < subsets.size(); isub++)
                          {
                            if (std::find(subsets_existing.begin(), subsets_existing.end(), subsets[isub]) == subsets_existing.end())
                              {
                                subsets_existing.push_back(subsets[isub]);
                              }
                          }
                      }
                  }
              }
            else
              {
                throw std::runtime_error("bad global_parts file");
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << "\n";
          throw std::runtime_error( e.what());
        }
        file.close();
      }

      /**
       *  loop over all files (m_iM...)
       *    loop over all my nodes
       *      add my proc rank to map value
       *
       *  cull node map (only nodes that are shared remain)
       *
       *  write nodeMap to file by procId: shared node with proc (map<procId, map<node,procOwner> > )
       *    loop over nodes, procs sharing node, NodeMapVector[procs][node] = owner
       *
       *  put NodeMapVector[proc] out in separate file
       *
       *  NodeMap read back in gives shared nodes mapped to owner - so, if not in map, not shared
       */

      void writeNodeMap(std::string fileName, int jM, bool is_local)
      {
        std::fstream file;
        file.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc);
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::writeNodeMap couldn't open file ")+fileName);
          }
        YAML::Emitter out;
        out << YAML::BeginMap; YAML_CHECK(out);
        NodeMap::iterator iter;
        for (iter = m_nodeMap->begin(); iter != m_nodeMap->end(); ++iter)
          {
            NodeMapValue& procs = iter->second;
            bool include_entry = true;
            if (is_local)
              {
                if (procs.size() <= 1) throw std::logic_error("SerializeNodeRegistry::writeNodeMap procs.size is <=1");
                NodeMapValue::iterator find_my_proc = std::find(procs.begin(), procs.end(), jM);
                include_entry = include_entry && (find_my_proc != procs.end());  // shared by me
              }
            if (include_entry)
              {
                if (procs[0] == jM)  // I own it
                  {
                  }

                out << YAML::Key << iter->first;     YAML_CHECK(out);
                out << YAML::Value;                  YAML_CHECK(out);
                // not strictly necessary, but we'll leave it as a sequence to be consistent with NodeMap type
                out << YAML::Flow;                YAML_CHECK(out);
                out << YAML::BeginSeq;            YAML_CHECK(out);
                if (is_local)
                  {
                    // put out only the owner of the node
                    out << procs[0];                  YAML_CHECK(out);
                  }
                else
                  {
                    // put out all procs shared by this node
                    for (unsigned ii = 0; ii < procs.size(); ii++)
                      {
                        out << procs[ii];                  YAML_CHECK(out);
                      }
                  }
                out << YAML::EndSeq;              YAML_CHECK(out);
              }
          }
        out << YAML::EndMap; YAML_CHECK(out);
        file << out.c_str();
        file.close();
      }

      /**
       * file with all node/proc pairs to determine shared/non-shared nodes and ownership
       */
      void createLocalNodeMapFiles(int iW)
      {
        VERIFY_OP_ON(m_nodeMap, !=, 0, "SerializeNodeRegistry::createLocalNodeMapFiles m_nodeMap is null");

        for (int jM = m_M_0; jM <= m_M_1; jM++)
          {
            std::string localNodeMapFile = std::string("local_nodeMap.yaml")+"."+toString(m_M)+"."+toString(jM)+"."+toString(m_W)+"."+toString(iW);

            writeNodeMap(localNodeMapFile, jM, true);
          }
      }

      /// Remove nodes that are only touched by one processor, and are thus not shared
      void cullNodeMap()
      {
        NodeMap::iterator iter;
        std::vector<NodeMapKey> non_shared;
        for (iter = m_nodeMap->begin(); iter != m_nodeMap->end(); ++iter)
          {
            const NodeMapValue& procs = iter->second;
            if (procs.size() == 1)
              {
                non_shared.push_back(iter->first);
              }
          }
        for (unsigned i = 0; i < non_shared.size(); i++)
          {
            m_nodeMap->erase(non_shared[i]);
          }
      }

      // node, proc owner
      void readNodeMap(std::string fileName, NodeMap& nodeMap, bool is_local)
      {
        std::fstream file;
        file.open(fileName.c_str(), std::ios_base::in);
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::readNodeMap couldn't open file ")+fileName);
          }

        //YAML::Parser parser(file);
        YAML::Node doc;

        try {
          //while(parser.GetNextDocument(doc)) {
          if (1) {
            doc = YAML::Load(file);
            std::cout << "tmp srk readNodeMap= " << fileName
                      << " doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (doc.Type() == YAML::NodeType::Map)
              {
                for(YAML::const_iterator iter=doc.begin();iter!=doc.end();++iter)
                  {
                    const YAML::Node key = iter->first;
                    stk::mesh::EntityId id = key.as<stk::mesh::EntityId>();
                    NodeMapValue procs;
                    const YAML::Node val = iter->second;
                    procs = val.as<NodeMapValue>();
                    //std::cout << "readNodeMap id= " << id << " procs= " << procs << std::endl;
                    if (is_local && procs.size() != 1)
                      throw std::logic_error(std::string("SerializeNodeRegistry::readNodeMap procs.size is != 1, = ")+toString(procs.size()));
                    NodeMapValue& procs_existing = nodeMap[id];
                    for (unsigned ip=0; ip < procs.size(); ip++)
                      {
                        if (std::find(procs_existing.begin(), procs_existing.end(), procs[ip]) == procs_existing.end())
                          {
                            procs_existing.push_back(procs[ip]);
                          }
                      }
                    std::sort(procs_existing.begin(), procs_existing.end());
                    //std::cout << "id, procs= " << id << " " << procs_existing << std::endl;
                    //nodeMap[id] = procs;
                  }
              }
            else
              {
                throw std::runtime_error("bad nodeMap file");
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << "\n";
          throw std::runtime_error( e.what());
        }
        file.close();
      }

      void readLocalNodeMapFile(NodeMap& nodeMap)
      {
        std::string localNodeMapFile = std::string("local_nodeMap.yaml")+"."+toString(m_M)+"."+toString(m_iM)+"."+toString(m_W)+"."+toString(m_iW);
        readNodeMap(localNodeMapFile, nodeMap, true);
      }

      // creates m_partMap if null, reads in global part map from file
      void getGlobalPartMap()
      {
        readGlobalPartsFile();
      }

      void setGlobalPartMap()
      {
        writeGlobalPartsFile();
      }

      void initializeGlobalPartMap()
      {
        if (m_partMap) delete m_partMap;
        m_partMap = new PartMap;
        setGlobalPartMap();
      }

      void getGlobalNodeMap()
      {
        std::cout << "tmp srk SerializeNodeRegistry::getGlobalNodeMap, m_nodeMap= " << m_nodeMap << std::endl;
        if (m_nodeMap) delete m_nodeMap;
        m_nodeMap = new NodeMap;
        std::string globalNodeMapFile = std::string("global_nodeMap.yaml")+".W."+toString(m_W)+"."+toString(m_iW);
        readNodeMap(globalNodeMapFile, *m_nodeMap, false);
      }

      void setGlobalNodeMap()
      {
        std::cout << "tmp srk SerializeNodeRegistry::setGlobalNodeMap\n";
        std::string globalNodeMapFile = std::string("global_nodeMap.yaml")+".W."+toString(m_W)+"."+toString(m_iW);
        writeNodeMap(globalNodeMapFile, 0, false);
      }

      void initializeGlobalNodeMap()
      {
        std::cout << "tmp srk SerializeNodeRegistry::initializeGlobalNodeMap\n";
        if (m_nodeMap) delete m_nodeMap;
        m_nodeMap = new NodeMap;
        setGlobalNodeMap();
      }


      void passM1_mergeGlobalParts()
      {
        const stk::mesh::PartVector& parts = m_eMesh.get_fem_meta_data()->get_parts();
        for (unsigned ipart=0; ipart < parts.size(); ipart++)
          {
            //const Part& part = *parts[ipart];
            stk::mesh::Part& part = *parts[ipart];
            if (stk::mesh::is_auto_declared_part(part))
              {
                continue;
              }

            // only io parts allowed
            if (!stk::io::is_part_io_part(part))
              {
                if (m_debug) std::cout << "part without io attr= " << part.name() << std::endl;
                continue;
              }

            stk::topology topo = m_eMesh.get_fem_meta_data()->get_topology(part);
            TopologyName topo_name = "null";
            if (topo != stk::topology::INVALID_TOPOLOGY) topo_name = topo.name();
            //std::cout << "part, topo= " << part.name() << " " << topo_name << std::endl;
            const stk::mesh::PartVector& part_subsets = part.subsets();
            PartSubsets subsets(part_subsets.size());
            for (unsigned isub=0; isub < subsets.size(); isub++)
              {
                subsets[isub] = part_subsets[isub]->name();
              }
            PartMapData pmd {part.primary_entity_rank(), topo_name, subsets};
            PartMap::iterator inMap = m_partMap->find(part.name());
            // if not yet in the map, or if already in map, choose to overwrite if this is the one with subset info
            if (inMap == m_partMap->end())
              {
                (*m_partMap)[part.name()] = pmd;
              }
            else
              {
                PartSubsets& inMapSubsets = std::get<2>(inMap->second);
                for (unsigned isub=0; isub < subsets.size(); isub++)
                  {
                    PartSubsets::iterator iv = find(inMapSubsets.begin(), inMapSubsets.end(), subsets[isub]);
                    if (iv == inMapSubsets.end())
                      {
                        inMapSubsets.push_back(subsets[isub]);
                      }
                  }
              }
          }
        // init/fini
        if (m_W == 1 && m_iM == m_M_1)
          {
            writeGlobalPartsFile();
          }
      }

      // open each global nodes file and add to the global node map
      void createGlobalNodeMap()
      {
        if (m_nodeMap) delete m_nodeMap;
        m_nodeMap = new NodeMap;
        for (m_iM = m_M_0; m_iM <= m_M_1; m_iM++)
          {
            readLocalNodeMapFile(*m_nodeMap);
          }
      }

      void passM1_createLocalNodeMap()
      {
        const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (selector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_entities_in_bucket = bucket.size();

              for (unsigned iEntity = 0; iEntity < num_entities_in_bucket; iEntity++)
                {
                  stk::mesh::Entity entity = bucket[iEntity];
                  stk::mesh::EntityId id = m_eMesh.identifier(entity);
                  NodeMapValue& procs = (*m_nodeMap)[id];
                  if (procs.size() == 0)
                    {
                      procs.push_back(m_iM);
                    }
                  else
                    {
                      NodeMapValue::iterator it = std::find(procs.begin(), procs.end(), m_iM);
                      if (it == procs.end())
                        {
                          procs.push_back(m_iM);
                          std::sort(procs.begin(), procs.end());
                        }
                    }
                }
            }
          }
        // init/fini
        if (m_W == 1 && m_iM == m_M_1)
          {
            cullNodeMap();
            createLocalNodeMapFiles(m_iW);
          }
      }

      static void getStreamingPiece(int M, int W, int iW, int& M_0, int& M_1)
      {
        VERIFY_OP_ON(W, >, 0, "AdaptMain::getStreamingPiece #0");
        VERIFY_OP_ON(M, >, 0, "AdaptMain::getStreamingPiece #1");
        VERIFY_OP_ON(M % W, ==, 0, "AdaptMain::getStreamingPiece: streaming_size M must be divisible by number of workers W.");
        int len = M / W;
        M_0 = iW*len;
        M_1 = (iW+1)*len - 1;   // we use these as start/end, as in for(int iM = M_0; iM <= M_1; iM++)
      }

      void pass0_init()
      {
        std::cout << "pass0_init\n";
        // create initial file, write 0's in it  for m_id_max (it should be initialized to 0, but just to be sure, we reset it here)
        for (unsigned irank=0; irank < m_id_max.size(); irank++)
          m_id_max[irank]=0u;
        for (int jW = 0; jW < m_W; jW++)
          {
            std::string globalIdFile = "global_id_data.yaml.W."+toString(m_W)+"."+toString(jW);
            std::cout << "pass0_init globalIdFile= " << globalIdFile << std::endl;
            setCurrentGlobalMaxId(globalIdFile, m_id_max);
          }
      }

      void pass0_final()
      {
        std::cout << "pass0_final, m_id_max= " << m_id_max << std::endl;
        for (unsigned irank=0; irank < m_id_max.size(); irank++)
          {
            m_id_max[irank]=0u;
          }

        // first, get the max of all
        IdVector id_max_local_0(m_id_max.size(), 0u);
        // now get how many each iW used
        IdVector id_used(m_id_max.size(), 0u);
        std::vector<IdVector> num_id_used(m_W);
        for (int jW = 0; jW < m_W; jW++)
          {
            std::string globalIdFile = "global_id_data.yaml.W."+toString(m_W)+"."+toString(jW);
            getCurrentGlobalMaxId(globalIdFile, id_max_local_0, true);
            getCurrentGlobalMaxId(globalIdFile, id_used, false);
            num_id_used[jW] = id_used;
            std::cout << "pass0_final: num_id_used = " << num_id_used[jW] << std::endl;
          }

        std::vector<IdVector> id_max_new(m_W);
        id_max_new[0] = id_max_local_0;
        std::cout << "pass0_final: id_max_new[0] = " << id_max_new[0] << std::endl;
        for (int jW = 1; jW < m_W; jW++)
          {
            id_max_new[jW].resize(m_id_max.size());
            for (unsigned irank=0; irank < id_max_local_0.size(); irank++)
              {
                id_max_new[jW][irank] = id_max_new[jW-1][irank] + (num_id_used[jW-1][irank] + 1);
              }
            std::cout << "pass0_final: id_max_new[jW] = " << id_max_new[jW] << std::endl;
          }

        for (int jW = 0; jW < m_W; jW++)
          {
            std::string globalIdFileWrite = "global_id_data.yaml.W."+toString(m_W)+"."+toString(jW);
            setCurrentGlobalMaxId(globalIdFileWrite, id_max_new[jW]);
          }
        //         for (int jW = 0; jW < m_W; jW++)
        //           {
        //             setCurrentGlobalMaxId(globalIdFile, m_id_max);
        //           }
        //printCurrentGlobalMaxId("pass0");
      }

      void passM1_init()
      {
        // init/fini
        VERIFY_OP_ON(m_iW, ==, -1, "passM1_init #0");
        for (m_iW = 0; m_iW < m_W; ++m_iW)
          {
            initializeGlobalPartMap();
            initializeGlobalNodeMap();
          }
        m_iW = -1;
      }

      void passM1_final()
      {
        PartMap partMap;
        for (int jW = 0; jW < m_W; jW++)
          {
            std::string globalPartMapFile = std::string("global_parts.yaml")+".W."+toString(m_W)+"."+toString(jW);
            // does a merge...
            readGlobalPartsFile(globalPartMapFile, partMap);
          }
        for (int jW = 0; jW < m_W; jW++)
          {
            std::string globalPartMapFile = std::string("global_parts.yaml")+".W."+toString(m_W)+"."+toString(jW);
            //std::string globalPartMapFile = std::string("global_parts.yaml")+".W."+toString(m_W)+"."+toString(m_W);
            writeGlobalPartsFile(globalPartMapFile, partMap);
          }

        if (m_nodeMap) delete m_nodeMap;
        m_nodeMap = new NodeMap;
        for (int jW = 0; jW < m_W; jW++)
          {
            std::string globalNodeMapFile = std::string("global_nodeMap.yaml")+".W."+toString(m_W)+"."+toString(jW);
            readNodeMap(globalNodeMapFile, *m_nodeMap, false);
          }
        cullNodeMap();
        std::cout << "m_nodeMap->size() = " << m_nodeMap->size() << " m_W= " << m_W << std::endl;
        for (int jW = 0; jW < m_W; jW++)
          {
            SerializeNodeRegistry::getStreamingPiece(m_M, m_W, jW, m_M_0, m_M_1);
            createLocalNodeMapFiles(jW);
          }
      }

      void passM1()
      {
        for (int m_iM = m_M_0; m_iM <= m_M_1; m_iM++)
          {
            //input_mesh = input_mesh_save+"."+toString(M)+"."+toString(m_iM);
            std::string input_mesh = Ioss::Utils::decode_filename(m_input_mesh_name, m_iM, m_M);
            PerceptMesh eMesh(m_spatialDim);
            eMesh.open_read_only(input_mesh);
            NodeRegistry *some_nr = 0;
            SerializeNodeRegistry snr(eMesh, some_nr, m_input_mesh_name, m_output_mesh_name, m_M, m_iM,  m_W, m_iW, m_M_0, m_M_1);

            // init/fini
            if (m_W == 1 && m_iM == m_M_0)
              {
                snr.initializeGlobalPartMap();
                snr.initializeGlobalNodeMap();
              }

            snr.getGlobalPartMap();
            snr.passM1_mergeGlobalParts();
            snr.setGlobalPartMap();

            snr.getGlobalNodeMap();
            snr.passM1_createLocalNodeMap();
            snr.setGlobalNodeMap();

          }
      }

      /**
       *   pass0: open unrefined mesh, refine it, find max id
       *   (m_iM = 0...M)
       *   1. if m_iM==0, setCurrentGlobalMaxId to [0,0,0,0]
       *   2. getCurrentGlobalMaxId()
       *   3. find new max id from current mesh
       *   4. setCurrentGlobalMaxId()
       */
      void pass0()
      {
        // init/fini
        if (m_W == 1 && m_iM == m_M_0 )
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
        //if (m_debug)
        std::cout << "M[" << m_iM << ", " << m_M << "] W[" << m_iW << ", " << m_W << "] SerializeNodeRegistry::printCurrentGlobalMaxId: = " << m_id_max << " for: " << msg << std::endl;
      }

      /**
       *   pass1: refine mesh, write local NodeRegistry, set new max id from refined mesh - assumes parent elements exist
       *
       *   (m_iM = 0...M)
       *   1. open file.e.M.m_iM, refine mesh
       *   2.
       *   3. getCurrentGlobalMaxId()
       *   4. resetNewElementIds() (resets new element ids by looking at child elements only)
       *   5. write NodeRegistry in name.yaml.M.m_iM
       *   6. setCurrentGlobalMaxId()
       *   7. save refined mesh
       */
      void pass1()
      {
        NodeRegistry& nodeRegistry = *m_nodeRegistry;
        getCurrentGlobalMaxId();
        resetNewElementIds(m_eMesh, nodeRegistry);
        resetNewNodeIds(nodeRegistry);
        setCurrentGlobalMaxId();
        printCurrentGlobalMaxId("pass1");

        if (m_nodeMap) delete m_nodeMap;
        m_nodeMap = new NodeMap;
        readLocalNodeMapFile(*m_nodeMap);
        writeNodeRegistry(nodeRegistry, m_localNodeRegistryFile);
      }

      /**
       *   pass2 - create global NodeRegistry from each local one by "last one wins"
       *   (single call, no loop over m_iM)
       *   1. create new (global) NodeRegistry
       *   2. getCurrentGlobalMaxId
       *   3. loop m_iM
       *      a. read NodeRegistry from name.yaml.M.m_iM -> input values into new NodeRegistry
       *   4. for each key/value pair, increment idserver, save new id in value
       *   5. write new global NodeRegistry
       *
       */
      void pass2()
      {
        if (m_iM != m_M_0) throw std::logic_error("SerializeNodeRegistry::pass2 logic error");
        if (m_W > 1) return;  // done in pass2_final

        getCurrentGlobalMaxId();
        printCurrentGlobalMaxId("pass2 initial");

        PerceptMesh eMeshNew(m_spatialDim);
        eMeshNew.openEmpty();
        NodeRegistry globalNR(eMeshNew);

        PerceptMesh eMeshLocal(m_spatialDim);
        eMeshLocal.openEmpty();

        for (m_iM = m_M_0; m_iM <= m_M_1; m_iM++)
          {
            NodeRegistry newLocalNR(eMeshLocal);
            m_localNodeRegistryFile = std::string("local_nodeRegistry.yaml.")+toString(m_M)+"."+toString(m_iM)+"."+toString(m_W)+"."+toString(m_iW);
            readNodeRegistry(newLocalNR, m_localNodeRegistryFile);
            addLocalNodeRegistryToGlobal(newLocalNR, globalNR);
          }

        createGlobalNodeMap();

        writeNodeRegistry(globalNR, m_globalNodeRegistryFile);
        setCurrentGlobalMaxId();
        printCurrentGlobalMaxId("pass2 done");
      }

      void pass2_final()
      {
//         getCurrentGlobalMaxId();
//         printCurrentGlobalMaxId("pass2 initial");

        PerceptMesh eMeshNew(m_spatialDim);
        eMeshNew.openEmpty();
        NodeRegistry globalNR(eMeshNew);

        PerceptMesh eMeshLocal(m_spatialDim);
        eMeshLocal.openEmpty();

        if (m_nodeMap) delete m_nodeMap;
        m_nodeMap = new NodeMap;

        int iWSave = m_iW;
        for (int jW = 0; jW < m_W; jW++)
          {
            SerializeNodeRegistry::getStreamingPiece(m_M, m_W, jW, m_M_0, m_M_1);
            m_iW = jW;

            for (m_iM = m_M_0; m_iM <= m_M_1; m_iM++)
              {
                NodeRegistry newLocalNR(eMeshLocal);
                std::string localNodeRegistryFile = std::string("local_nodeRegistry.yaml.")+toString(m_M)+"."+toString(m_iM)+"."+toString(m_W)+"."+toString(jW);
                readNodeRegistry(newLocalNR, localNodeRegistryFile);
                addLocalNodeRegistryToGlobal(newLocalNR, globalNR);

                // replacement for: createGlobalNodeMap();
                readLocalNodeMapFile(*m_nodeMap);
              }
          }
        m_iW = iWSave;

        for (int jW = 0; jW < m_W; jW++)
          {
            std::string globalNodeRegistryFile = "global_nodeRegistry.yaml.W."+toString(m_W)+"."+toString(jW);
            writeNodeRegistry(globalNR, globalNodeRegistryFile);
          }
        //setCurrentGlobalMaxId();
        //printCurrentGlobalMaxId("pass2 done");
      }

      /**
       *   pass 3:
       *   1. read global NodeRegistry
       *   2. read each refined file-ref.M.m_iM
       *   3. lookup edge/face/elem in global NodeRegistry, reset id to that found in NR
       *   4. write refined file-ref-id.M.m_iM
       */
    private:
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
        readNodeRegistry(localNR, m_localNodeRegistryFile);
        readNodeRegistry(globalNR, m_globalNodeRegistryFile);
        lookupAndSetNewNodeIds(localNR, globalNR);
      }

    public:
      // single call instead of multiple loops on iM from caller - sets global parts, sets new node Ids from global nodeRegistry
      void pass3_new()
      {
        for (int jM = m_M_0; jM <= m_M_1; jM++)
          {
            std::string input_mesh_new = m_output_mesh_name+"-pass1";
            input_mesh_new = Ioss::Utils::decode_filename(input_mesh_new, jM, m_M);
            std::string output_mesh_new = Ioss::Utils::decode_filename(m_output_mesh_name, jM, m_M);
            PerceptMesh eMesh;
            eMesh.open(input_mesh_new);
            m_spatialDim = eMesh.get_spatial_dim();
            //eMesh.setStreamingSize(m_M);
            NodeRegistry *some_nr0 = 0;
            SerializeNodeRegistry snr0(eMesh, some_nr0, m_input_mesh_name, m_output_mesh_name, m_M, jM,  m_W, m_iW, m_M_0, m_M_1);
            snr0.getGlobalPartMap();
            snr0.declareGlobalParts();
            eMesh.commit();

            NodeRegistry *some_nr = 0;
            SerializeNodeRegistry snr(eMesh, some_nr, m_input_mesh_name, m_output_mesh_name, m_M, jM,  m_W, m_iW, m_M_0, m_M_1);
            snr.pass(3);
            std::cout << "tmp srk debug SerializeNodeRegistry: " << PERCEPT_OUT(m_M) << PERCEPT_OUT(jM) << PERCEPT_OUT(m_W) << PERCEPT_OUT(m_iW) << PERCEPT_OUT(m_M_0) << PERCEPT_OUT(m_M_1) << " output_mesh_new= " << output_mesh_new << std::endl;
            if (m_geomFile != "")
              eMesh.remove_geometry_blocks_on_output(m_geomFile);
            eMesh.save_as(output_mesh_new);
          }
      }

      // for each subDimEntity in local NodeRegistry, find it in the global NodeRegistry and reset ID's to the global values
      void lookupAndSetNewNodeIds(NodeRegistry& localNR, NodeRegistry& globalNR)
      {
        m_eMesh.get_bulk_data()->modification_begin();
        SubDimCellToDataMap& localMap = localNR.getMap();
        //SubDimCellToDataMap& globalMap = globalNR.getMap();
        //std::cout << " tmp SerializeNodeRegistry::lookupAndSetNewNodeIds localMap size: " << localMap.size() << std::endl;

        SubDimCellToDataMap::iterator iter;
        for (iter = localMap.begin(); iter != localMap.end(); ++iter)
          {
            const SubDimCell_SDCEntityType& subDimEntity = iter->first;
            SubDimCellData& nodeId_elementOwnerId = iter->second;

            // special case for "interior" subDimEntity's which are centroid nodes for quad or hex elements -
            //   by definition they aren't shared
            if (subDimEntity.size() == 1)
              continue;

            // lookup from global...
            SubDimCellData* global_nodeId_elementOwnerId_ptr = globalNR.getFromMapPtr(subDimEntity);
            if (!global_nodeId_elementOwnerId_ptr)
              {
                std::cout << "M[" << m_iM << "] SerializeNodeRegistry::lookupAndSetNewNodeIds couldn't find subDimEntity= " << subDimEntity;
                for (unsigned kk=0; kk < subDimEntity.size(); kk++)
                  {
                    std::cout << " [" << localNR.getMesh().get_bulk_data()->identifier(subDimEntity[kk]) << "] ";
                  }
                std::cout << std::endl;
              }
            else
              {
              }
            VERIFY_OP_ON(global_nodeId_elementOwnerId_ptr, !=, 0, "SerializeNodeRegistry::lookupAndSetNewNodeIds couldn't find subDimEntity");

            NodeIdsOnSubDimEntityType& global_nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(*global_nodeId_elementOwnerId_ptr);
            NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
            unsigned global_nnodes = global_nodeIds_onSE.size();
            unsigned nnodes = nodeIds_onSE.size();
            VERIFY_OP_ON(global_nnodes, ==, nnodes, "SerializeNodeRegistry::lookupAndSetNewNodeIds: mismatch in nnodes");

            for (unsigned inode=0; inode < nnodes; inode++)
              {
                stk::mesh::EntityId id_new = global_nodeIds_onSE.m_entity_id_vector[inode];
                stk::mesh::EntityId id_old = nodeIds_onSE.m_entity_id_vector[inode];
                stk::mesh::Entity tmp_global_node = globalNR.getMesh().get_bulk_data()->get_entity(stk::topology::NODE_RANK, id_new);
                stk::mesh::Entity local_node_to_change = m_eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, id_old);
                if (m_debug) std::cout << "m_iM= " << m_iM << " id_new= " << id_new << " id_old= " << id_old << std::endl;
                VERIFY_OP_ON(local_node_to_change, !=, stk::mesh::Entity(), "SerializeNodeRegistry::lookupAndSetNewNodeIds null local_node_to_change");
                VERIFY_OP_ON(tmp_global_node, !=, stk::mesh::Entity(), "SerializeNodeRegistry::lookupAndSetNewNodeIds null tmp_global_node");
                VERIFY_OP_ON(tmp_global_node, ==, global_nodeIds_onSE[inode], "SerializeNodeRegistry::lookupAndSetNewNodeIds new node mistmatch");
                if (id_new != id_old)
                  {
                    if (m_debug) std::cout << "DIFF m_iM= " << m_iM << " id_new= " << id_new << " id_old= " << id_old << std::endl;
                    stk::mesh::Entity local_node_new_id_check = m_eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, id_new);
                    if (!m_eMesh.is_valid(local_node_new_id_check))
                      {
                        m_eMesh.get_bulk_data()->change_entity_id(id_new, local_node_to_change);
                      }
                  }
              }
          }
        stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
        m_eMesh.get_bulk_data()->modification_end();
      }

      void writeNodeRegistry(NodeRegistry& nodeRegistry, std::string filename)
      {
        YAML::Emitter yaml;
        if (m_debug) std::cout << "\nnodeRegistry.serialize_write(yaml) to file= " << filename << std::endl;
        // old way was to not use the node map to filter and so we wrote the whole interior and boundary map
        const bool do_filter_for_shared_nodes = true;
        if (do_filter_for_shared_nodes)
          serialize_write(nodeRegistry, yaml, m_nodeMap);
        else
          serialize_write(nodeRegistry, yaml, 0);
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
        if (m_debug) std::cout << "\nnodeRegistry.serialize_read() from file= " << filename << std::endl;
        std::ifstream file(filename.c_str());
        if (!file.is_open() || !file.good())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::readNodeRegistry couldn't open file ")+filename);
          }

        serialize_read(nodeRegistry, file);
      }

      /// clones {KEY,VALUE} into {key,value} and puts it in globalNR
      void addKeyValuePair(NodeRegistry& globalNR, NodeRegistry& localNR, const SubDimCell_SDCEntityType& KEY,
                           SubDimCellData& VALUE )
      {
        bool force_have_node=false;
        PerceptMesh& eMesh = globalNR.getMesh();

        typedef stk::mesh::EntityId SDCEntityType_ID;
        SDCEntityType_ID key_nodeId;
        SubDimCell_SDCEntityType key(&eMesh); // subDimEntity = (*iter).first;
        stk::mesh::EntityId value_nodeId;

        SubDimCellData value; // nodeId_elementOwnerId = (*iter).second;
        NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(value);
        nodeIds_onSE.resize(0);
        stk::mesh::EntityKey& value_entity_key = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(value);
        for(unsigned ikeyd=0; ikeyd < KEY.size(); ++ikeyd) 
          {
            key_nodeId = localNR.getMesh().identifier(KEY[ikeyd]);
            SDCEntityType node = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, key_nodeId);
            if (!eMesh.is_valid(node))
              {
                if (force_have_node)
                  throw std::runtime_error("SerializeNodeRegistry::addKeyValuePair: null node returned from get_entity");
                else
                  {
                    stk::mesh::PartVector parts(1, &eMesh.get_fem_meta_data()->universal_part());
                    node = eMesh.get_bulk_data()->declare_node(static_cast<stk::mesh::EntityId>(key_nodeId), parts);
                  }
              }
            key.insert( node );
          }

        stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>((UInt)std::get<SDC_DATA_OWNING_ELEMENT_KEY>(VALUE).rank());
        size_t id = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(VALUE).id();

        value_entity_key = stk::mesh::EntityKey(rank,id);
        NodeIdsOnSubDimEntityType& NODEIDS_ONSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(VALUE);
        value_nodeId = NODEIDS_ONSE.m_entity_id_vector[0];
        nodeIds_onSE.m_entity_id_vector.push_back(value_nodeId);
        stk::mesh::Entity entity = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, value_nodeId);
        if (!eMesh.is_valid(entity))
          {
            if (force_have_node)
              throw std::runtime_error("SerializeNodeRegistry::addKeyValuePair_read: null node returned from get_entity 2");
            else
              {
                stk::mesh::PartVector parts(1, &eMesh.get_fem_meta_data()->universal_part());
                entity = eMesh.get_bulk_data()->declare_node(static_cast<stk::mesh::EntityId>(value_nodeId),
                                                               parts);
              }
          }

        nodeIds_onSE.push_back(stk::mesh::Entity(entity));

        int sz = nodeIds_onSE.size();
        int sz1 = nodeIds_onSE.m_entity_id_vector.size();
        VERIFY_OP_ON(sz, ==, sz1, "SerializeNodeRegistry::addKeyValuePair sz");
        stk::mesh::EntityId id0 = nodeIds_onSE.m_entity_id_vector[sz-1];
        stk::mesh::EntityId id_check = globalNR.getMesh().identifier(nodeIds_onSE[sz-1]);
        VERIFY_OP_ON(id_check, ==, id0, "SerializeNodeRegistry::addKeyValuePair_read id");

        //if (DEBUG_YAML) std::cout << "s_r value_nodeId= " << value_nodeId << " entity= " << entity
        //                                                <<  " len0= " << nodeIds_onSE.size() << " len1= " << nodeIds_onSE.m_entity_id_vector.size() << std::endl;
        globalNR.getMap()[key] = value;
      }

      void addLocalNodeRegistryToGlobal(NodeRegistry& newLocalNR, NodeRegistry& globalNR)
      {
        globalNR.getMesh().get_bulk_data()->modification_begin();

        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& localMap = newLocalNR.getMap();
        SubDimCellToDataMap& globalMap = globalNR.getMap();
        //std::cout << " tmp SerializeNodeRegistry::processNodeRegistry localMap size: " << localMap.size() << std::endl;

        // key.serialized = { nodeid_0,... : set<EntityId> }
        // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }

        for (iter = localMap.begin(); iter != localMap.end(); ++iter)
          {
            const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnerId = (*iter).second;
            if (m_debug)
              {
                std::cout << "SerializeNodeRegistry::processNodeRegistry inserting localMap entry = " << subDimEntity ;
                for (unsigned kk=0; kk < subDimEntity.size(); kk++)
                  {
                    std::cout << " [" << newLocalNR.getMesh().get_bulk_data()->identifier(subDimEntity[kk]) << "] ";
                  }
                std::cout << " data= " << nodeId_elementOwnerId << " nid=" << std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId).m_entity_id_vector[0] << std::endl;
              }
            /// clone subDimEntity...
            //globalMap[subDimEntity] = nodeId_elementOwnerId;
            addKeyValuePair(globalNR, newLocalNR, subDimEntity, nodeId_elementOwnerId);
          }
        if (m_debug) std::cout << "SerializeNodeRegistry::processNodeRegistry globalMap size= " << globalMap.size() << std::endl;
        stk::mesh::fixup_ghosted_to_shared_nodes(*globalNR.getMesh().get_bulk_data());
        globalNR.getMesh().get_bulk_data()->modification_end();
      }

      ///Note: we read and write to a file instead of storing in memory to allow adapt_exe to be called on a
      //   single piece of the partition - useful for testing and for potentially parallelizing the global outer loop
      void getCurrentGlobalMaxId()
      {
        getCurrentGlobalMaxId(m_globalIdFile, m_id_max, false);
      }

      void getCurrentGlobalMaxId(std::string fileName, IdVector& id_max, bool merge=false)
      {
        // get max node, element id, dump to global file, if global file exists, read from it, use those as start values
        std::fstream file;
        file.open(fileName.c_str(), std::ios_base::in);
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::getCurrentGlobalMaxId couldn't open file ")+fileName);
          }

        //YAML::Parser parser(file);
        YAML::Node doc;
        stk::mesh::EntityId id_max_local;
        try {
          //while(parser.GetNextDocument(doc)) {
          if (1) {
            doc = YAML::Load(file);
            if (m_debug) std::cout << "\n read doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (doc.Type() == YAML::NodeType::Map)
              {
                for (unsigned irank=0; irank < m_id_max.size(); irank++)
                  {
                    id_max_local = doc[m_entity_rank_names[irank]].as<stk::mesh::EntityId>();
                    if (merge)
                      {
                        id_max[irank] = std::max(id_max_local, id_max[irank]);
                      }
                    else
                      {
                        id_max[irank] = id_max_local;
                      }
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
        setCurrentGlobalMaxId(m_globalIdFile, m_id_max);
      }

      void setCurrentGlobalMaxId(std::string fileName, IdVector& id_max)
      {
        std::fstream file;
        file.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc);
        if (!file.is_open())
          {
            throw std::runtime_error(std::string("SerializeNodeRegistry::setCurrentGlobalMaxId couldn't open file ")+fileName);
          }

        YAML::Emitter out;
        out << YAML::BeginMap;
        for (unsigned irank=0; irank < id_max.size(); irank++)
          {
            out << YAML::Key << m_entity_rank_names[irank] << YAML::Value << id_max[irank];
          }
        out << YAML::EndMap;
        file << out.c_str();
        file.close();
      }

      void findCurrentMaxId(PerceptMesh& eMesh)
      {
        for (stk::mesh::EntityRank irank=stk::topology::NODE_RANK; irank < static_cast<stk::mesh::EntityRank>(m_id_max.size()); irank++)
          {
            // bucket loop
            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( irank );

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (selector(**k))
                {
                  stk::mesh::Bucket & bucket = **k ;
                  const unsigned num_entities_in_bucket = bucket.size();

                  for (unsigned iEntity = 0; iEntity < num_entities_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity entity = bucket[iEntity];
                      stk::mesh::EntityId id = eMesh.identifier(entity);
                      m_id_max[irank] = std::max(m_id_max[irank], id);
                    }
                }
              }
          }
      }

      // use id_max values as a simple id server to reset new element ids
      void resetNewElementIds(PerceptMesh& eMesh, NodeRegistry& /*nodeRegistry*/)
      {

        typedef std::pair<stk::mesh::EntityId, stk::mesh::EntityId> EntityPair;
        for (stk::mesh::EntityRank irank=stk::topology::EDGE_RANK; irank < static_cast<stk::mesh::EntityRank>(m_id_max.size()); irank++)
          {
            if (irank == FAMILY_TREE_RANK) continue;
            // bucket loop
            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( irank );
            std::vector<EntityPair > id_change;

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (selector(**k))
                {
                  stk::mesh::Bucket & bucket = **k ;
                  const unsigned num_entities_in_bucket = bucket.size();

                  for (unsigned iEntity = 0; iEntity < num_entities_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity entity = bucket[iEntity];
                      //! FIXME
                      //if (eMesh.hasFamilyTree(entity) && eMesh.isChildElement(entity))
                        {
                          stk::mesh::EntityId id = m_id_max[irank] + 1;
                          m_id_max[irank] = id;

                          stk::mesh::EntityId id_old = eMesh.identifier(entity);
                          id_change.push_back(EntityPair(id_old, id));
                        }
                    }
                }
              }

            eMesh.get_bulk_data()->modification_begin();

            for (unsigned ii=0; ii< id_change.size(); ii++)
              {
                stk::mesh::Entity entity = eMesh.get_bulk_data()->get_entity(irank, id_change[ii].first);
                VERIFY_OP_ON(entity, !=, stk::mesh::Entity(), "SerializeNodeRegistry::resetNewElementIds");
                VERIFY_OP_ON(eMesh.identifier(entity), ==, id_change[ii].first, "SerializeNodeRegistry::resetNewElementIds bad00");
                stk::mesh::EntityId id_new = id_change[ii].second;
                //stk::mesh::EntityId id_old = entity->identifier();
                eMesh.get_bulk_data()->change_entity_id(id_new, entity);
              }
          }
      }

      void checkNR(NodeRegistry& nodeRegistry, std::string msg="")
      {
        nodeRegistry.getMesh().get_bulk_data()->modification_begin();

        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = nodeRegistry.getMap();
        std::cout << " tmp SerializeNodeRegistry::checkNR msg= " << msg << " map size: " << map.size() << std::endl;

        if (1)
          {
            for (iter = map.begin(); iter != map.end(); ++iter)
              {
                //const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
                SubDimCellData& nodeId_elementOwnerId = (*iter).second;
                NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
                unsigned nnodes = nodeIds_onSE.size();
                for (unsigned inode=0; inode < nnodes; inode++)
                  {
                    VERIFY_OP_ON(nodeIds_onSE[inode], !=, stk::mesh::Entity(), "SerializeNodeRegistry::checkNR node is null 0");
                    stk::mesh::EntityId id_old = nodeIds_onSE.m_entity_id_vector[inode];
                    VERIFY_OP_ON(id_old, !=, 0, "SerializeNodeRegistry::checkNR node id is  0");
                    stk::mesh::EntityId id_old_check = nodeRegistry.getMesh().get_bulk_data()->identifier(nodeIds_onSE[inode]);
                    if (nodeRegistry.getMesh().get_bulk_data()->entity_rank(nodeIds_onSE[inode]) != 0)
                      std::cout <<  "SerializeNodeRegistry checkNR 1" << nodeRegistry.getMesh().get_bulk_data()->identifier(nodeIds_onSE[inode]) << " " << msg << std::endl;

                    VERIFY_OP_ON(nodeRegistry.getMesh().get_bulk_data()->entity_rank(nodeIds_onSE[inode]), ==, 0, "SerializeNodeRegistry checkNR 1");
                    if (id_old_check != id_old)
                      std::cout <<  "SerializeNodeRegistry checkNR 2" << nodeRegistry.getMesh().get_bulk_data()->identifier(nodeIds_onSE[inode]) << " " << msg << std::endl;

                    VERIFY_OP_ON(id_old_check, ==, id_old, "SerializeNodeRegistry::checkNR id_old");

                    stk::mesh::Entity node = nodeRegistry.getMesh().get_bulk_data()->get_entity(stk::topology::NODE_RANK, id_old);
                    VERIFY_OP_ON(node, !=, stk::mesh::Entity(), "SerializeNodeRegistry::checkNR node is null");
                    VERIFY_OP_ON(node, ==, nodeIds_onSE[inode], "SerializeNodeRegistry::checkNR node is not same");
                  }
              }
          }
        std::cout << " tmp SerializeNodeRegistry::checkNR done msg= " << msg << " map size: " << map.size() << std::endl;
        stk::mesh::fixup_ghosted_to_shared_nodes(*nodeRegistry.getMesh().get_bulk_data());
        nodeRegistry.getMesh().get_bulk_data()->modification_end();
      }

      void fixNR(NodeRegistry& nodeRegistry, std::string msg="")
      {
        nodeRegistry.getMesh().get_bulk_data()->modification_begin();
        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = nodeRegistry.getMap();
        std::cout << " tmp SerializeNodeRegistry::fixNR msg= " << msg << " map size: " << map.size() << std::endl;

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            //const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnerId = (*iter).second;
            NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
            unsigned nnodes = nodeIds_onSE.size();
            for (unsigned inode=0; inode < nnodes; inode++)
              {
                stk::mesh::EntityId id = nodeIds_onSE.m_entity_id_vector[inode];
                stk::mesh::Entity node = nodeRegistry.getMesh().get_bulk_data()->get_entity(stk::topology::NODE_RANK, id);
                VERIFY_OP_ON(node, !=, stk::mesh::Entity(), "SerializeNodeRegistry::fixNR node is null "+toString(id));
                nodeIds_onSE[inode] = node;
              }
          }
        stk::mesh::fixup_ghosted_to_shared_nodes(*nodeRegistry.getMesh().get_bulk_data());
        nodeRegistry.getMesh().get_bulk_data()->modification_end();
        std::cout << " tmp SerializeNodeRegistry::fixNR done msg= " << msg << " map size: " << map.size() << std::endl;
      }

      /// Using the current global max id, in a simple id-server manner, generate new id's and assign to
      ///   the shared nodes.
      void resetNewNodeIds(NodeRegistry& nodeRegistry)
      {
        nodeRegistry.getMesh().get_bulk_data()->modification_begin();

        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = nodeRegistry.getMap();
        std::cout << " tmp SerializeNodeRegistry::resetNewNodeIds map size: " << map.size() << std::endl;

        if (1)
          {
            for (iter = map.begin(); iter != map.end(); ++iter)
              {
                //const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
                SubDimCellData& nodeId_elementOwnerId = (*iter).second;
                NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
                unsigned nnodes = nodeIds_onSE.size();
                for (unsigned inode=0; inode < nnodes; inode++)
                  {
                    stk::mesh::EntityId id_old = nodeIds_onSE.m_entity_id_vector[inode];
                    stk::mesh::EntityId id_old_check = nodeRegistry.getMesh().get_bulk_data()->identifier(nodeIds_onSE[inode]);
                    VERIFY_OP_ON(id_old_check, ==, id_old, "SerializeNodeRegistry::resetNewNodeIds id_old");
                  }
              }
          }
        // key.serialized = { nodeid_0,... : set<EntityId> }
        // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }

        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            //const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnerId = (*iter).second;
            NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnerId);
            unsigned nnodes = nodeIds_onSE.size();
            for (unsigned inode=0; inode < nnodes; inode++)
              {
                stk::mesh::EntityId id_old = nodeIds_onSE.m_entity_id_vector[inode];
                stk::mesh::EntityId id_old_check = nodeRegistry.getMesh().get_bulk_data()->identifier(nodeIds_onSE[inode]);
                VERIFY_OP_ON(id_old_check, ==, id_old, "SerializeNodeRegistry::resetNewNodeIds id_old");

                stk::mesh::EntityId id_new = m_id_max[0]+1;
                m_id_max[0] = id_new;
                nodeIds_onSE.m_entity_id_vector[inode] = id_new;

                nodeRegistry.getMesh().get_bulk_data()->change_entity_id(id_new, nodeIds_onSE[inode]);

                stk::mesh::EntityId id_new_check = nodeRegistry.getMesh().get_bulk_data()->identifier(nodeIds_onSE[inode]);
                VERIFY_OP_ON(id_new_check, ==, id_new, "SerializeNodeRegistry::resetNewNodeIds id_new");
                id_new = nodeIds_onSE.m_entity_id_vector[inode];

                VERIFY_OP_ON(id_new_check, ==, id_new, "SerializeNodeRegistry::resetNewNodeIds id_new 2");
              }
          }
        stk::mesh::fixup_ghosted_to_shared_nodes(*nodeRegistry.getMesh().get_bulk_data());
        nodeRegistry.getMesh().get_bulk_data()->modification_end();
      }


      static void serialize_write(NodeRegistry& nodeRegistry, YAML::Emitter& emitter, NodeMap *nodeMapFilter = 0, std::string /*msg*/="")
      {
        SubDimCellToDataMap::iterator iter;
        SubDimCellToDataMap& map = nodeRegistry.getMap();
        //std::cout << msg << " tmp serialize_write map size: " << map.size() << std::endl;

        //if (0) emitter << YAML::Anchor("NodeRegistry::map");   YAML_ERRCHECK;
        //emitter << YAML::Flow;      YAML_ERRCHECK;
        emitter << YAML::BeginMap;      YAML_ERRCHECK;

        // key.serialized = { nodeid_0,... : set<EntityId> }
        // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }

        //int jj=0;
        for (iter = map.begin(); iter != map.end(); ++iter)
          {
            const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnerId = (*iter).second;

            // check if all nodes are on the boundary (defined by shared nodes in the NodeMap)
            if (nodeMapFilter)
              {
                // special case for "interior" subDimEntity's which are centroid nodes for quad or hex elements -
                //   by definition they aren't shared
                if (subDimEntity.size() == 1)
                  continue;

                bool notFound = false;
                for (unsigned k=0; k < subDimEntity.size(); k++)
                  {
                    //std::cout << " " << subDimEntity[k]->identifier() << " ";
                    NodeMap::iterator filter_iter = nodeMapFilter->find(nodeRegistry.getMesh().identifier(subDimEntity[k]));
                    if (filter_iter == nodeMapFilter->end())
                      {
                        notFound = true;
                        break;
                      }
                  }
                if (notFound)
                  continue;
              }

            //emitter << YAML::Key << subDimEntity;
            emitter << YAML::Key;       YAML_ERRCHECK;
            emitter << YAML::Flow;      YAML_ERRCHECK;
            emitter << YAML::BeginSeq;  YAML_ERRCHECK;
            //if (0) emitter << YAML::Anchor(std::string("seq")+std::to_string(jj++));   YAML_ERRCHECK;
            for (unsigned k=0; k < subDimEntity.size(); k++)
              {
                //std::cout << " " << subDimEntity[k]->identifier() << " ";
                emitter << nodeRegistry.getMesh().identifier(subDimEntity[k]);          YAML_ERRCHECK;
              }
            emitter << YAML::EndSeq;      YAML_ERRCHECK;

            NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>( nodeId_elementOwnerId);
            stk::mesh::EntityKey& value_entity_key = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnerId);

            //emitter << YAML::Scalar << nodeIds_onSE.size()
            emitter << YAML::Value;          YAML_ERRCHECK;
            emitter << YAML::Flow;           YAML_ERRCHECK;
            emitter << YAML::BeginSeq;       YAML_ERRCHECK;
            emitter << value_entity_key.rank();
            emitter << value_entity_key.id();
            for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
              {
                //emitter << (int)nodeIds_onSE[ii]->identifier();      YAML_ERRCHECK;
                stk::mesh::EntityId id = nodeIds_onSE.m_entity_id_vector[ii];
                stk::mesh::EntityId id_check = nodeRegistry.getMesh().identifier(nodeIds_onSE[ii]);
                VERIFY_OP_ON(id_check, ==, id, "SerializeNodeRegistry::serialize_write id");

                emitter << (int)nodeIds_onSE.m_entity_id_vector[ii];      YAML_ERRCHECK;
              }
            emitter << YAML::EndSeq;     YAML_ERRCHECK;
          }

        emitter << YAML::EndMap;    YAML_ERRCHECK;

      }

      static void serialize_read(NodeRegistry& nodeRegistry, std::ifstream& file_in,  std::string /*msg*/="", bool force_have_node=false)
      {
        PerceptMesh& eMesh = nodeRegistry.getMesh();
        eMesh.get_bulk_data()->modification_begin();

        //YAML::Parser parser(file_in);
        YAML::Node doc;
        if (DEBUG_YAML)
          std::cout
            << "\n serialize_read..."
            << " YAML::NodeType::Null= " <<  YAML::NodeType::Null
            << " YAML::NodeType::Scalar= " <<  YAML::NodeType::Scalar
            << " YAML::NodeType::Sequence= " <<  YAML::NodeType::Sequence
            << " YAML::NodeType::Map= " <<  YAML::NodeType::Map
            << std::endl;

        SubDimCellToDataMap& map = nodeRegistry.getMap();
        //std::cout << msg << " tmp serialize_read map size: " << map.size() << std::endl;

        try {
          //while(parser.GetNextDocument(doc)) {
          if (1) {
            doc = YAML::Load(file_in);
            if (DEBUG_YAML) std::cout << "s_r doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (doc.Type() == YAML::NodeType::Map)
              {
                for(YAML::const_iterator it=doc.begin();it!=doc.end();++it) {
                  typedef stk::mesh::EntityId SDCEntityType_ID;
                  //typedef stk::mesh::Entity SDCEntityType;
                  SDCEntityType_ID key_quantum;
                  //typedef SubDimCell<SDCEntityType> SubDimCell_SDCEntityType;
                  SubDimCell_SDCEntityType key(&eMesh); // subDimEntity = (*iter).first;

                  //struct NodeIdsOnSubDimEntityType : public std::vector<stk::mesh::Entity>
                  // {
                  //     typedef IdVector entity_id_vector_type;

                  stk::mesh::EntityId value_tuple_0_quantum;
                  NodeIdsOnSubDimEntityType value_tuple_0;

                  SubDimCellData value; // nodeId_elementOwnerId = (*iter).second;
                  NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(value);
                  nodeIds_onSE.resize(0);
                  stk::mesh::EntityKey& value_entity_key = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(value);
                  // value = { {new_node0, new_node1,...}:[vector<Entity>,vector<EntityId>], {elem_own[rank, ele_id]:EntityKey} }
                  // key = { nodePtr_0,... : set<Entity> }
                  // value.serialized = { {new_nid0, new_nid1,...}:vector<EntityId>, {elem_own[rank, ele_id]:EntityKey} }
                  // key.serialized = { nodeid_0,... : set<EntityId> }


                  //if (DEBUG_YAML) std::cout << "it.first().Type() = " << it.first().Type() << " it.first().Tag()= " << it.first().Tag() << std::endl;
                  //if (DEBUG_YAML) std::cout << "it.second().Type() = " << it.second().Type() << " it.second().Tag()= " << it.second().Tag() << std::endl;
                  const YAML::Node keySeq = it->first;
                  for(YAML::const_iterator itk=keySeq.begin();itk!=keySeq.end();++itk) {
                    key_quantum = itk->as<SDCEntityType_ID>();
                    if (DEBUG_YAML) std::cout << "s_r key_quantum= " << key_quantum << std::endl;
                    SDCEntityType node = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, key_quantum);
                    //key.insert(const_cast<stk::mesh::Entity>(&element) );
                    if (!eMesh.is_valid(node))
                      {
                        if (force_have_node)
                          throw std::runtime_error("NodeRegistry::serialize_read: null node returned from get_entity");
                        else
                          {
                            stk::mesh::PartVector parts(1, &eMesh.get_fem_meta_data()->universal_part());
                            node = eMesh.get_bulk_data()->declare_node(static_cast<stk::mesh::EntityId>(key_quantum), parts);
                          }
                      }

                    key.insert( node );

                  }

                  int iseq=0;
                  const YAML::Node valSeq = it->second;
                  stk::mesh::EntityRank rank = stk::topology::INVALID_RANK;
                  size_t id;
                  for(YAML::const_iterator itv=valSeq.begin();itv!=valSeq.end();++itv,++iseq) {
                    if (iseq == 0)
                      {
                        UInt rank_input = itv->as<UInt>();
                        rank = static_cast<stk::mesh::EntityRank>(rank_input);
                      }
                    else if (iseq == 1)
                      {
                        id = itv->as<size_t>();
                        stk::mesh::EntityKey entityKey(rank,id);
                        if (DEBUG_YAML) std::cout << "s_r value_tuple_1= " << rank << " " << id << std::endl;
                        value_entity_key = stk::mesh::EntityKey(rank,id);
                        if (DEBUG_YAML) std::cout << "s_r owning element rank= " << std::get<SDC_DATA_OWNING_ELEMENT_KEY>(value).rank()
                                                  << " owning element id= " << std::get<SDC_DATA_OWNING_ELEMENT_KEY>(value).id()
                                                  << std::endl;
                      }
                    else
                      {
                        value_tuple_0_quantum = itv->as<stk::mesh::EntityId>();

                        nodeIds_onSE.m_entity_id_vector.push_back(value_tuple_0_quantum);
                        stk::mesh::Entity entity = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, value_tuple_0_quantum);
                        if (!eMesh.is_valid(entity))
                          {
                            if (force_have_node)
                              throw std::runtime_error("NodeRegistry::serialize_read: null node returned from get_entity 2");
                            else
                              {
                                stk::mesh::PartVector parts(1, &eMesh.get_fem_meta_data()->universal_part());
                                entity = eMesh.get_bulk_data()->declare_node(static_cast<stk::mesh::EntityId>(value_tuple_0_quantum),
                                                                              parts);
                              }
                          }

                        nodeIds_onSE.push_back(entity);

                        int sz = nodeIds_onSE.size();
                        int sz1 = nodeIds_onSE.m_entity_id_vector.size();
                        VERIFY_OP_ON(sz, ==, sz1, "SerializeNodeRegistry::serialize_read sz");
                        stk::mesh::EntityId id0 = nodeIds_onSE.m_entity_id_vector[sz-1];
                        stk::mesh::EntityId id_check = nodeRegistry.getMesh().identifier(nodeIds_onSE[sz-1]);
                        VERIFY_OP_ON(id_check, ==, id0, "SerializeNodeRegistry::serialize_read id");

                        if (DEBUG_YAML) std::cout << "s_r value_tuple_0_quantum= " << value_tuple_0_quantum << " entity= " << entity
                                                  <<  " len0= " << nodeIds_onSE.size() << " len1= " << nodeIds_onSE.m_entity_id_vector.size() << std::endl;
                      }
                  }

                  map[key] = value;
                }
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << "\n";
          throw std::runtime_error(std::string("yaml parsing error: ")+e.what());
        }

        stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
        eMesh.get_bulk_data()->modification_end();
      }


    };

  }

#endif // HAVE_YAML
